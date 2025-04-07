from functools import partial
from jinja2 import Environment
from json import dumps
import logging
from os import listdir, makedirs
from os.path import join, basename
from re import sub
from sequence_processing_pipeline.Job import Job, KISSLoader
from sequence_processing_pipeline.PipelineError import (PipelineError,
                                                        JobFailedError)


class FastQCJob(Job):
    def __init__(self, run_dir, output_path, raw_fastq_files_path,
                 processed_fastq_files_path, nprocs, nthreads, fastqc_path,
                 modules_to_load, qiita_job_id, queue_name, node_count,
                 wall_time_limit, jmem, pool_size,
                 max_array_length, is_amplicon):
        super().__init__(run_dir,
                         output_path,
                         'FastQCJob',
                         [fastqc_path],
                         max_array_length,
                         modules_to_load=modules_to_load)

        self.nprocs = nprocs
        self.nthreads = nthreads
        self.fastqc_path = fastqc_path
        self.queue_name = queue_name
        self.node_count = node_count
        self.wall_time_limit = wall_time_limit
        self.jmem = jmem
        self.qiita_job_id = qiita_job_id
        self.pool_size = pool_size
        self.raw_fastq_files_path = raw_fastq_files_path
        self.processed_fastq_files_path = processed_fastq_files_path
        self.is_amplicon = is_amplicon

        self.job_script_path = join(self.output_path, f"{self.job_name}.sh")

        self.commands, self.project_names = self._get_commands()
        # for lists greater than n commands, chain the extra commands,
        # distributing them evenly throughout the first n commands.
        self.commands = self._group_commands(self.commands)
        self.suffix = 'fastqc.html'

        # for projects that use sequence_processing_pipeline as a dependency,
        # jinja_env must be set to sequence_processing_pipeline's root path,
        # rather than the project's root path.
        self.jinja_env = Environment(loader=KISSLoader('templates'))

        self._generate_job_script()

    def _get_commands(self):
        """
        Generate a set of commands to execute, based on the input metadata.
        :return: A list of commands to execute w/in a job script
        """
        results = []

        # gather the parameters for processing all relevant raw fastq
        # files.
        params, project_names = self._scan_fastq_files(True)
        for fwd_file_path, rev_file_path, output_path in params:
            command = ['fastqc', '--noextract', '-t', str(self.nthreads),
                       fwd_file_path, rev_file_path, '-o', output_path]
            results.append(' '.join(command))

        if not self.is_amplicon:
            # next, do the same for the trimmed/filtered fastq files.
            params, additional_project_names = self._scan_fastq_files(False)
            for fwd_file_path, rev_file_path, output_path in params:
                command = ['fastqc', '--noextract', '-t', str(self.nthreads),
                           fwd_file_path, rev_file_path, '-o', output_path]
                results.append(' '.join(command))
            # remove duplicate project names from the list
            project_names = list(set(project_names + additional_project_names))

        return results, project_names

    def _find_projects(self, path_to_run_id_data_fastq_dir, is_raw_input):
        results = []
        for directory in listdir(path_to_run_id_data_fastq_dir):
            project_dir = join(path_to_run_id_data_fastq_dir, directory)
            files = self._find_files(project_dir)

            # extract only fastq files from the list and remove all files
            # in 'zero_files' sub-folder.
            files = [x for x in files if x.endswith('.fastq.gz') and
                     'zero_files' not in x]

            # remove fastq files in the only-adapter-filtered
            # folder from consideration if they are present.
            files = [x for x in files if 'only-adapter-filtered' not in x]

            # break files up into R1, R2, I1, I2
            # assume _R1_ does not occur in the path as well.
            r1_only = [x for x in files if '_R1_' in x]
            r2_only = [x for x in files if '_R2_' in x]

            # amplicon runs may or may not have an i2. this is okay.
            i1_only = [x for x in files if '_I1_' in x]
            i2_only = [x for x in files if '_I2_' in x]

            if r1_only:
                tmp = ' '.join(r1_only)
                if 'trimmed_sequences' in tmp:
                    # a_trim = True, h_filter= = False
                    filter_type = 'trimmed_sequences'
                elif 'filtered_sequences' in tmp:
                    # a_trim = True, h_filter= = True
                    # (trimmed AND filtered)
                    filter_type = 'filtered_sequences'
                elif 'amplicon' in tmp:
                    # a_trim = False, h_filter= = False
                    filter_type = 'amplicon'
                else:
                    if is_raw_input:
                        filter_type = 'raw'
                    else:
                        raise ValueError("indeterminate type")

                r1_only.sort()
                r2_only.sort()

                if filter_type == 'amplicon':
                    i1_only.sort()
                    i2_only.sort()
                    results.append((directory, filter_type, project_dir,
                                    r1_only + i1_only, r2_only + i2_only))

                results.append(
                    (directory, filter_type, project_dir, r1_only, r2_only))

        return results if results else None

    def _scan_fastq_files(self, is_raw_input=False):
        find_path = (self.raw_fastq_files_path if is_raw_input else
                     self.processed_fastq_files_path)

        projects = self._find_projects(find_path, is_raw_input)

        if projects is None:
            raise PipelineError("There are no fastq files for FastQCJob to "
                                "process in %s." % find_path)

        fastqc_results = []

        # rather than get a list of project_names from a sample-sheet,
        # gather a list of names from the outputs of previous jobs.
        project_names = []

        for proj_name, fltr_type, fastq_fp, fwd_files, rev_files in projects:
            project_names.append(proj_name)
            p_path = partial(join, self.output_path, 'fastqc', proj_name)
            output_dir = p_path('bclconvert' if is_raw_input else fltr_type)
            makedirs(output_dir, exist_ok=True)

            for some_fwd_file, some_rev_file in zip(fwd_files, rev_files):
                fastqc_results.append((some_fwd_file, some_rev_file,
                                       output_dir))
        # remove duplicates
        project_names = list(set(project_names))

        return fastqc_results, project_names

    def _get_failed_indexes(self, job_id):
        completed_files = self._find_files(self.output_path)

        # remove path and .completed extension from file-name. e.g.:
        # 'project1_0', 'project1_1', ..., 'project1_n'
        completed_files = [sub(r'\.completed$', '', basename(fp)) for fp in
                           completed_files if fp.endswith('.completed')]

        # extract the line number in the .detailed file corresponding to
        # the command used for this job
        completed_indexes = [int(cf.split('_')[-1]) for cf in completed_files]

        # a successfully completed job array should have a list of array
        # numbers from 0 - len(self.commands).
        all_indexes = [x for x in range(1, len(self.commands) + 1)]
        failed_indexes = list(set(all_indexes) - set(completed_indexes))
        failed_indexes.sort()

        # generate log-file here instead of in run() where it can be
        # unittested more easily.
        log_fp = join(self.output_path,
                      'logs',
                      f'failed_indexes_{job_id}.json')

        if failed_indexes:
            with open(log_fp, 'w') as f:
                f.write(dumps({'job_id': job_id,
                               'failed_indexes': failed_indexes}, indent=2))

        return failed_indexes

    def run(self, callback=None):
        try:
            job_info = self.submit_job(self.job_script_path,
                                       exec_from=self.log_path,
                                       callback=callback)
        except JobFailedError as e:
            # When a job has failed, parse the logs generated by this specific
            # job to return a more descriptive message to the user.
            info = self.parse_logs()
            # prepend just the message component of the Error.
            info.insert(0, str(e))
            raise JobFailedError('\n'.join(info)) from None

        logging.debug(job_info)

        if self._get_failed_indexes(job_info['job_id']):
            # raise error if list isn't empty.
            raise PipelineError("FastQCJob did not complete successfully.")

        self.mark_job_completed()

    def _generate_job_script(self):
        # bypass generating job script for a force-fail job, since it is
        # not needed.
        if self.force_job_fail:
            return None

        template = self.jinja_env.get_template("fastqc_job.sh")

        job_name = f'{self.qiita_job_id}_{self.job_name}'
        details_file_name = f'{self.job_name}.array-details'
        array_details = join(self.output_path, details_file_name)
        array_params = "1-%d%%%d" % (len(self.commands), self.pool_size)
        modules_to_load = ' '.join(self.modules_to_load)

        with open(self.job_script_path, mode="w", encoding="utf-8") as f:
            f.write(template.render(job_name=job_name,
                                    array_details=array_details,
                                    queue_name=self.queue_name,
                                    node_count=self.node_count,
                                    nprocs=self.nprocs,
                                    wall_time_limit=self.wall_time_limit,
                                    mem_in_gb=self.jmem,
                                    array_params=array_params,
                                    output_path=self.output_path,
                                    modules_to_load=modules_to_load))

        # save the .details file as well
        with open(array_details, 'w') as f:
            f.write('\n'.join(self.commands))

        return self.job_script_path
