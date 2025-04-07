from functools import partial
from jinja2 import Environment
from json import dumps
import logging
from os import listdir
from os.path import join, basename, exists, sep, split
from sequence_processing_pipeline.Job import Job, KISSLoader
from sequence_processing_pipeline.PipelineError import (PipelineError,
                                                        JobFailedError)
from sequence_processing_pipeline.util import determine_orientation
from re import sub


class MultiQCJob(Job):
    def __init__(self, run_dir, output_path, raw_fastq_files_path,
                 processed_fastq_files_path, nprocs, nthreads, multiqc_path,
                 modules_to_load, qiita_job_id, queue_name, node_count,
                 wall_time_limit, jmem, pool_size, fastqc_root_path,
                 max_array_length, multiqc_config_file_path, is_amplicon):
        super().__init__(run_dir,
                         output_path,
                         'MultiQCJob',
                         [multiqc_path],
                         max_array_length,
                         modules_to_load=modules_to_load)

        self.nprocs = nprocs
        self.nthreads = nthreads
        self.multiqc_path = multiqc_path
        self.queue_name = queue_name
        self.node_count = node_count
        self.wall_time_limit = wall_time_limit
        self.jmem = jmem
        self.qiita_job_id = qiita_job_id
        self.pool_size = pool_size
        self.raw_fastq_files_path = raw_fastq_files_path
        self.processed_fastq_files_path = processed_fastq_files_path
        self.multiqc_config_file_path = multiqc_config_file_path
        self.is_amplicon = is_amplicon
        self.fastqc_root_path = fastqc_root_path

        self.job_script_path = join(self.output_path, f"{self.job_name}.sh")

        # for projects that use sequence_processing_pipeline as a dependency,
        # jinja_env must be set to sequence_processing_pipeline's root path,
        # rather than the project's root path.
        self.jinja_env = Environment(loader=KISSLoader('templates'))

        # bypass generating job script for a force-fail job, since it is
        # not needed.
        if not self.force_job_fail:
            self._generate_job_script()

    def _find_projects(self):
        find_paths = [self.processed_fastq_files_path]

        if not self.is_amplicon:
            # avoid processing the raw fastq files for amplicon runs because
            # they are identical to files in self.processed_fastq_files_path.
            find_paths += [self.raw_fastq_files_path]

        projects = []

        for fastq_files_path in find_paths:
            for directory in listdir(fastq_files_path):
                # confirm that this directory has data we want to show to
                # multiqc.

                # generate a list of all files in this directory.
                files = self._find_files(join(fastq_files_path, directory))

                # filter out all files that aren't fastq.gz files.
                files = [x for x in files if x.endswith('.fastq.gz')]

                for _file in files:
                    # split path into a list of folder names and the filename.
                    # filter out the contents of any folders that we don't
                    # want included in the report.
                    file_path, file_name = split(_file)

                    ignore_this_file = [x for x in file_path.split(sep)
                                        if x in ['zero_files',
                                                 'only-adapter-filtered']]

                    if ignore_this_file:
                        # if one or more of the folders are present in _file's
                        # path, then do not consider this file.
                        continue

                    # lastly, only consider folders that contain at least one
                    # R1 file.
                    if determine_orientation(file_name) != 'R1':
                        continue

                    # according to legacy behavior, if _file has met the above
                    # criteria, then add the value of directory as a project
                    # name.
                    projects.append(directory)

        if projects:
            # remove duplicates
            return sorted(set(projects))

        raise PipelineError("There are no fastq files for MultiQCJob to "
                            "process")

    def _get_failed_indexes(self, job_id):
        completed_files = self._find_files(self.output_path)

        # remove path and .completed extension from file-name. e.g.:
        # 'project1_0', 'project1_1', ..., 'project1_n'
        completed_files = [sub(r'\.completed$', '', basename(fp)) for fp in
                           completed_files if fp.endswith('.completed')]

        # extract the line number in the .detailed file corresponding to
        # the command used for this job
        completed_indexes = [int(cf.split('_')[-1]) for cf in completed_files]

        all_indexes = list(range(1, len(self.array_cmds) + 1))
        failed_indexes = sorted(set(all_indexes) - set(completed_indexes))

        # generate log-file here instead of in run() where it can be
        # unittested more easily.
        if failed_indexes:
            with open(join(self.output_path, 'logs',
                           f'failed_indexes_{job_id}.json'), 'w') as f:
                f.write(dumps({'job_id': job_id,
                               'failed_indexes': failed_indexes}, indent=2))

        return failed_indexes

    def _get_commands(self):
        # If project-level reports were not needed, MultiQC could simply be
        # given the path to the run-directory itself and it will discover all
        # of the relevant data files. Confirmed that for a one-project sample-
        # sheet, this produces an equivalent report.

        array_cmds = []

        for project in self._find_projects():
            # MultiQC doesn't like input paths that don't exist. Simply add
            # all paths that do exist as input.
            input_path_list = []

            p_path = partial(join, self.fastqc_root_path, 'fastqc')

            for filter_type in ['bclconvert', 'trimmed_sequences',
                                'filtered_sequences', 'amplicon']:
                input_path_list.append(p_path(project, filter_type))

            input_path_list.append(p_path(project, 'Reports'))

            p_path = partial(join, self.processed_fastq_files_path, project)
            input_path_list.append(p_path('fastp_reports_dir', 'json'))

            # I don't usually see a json directory associated with raw data.
            # It looks to be metadata coming directly off the machine, in the
            # few instances I've seen it in /sequencing...
            p_path = partial(join, self.raw_fastq_files_path, project)
            input_path_list.append(p_path('json'))

            input_path_list = [x for x in input_path_list if exists(x)]

            cmd_head = ['multiqc', '-c', self.multiqc_config_file_path,
                        '--fullnames', '--force']

            # --interactive graphs is set to True in MultiQC configuration
            # file and hence this switch was redunant and now removed.
            cmd_tail = ['-o', join(self.output_path, 'multiqc', project)]

            array_cmds.append(' '.join(cmd_head + input_path_list + cmd_tail))

        # These commands are okay to execute in parallel because each command
        # is limited to a specific project and each invocation creates its own
        # multiqc/project output directory so there will not be collisions.
        # These commands must be executed after FastQCJob has completed for
        # FastQC report results to be included, however.
        return array_cmds

    def _generate_job_script(self):
        template = self.jinja_env.get_template("multiqc_job.sh")

        self.array_cmds = self._get_commands()

        job_name = f'{self.qiita_job_id}_{self.job_name}'
        details_file_name = f'{self.job_name}.array-details'
        array_details = join(self.output_path, details_file_name)
        array_params = "1-%d%%%d" % (len(self.array_cmds), self.pool_size)
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
            f.write('\n'.join(self.array_cmds) + '\n')

        return self.job_script_path

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
            raise PipelineError("MultiQCJob did not complete successfully.")
