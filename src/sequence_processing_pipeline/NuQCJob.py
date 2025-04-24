from metapool import load_sample_sheet
from os import stat, makedirs, rename
from os.path import join, basename, dirname, exists, abspath, split
from sequence_processing_pipeline.Job import Job, KISSLoader
from sequence_processing_pipeline.PipelineError import (PipelineError,
                                                        JobFailedError)
from sequence_processing_pipeline.Pipeline import Pipeline
from shutil import move
import logging
from sequence_processing_pipeline.Commands import split_similar_size_bins
from sequence_processing_pipeline.util import iter_paired_files
from jinja2 import Environment
from glob import glob
import re
from sys import executable


logging.basicConfig(level=logging.DEBUG)


class NuQCJob(Job):
    def __init__(self, fastq_root_dir, output_path, sample_sheet_path,
                 minimap_database_paths, queue_name, node_count,
                 wall_time_limit, jmem, fastp_path, minimap2_path,
                 samtools_path, modules_to_load, qiita_job_id,
                 max_array_length, known_adapters_path, movi_path, gres_value,
                 pmls_path, additional_fastq_tags, bucket_size=8,
                 length_limit=100, cores_per_task=4):
        """
        Submit a slurm job where the contents of fastq_root_dir are processed
        using fastp, minimap2, and samtools. Human-genome sequences will be
        filtered out, if needed.
        :param fastq_root_dir: Path to a dir of Fastq files, org. by project.
        :param output_path: Path where all pipeline-generated files live.
        :param sample_sheet_path: Path to a sample sheet file.
        :param minimap_database_paths: Path to human genome databases in env.
        :param queue_name: Torque queue name to use in running env.
        :param node_count: Number of nodes to use in running env.
        :param wall_time_limit: Hard wall-clock-time limit (in min) for procs.
        :param jmem: String representing total memory limit for entire job.
        :param fastp_path: The path to the fastp executable
        :param minimap2_path: The path to the minimap2 executable
        :param samtools_path: The path to the samtools executable
        :param modules_to_load: A list of Linux module names to load
        :param qiita_job_id: identify Torque jobs using qiita_job_id
        :param known_adapters_path: The path to an .fna file of known adapters.
        :param movi_path: The path to the Movi executable.
        :param gres_value: Cluster-related parameter.
        :param pmls_path: The path to the pmls script.
        :param bucket_size: the size in GB of each bucket to process
        :param length_limit: reads shorter than this will be discarded.
        :param cores_per_task: Number of CPU cores per node to request.
        :param additional_fastq_tags: A list of fastq tags to preserve during
        filtering.
        """
        super().__init__(fastq_root_dir,
                         output_path,
                         'NuQCJob',
                         [fastp_path, minimap2_path, samtools_path],
                         max_array_length,
                         modules_to_load=modules_to_load)
        self.sample_sheet_path = sample_sheet_path
        self._file_check(self.sample_sheet_path)
        metadata = self._process_sample_sheet()
        self.sample_ids = metadata['sample_ids']
        self.project_data = metadata['projects']
        self.chemistry = metadata['chemistry']
        self.mmi_file_paths = minimap_database_paths
        self.queue_name = queue_name
        self.node_count = node_count
        self.wall_time_limit = wall_time_limit
        # raise an Error if jmem is not a valid floating point value.
        self.jmem = str(int(jmem))
        self.fastp_path = fastp_path
        self.minimap2_path = minimap2_path
        self.samtools_path = samtools_path
        self.qiita_job_id = qiita_job_id
        self.suffix = 'fastq.gz'
        self.movi_path = movi_path
        self.gres_value = gres_value
        self.pmls_path = pmls_path
        self.additional_fastq_tags = additional_fastq_tags
        self.audit_folders = ['filtered_sequences']

        # for projects that use sequence_processing_pipeline as a dependency,
        # jinja_env must be set to sequence_processing_pipeline's root path,
        # rather than the project's root path.
        self.jinja_env = Environment(loader=KISSLoader('templates'))

        self.counts = {}
        self.known_adapters_path = known_adapters_path
        self.bucket_size = bucket_size
        self.length_limit = length_limit

        # NuQCJob() impl uses -c (--cores-per-task) switch instead of
        # -n (--tasks-per-node). --cores-per-task requests the number of cpus
        # per process. This is to support multithreaded jobs that require more
        # than one cpu per task. All cores will be allocated on a single node.
        #
        # This is different than using -n + -N (number of nodes to request)
        # because it's allowable to request more cores than are available on
        # one node using this pair of switches (N nodes * n tasks per node).
        self.cores_per_task = cores_per_task

        self.temp_dir = join(self.output_path, 'tmp')
        makedirs(self.temp_dir, exist_ok=True)

        self.batch_prefix = f"hds-{self.qiita_job_id}"
        self.minimum_bytes = 3100
        self.fastq_regex = re.compile(r'^(.*)_S\d{1,4}_L\d{3}_R\d_\d{3}'
                                      r'\.fastq\.gz$')
        self.interleave_fastq_regex = re.compile(r'^(.*)_S\d{1,4}_L\d{3}_R\d'
                                                 r'_\d{3}\.interleave\.fastq'
                                                 r'\.gz$')
        self.html_regex = re.compile(r'^(.*)_S\d{1,4}_L\d{3}_R\d_\d{3}\.html$')
        self.json_regex = re.compile(r'^(.*)_S\d{1,4}_L\d{3}_R\d_\d{3}\.json$')

        self._validate_project_data()

    def _validate_project_data(self):
        # Validate project settings in [Bioinformatics]
        for project in self.project_data:
            if project['ForwardAdapter'] == 'NA':
                project['ForwardAdapter'] = None

            if project['ReverseAdapter'] == 'NA':
                project['ReverseAdapter'] = None

            if project['ForwardAdapter'] is None:
                if project['ReverseAdapter'] is not None:
                    raise ValueError(("ForwardAdapter is declared but not "
                                      "ReverseAdapter."))

            if project['ReverseAdapter'] is None:
                if project['ForwardAdapter'] is not None:
                    raise ValueError(("ReverseAdapter is declared but not "
                                      "ForwardAdapter."))

            if not isinstance(project['HumanFiltering'], bool):
                raise ValueError("needs_adapter_trimming must be boolean.")

    def _filter_empty_fastq_files(self, filtered_directory,
                                  empty_files_directory,
                                  minimum_bytes):
        '''
        Filters out and moves fastq files that are below threshold.
        :param filtered_directory:
        :param empty_files_directory:
        :param minimum_bytes:
        :return:
        '''
        empty_list = []

        files = glob(join(filtered_directory, f'*.{self.suffix}'))

        for r1, r2 in iter_paired_files(files):
            full_path = join(filtered_directory, r1)
            full_path_reverse = join(filtered_directory, r2)
            if stat(full_path).st_size <= minimum_bytes or stat(
                    full_path_reverse).st_size <= minimum_bytes:
                logging.debug(f'moving {full_path} and {full_path_reverse}'
                              f' to empty list.')
                empty_list.append(full_path)
                empty_list.append(full_path_reverse)

        if empty_list:
            logging.debug(f'making directory {empty_files_directory}')
            makedirs(empty_files_directory, exist_ok=True)

        for item in empty_list:
            logging.debug(f'moving {item}')
            move(item, empty_files_directory)

    def _move_helper(self, completed_files, regex, samples_in_project, dst):
        files_to_move = []
        for fp in completed_files:
            file_name = basename(fp)
            substr = regex.search(file_name)
            if substr is None:
                raise ValueError(f"{file_name} does not follow naming "
                                 "pattern.")
            else:
                # check if found substring is a member of this
                # project. Note sample-name != sample-id
                if substr[1] in samples_in_project:
                    if fp.endswith('.fastq.gz'):
                        # legacy QC'ed files were always denoted with
                        # 'trimmed' to distinguish them from raw files.
                        renamed_fp = fp.replace('.fastq.gz',
                                                '.trimmed.fastq.gz')
                        rename(fp, renamed_fp)
                        # move file into destination w/new filename
                        files_to_move.append(renamed_fp)
                    else:
                        # move file into destination folder w/no namechange.
                        files_to_move.append(fp)

        for fp in files_to_move:
            move(fp, dst)

    def _move_trimmed_files(self, project_name, output_path):
        '''
        Given output_path, move all fastqs to a new subdir named project_name.
        :param project_name: The name of the new folder to be created.
        :param output_path: The path to scan for fastq files.
        :return: None
        '''

        if exists(output_path):
            pattern = f"{output_path}/*.fastq.gz"

            # this directory shouldn't already exist.
            makedirs(join(output_path, project_name), exist_ok=False)

            sample_ids = [x[0] for x in self.sample_ids
                          if x[1] == project_name]

            for trimmed_file in list(glob(pattern)):
                file_name = split(trimmed_file)[1]
                substr = self.interleave_fastq_regex.search(file_name)
                if substr is not None:
                    # only move the sample_ids in this project.
                    if substr[1] in sample_ids:
                        move(trimmed_file, join(output_path, project_name))
        else:
            raise ValueError(f"'{output_path}' does not exist")

    def run(self, callback=None):
        # now a single job-script will be created to process all projects at
        # the same time, and intelligently handle adapter-trimming as needed
        # as well as human-filtering.

        batch_location = join(self.temp_dir, self.batch_prefix)

        if self.force_job_fail:
            batch_count = 0
            max_size = 0
        else:
            batch_count, max_size = split_similar_size_bins(self.root_dir,
                                                            self.bucket_size,
                                                            batch_location)

        job_script_path = self._generate_job_script(max_size)

        self.counts[self.batch_prefix] = batch_count

        export_params = [f"PREFIX={batch_location}",
                         f"OUTPUT={self.output_path}",
                         f"TMPDIR={self.temp_dir}"]

        job_params = ['-J', self.batch_prefix, f'--array 1-{batch_count}',
                      '--export', ','.join(export_params)]

        # job_script_path formerly known as:
        #  process.multiprep.pangenome.adapter-filter.pe.sbatch

        try:
            job_info = self.submit_job(job_script_path,
                                       job_parameters=' '.join(job_params),
                                       exec_from=self.log_path,
                                       callback=callback)
        except JobFailedError as e:
            # When a job has failed, parse the logs generated by this specific
            # job to return a more descriptive message to the user.
            info = self.parse_logs()
            # prepend just the message component of the Error.
            info.insert(0, str(e))
            raise JobFailedError('\n'.join(info))

        job_id = job_info['job_id']

        self.mark_job_completed()

        logging.debug(f'NuQCJob {job_id} completed')

        for project in self.project_data:
            project_name = project['Sample_Project']
            needs_human_filtering = project['HumanFiltering']
            source_dir = join(self.output_path, project_name)
            pattern = f"{source_dir}/*.fastq.gz"
            completed_files = list(glob(pattern))

            # if the 'only-adapter-filtered' directory exists, move the files
            # into a unique location so that files from multiple projects
            # don't overwrite each other.
            trimmed_only_path = join(self.output_path,
                                     'only-adapter-filtered')

            if exists(trimmed_only_path):
                self._move_trimmed_files(project_name, trimmed_only_path)

            if needs_human_filtering is True:
                filtered_directory = join(source_dir, 'filtered_sequences')
            else:
                filtered_directory = join(source_dir, 'trimmed_sequences')

            # create the properly named directory to move files to in
            # in order to preserve legacy behavior.
            makedirs(filtered_directory, exist_ok=True)

            # get the list of sample-names in this project.
            samples_in_project = [x[0] for x in self.sample_ids
                                  if x[1] == project_name]

            # Tissue_1_Mag_Hom_DNASe_RIBO_S16_L001_R2_001.fastq.gz
            # Nislux_SLC_Trizol_DNASe_S7_L001_R2_001.fastq.gz
            self._move_helper(completed_files,
                              self.fastq_regex,
                              samples_in_project,
                              filtered_directory)

            # once fastq.gz files have been moved into the right project,
            # we now need to consider the html and json fastp_reports
            # files.
            old_html_path = join(self.output_path, 'fastp_reports_dir', 'html')
            old_json_path = join(self.output_path, 'fastp_reports_dir', 'json')

            new_html_path = join(source_dir, 'fastp_reports_dir', 'html')
            new_json_path = join(source_dir, 'fastp_reports_dir', 'json')

            makedirs(new_html_path, exist_ok=True)
            makedirs(new_json_path, exist_ok=True)

            # move all html files underneath the subdirectory for this project.
            pattern = f"{old_html_path}/*.html"
            completed_htmls = list(glob(pattern))
            self._move_helper(completed_htmls,
                              # Tissue_1_Super_Trizol_S19_L001_R1_001.html
                              self.html_regex,
                              samples_in_project,
                              new_html_path)

            # move all json files underneath the subdirectory for this project.
            pattern = f"{old_json_path}/*.json"
            completed_jsons = list(glob(pattern))
            self._move_helper(completed_jsons,
                              # Tissue_1_Super_Trizol_S19_L001_R1_001.json
                              self.json_regex,
                              samples_in_project,
                              new_json_path)

            # now that files are separated by project as per legacy
            # operation, continue normal processing.
            empty_files_directory = join(source_dir, 'zero_files')
            self._filter_empty_fastq_files(filtered_directory,
                                           empty_files_directory,
                                           self.minimum_bytes)

        self.mark_post_processing_completed()

    def _confirm_job_completed(self):
        # since NuQCJob processes across all projects in a run, there isn't
        # a need to iterate by project_name and job_id.
        pattern = f"{self.output_path}/hds-{self.qiita_job_id}.*.completed"
        completed_files = list(glob(pattern))
        if completed_files:
            return True

        return False

    def _process_sample_sheet(self):
        sheet = load_sample_sheet(self.sample_sheet_path)

        if not sheet.validate_and_scrub_sample_sheet():
            s = "Sample sheet %s is not valid." % self.sample_sheet_path
            raise PipelineError(s)

        header = sheet.Header
        chemistry = header['chemistry']

        if header['Assay'] not in Pipeline.assay_types:
            s = "Assay value '%s' is not recognized." % header['Assay']
            raise PipelineError(s)

        sample_ids = []
        for sample in sheet.samples:
            sample_ids.append((sample['Sample_ID'], sample['Sample_Project']))

        bioinformatics = sheet.Bioinformatics

        # reorganize the data into a list of dictionaries, one for each row.
        # the ordering of the rows will be preserved in the order of the list.
        lst = bioinformatics.to_dict('records')

        # human-filtering jobs are scoped by project. Each job requires
        # particular knowledge of the project.
        return {'chemistry': chemistry,
                'projects': lst,
                'sample_ids': sample_ids}

    def _generate_mmi_filter_cmds(self, working_dir):
        initial_input = join(working_dir, "seqs.interleaved.fastq")
        final_output = join(working_dir, "seqs.interleaved.filter_"
                                         "alignment.fastq")

        cmds = []

        tmp_file1 = join(working_dir, "foo")
        tmp_file2 = join(working_dir, "bar")

        cores_to_allocate = int(self.cores_per_task / 2)

        # the default setting.
        tags = ""
        t_switch = ""

        if self.additional_fastq_tags:
            # add tags for known metadata types that fastq files may have
            # been annotated with. Samtools will safely ignore tags that
            # are not present.
            # NB: This doesn't appear to be true, actually. if there is
            # a metadata element but it does not begin with 'BX', supplying
            # '-T BX' will cause an error writing output to disk.
            tags = " -T %s" % ','.join(self.additional_fastq_tags)
            t_switch = " -y"

        for count, mmi_db_path in enumerate(self.mmi_file_paths):
            if count == 0:
                # prime initial state with unfiltered file and create first of
                # two tmp files.
                input = initial_input
                output = tmp_file1
            elif count % 2 == 0:
                # alternate between two tmp file names so that the input and
                # output are never the same file.
                input = tmp_file2
                output = tmp_file1
            else:
                input = tmp_file1
                output = tmp_file2

            cmds.append(f"minimap2 -2 -ax sr{t_switch} -t {cores_to_allocate} "
                        f"{mmi_db_path} {input} -a | samtools fastq -@ "
                        f"{cores_to_allocate} -f 12 -F 256{tags} > "
                        f"{output}")

        # rename the latest tmp file to the final output filename.
        cmds.append(f"mv {output} {final_output}")

        # simple cleanup erases tmp files if they exist.
        cmds.append(f"[ -e {tmp_file1} ] && rm {tmp_file1}")
        cmds.append(f"[ -e {tmp_file2} ] && rm {tmp_file2}")

        return "\n".join(cmds)

    def _generate_job_script(self, max_bucket_size):
        # bypass generating job script for a force-fail job, since it is
        # not needed.
        if self.force_job_fail:
            return None

        job_script_path = join(self.output_path, 'process_all_fastq_files.sh')
        template = self.jinja_env.get_template("nuqc_job.sh")

        job_name = f'{self.qiita_job_id}_{self.job_name}'

        html_path = join(self.output_path, 'fastp_reports_dir', 'html')
        json_path = join(self.output_path, 'fastp_reports_dir', 'json')

        # get location of python executable in this environment.
        # demux script should be present in the same location.
        demux_path = join(dirname(executable), 'demux')

        if not exists(demux_path):
            raise ValueError(f"{demux_path} does not exist.")

        # get this file location and add splitter as it should live there
        splitter_binary = join(
            dirname(abspath(__file__)), 'scripts', 'splitter')
        if not exists(splitter_binary):
            raise ValueError(f'{splitter_binary} does not exist.')

        # this method relies on an environment variable defined in nu_qc.sh
        # used to define where unfiltered fastq files are and where temp
        # files can be created. (${jobd})
        mmi_filter_cmds = self._generate_mmi_filter_cmds("${jobd}")

        with open(job_script_path, mode="w", encoding="utf-8") as f:
            # the job resources should come from a configuration file

            # generate a string of linux system modules to load before
            # processing begins.
            mtl = ' '.join(self.modules_to_load)

            f.write(template.render(job_name=job_name,
                                    queue_name=self.queue_name,
                                    # should be 4 * 24 * 60 = 4 days
                                    wall_time_limit=self.wall_time_limit,
                                    mem_in_gb=self.jmem,
                                    # Note NuQCJob now maps node_count to
                                    # SLURM -N parameter to act like other
                                    # Job classes.
                                    # self.node_count should be 1
                                    node_count=self.node_count,
                                    # cores-per-task (-c) should be 4
                                    cores_per_task=self.cores_per_task,
                                    knwn_adpt_path=self.known_adapters_path,
                                    output_path=self.output_path,
                                    html_path=html_path,
                                    json_path=json_path,
                                    demux_path=demux_path,
                                    temp_dir=self.temp_dir,
                                    splitter_binary=splitter_binary,
                                    modules_to_load=mtl,
                                    length_limit=self.length_limit,
                                    gres_value=self.gres_value,
                                    movi_path=self.movi_path,
                                    mmi_filter_cmds=mmi_filter_cmds,
                                    pmls_path=self.pmls_path))

        return job_script_path

    def parse_logs(self):
        log_path = join(self.output_path, 'logs')
        files = sorted(glob(join(log_path, '*')))
        msgs = []

        # assume that the only possible files in logs directory are '.out'
        # files and zero, one, or many 'seqs.movi.n.txt.gz' files.
        # the latter may be present because the last step of a successful
        # job is to rename and move this file into its final location while
        # the logs directory is the default 'working' directory for this job
        # as this ensures slurm.err and slurm.out files will always be in
        # a known location.

        # for now, construct lists of both of these types of files.
        output_logs = [x for x in files if x.endswith('.out')]

        for some_file in output_logs:
            with open(some_file, 'r') as f:
                msgs += [line.strip() for line in f.readlines()
                         if 'error:' in line.lower()]

        return msgs
