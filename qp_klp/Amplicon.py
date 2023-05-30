from os import listdir, makedirs
from os.path import join, isfile, basename
from shutil import copyfile
from qp_klp.Step import Step


class Amplicon(Step):
    def __init__(self, pipeline, master_qiita_job_id,
                 status_update_callback=None):
        super().__init__(pipeline,
                         master_qiita_job_id,
                         status_update_callback)

        if pipeline.pipeline_type != Step.AMPLICON_TYPE:
            raise ValueError("Cannot create an Amplicon run using a "
                             f"{pipeline.pipeline_type}-configured Pipeline.")

    def convert_bcl_to_fastq(self):
        # The 'bcl2fastq' key is a convention hard-coded into mg-scripts and
        # qp-klp projects. By design and history, amplicon jobs use EMP primers
        # and are demuxed downstream of qp-klp by Qiita. This necessitates the
        # use of bcl2fastq and a 'dummy' sample-sheet to avoid the
        # demultiplexing that otherwise occurs at this stage. The name and
        # path of the executable, the resource requirements to instantiate a
        # SLURM job with, etc. are stored in configuration['bcl2fastq'].
        config = self.pipeline.configuration['bcl2fastq']
        super()._convert_bcl_to_fastq(config, self.pipeline.sample_sheet)

    def quality_control(self):
        # Quality control for Amplicon runs occurs downstream.
        # Do not perform QC at this time.

        # Simulate QCJob's output directory for use as input into FastQCJob.
        projects = self.pipeline.get_project_info()
        projects = [x['project_name'] for x in projects]

        for project_name in projects:
            # copy the files from ConvertJob output to faked QCJob output
            # folder: $WKDIR/$RUN_ID/QCJob/$PROJ_NAME/amplicon
            output_folder = join(self.pipeline.output_path,
                                 'QCJob',
                                 project_name,
                                 # for legacy purposes, output folders are
                                 # either 'trimmed_sequences', 'amplicon', or
                                 # 'filtered_sequences'. Hence, this folder
                                 # is not defined using AMPLICON_TYPE as that
                                 # value may or may not equal the needed value.
                                 'amplicon')

            makedirs(output_folder)

            raw_fastq_files_path = join(self.pipeline.output_path,
                                        'ConvertJob')

            # get list of all raw output files to be copied.
            job_output = [join(raw_fastq_files_path, x) for x in
                          listdir(raw_fastq_files_path)]
            job_output = [x for x in job_output if isfile(x)]
            job_output = [x for x in job_output if x.endswith('fastq.gz')]
            # Undetermined files are very small and should be filtered from
            # results.
            job_output = [x for x in job_output if
                          not basename(x).startswith('Undetermined')]

            # copy the file
            for fastq_file in job_output:
                new_path = join(output_folder, basename(fastq_file))
                copyfile(fastq_file, new_path)

            # FastQC expects the ConvertJob output to also be organized by
            # project. Since this would entail running the same ConvertJob
            # multiple times on the same input with just a name-change in
            # the dummy sample-sheet, we instead perform ConvertJob once
            # and copy the output from ConvertJob's output folder into
            # ConvertJob's output folder/project1, project2 ... projectN.
            output_folder = join(raw_fastq_files_path, project_name)
            makedirs(output_folder)

            job_output = [join(raw_fastq_files_path, x) for x in
                          listdir(raw_fastq_files_path)]
            job_output = [x for x in job_output if isfile(x) and x.endswith(
                'fastq.gz') and not basename(x).startswith('Undetermined')]

            for raw_fastq_file in job_output:
                new_path = join(output_folder, basename(raw_fastq_file))
                copyfile(raw_fastq_file, new_path)

    def generate_reports(self):
        super()._generate_reports()
        return None  # amplicon doesn't need project names

    def _get_data_type(self, prep_file_path):
        metadata = Step.parse_prep_file(prep_file_path)
        if 'target_gene' in metadata.columns:
            # remove duplicate values, then convert back to list for
            # accession.
            tg = list(set(metadata['sample target_gene']))
            if len(tg) != 1:
                raise ValueError("More than one value for target_gene")

            if tg[0] in Step.AMPLICON_SUB_TYPES:
                return tg[0]

            raise ValueError(f"'{tg[0]}' is not a valid type - valid data-"
                             f"types are {Step.AMPLICON_SUB_TYPES}")
        else:
            raise ValueError("'target_gene' column not present in "
                             "generated prep-files")

    def generate_touched_studies(self, qclient):
        results = {}
        for study_id, pf_paths in self.prep_file_paths.items():
            for pf_path in pf_paths:
                results[pf_path] = self._get_data_type(pf_path)

    def generate_prep_file(self):
        config = self.pipeline.configuration['seqpro']
        seqpro_path = config['seqpro_path'].replace('seqpro', 'seqpro_mf')
        project_names = [x['project_name'] for x in
                         self.pipeline.get_project_info()]

        job = super()._generate_prep_file(config,
                                          self.pipeline.mapping_file_path,
                                          seqpro_path,
                                          project_names)

        self.prep_file_paths = job.prep_file_paths

    def generate_commands(self, qclient):
        super()._generate_commands()
        self.cmds.append(f'cd {self.pipeline.output_path}; '
                         'tar zcvf reports-ConvertJob.tgz ConvertJob/Reports')
        self.write_commands_to_output_path()
