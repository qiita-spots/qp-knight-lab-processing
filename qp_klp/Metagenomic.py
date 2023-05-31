from qp_klp.Step import FailedSamplesRecord, Step


class Metagenomic(Step):
    def __init__(self, pipeline, master_qiita_job_id,
                 status_update_callback=None, lane_number=None):
        super().__init__(pipeline,
                         master_qiita_job_id,
                         status_update_callback,
                         lane_number)

        if pipeline.pipeline_type not in Step.META_TYPES:
            raise ValueError("Cannot instantiate Metagenomic object from "
                             f"pipeline of type '{pipeline.pipeline_type}'")

        # Note: FailedSamplesRecord is not used when processing amplicon as the
        # samples are processed as a single fastq file and hence that info
        # is not available.
        self.fsr = FailedSamplesRecord(self.pipeline.output_path,
                                       pipeline.sample_sheet.samples)

    def convert_bcl_to_fastq(self):
        # The 'bcl-convert' key is a convention hard-coded into mg-scripts and
        # qp-klp projects. Currently meta*omic jobs use bcl-convert for its
        # improved performance over bcl2fastq. The name and path of the
        # executable, the resource requirements to instantiate a SLURM job
        # with, etc. are stored in configuration['bcl-convert''].
        config = self.pipeline.configuration['bcl-convert']
        job = super()._convert_bcl_to_fastq(config,
                                            self.pipeline.sample_sheet.path)
        self.fsr.write(job.audit(self.pipeline.get_sample_ids()), 'ConvertJob')

    def quality_control(self):
        config = self.pipeline.configuration['qc']
        job = super()._quality_control(config, self.pipeline.sample_sheet.path)
        self.fsr.write(job.audit(self.pipeline.get_sample_ids()), 'QCJob')

    def generate_reports(self):
        job = super()._generate_reports()
        self.fsr.write(job.audit(self.pipeline.get_sample_ids()), 'FastQCJob')

        self.project_names = job.project_names

    def generate_prep_file(self):
        config = self.pipeline.configuration['seqpro']

        if self.project_names is None:
            raise ValueError("reports not yet generated")

        job = super()._generate_prep_file(config,
                                          self.pipeline.sample_sheet.path,
                                          config['seqpro_path'],
                                          self.project_names)

        self.prep_file_paths = job.prep_file_paths

    def generate_touched_studies(self, qclient):
        results = {}

        for study_id, pf_paths in self.prep_file_paths.items():
            for pf_path in pf_paths:
                # record the data-type as either metagenomic or
                # metatranscriptomic, according to what's stored in the
                # pipeline.
                results[pf_path] = self.pipeline.pipeline_type

    def generate_commands(self, qclient):
        self.qclient = qclient
        super()._generate_commands()
        out_dir = self.pipeline.output_path
        self.cmds.append(f'cd {self.pipeline.output_path}; '
                         'tar zcvf logs-QCJob.tgz QCJob/logs')

        # copy all tgz files, including sample-files.tgz, to final_results.
        self.cmds.append(f'cd {out_dir}; mv *.tgz final_results')
        self.cmds.append(f'cd {out_dir}; mv FastQCJob/multiqc final_results')

        self.cmds.append(f'cd {self.pipeline.output_path}; '
                         'tar zcvf reports-ConvertJob.tgz ConvertJob/Reports '
                         'ConvertJob/Logs')

        self.write_commands_to_output_path()
