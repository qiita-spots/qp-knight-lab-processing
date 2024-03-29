from qp_klp.Step import FailedSamplesRecord, Step


class Metagenomic(Step):
    def __init__(self, pipeline, master_qiita_job_id,
                 status_update_callback=None, lane_number=None,
                 is_restart=False):
        super().__init__(pipeline,
                         master_qiita_job_id,
                         status_update_callback,
                         lane_number,
                         is_restart)

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
        config = self.pipeline.config_profile['profile']['configuration']
        job = super()._convert_bcl_to_fastq(config['bcl-convert'],
                                            self.pipeline.sample_sheet.path)
        self.fsr.write(job.audit(self.pipeline.get_sample_ids()), 'ConvertJob')

    def quality_control(self):
        config = self.pipeline.config_profile['profile']['configuration']
        job = super()._quality_control(config['nu-qc'],
                                       self.pipeline.sample_sheet.path)
        self.fsr.write(job.audit(self.pipeline.get_sample_ids()), 'NuQCJob')

    def generate_reports(self):
        job = super()._generate_reports()
        self.fsr.write(job.audit(self.pipeline.get_sample_ids()), 'FastQCJob')

        # as generate_reports is the final step that updates self.fsr,
        # we can generate the final human-readable report following this step.
        self.fsr.generate_report()

    def get_data_type(self, prep_file_path):
        # prep_file_path is unused. It's kept for compatability with Amplicon
        # and Step.
        return self.pipeline.pipeline_type

    def generate_prep_file(self):
        config = self.pipeline.config_profile['profile']['configuration']

        job = super()._generate_prep_file(config['seqpro'],
                                          self.pipeline.sample_sheet.path,
                                          config['seqpro']['seqpro_path'])

        self.prep_file_paths = job.prep_file_paths
        self.has_replicates = job.has_replicates
