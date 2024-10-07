from sequence_processing_pipeline.ConvertJob import ConvertJob
from sequence_processing_pipeline.TellReadJob import TellReadJob


INSTRUMENT_NAME_NONE = "Instrument"
INSTRUMENT_NAME_ILLUMINA = "Illumina"
INSTRUMENT_NAME_TELLSEQ = "TellSeq"


class Instrument():
    """
    Instruments encapsulate Job()s and other functionality that vary on the
     nature of the Instrument used to create the raw data. All Instruments are
     mixins for Workflow() classes and shouldn't define their own
     initialization.
    """
    @property
    def instrument_name(self):
        return INSTRUMENT_NAME_NONE


class Illumina(Instrument):
    @property
    def instrument_name(self):
        return INSTRUMENT_NAME_ILLUMINA

    def convert_raw_to_fastq(self):
        config = self.pipeline.get_software_configuration('bcl-convert')

        job = ConvertJob(self.pipeline.run_dir,
                         self.pipeline.output_path,
                         self.pipeline.input_file_path,
                         config['queue'],
                         config['nodes'],
                         config['nprocs'],
                         config['wallclock_time_in_minutes'],
                         config['per_process_memory_limit'],
                         config['executable_path'],
                         config['modules_to_load'],
                         self.master_qiita_job_id)

        job.run(callback=self.update_callback)

        # NB: This isn't currently needed for Amplicon runs.
        return job.audit(self.pipeline.get_sample_ids())


class TellSeq(Instrument):
    @property
    def instrument_name(self):
        # TODO: Not certain TellSeq is the proper name of the Instrument but
        #  for now let's run with it.
        return INSTRUMENT_NAME_TELLSEQ

    def convert_raw_to_fastq(self):
        config = self.pipeline.get_software_configuration('tell-read')

        job = TellReadJob(self.pipeline.run_dir,
                          self.pipeline.output_path,
                          self.pipeline.input_file_path,
                          config['queue'],
                          config['nodes'],
                          config['nprocs'],
                          config['wallclock_time_in_minutes'],
                          config['per_process_memory_limit'],
                          config['executable_path'],
                          config['modules_to_load'],
                          self.master_qiita_job_id)

        job.run(callback=self.update_callback)

        # NB: Potentially needed for tell-read.
        return job.audit(self.pipeline.get_sample_ids())
