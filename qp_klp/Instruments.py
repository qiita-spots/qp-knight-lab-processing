from sequence_processing_pipeline.ConvertJob import ConvertJob
from sequence_processing_pipeline.TellReadJob import TellReadJob
from sequence_processing_pipeline.TRNormCountsJob import TRNormCountsJob
from sequence_processing_pipeline.TRIntegrateJob import TRIntegrateJob
from qp_klp.FailedSamplesRecord import FailedSamplesRecord


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
    instrument_type = INSTRUMENT_NAME_NONE


class Illumina(Instrument, FailedSamplesRecord):
    instrument_type = INSTRUMENT_NAME_ILLUMINA

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

        job.run(callback=self.status_update_callback)

        # audit the results to determine which samples failed to convert
        # properly. Append these to the failed-samples report and also
        # return the list directly to the caller.
        failed_samples = job.audit(self.pipeline.get_sample_ids())
        self.fsr_write(failed_samples, job.__class__.__name__)
        return failed_samples


class TellSeq(Instrument):
    instrument_type = INSTRUMENT_NAME_TELLSEQ

    def convert_raw_to_fastq(self):
        config = self.pipeline.get_software_configuration('tell-read')

        tr_job = TellReadJob(self.pipeline.run_dir,
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

        tr_job.run(callback=self.status_update_callback)

        # TODO: determine these appropriately
        max_array_length = "foo"
        label = "bar"

        if self.iseq_run:
            # NB: the original master script performed this job prior to
            # integration and after the main job completed, but only
            # for iSeq jobs. This may be useful for additional
            # situations as well.
            nc_job = TRNormCountsJob(self.pipeline.run_dir,
                                     self.pipeline.output_path,
                                     self.pipeline.input_file_path,
                                     config['queue'],
                                     config['nodes'],
                                     config['wallclock_time_in_minutes'],
                                     config['per_process_memory_limit'],
                                     config['modules_to_load'],
                                     self.master_qiita_job_id,
                                     max_array_length,
                                     config['indicies_script_path'],
                                     label)

            nc_job.run(callback=self.status_update_callback)

        # after the primary job and the optional counts job is completed,
        # the job to integrate results and add metadata to the fastq files
        # is performed.
        i_job = TRIntegrateJob(self.pipeline.run_dir,
                               self.pipeline.output_path,
                               self.pipeline.input_file_path,
                               config['queue'],
                               config['nodes'],
                               config['wallclock_time_in_minutes'],
                               config['per_process_memory_limit'],
                               config['modules_to_load'],
                               self.master_qiita_job_id,
                               max_array_length,
                               config['indicies_script_path'],
                               label)

        i_job.run(callback=self.status_update_callback)

        # TODO: after i_job is completed, there are two optional jobs that
        # can be performed in parallel using the new functionality in Job()
        # class.

        # TODO: take post-processing code from Version 1 impl and tack on
        # the cleanup line from the original cleanup script in order to
        # organize the fastq files into the form downstream operations
        # expect in SPP.

        # NB: Note that unlike convert_raw_to_fastq() in other Workflows(),
        # this implementation uses 3-5 separate Job()s in order to have more
        # flexibility in implementing workflows. The integrate job is
        # performed here because we identified that we want the integrated
        # results. Thhe two optional jobs could instead be called from
        # a separate method in this Instrument class if they are valuable
        # but not part of generating raw (that is unfiltered) fastq files.

        # Potentially needed for tell-seq jobs but perhaps not.
        # return job.audit(self.pipeline.get_sample_ids())
