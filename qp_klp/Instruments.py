from sequence_processing_pipeline.ConvertJob import ConvertJob
from sequence_processing_pipeline.TellReadJob import TellReadJob
from sequence_processing_pipeline.TRNormCountsJob import TRNormCountsJob
from sequence_processing_pipeline.TRIntegrateJob import TRIntegrateJob
from sequence_processing_pipeline.PipelineError import PipelineError


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


class Illumina(Instrument):
    instrument_type = INSTRUMENT_NAME_ILLUMINA

    def _get_configuration(self, command):
        # NB: This helper method is to change the behavior of Pipeline.
        # get_software_command. In most cases we know and expect a specific
        # configuration dictionary and it's beneficial to raise an Error if
        # it's not found. In this particular case however, the user expects
        # that one or both of two different config dictionaries can be
        # present.
        try:
            return self.pipeline.get_software_configuration(command)
        except PipelineError as e:
            if str(e) != f"'{command}' is not defined in configuration":
                # re-raise the Error if it's not the one we are expecting and
                # can properly handle here.
                raise PipelineError(e)

    def convert_raw_to_fastq(self):
        # NB: it is preferable to use bcl-convert over bcl2fastq to generate
        # fastq files in most situations. bcl2fastq is used for processing
        # 16S runs, as it handles the need to keep the samples muxed better
        # than bcl-convert.
        #
        # ConvertJob knows how to generate the proper job-script for each
        # utility. It determines the software to use by processing the
        # executable_path for a file_name. Set the software to use by
        # defining either a bcl-convert or bcl2fastq section in the config
        # profile. If both are defined, this class will perfer bcl-convert.
        bcl2fastq_conf = self._get_configuration('bcl2fastq')
        bclcnvrt_conf = self._get_configuration('bcl-convert')

        if bclcnvrt_conf is None and bcl2fastq_conf is None:
            raise PipelineError("bcl-convert and bcl2fastq sections not "
                                "defined in configuation profile.")

        config = bcl2fastq_conf if bclcnvrt_conf is None else bclcnvrt_conf

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
        if hasattr(self, 'fsr'):
            # NB 16S does not require a failed samples report and
            # it is not performed by SPP.
            self.fsr.write(failed_samples, job.__class__.__name__)
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
