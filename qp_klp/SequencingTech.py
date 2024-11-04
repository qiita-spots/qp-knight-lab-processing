from sequence_processing_pipeline.ConvertJob import ConvertJob
from sequence_processing_pipeline.TellReadJob import TellReadJob
from sequence_processing_pipeline.TRNormCountsJob import TRNormCountsJob
from sequence_processing_pipeline.TRIntegrateJob import TRIntegrateJob
from sequence_processing_pipeline.PipelineError import PipelineError
from os.path import join


SEQTECH_NAME_NONE = "None"
SEQTECH_NAME_ILLUMINA = "Illumina"
SEQTECH_NAME_TELLSEQ = "TellSeq"


class SequencingTech():
    """
    SequencingTechs encapsulate Job()s and other functionality that vary on
     the nature of the Instrument used to create the raw data. All Instruments
     are mixins for Workflow() classes and shouldn't define their own
     initialization.
    """
    seqtech_type = SEQTECH_NAME_NONE


class Illumina(SequencingTech):
    seqtech_type = SEQTECH_NAME_ILLUMINA

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


class TellSeq(SequencingTech):
    seqtech_type = SEQTECH_NAME_TELLSEQ

    def convert_raw_to_fastq(self):
        config = self.pipeline.get_software_configuration('tell-read')

        print("RUN DIR: %s" % self.pipeline.run_dir)
        print("OUTPUT PATH: %s" % self.pipeline.output_path)
        print("INPUT FILE PATH: %s" % self.pipeline.input_file_path)

        tr_job = TellReadJob(self.pipeline.run_dir,
                             self.pipeline.output_path,
                             self.pipeline.input_file_path,
                             config['queue'],
                             config['nodes'],
                             config['wallclock_time_in_minutes'],
                             config['per_process_memory_limit'],
                             config['modules_to_load'],
                             self.master_qiita_job_id,
                             config['label'],
                             config['reference_base'],
                             config['reference_map'],
                             config['tmp1_path'],
                             config['sing_script_path'],
                             config['cores_per_task'])

        tr_job.run(callback=self.status_update_callback)

        '''
        when run is run, we're going to create a job script which we can test to see if it checks out.
        then we're going to submit the job, and we need to fake the sbatch command.
        then we need to fake the results directory.
        all we're testing is that the job script is what we expect.

        '''

        """


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

        # NB: after i_job is completed, there are two optional jobs that
        # can be performed in parallel using the new functionality in Job()
        # class. However we are not using the output from this step right now
        # so we will leave it unimplemented temporarily.

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


        # post-process working directory to make it appear like results
        # generated by ConvertJob

        integrated_files_path = join(self.pipeline.output_path, 'output', "integrated")

        if not exists(integrated_files_path):
            raise ValueError(f"{integrated_files_path} does not exist")

        # move integrated directory to TRConvertJob directory, co-level with
        # output directory. This makes it easier to delete the rest of the
        # output that we don't need.

        # move err and out logs into logs subdirectory.
        for root, dirs, files in walk(self.output_path):
            for _file in files:
                _path = join(root, _file)
                if _path.endswith('.err'):
                    move(_path, join(self.output_path, 'logs'))
                elif _path.endswith('.out'):
                    move(_path, join(self.output_path, 'logs'))
            # don't go below one level.
            break

        # save two logs and move them into standard Job logs directory.
        move(join(self.output_path, 'output', 'log'),
             join(self.output_path, 'logs'))
        move(join(self.output_path, 'output', 'output.log'),
             join(self.output_path, 'logs'))

        # rename the files and move them into project directories.
        for root, dirs, files in walk(integrated_files_path):
            for _file in files:
                fastq_file = join(root, _file)
                self._post_process_file(fastq_file, self.mapping)

        # move project folders from integrated directory to working_dir.
        contents = listdir(integrated_files_path)
        for name in contents:
            move(join(integrated_files_path, name),
                 self.output_path)

        # delete the original output directory.
        rmtree(join(self.output_path, 'output'))




\

        def _post_process_file(self, fastq_file, mapping):
        # generate names of the form generated by bcl-convert/bcl2fastq:
        # <Sample_ID>_S#_L00#_<R# or I#>_001.fastq.gz
        # see:
        # https://help.basespace.illumina.com/files-used-by-basespace/
        # fastq-files
        _dir, _file = split(fastq_file)

        # ex: integrated/C544.R2.fastq.gz
        m = match(r"(C5\d\d)\.([R,I]\d)\.fastq.gz", _file)

        if m is None:
            raise ValueError(f"The filename '{_file}' is not of a "
                             "recognizable form")

        adapter_id = m[1]
        read_type = m[2]

        if adapter_id not in mapping:
            raise ValueError(f"{adapter_id} is not present in mapping")

        sample_name, sample_index, project_name = mapping[adapter_id]

        # generate the new filename for the fastq file, and reorganize the
        # files by project.
        new_name = "%s_S%d_%s_%s_001.fastq.gz" % (sample_name,
                                                  sample_index,
                                                  self.lane,
                                                  read_type)

        # ensure that the project directory exists before we rename and move
        # the file to that location.
        makedirs(join(_dir, project_name), exist_ok=True)

        # if there's an error renaming and moving the file, let it pass up to
        # the user.
        final_path = join(_dir, project_name, new_name)
        rename(fastq_file, final_path)
        return final_path

    def _generate_sample_mapping(self):
        # this generates a sample mapping for the C501-C596 adapters used by
        # the vendor to a sample-name and project. In production use this
        # mapping would need to be created from the future sample-sheet.
        project_names = ['Project1', 'Project2', 'Project3']
        sample_mapping = {}

        for sample_index in range(1, 97):
            adapter_id = "C%s" % str(sample_index + 500)
            sample_name = "MySample%d" % sample_index
            project_name = project_names[sample_index % 3]
            sample_mapping[adapter_id] = (sample_name, sample_index,
                                          project_name)

        return sample_mapping
        """


