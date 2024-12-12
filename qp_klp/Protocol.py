from sequence_processing_pipeline.ConvertJob import ConvertJob
from sequence_processing_pipeline.TellReadJob import TellReadJob
from sequence_processing_pipeline.SeqCountsJob import SeqCountsJob
from sequence_processing_pipeline.TRIntegrateJob import TRIntegrateJob
from sequence_processing_pipeline.PipelineError import PipelineError
from os.path import join, split
from re import match
from os import makedirs, rename, walk
from metapool import load_sample_sheet


PROTOCOL_NAME_NONE = "None"
PROTOCOL_NAME_ILLUMINA = "Illumina"
PROTOCOL_NAME_TELLSEQ = "TellSeq"


class Protocol():
    """
    Protocols encapsulate Job()s and other functionality that vary on
     the nature of the Instrument used to create the raw data. All Instruments
     are mixins for Workflow() classes and shouldn't define their own
     initialization.
    """
    protocol_type = PROTOCOL_NAME_NONE


class Illumina(Protocol):
    protocol_type = PROTOCOL_NAME_ILLUMINA

    def convert_raw_to_fastq(self):
        def get_config(command):
            try:
                return self.pipeline.get_software_configuration(command)
            except PipelineError as e:
                if str(e) != f"'{command}' is not defined in configuration":
                    # re-raise the Error if it's not the one we are expecting
                    # and can properly handle here.
                    raise PipelineError(e)

        # The choice to use either bcl2fastq or bcl-convert is set in the
        # configuration file that matches the parameters of the sequencing
        # run found in the user input file. The Pipeline object is smart
        # enough to select the right configuration file based on the
        # sample sheet passed to it.

        # However in this mixin, we don't know which configuration file we
        # will get and thus we must check for the presence or absence of
        # both dictionaries in the configuration file at run-time and make
        # a determination based on that.
        bcl2fastq_conf = get_config('bcl2fastq')
        bclcnvrt_conf = get_config('bcl-convert')

        if bclcnvrt_conf is None and bcl2fastq_conf is None:
            raise PipelineError("bcl-convert and bcl2fastq sections not "
                                "defined in configuation profile.")

        # if both are defined, we will use bcl-convert by default.
        config = bcl2fastq_conf if bclcnvrt_conf is None else bclcnvrt_conf

        job = ConvertJob(self.pipeline.run_dir,
                         self.pipeline.output_path,
                         # for amplicon runs, get_sample_sheet_path() will
                         # return the path for a generated dummy sheet.
                         self.pipeline.get_sample_sheet_path(),
                         config['queue'],
                         config['nodes'],
                         config['nprocs'],
                         config['wallclock_time_in_minutes'],
                         config['per_process_memory_limit'],
                         config['executable_path'],
                         config['modules_to_load'],
                         self.master_qiita_job_id)

        self.raw_fastq_files_path = join(self.pipeline.output_path,
                                         'ConvertJob')

        if 'ConvertJob' not in self.skip_steps:
            job.run(callback=self.job_callback)

        # audit the results to determine which samples failed to convert
        # properly. Append these to the failed-samples report and also
        # return the list directly to the caller.
        failed_samples = job.audit(self.pipeline.get_sample_ids())
        if hasattr(self, 'fsr'):
            self.fsr.write(failed_samples, job.__class__.__name__)

        return failed_samples


class TellSeq(Protocol):
    protocol_type = PROTOCOL_NAME_TELLSEQ

    def convert_raw_to_fastq(self):
        config = self.pipeline.get_software_configuration('tell-seq')

        job = TellReadJob(self.pipeline.run_dir,
                          self.pipeline.output_path,
                          self.pipeline.input_file_path,
                          config['queue'],
                          config['nodes'],
                          config['wallclock_time_in_minutes'],
                          config['tellread_mem_limit'],
                          config['modules_to_load'],
                          self.master_qiita_job_id,
                          config['reference_base'],
                          config['reference_map'],
                          config['sing_script_path'],
                          config['tellread_cores'])

        self.raw_fastq_files_path = join(self.pipeline.output_path,
                                         'TellReadJob', 'Full')
        self.my_sil_path = join(self.pipeline.output_path,
                                'TellReadJob',
                                'sample_index_list_TellReadJob.txt')

        # if TellReadJob already completed, then skip the over the time-
        # consuming portion but populate the needed member variables.
        if 'TellReadJob' not in self.skip_steps:
            job.run(callback=self.job_callback)

        self.pipeline.get_sample_ids()
        failed_samples = []
        if hasattr(self, 'fsr'):
            # NB 16S does not require a failed samples report and
            # it is not performed by SPP.
            self.fsr.write(failed_samples, job.__class__.__name__)

        return failed_samples

    def generate_sequence_counts(self):
        config = self.pipeline.get_software_configuration('tell-seq')

        job = SeqCountsJob(self.pipeline.run_dir,
                           self.pipeline.output_path,
                           self.pipeline.input_file_path,
                           config['queue'],
                           config['nodes'],
                           config['wallclock_time_in_minutes'],
                           config['normcount_mem_limit'],
                           config['modules_to_load'],
                           self.master_qiita_job_id,
                           '',
                           config['integrate_script_path'],
                           self.pipeline.qiita_job_id)

        if 'SeqCountsJob' not in self.skip_steps:
            job.run(callback=self.job_callback)

        # audit the results to determine which samples failed to convert
        # properly. Append these to the failed-samples report and also
        # return the list directly to the caller.
        failed_samples = job.audit_me(self.pipeline.get_sample_ids())
        if hasattr(self, 'fsr'):
            # NB 16S does not require a failed samples report and
            # it is not performed by SPP.
            self.fsr.write(failed_samples, job.__class__.__name__)

        return failed_samples

    def integrate_results(self):
        config = self.pipeline.get_software_configuration('tell-seq')

        # after the primary job and the optional counts job is completed,
        # the job to integrate results and add metadata to the fastq files
        # is performed.

        job = TRIntegrateJob(self.pipeline.run_dir,
                             self.pipeline.output_path,
                             self.pipeline.input_file_path,
                             config['queue'],
                             config['nodes'],
                             config['wallclock_time_in_minutes'],
                             config['integrate_mem_limit'],
                             config['modules_to_load'],
                             self.master_qiita_job_id,
                             "foo",
                             config['integrate_script_path'],
                             # NB: sample_index_list used may vary
                             # from project to project in the future.
                             # If so replace config entry with a user
                             # supplied entry or an entry in the sample
                             # sheet.

                             # NB: the version of the sil needed is the one
                             # generated by TellReadJob(), not the master
                             # sil. The former will account for a set of
                             # barcode_ids that don't begin at C501 and/or
                             # skip over some values like C509.
                             self.my_sil_path,
                             self.raw_fastq_files_path,
                             "",
                             "",
                             config['integrate_cores'])

        if 'TRIntegrateJob' not in self.skip_steps:
            job.run(callback=self.job_callback)

        # raw_fastq_files_path is used by downstream processes to know
        # where to locate 'raw' fastq files. Before we could assume that
        # it would always be in ConvertJob's working directory but now
        # this is no longer the case. Currently used by NuQCJob.
        self.raw_fastq_files_path = join(self.pipeline.output_path,
                                         'TRIntegrateJob',
                                         'integrated')

        if 'TRIJ_Post_Processing' in self.skip_steps:
            # if 'post_processing_completed' is found in the same location,
            # this means the steps below were already completed.
            return

        # post-processing the results is a relatively trivial task in terms
        # of time. It also relies on metadata stored in the job object.
        # Hence, it is performed here.
        mapping = self._generate_mapping()

        # rename the files and move them into project directories.
        for root, dirs, files in walk(self.raw_fastq_files_path):
            for _file in files:
                fastq_file = join(root, _file)
                self._post_process_file(fastq_file,
                                        mapping,
                                        self.lane_number)

        # audit the results to determine which samples failed to convert
        # properly. Append these to the failed-samples report and also
        # return the list directly to the caller.
        failed_samples = job.audit_me(self.pipeline.get_sample_ids())

        if hasattr(self, 'fsr'):
            # NB 16S does not require a failed samples report and
            # it is not performed by SPP.
            self.fsr.write(failed_samples, job.__class__.__name__)

        return failed_samples

    def _generate_mapping(self):
        # NB: This method should only be implemented for tellseq-related
        # SequencingTech() objects where a barcode_id column can be expected
        # to exist and a mapping is needed.
        sheet = load_sample_sheet(self.pipeline.get_sample_sheet_path())

        results = {}

        count = 1
        for sample in sheet.samples:
            barcode_id = sample['barcode_id']
            results[barcode_id] = {'sample_name': sample['Sample_Name'],
                                   'sample_id': sample['Sample_ID'],
                                   'sample_index': count,
                                   'project_name': sample['Sample_Project']}
            count += 1
        return results

    def _post_process_file(self, fastq_file, mapping, lane):
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

        barcode_id = m[1]
        read_type = m[2]

        if barcode_id not in mapping:
            raise ValueError(f"{barcode_id} is not present in sample-sheet")

        sample_id = mapping[barcode_id]['sample_id']
        project_name = mapping[barcode_id]['project_name']
        sample_index = mapping[barcode_id]['sample_index']

        # generate the new filename for the fastq file, and reorganize the
        # files by project.
        new_name = "%s_S%d_%s_%s_001.fastq.gz" % (sample_id,
                                                  sample_index,
                                                  "L%s" % str(lane).zfill(3),
                                                  read_type)

        # ensure that the project directory exists before we rename and move
        # the file to that location.
        makedirs(join(_dir, project_name), exist_ok=True)

        # if there's an error renaming and moving the file, let it pass up to
        # the user.
        final_path = join(_dir, project_name, new_name)
        rename(fastq_file, final_path)

        return final_path
