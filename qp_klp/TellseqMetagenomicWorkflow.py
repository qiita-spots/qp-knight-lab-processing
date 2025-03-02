from .Protocol import TellSeq
from os.path import join, abspath, exists
from os import walk
from sequence_processing_pipeline.Pipeline import Pipeline, InstrumentUtils
from .Assays import Metagenomic
from .Assays import ASSAY_NAME_METAGENOMIC
from .Workflows import Workflow
from .FailedSamplesRecord import FailedSamplesRecord
from collections import defaultdict


class TellSeqMetagenomicWorkflow(Workflow, Metagenomic, TellSeq):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.mandatory_attributes = ['qclient', 'uif_path', 'config_fp',
                                     'run_identifier', 'output_dir', 'job_id',
                                     'is_restart']

        self.confirm_mandatory_attributes()

        # second stage initializer that could conceivably be pushed down into
        # specific children requiring specific parameters.
        self.qclient = self.kwargs['qclient']

        run_id = self.kwargs['run_identifier']

        self.pipeline = Pipeline(self.kwargs['config_fp'],
                                 run_id,
                                 self.kwargs['uif_path'],
                                 self.kwargs['output_dir'],
                                 self.kwargs['job_id'],
                                 ASSAY_NAME_METAGENOMIC,
                                 lane_number=self.kwargs['lane_number'])

        self.fsr = FailedSamplesRecord(self.kwargs['output_dir'],
                                       self.pipeline.sample_sheet.samples)

        # given run_id, Pipeline should have found the appropriate run_dir.
        type = InstrumentUtils.get_instrument_type(self.pipeline.run_dir)

        self.iseq_run = True if type == 'iSeq' else False

        self.master_qiita_job_id = None

        self.lane_number = self.kwargs['lane_number']
        self.is_restart = bool(self.kwargs['is_restart'])

        if self.is_restart is True:
            self.determine_steps_to_skip()

        self.update = True

        if 'update_qiita' in kwargs:
            if not isinstance(kwargs['update_qiita'], bool):
                raise ValueError("value for 'update_qiita' must be of "
                                 "type bool")

            self.update = kwargs['update_qiita']

    def determine_steps_to_skip(self):
        out_dir = self.pipeline.output_path

        directories_to_check = ['TellReadJob', 'TRIntegrateJob', 'NuQCJob',
                                'FastQCJob', 'GenPrepFileJob']

        for directory in directories_to_check:
            if exists(join(out_dir, directory)):
                if exists(join(out_dir, directory, 'job_completed')):
                    # this step completed successfully.
                    self.skip_steps.append(directory)
                    if exists(join(out_dir, directory,
                                   'post_processing_completed')):
                        self.skip_steps.append('TRIJ_Post_Processing')
                else:
                    # work stopped before this job could be completed.
                    msg = "%s doesn't have job completed" % join(out_dir,
                                                                 directory)
                    raise ValueError(msg)

    def execute_pipeline(self):
        '''
        Executes steps of pipeline in proper sequence.
        :return: None
        '''

        # perform some (re)initialization steps on (re)startup.
        self.pre_check()

        # this is performed even in the event of a restart.
        self.generate_special_map()

        # even if a job is being skipped, it's being skipped because it was
        # determined that it already completed successfully. Hence,
        # increment the status because we are still iterating through them.

        self.update_status("Converting data", 1, 9)

        # convert_raw_to_fastq() now performs its own checking of skip_steps.
        # convert_raw_to_fastq() now performs its own write to fsr reports.
        # This means fsr reports will be accurate even on restarts.
        self.convert_raw_to_fastq()

        self.integrate_results()

        self.generate_sequence_counts()

        self.update_status("Performing quality control", 2, 9)
        self.quality_control()

        self.update_status("Generating reports", 3, 9)
        self.generate_reports()

        self.update_status("Generating preps", 4, 9)
        self.generate_prep_file()

        # moved final component of genprepfilejob outside of object.
        # obtain the paths to the prep-files generated by GenPrepFileJob
        # w/out having to recover full state.
        tmp = join(self.pipeline.output_path, 'GenPrepFileJob', 'PrepFiles')

        self.has_replicates = False

        prep_paths = []
        self.prep_file_paths = defaultdict(list)

        for root, dirs, files in walk(tmp):
            for _file in files:
                # breakup the prep-info-file into segments
                # (run-id, project_qid, other) and cleave
                # the qiita-id from the project_name.
                qid = _file.split('.')[1].split('_')[-1]

                if _file.endswith('.tsv'):
                    _path = abspath(join(root, _file))
                    prep_paths.append(_path)
                    self.prep_file_paths[qid].append(_path)

            for _dir in dirs:
                if _dir == '1':
                    # if PrepFiles contains the '1' directory, then it's a
                    # given that this sample-sheet contains replicates.
                    self.has_replicates = True

        # currently imported from Assay although it is a base method. it
        # could be imported into Workflows potentially, since it is a post-
        # processing step. All pairings of assay and instrument type need to
        # generate prep-info files in the same format.
        self.overwrite_prep_files(prep_paths)

        # for now, simply re-run any line below as if it was a new job, even
        # for a restart. functionality is idempotent, except for the
        # registration of new preps in Qiita. These will simply be removed
        # manually.

        # post-processing steps are by default associated with the Workflow
        # class, since they deal with fastq files and Qiita, and don't depend
        # on assay or instrument type.
        self.update_status("Generating sample information", 5, 9)
        self.sifs = self.generate_sifs()

        # post-processing step.
        self.update_status("Registering blanks in Qiita", 6, 9)
        if self.update:
            self.update_blanks_in_qiita()

        self.update_status("Loading preps into Qiita", 7, 9)
        if self.update:
            self.update_prep_templates()

        # before we load preps into Qiita we need to copy the fastq
        # files n times for n preps and correct the file-paths each
        # prep is pointing to.
        self.load_preps_into_qiita()

        self.fsr.generate_report()

        self.update_status("Generating packaging commands", 8, 9)
        self.generate_commands()

        self.update_status("Packaging results", 9, 9)
        if self.update:
            self.execute_commands()
