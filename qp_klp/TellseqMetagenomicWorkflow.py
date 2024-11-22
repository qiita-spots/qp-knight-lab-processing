
from .SequencingTech import Illumina, TellSeq
from os.path import join, abspath, exists, split
from os import walk, makedirs, listdir
import pandas as pd
from json import dumps
from subprocess import Popen, PIPE
from glob import glob
from shutil import copyfile, rmtree
from sequence_processing_pipeline.Pipeline import Pipeline
from metapool import load_sample_sheet
import logging
from .Assays import Amplicon, Metagenomic, Metatranscriptomic
from .Assays import (METAOMIC_ASSAY_NAMES, ASSAY_NAME_AMPLICON,
                     ASSAY_NAME_METAGENOMIC, ASSAY_NAME_METATRANSCRIPTOMIC)
from .SequencingTech import SEQTECH_NAME_ILLUMINA, SEQTECH_NAME_TELLSEQ
from .FailedSamplesRecord import FailedSamplesRecord
from .Workflows import Workflow, WorkflowError




class TellSeqMetagenomicWorkflow(Workflow, Metagenomic, TellSeq):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # TODO: For now set to False
        self.iseq_run = False

        # TODO: Replace these with frozen set() or similar.
        self.mandatory_attributes = ['qclient', 'uif_path', 'config_fp',
                                     'run_identifier', 'output_dir', 'job_id',
                                     'is_restart']

        self.confirm_mandatory_attributes()

        # second stage initializer that could conceivably be pushed down into
        # specific children requiring specific parameters.
        self.qclient = self.kwargs['qclient']

        self.pipeline = Pipeline(self.kwargs['config_fp'],
                                 self.kwargs['run_identifier'],
                                 self.kwargs['uif_path'],
                                 self.kwargs['output_dir'],
                                 self.kwargs['job_id'],
                                 ASSAY_NAME_METAGENOMIC,
                                 lane_number=self.kwargs['lane_number'])

        self.master_qiita_job_id = None

        self.lane_number = self.kwargs['lane_number']
        self.is_restart = bool(self.kwargs['is_restart'])

        if self.is_restart is True:
            self.determine_steps_to_skip()

        self.update = True

        if 'update_qiita' in kwargs:
            if bool(kwargs['update_qiita']) is False:
                self.update = False


    def determine_steps_to_skip(self):
        out_dir = self.pipeline.output_path

        directories_to_check = ['TellReadJob', 'TRIntegrateJob', 'NuQCJob', 'FastQCJob', 'GenPrepFileJob']

        for directory in directories_to_check:
            if exists(join(out_dir, directory)):
                if exists(join(out_dir, directory, 'job_completed')):
                    # this step completed successfully.
                    self.skip_steps.append(directory)
                else:
                    # work stopped before this job could be completed.
                    rmtree(join(out_dir, directory))

    def execute_pipeline(self):
        '''
        Executes steps of pipeline in proper sequence.
        :param qclient: Qiita client library or equivalent.
        :param update_status: callback function to update status.
        :param update: Set False to prevent updates to Qiita.
        :return: None
        '''
        if not self.is_restart:
            self.pre_check()

        # this is performed even in the event of a restart.
        self.generate_special_map()

        # even if a job is being skipped, it's being skipped because it was
        # determined that it already completed successfully. Hence,
        # increment the status because we are still iterating through them.

        self.update_status("Step 1 of 9: Converting data")

        # this particular convert_raw_to_fastq pipeline can read self.skip_steps
        # itself and determine whether or not to redo TRIntegrate and so on.
        # no need to make that decision here.


        # converting raw data to fastq depends heavily on the instrument
        # used to generate the run_directory. Hence this method is
        # supplied by the instrument mixin.
        if 'TellReadJob' not in self.skip_steps:
            results = self.convert_raw_to_fastq()
            self.fsr_write(results, 'TellReadJob')

        if 'TRNormCountsJob' not in self.skip_steps:
            if self.iseq_run:
                results = self.normalize_counts()
                self.fsr_write(results, 'TRNormCountsJob')

        if 'TRIntegrateJob' not in self.skip_steps:
            results = self.integrate_results()
            self.fsr_write(results, 'TRIntegrateJob')

        # NB: after i_job is completed, there are two optional jobs that
        # can be performed in parallel using the new functionality in Job()
        # class. However we are not using the output from this step right now
        # so we will leave it unimplemented temporarily.

        self.reorganize_results



        self.update_status("Step 2 of 9: Performing quality control")
        self.quality_control()
        raise ValueError("Aborted early after host filtering")

