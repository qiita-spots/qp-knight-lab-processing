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


class StandardAmpliconWorkflow(Workflow, Amplicon, Illumina):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

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
                                 ASSAY_NAME_AMPLICON,
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

        directories_to_check = ['ConvertJob', 'NuQCJob', 'FastQCJob', 'GenPrepFileJob']

        for directory in directories_to_check:
            if exists(join(out_dir, directory)):
                if exists(join(out_dir, directory, 'job_completed')):
                    # this step completed successfully.
                    self.skip_steps.append(directory)
                else:
                    # work stopped before this job could be completed.
                    rmtree(join(out_dir, directory))

    def execute_pipeline(self):
        """
        Executes steps of pipeline in proper sequence.
        :param qclient: Qiita client library or equivalent.
        :param update_status: callback function to update status.
        :param update: Set False to prevent updates to Qiita.
        :return: None
        """
        # this is performed even in the event of a restart.
        self.generate_special_map()

        # even if a job is being skipped, it's being skipped because it was
        # determined that it already completed successfully. Hence,
        # increment the status because we are still iterating through them.

        self.update_status("Step 1 of 9: Converting data")

        if 'ConvertJob' not in self.skip_steps:
            results = self.convert_raw_to_fastq()
            self.fsr_write(results, 'ConvertJob')

        self.update_status("Step 2 of 9: Performing quality control")

        # Amplicon workflows do not perform human filtering (quality control)
        # as the fastq files will not be demuxed. Human filtering will instead
        # be performed downstream by Qiita. In this case, we don't seek to
        # perform quality control. Instead we want to simulate the output of
        # NuQCJob as needed for downstream steps like FastQCJob to be able
        # to process w/out modification.

        if 'NuQCJob' not in self.skip_steps:
            results = self.post_process_raw_fastq_output()
            self.fsr_write(results, 'NuQCJob')
