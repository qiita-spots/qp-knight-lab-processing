from .Protocol import PacBio
from sequence_processing_pipeline.Pipeline import Pipeline
from .Assays import Metagenomic
from .Assays import ASSAY_NAME_METAGENOMIC
from .FailedSamplesRecord import FailedSamplesRecord
from .Workflows import Workflow
import pandas as pd


class PacBioMetagenomicWorkflow(Workflow, Metagenomic, PacBio):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.mandatory_attributes = ['qclient', 'uif_path',
                                     'lane_number', 'config_fp',
                                     'run_identifier', 'output_dir', 'job_id',
                                     'lane_number', 'is_restart']

        self.confirm_mandatory_attributes()

        # second stage initializer that could conceivably be pushed down into
        # specific children requiring specific parameters.
        self.qclient = self.kwargs['qclient']

        self.overwrite_prep_with_original = False
        if 'overwrite_prep_with_original' in self.kwargs:
            self.overwrite_prep_with_original = \
                self.kwargs['overwrite_prep_with_original']
        self.files_regex = 'long_read'
        self.pipeline = Pipeline(self.kwargs['config_fp'],
                                 self.kwargs['run_identifier'],
                                 self.kwargs['uif_path'],
                                 self.kwargs['output_dir'],
                                 self.kwargs['job_id'],
                                 ASSAY_NAME_METAGENOMIC,
                                 lane_number=self.kwargs['lane_number'])

        self.fsr = FailedSamplesRecord(self.kwargs['output_dir'],
                                       self.pipeline.sample_sheet.samples)

        samples = [
            {'barcode': sample['Sample_ID'],
             'sample_name': sample['Sample_Name'],
             'project_name': sample['Sample_Project']}
            for sample in self.pipeline.sample_sheet.samples]
        df = pd.DataFrame(samples)
        sample_list_fp = f"{self.kwargs['output_dir']}/sample_list.tsv"
        df.to_csv(sample_list_fp, sep='\t', index=False)

        self.master_qiita_job_id = self.kwargs['job_id']

        self.lane_number = self.kwargs['lane_number']
        self.is_restart = bool(self.kwargs['is_restart'])

        if self.is_restart is True:
            self.determine_steps_to_skip()

        # this is a convenience member to allow testing w/out updating Qiita.
        self.update = True

        if 'update_qiita' in kwargs:
            if not isinstance(kwargs['update_qiita'], bool):
                raise ValueError("value for 'update_qiita' must be of "
                                 "type bool")

            self.update = kwargs['update_qiita']
