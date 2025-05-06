from .Protocol import Illumina
from sequence_processing_pipeline.Pipeline import Pipeline
from .Assays import Amplicon
from .Assays import ASSAY_NAME_AMPLICON
from .Workflows import Workflow


class StandardAmpliconWorkflow(Workflow, Amplicon, Illumina):
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

        self.pipeline = Pipeline(self.kwargs['config_fp'],
                                 self.kwargs['run_identifier'],
                                 self.kwargs['uif_path'],
                                 self.kwargs['output_dir'],
                                 self.kwargs['job_id'],
                                 ASSAY_NAME_AMPLICON,
                                 # amplicon runs always use lane 1.
                                 lane_number=1)

        # NB: Amplicon workflows don't have failed samples records because
        # the fastq files are not demultiplexed.

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
