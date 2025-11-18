from sequence_processing_pipeline.Pipeline import InstrumentUtils, Pipeline

from .Assays import ASSAY_NAME_METAGENOMIC, Metagenomic
from .FailedSamplesRecord import FailedSamplesRecord
from .Protocol import TellSeq
from .Workflows import Workflow


class TellSeqMetagenomicWorkflow(Workflow, Metagenomic, TellSeq):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.mandatory_attributes = [
            "qclient",
            "uif_path",
            "config_fp",
            "run_identifier",
            "output_dir",
            "job_id",
            "is_restart",
        ]

        self.confirm_mandatory_attributes()

        # second stage initializer that could conceivably be pushed down into
        # specific children requiring specific parameters.
        self.qclient = self.kwargs["qclient"]

        run_id = self.kwargs["run_identifier"]

        self.pipeline = Pipeline(
            self.kwargs["config_fp"],
            run_id,
            self.kwargs["uif_path"],
            self.kwargs["output_dir"],
            self.kwargs["job_id"],
            ASSAY_NAME_METAGENOMIC,
            lane_number=self.kwargs["lane_number"],
        )

        self.fsr = FailedSamplesRecord(
            self.kwargs["output_dir"], self.pipeline.sample_sheet.samples
        )

        # given run_id, Pipeline should have found the appropriate run_dir.
        type = InstrumentUtils.get_instrument_type(self.pipeline.run_dir)

        self.iseq_run = True if type == "iSeq" else False

        self.master_qiita_job_id = self.kwargs["job_id"]

        self.lane_number = self.kwargs["lane_number"]
        self.is_restart = bool(self.kwargs["is_restart"])

        self.directories_to_check = [
            "TellReadJob",
            "TRIntegrateJob",
            "NuQCJob",
            "FastQCJob",
            "SeqCountsJob",
            "GenPrepFileJob",
        ]

        if self.is_restart is True:
            self.determine_steps_to_skip()

        self.update = True

        if "update_qiita" in kwargs:
            if not isinstance(kwargs["update_qiita"], bool):
                raise ValueError("value for 'update_qiita' must be of type bool")

            self.update = kwargs["update_qiita"]
