from os.path import join

import pandas as pd
from metapool.sample_sheet import (
    BARCODE_ID_KEY,
    LANE_KEY,
    SS_SAMPLE_ID_KEY,
    SS_SAMPLE_PROJECT_KEY,
    SYNDNA_IS_TWISTED_KEY,
    TWIST_ADAPTOR_ID_KEY,
)

from sequence_processing_pipeline.Pipeline import Pipeline

from .Assays import ASSAY_NAME_METAGENOMIC, Metagenomic
from .FailedSamplesRecord import FailedSamplesRecord
from .Protocol import PacBio
from .Workflows import Workflow


class PacBioMetagenomicWorkflow(Workflow, Metagenomic, PacBio):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.mandatory_attributes = [
            "qclient",
            "uif_path",
            "lane_number",
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

        self.overwrite_prep_with_original = False
        if "overwrite_prep_with_original" in self.kwargs:
            self.overwrite_prep_with_original = self.kwargs[
                "overwrite_prep_with_original"
            ]
        self.pipeline = Pipeline(
            self.kwargs["config_fp"],
            self.kwargs["run_identifier"],
            self.kwargs["uif_path"],
            self.kwargs["output_dir"],
            self.kwargs["job_id"],
            ASSAY_NAME_METAGENOMIC,
            lane_number=self.kwargs["lane_number"],
        )

        self.fsr = FailedSamplesRecord(
            self.kwargs["output_dir"], self.pipeline.sample_sheet.samples
        )

        samples = []
        for sample in self.pipeline.sample_sheet.samples:
            tai, sitk = None, None
            if TWIST_ADAPTOR_ID_KEY in sample:
                tai = sample[TWIST_ADAPTOR_ID_KEY]
                if SYNDNA_IS_TWISTED_KEY in sample:
                    sitk = sample[SYNDNA_IS_TWISTED_KEY]
            samples.append(
                {
                    "barcode": sample[BARCODE_ID_KEY],
                    "sample_name": sample[SS_SAMPLE_ID_KEY],
                    "project_name": sample[SS_SAMPLE_PROJECT_KEY],
                    "lane": sample[LANE_KEY],
                    "twist_adaptor_id": tai,
                    "syndna_is_twisted_key": sitk,
                }
            )
        df = pd.DataFrame(samples)
        sample_list_fp = f"{self.kwargs['output_dir']}/sample_list.tsv"
        df.to_csv(sample_list_fp, sep="\t", index=False)

        self.master_qiita_job_id = self.kwargs["job_id"]

        self.lane_number = self.kwargs["lane_number"]
        self.is_restart = bool(self.kwargs["is_restart"])

        if self.is_restart is True:
            self.raw_fastq_files_path = join(self.pipeline.output_path, "ConvertJob")
            self.reports_path = join(self.raw_fastq_files_path, "SeqCounts.csv")
            self.determine_steps_to_skip()

        # this is a convenience member to allow testing w/out updating Qiita.
        self.update = True

        if "update_qiita" in kwargs:
            if not isinstance(kwargs["update_qiita"], bool):
                raise ValueError("value for 'update_qiita' must be of type bool")

            self.update = kwargs["update_qiita"]
