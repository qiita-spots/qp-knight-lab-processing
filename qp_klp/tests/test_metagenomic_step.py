# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from qp_klp.tests.test_step import BaseStepTests
from qp_klp.Metagenomic import Metagenomic
from sequence_processing_pipeline.Pipeline import Pipeline
from os import makedirs


class MetagenomicTests(BaseStepTests):
    def setUp(self):
        super().setUp()

    def test_creation(self):
        # Test base-class creation method, even though base-class will never
        # be instantiated by itself in normal usage.
        self._delete_test_output()

        with self.assertRaisesRegex(ValueError, "A pipeline object is needed"
                                                " to initialize Step"):
            Metagenomic(None, self.qiita_id, None)

        with self.assertRaisesRegex(ValueError, "A Qiita job-id is needed to "
                                                "initialize Step"):
            Metagenomic(self.pipeline, None, None)

        step = Metagenomic(self.pipeline, self.qiita_id, None)

        self.assertIsNotNone(step)

        makedirs(self.output_file_path, exist_ok=True)

        trans_pipeline = Pipeline(None,
                                  self.good_run_id,
                                  self.good_transcript_sheet_path,
                                  None,
                                  self.output_file_path,
                                  self.qiita_id,
                                  Metagenomic.METATRANSCRIPTOMIC_TYPE,
                                  BaseStepTests.CONFIGURATION)

        step = Metagenomic(trans_pipeline, self.qiita_id, None)

        self.assertIsNotNone(step)

    def test_convert_bcl_to_fastq(self):
        self._delete_test_output()
        self._create_test_input(1)

        sn_tid_map_by_project = {}
        step = Metagenomic(self.pipeline, self.qiita_id, None)

        step.convert_bcl_to_fastq()

    def test_quality_control(self):
        self._delete_test_output()
        self._create_test_input(2)

        sn_tid_map_by_project = {}

        step = Metagenomic(self.pipeline, self.qiita_id, None)
        step.quality_control()
