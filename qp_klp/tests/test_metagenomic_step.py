# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from qp_klp.tests.test_step import BaseStepTests
from qp_klp.Metagenomic import Metagenomic


class MetagenomicTests(BaseStepTests):
    def setUp(self):
        super().setUp()

    def test_creation(self):
        # Test base-class creation method, even though base-class will never
        # be instantiated by itself in normal usage.
        self._delete_test_output()

        # TODO: Note we don't do much with this variable yet.
        sn_tid_map_by_project = {}

        with self.assertRaisesRegex(ValueError, "A pipeline object is needed"
                                                " to initialize Step"):
            Metagenomic(None, self.qiita_id, sn_tid_map_by_project, None)

        with self.assertRaisesRegex(ValueError, "A Qiita job-id is needed to "
                                                "initialize Step"):
            Metagenomic(self.pipeline, None, sn_tid_map_by_project, None)

        with self.assertRaisesRegex(ValueError, "sn_tid_map_by_project is "
                                                "needed to initialize Step"):
            Metagenomic(self.pipeline, self.qiita_id, None, None)

        step = Metagenomic(self.pipeline, self.qiita_id,
                           sn_tid_map_by_project, None)

        self.assertIsNotNone(step)

    def test_convert_bcl_to_fastq(self):
        self._delete_test_output()
        self._create_test_input(1)

        sn_tid_map_by_project = {}
        step = Metagenomic(self.pipeline, self.qiita_id,
                           sn_tid_map_by_project, None)

        step.convert_bcl_to_fastq()

    def test_quality_control(self):
        self._delete_test_output()
        self._create_test_input(2)

        sn_tid_map_by_project = {}

        step = Metagenomic(self.pipeline, self.qiita_id,
                           sn_tid_map_by_project, None)
        step.quality_control()
