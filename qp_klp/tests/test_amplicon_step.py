# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from qp_klp.Amplicon import Amplicon
from sequence_processing_pipeline.Pipeline import Pipeline
from os.path import join, abspath
from qp_klp.tests.test_step import BaseStepTests


class AmpliconTests(BaseStepTests):
    def setUp(self):
        super().setUp()
        self.good_mapping_file = join(abspath('./qp_klp'), 'tests', 'data',
                                      'good-mapping-file.txt')

    def test_creation(self):
        # Test base-class creation method, even though base-class will never
        # be instantiated by itself in normal usage.

        # dummy used for constructor testing.
        sn_tid_map_by_project = {}

        with self.assertRaisesRegex(ValueError, "A pipeline object is needed"
                                                " to initialize Step"):
            Amplicon(None, self.qiita_id, sn_tid_map_by_project, None)

        # create amplicon pipeline for failure tests.
        amplicon_pipeline = Pipeline(None,
                                     self.good_run_id,
                                     None,
                                     self.good_mapping_file,
                                     self.output_file_path,
                                     self.qiita_id,
                                     'amplicon',
                                     AmpliconTests.CONFIGURATION)

        with self.assertRaisesRegex(ValueError, "A Qiita job-id is needed to "
                                                "initialize Step"):
            Amplicon(amplicon_pipeline, None,
                     sn_tid_map_by_project, None)

        with self.assertRaisesRegex(ValueError, "sn_tid_map_by_project is "
                                                "needed to initialize Step"):
            Amplicon(amplicon_pipeline, self.qiita_id, None, None)

        # create metagenomic pipeline for final test.
        metagenomic_pipeline = Pipeline(None,
                                        self.good_run_id,
                                        self.good_sample_sheet_path,
                                        None,
                                        self.output_file_path,
                                        self.qiita_id,
                                        'metagenomic',
                                        AmpliconTests.CONFIGURATION)

        with self.assertRaisesRegex(ValueError, "Cannot create an Amplicon run"
                                                " using a metagenomic-"
                                                "configured Pipeline."):
            Amplicon(metagenomic_pipeline, self.qiita_id,
                     sn_tid_map_by_project, None)
