# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from unittest import TestCase
from qp_klp.AmpliconStep import AmpliconStep
from sequence_processing_pipeline.Pipeline import Pipeline
from os.path import join, abspath
from functools import partial
from os import makedirs, chmod
import json
from shutil import rmtree


class AmpliconStepTests(TestCase):
    def setUp(self):
        package_root = abspath('./qp_klp')
        self.path = partial(join, package_root, 'tests', 'data')
        self.good_config_file = join(package_root, 'configuration.json')
        self.good_run_id = '211021_A00000_0000_SAMPLE'
        self.good_sample_sheet_path = self.path('good-sample-sheet.csv')
        self.output_file_path = self.path('output_dir')
        self.qiita_id = '077c4da8-74eb-4184-8860-0207f53623be'
        makedirs(self.output_file_path, exist_ok=True)

        tmp = json.load(open(self.good_config_file, 'r'))['configuration']
        self.config = tmp

    def test_creation(self):
        # Test base-class creation method, even though base-class will never
        # be instantiated by itself in normal usage.

        # dummy used for constructor testing.
        sn_tid_map_by_project = {}

        # create metagenomic pipeline for failure tests.
        metagenomic_pipeline = Pipeline(self.good_config_file,
                                        self.good_run_id,
                                        self.good_sample_sheet_path,
                                        None,
                                        self.output_file_path,
                                        self.qiita_id,
                                        'metagenomic',
                                        None)

        with self.assertRaisesRegex(ValueError, "A pipeline object is needed"
                                                " to initialize Step"):
            AmpliconStep(None, self.qiita_id, sn_tid_map_by_project, None)

        with self.assertRaisesRegex(ValueError, "A Qiita job-id is needed to "
                                                "initialize Step"):
            AmpliconStep(metagenomic_pipeline, None,
                         sn_tid_map_by_project, None)

        with self.assertRaisesRegex(ValueError, "sn_tid_map_by_project is "
                                                "needed to initialize Step"):
            AmpliconStep(metagenomic_pipeline, self.qiita_id, None, None)

        with self.assertRaisesRegex(ValueError, "Cannot instantiate Amplicon"
                                                "Step object from pipeline of"
                                                " type 'metagenomic'"):
            AmpliconStep(metagenomic_pipeline, self.qiita_id,
                         sn_tid_map_by_project, None)
