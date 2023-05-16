# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from unittest import TestCase
from qp_klp.Step import Step
from sequence_processing_pipeline.Pipeline import Pipeline
from os.path import join, abspath
from functools import partial
from os import makedirs


class BaseStepTests(TestCase):
    def test_creation(self):
        # Test base-class creation method, even though base-class will never
        # be instantiated by itself in normal usage.
        package_root = abspath('./qp_klp')
        self.path = partial(join, package_root, 'tests', 'data')
        self.good_config_file = join(package_root, 'configuration.json')
        self.good_run_id = '211021_A00000_0000_SAMPLE'
        self.good_sample_sheet_path = self.path('good-sample-sheet.csv')
        self.output_file_path = self.path('output_dir')
        self.qiita_id = '077c4da8-74eb-4184-8860-0207f53623be'
        makedirs(self.output_file_path, exist_ok=True)

        pipeline = Pipeline(self.good_config_file, self.good_run_id,
                            self.good_sample_sheet_path, None,
                            self.output_file_path, self.qiita_id,
                            'metagenomic', None)

        # TODO: Note we don't do much with this variable yet.
        sn_tid_map_by_project = {}

        with self.assertRaisesRegex(ValueError, "A pipeline object is needed"
                                                " to initialize Step"):
            Step(None, self.qiita_id, sn_tid_map_by_project, None)

        with self.assertRaisesRegex(ValueError, "A Qiita job-id is needed to "
                                                "initialize Step"):
            Step(pipeline, None, sn_tid_map_by_project, None)

        with self.assertRaisesRegex(ValueError, "sn_tid_map_by_project is "
                                                "needed to initialize Step"):
            Step(pipeline, self.qiita_id, None, None)

    def test_convert_bcl_to_fastq(self):
        pass

    def test_quality_control(self):
        pass

    def test_generate_reports(self):
        pass

    def test_generate_prep_file(self):
        pass

    def test_generate_commands(self):
        pass

    def test_write_commands_to_output_path(self):
        pass

    def test_execute_commands(self):
        pass

    def test_generate_sifs(self):
        pass

    def test_get_prep_file_paths(self):
        pass

