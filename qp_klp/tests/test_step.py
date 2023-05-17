# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
import os
from unittest import TestCase
from qp_klp.Step import Step
from sequence_processing_pipeline.Pipeline import Pipeline
from os.path import join, abspath
from functools import partial
from os import makedirs, chmod
import json


class BaseStepTests(TestCase):
    def setUp(self):
        package_root = abspath('./qp_klp')
        print("ROOT: %s" % package_root)

        self.path = partial(join, package_root, 'tests', 'data')

        print("SELF.PATH: %s" % self.path())

        self.good_config_file = join(package_root, 'configuration.json')
        self.good_run_id = '211021_A00000_0000_SAMPLE'
        self.good_sample_sheet_path = self.path('good-sample-sheet.csv')
        self.output_file_path = self.path('output_dir')

        print("OUTPUT PATH: %s" % self.output_file_path)
        self.qiita_id = '077c4da8-74eb-4184-8860-0207f53623be'
        makedirs(self.output_file_path, exist_ok=True)

        self.pipeline = Pipeline(self.good_config_file, self.good_run_id,
                                 self.good_sample_sheet_path, None,
                                 self.output_file_path, self.qiita_id,
                                 'metagenomic', None)

        tmp = json.load(open(self.good_config_file, 'r'))['configuration']
        self.config = tmp

    def test_creation(self):
        # Test base-class creation method, even though base-class will never
        # be instantiated by itself in normal usage.

        # TODO: Note we don't do much with this variable yet.
        sn_tid_map_by_project = {}

        with self.assertRaisesRegex(ValueError, "A pipeline object is needed"
                                                " to initialize Step"):
            Step(None, self.qiita_id, sn_tid_map_by_project, None)

        with self.assertRaisesRegex(ValueError, "A Qiita job-id is needed to "
                                                "initialize Step"):
            Step(self.pipeline, None, sn_tid_map_by_project, None)

        with self.assertRaisesRegex(ValueError, "sn_tid_map_by_project is "
                                                "needed to initialize Step"):
            Step(self.pipeline, self.qiita_id, None, None)

        Step(self.pipeline, self.qiita_id, sn_tid_map_by_project, None)

    def test_convert_bcl_to_fastq(self):
        sn_tid_map_by_project = {}
        step = Step(self.pipeline, self.qiita_id, sn_tid_map_by_project, None)

        fake_path = join(self.output_file_path, 'ConvertJob', 'logs')
        print("OUTPUT PATH2: %s" % self.output_file_path)

        makedirs(fake_path, exist_ok=True)

        with open(join(fake_path, 'sbatch'), 'w') as f:
            f.write("#!/bin/sh\necho 'Submitted batch job 9999999'\n")

        chmod(join(fake_path, 'sbatch'), 0o777)

        print(os.listdir(fake_path))

        fake_path = join(abspath('.'), 'sacct')

        with open(fake_path, 'w') as f:
            f.write("echo '9999999|99999999-9999-9999-9999-999999999999.txt|"
                    "COMPLETED|09:53:41|0:0'")

        chmod(fake_path, 0o777)

        fake_path = join(abspath('.'), 'sbatch')

        with open(fake_path, 'w') as f:
            f.write("echo 'Submitted batch job 9999998\n'")

        chmod(fake_path, 0o777)

        step._convert_bcl_to_fastq(self.config['bcl-convert'],
                                   self.good_sample_sheet_path)

    def test_quality_control(self):
        # sn_tid_map_by_project = {}
        # step = Step(self.pipeline, self.qiita_id, sn_tid_map_by_project,
        # None)
        # step._quality_control(self.config['qc'], self.good_sample_sheet_path)
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
