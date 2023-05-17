# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from unittest import TestCase
from qp_klp.MetagenomicStep import MetagenomicStep
from sequence_processing_pipeline.Pipeline import Pipeline
from os.path import join, abspath
from functools import partial
from os import makedirs, chmod
import json
from shutil import rmtree


class MetagenomicStepTests(TestCase):
    def setUp(self):
        package_root = abspath('./qp_klp')
        self.path = partial(join, package_root, 'tests', 'data')
        self.good_config_file = join(package_root, 'configuration.json')
        self.good_run_id = '211021_A00000_0000_SAMPLE'
        self.good_sample_sheet_path = self.path('good-sample-sheet.csv')
        self.output_file_path = self.path('output_dir')
        self.qiita_id = '077c4da8-74eb-4184-8860-0207f53623be'
        makedirs(self.output_file_path, exist_ok=True)

        self.pipeline = Pipeline(self.good_config_file, self.good_run_id,
                                 self.good_sample_sheet_path, None,
                                 self.output_file_path, self.qiita_id,
                                 'metagenomic', None)

        tmp = json.load(open(self.good_config_file, 'r'))['configuration']
        self.config = tmp

    def _create_test_input(self, stage):
        if stage >= 1:
            fake_path = join(self.output_file_path, 'ConvertJob', 'logs')
            makedirs(fake_path, exist_ok=True)
            with open(join(fake_path, 'sbatch'), 'w') as f:
                f.write("#!/bin/sh\necho 'Submitted batch job 9999999'\n")
            chmod(join(fake_path, 'sbatch'), 0o777)

            fake_path = join(abspath('.'), 'sacct')
            with open(fake_path, 'w') as f:
                f.write(
                    "echo '9999999|99999999-9999-9999-9999-999999999999.txt|"
                    "COMPLETED|09:53:41|0:0'")
            chmod(fake_path, 0o777)

            fake_path = join(abspath('.'), 'sbatch')
            with open(fake_path, 'w') as f:
                f.write("echo 'Submitted batch job 9999998\n'")
            chmod(fake_path, 0o777)

        if stage >= 2:
            fake_path = join(self.output_file_path, 'QCJob', 'logs')
            makedirs(fake_path, exist_ok=True)
            with open(join(fake_path, 'sbatch'), 'w') as f:
                f.write("#!/bin/sh\necho 'Submitted batch job 9999999'\n")
            chmod(join(fake_path, 'sbatch'), 0o777)

            exp = {'Feist_11661': ['CDPH-SAL_Salmonella_Typhi_MDL-143',
                                   'CDPH-SAL_Salmonella_Typhi_MDL-144',
                                   'CDPH-SAL_Salmonella_Typhi_MDL-145',
                                   'CDPH-SAL_Salmonella_Typhi_MDL-146',
                                   'CDPH-SAL_Salmonella_Typhi_MDL-147'],
                   'Gerwick_6123': ['3A', '4A', '5B', '6A', '7A'],
                   'NYU_BMS_Melanoma_13059': ['AP581451B02', 'EP256645B01',
                                              'EP112567B02', 'EP337425B01',
                                              'LP127890A01']}
            for project in exp:
                fake_path = join(self.output_file_path, 'ConvertJob', project)
                makedirs(fake_path, exist_ok=True)

                for sample in exp[project]:
                    r1 = join(fake_path, f'{sample}_SXXX_L001_R1_001.fastq.gz')
                    r2 = join(fake_path, f'{sample}_SXXX_L001_R2_001.fastq.gz')

                    for file_path in [r1, r2]:
                        with open(file_path, 'w') as f:
                            f.write("This is a file.")

    def _delete_test_output(self):
        rmtree(self.output_file_path)

    def test_creation(self):
        # Test base-class creation method, even though base-class will never
        # be instantiated by itself in normal usage.
        self._delete_test_output()

        # TODO: Note we don't do much with this variable yet.
        sn_tid_map_by_project = {}

        with self.assertRaisesRegex(ValueError, "A pipeline object is needed"
                                                " to initialize Step"):
            MetagenomicStep(None, self.qiita_id, sn_tid_map_by_project, None)

        with self.assertRaisesRegex(ValueError, "A Qiita job-id is needed to "
                                                "initialize Step"):
            MetagenomicStep(self.pipeline, None, sn_tid_map_by_project, None)

        with self.assertRaisesRegex(ValueError, "sn_tid_map_by_project is "
                                                "needed to initialize Step"):
            MetagenomicStep(self.pipeline, self.qiita_id, None, None)

        step = MetagenomicStep(self.pipeline, self.qiita_id,
                               sn_tid_map_by_project, None)

        self.assertIsNotNone(step)

    def test_convert_bcl_to_fastq(self):
        self._delete_test_output()
        self._create_test_input(1)

        sn_tid_map_by_project = {}
        step = MetagenomicStep(self.pipeline, self.qiita_id,
                               sn_tid_map_by_project, None)

        step.convert_bcl_to_fastq()

    def test_quality_control(self):
        self._delete_test_output()
        self._create_test_input(2)

        sn_tid_map_by_project = {}

        step = MetagenomicStep(self.pipeline, self.qiita_id,
                               sn_tid_map_by_project, None)
        step.quality_control()
