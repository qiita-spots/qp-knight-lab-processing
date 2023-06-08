# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from qp_klp.tests.test_step import BaseStepTests, FakeClient
from qp_klp.Metagenomic import Metagenomic
from sequence_processing_pipeline.Pipeline import Pipeline
from os import makedirs
from os.path import join
from os.path import exists


class MetagenomicTests(BaseStepTests):
    def setUp(self):
        super().setUp()

    def test_metagenomic_creation(self):
        # Test base-class creation method, even though base-class will never
        # be instantiated by itself in normal usage.

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

    def test_metagenomic_convert_bcl_to_fastq(self):
        self._create_test_input(1)

        step = Metagenomic(self.pipeline, self.qiita_id, None)
        step.convert_bcl_to_fastq()

    def test_metagenomic_quality_control(self):
        self._create_test_input(2)

        step = Metagenomic(self.pipeline, self.qiita_id, None)
        step.quality_control()

    def atest_generate_commands(self):
        self._create_test_input(4)
        fake_client = FakeClient()

        # confirm metagenomic version is correct. Output should be written to
        # cmds.log file.
        mg_step = Metagenomic(self.pipeline, self.qiita_id, None, 1)

        # need to generate some metadata in order to generate commands.
        mg_step.generate_special_map(fake_client)

        mg_step.generate_commands()

        exp = [(f'cd {self.output_file_path}; '
                'tar zcvf logs-ConvertJob.tgz ConvertJob/logs'),
               (f'cd {self.output_file_path}; '
                'tar zcvf logs-FastQCJob.tgz FastQCJob/logs'),
               (f'cd {self.output_file_path}; '
                'tar zcvf reports-FastQCJob.tgz FastQCJob/fastqc'),
               (f'cd {self.output_file_path}; '
                'tar zcvf logs-GenPrepFileJob.tgz GenPrepFileJob/logs'),
               (f'cd {self.output_file_path}; '
                'tar zcvf prep-files.tgz GenPrepFileJob/PrepFiles'),
               (f'cd {self.output_file_path}; '
                'tar zcvf reports-QCJob.tgz QCJob/NYU_BMS_Melanoma_13059/'
                'fastp_reports_dir'),
               (f'cd {self.output_file_path}; '
                'tar zcvf reports-QCJob.tgz QCJob/Feist_11661/'
                'fastp_reports_dir'),
               (f'cd {self.output_file_path}; '
                'tar zcvf reports-QCJob.tgz QCJob/Gerwick_6123/'
                'fastp_reports_dir'),
               (f'cd {self.output_file_path}; '
                'tar zcvf logs-QCJob.tgz QCJob/logs'),
               (f'cd {self.output_file_path}; '
                'mv *.tgz final_results'),
               (f'cd {self.output_file_path}; '
                'mv FastQCJob/multiqc final_results'),
               (f'cd {self.output_file_path}; '
                'tar zcvf reports-ConvertJob.tgz ConvertJob/'
                'Reports ConvertJob/Logs')]

        tmp = join(self.output_file_path, 'cmds.log')
        self.assertTrue(exists(tmp))

        with open(tmp, 'r') as f:
            obs = f.readlines()
            obs = [x.strip() for x in obs]

        self.assertEqual(obs, exp)
