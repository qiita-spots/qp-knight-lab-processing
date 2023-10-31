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
from os.path import join, split, exists, basename, dirname
from random import randrange
from glob import glob
from shutil import copyfile


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

    def _generate_fake_fastq_files(self, output_path, needs_filtering, max_fastq_size_in_mb):
        if not exists(output_path):
            raise ValueError("%s doesn't exist" % output_path)

        makedirs(output_path, exist_ok=True)

        for i in range(1, 15):
            forward_read_fp = join(output_path, "SAMPLE%d_L001_R1_001.fastq.gz" % i)
            reverse_read_fp = join(output_path, "SAMPLE%d_L001_R2_001.fastq.gz" % i)

            # where 1048576 equals 1MB in bytes
            file_size = randrange(1048576, max_fastq_size_in_mb * 1048576)

            with open(forward_read_fp, 'w') as f:
                # assume 'A' is one byte in size on disk.
                f.write("A" * file_size)

            # for convenience, let reverse read be equal in size to forward read.
            with open(reverse_read_fp, 'w') as f:
                f.write("A" * file_size)

    def test_metagenomic_quality_control(self):
        self._create_test_input(2)

        # dp = '/Users/ccowart/NEW_QC/qp-knight-lab-processing/qp_klp/tests/data/output_dir/ConvertJob'
        dp = 'qp_klp/tests/data/output_dir/ConvertJob'

        d = {'NYU_BMS_Melanoma_13059': {'needs_filtering': False, 'samples': []},
             'Feist_11661': {'needs_filtering': False, 'samples': []},
             'Gerwick_6123': {'needs_filtering': True, 'samples': []}}

        for project_name in d:
            self._generate_fake_fastq_files(join(dp, project_name), d[project_name]['needs_filtering'], 10)

        step = Metagenomic(self.pipeline, self.qiita_id, None)

        # dp2 = '/Users/ccowart/NEW_QC/qp-knight-lab-processing/qp_klp/tests/data/output_dir/NuQCJob'
        dp2 = 'qp_klp/tests/data/output_dir/NuQCJob'

        for project_name in d:
            # after initialization of step object but before run() is called,
            # copy the raw files and rename them into output files to simulate
            # the job script running in SLURM.
            src = join(dp, project_name)
            sub_dir = 'filtered_sequences' if d[project_name]['needs_filtering'] else 'trimmed_sequences'
            dst = join(dp2, project_name, sub_dir)

            makedirs(dst, exist_ok=True)
           
            # TODO: Figure out why other samples are being generated in
            # '../qp_klp/tests/data/output_dir/ConvertJob/NYU_BMS_Melanoma_13059/EP337425B01_SXXX_L001_R2_001.fastq.gz
            foo = glob(src + '/SAMPLE*.fastq.gz')
            for fp in foo:
                _, file_name = split(fp)
                file_name = file_name.replace('.fastq.gz', '.trimmed.fastq.gz')
                file_path = join(dst, file_name)
                copyfile(fp, file_path)

        # execute the quality_control() method, which will call NuQC's run() method in turn.
        step.quality_control()





