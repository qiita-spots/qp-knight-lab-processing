import unittest
from sequence_processing_pipeline.util import (
    iter_paired_files, determine_orientation)


class TestUtil(unittest.TestCase):
    def test_iter_paired_files(self):
        # tuples of randomly ordered fastq files and thier expected
        # sorted and organized output from iter_paired_files().

        # underscore filenames updated to require '_001.fastq.gz'.
        # legacy dot filenames test remains as-is.
        tests = [(['b_R2_001.fastq.gz', 'a_R1_001.fastq.gz',
                   'a_R2_001.fastq.gz', 'b_R1_001.fastq.gz'],
                  [('a_R1_001.fastq.gz', 'a_R2_001.fastq.gz'),
                   ('b_R1_001.fastq.gz', 'b_R2_001.fastq.gz')]),
                 (['a.R1.foo', 'b.R2.bar', 'a.R2.baz', 'b.R1.bing'],
                  [('a.R1.foo', 'a.R2.baz'), ('b.R1.bing', 'b.R2.bar')])]

        for files, exp in tests:
            obs = list(iter_paired_files(files))
            self.assertEqual(obs, exp)

    def test_iter_paired_files_uneven(self):
        files = ['a', 'b', 'c']
        with self.assertRaisesRegex(ValueError, "not paired"):
            list(iter_paired_files(files))

    def test_iter_paired_files_no_match(self):
        files = ['a-R1-foo', 'a-R2-foo']
        with self.assertRaisesRegex(ValueError, "Unable to"):
            list(iter_paired_files(files))

    def test_iter_paired_files_bad_pair(self):
        files = ['a.R1.foo', 'a_R2_bar']
        with self.assertRaisesRegex(ValueError, "Cannot find"):
            list(iter_paired_files(files))

    def test_iter_paired_files_mismatch_prefix(self):
        files = ['a_R1_001.fastq.gz', 'ab_R2_001.fastq.gz']
        with self.assertRaisesRegex(ValueError, "Mismatch prefixes"):
            list(iter_paired_files(files))

        files = ['/foo/bar/a_R1_001.fastq.gz', '/foo/bar/ab_R2_001.fastq.gz']
        with self.assertRaisesRegex(ValueError, "Mismatch prefixes"):
            list(iter_paired_files(files))

    def test_determine_orientation(self):
        test_names = [
            # single additional occurrence: R1
            ("ABC_7_04_1776_R1_SRE_S3_L007_R1_001.trimmed.fastq.gz", "R1"),
            ("ABC_7_04_1776_R1_SRE_S3_L007_R2_001.trimmed.fastq.gz", "R2"),
            ("ABC_7_04_1776_R1_SRE_S3_L007_I1_001.trimmed.fastq.gz", "I1"),
            ("ABC_7_04_1776_R1_SRE_S3_L007_I2_001.trimmed.fastq.gz", "I2"),

            # test w/dots.
            ("ABC_7_04_1776.R1.SRE_S3_L007.R1.001.trimmed.fastq.gz", "R1"),
            ("ABC_7_04_1776.R1.SRE_S3_L007.R2.001.trimmed.fastq.gz", "R2"),
            ("ABC_7_04_1776.R1.SRE_S3_L007.I1.001.trimmed.fastq.gz", "I1"),
            ("ABC_7_04_1776.R1.SRE_S3_L007.I2.001.trimmed.fastq.gz", "I2"),

            # single additional occurrence: R2
            ("ABC_7_04_1776_R2_SRE_S3_L007_R1_001.trimmed.fastq.gz", "R1"),
            ("ABC_7_04_1776_R2_SRE_S3_L007_R2_001.trimmed.fastq.gz", "R2"),
            ("ABC_7_04_1776_R2_SRE_S3_L007_I1_001.trimmed.fastq.gz", "I1"),
            ("ABC_7_04_1776_R2_SRE_S3_L007_I2_001.trimmed.fastq.gz", "I2"),

            # single additional occurrence: In
            ("ABC_7_04_1776_I2_SRE_S3_L007_R1_001.trimmed.fastq.gz", "R1"),
            ("ABC_7_04_1776_I1_SRE_S3_L007_R2_001.trimmed.fastq.gz", "R2"),
            ("ABC_7_04_1776_I2_SRE_S3_L007_I1_001.trimmed.fastq.gz", "I1"),
            ("ABC_7_04_1776_I1_SRE_S3_L007_I2_001.trimmed.fastq.gz", "I2"),

            # no additional occurrences
            ("ABC_7_04_1776_SRE_S3_L007_R1_001.trimmed.fastq.gz", "R1"),
            ("ABC_7_04_1776_SRE_S3_L007_R2_001.trimmed.fastq.gz", "R2"),
            ("ABC_7_04_1776_SRE_S3_L007_I1_001.trimmed.fastq.gz", "I1"),
            ("ABC_7_04_1776_SRE_S3_L007_I2_001.trimmed.fastq.gz", "I2"),

            # two additional occurrences
            ("ABC_7_04_1776_I2_SRE.R1.S3_L007_R1_001.trimmed.fastq.gz", "R1"),
            ("ABC_7_04_1776_I1_SRE.R1.S3_L007_R2_001.trimmed.fastq.gz", "R2"),
            ("ABC_7_04_1776_I2_SRE.R1.S3_L007_I1_001.trimmed.fastq.gz", "I1"),
            ("ABC_7_04_1776_I1_SRE.R1.S3_L007_I2_001.trimmed.fastq.gz", "I2"),
        ]

        for file_name, exp in test_names:
            self.assertEqual(determine_orientation(file_name), exp)


if __name__ == '__main__':
    unittest.main()
