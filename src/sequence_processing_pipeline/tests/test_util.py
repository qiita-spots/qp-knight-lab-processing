import unittest
from sequence_processing_pipeline.util import iter_paired_files


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


if __name__ == '__main__':
    unittest.main()
