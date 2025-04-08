from os.path import join
from sequence_processing_pipeline.SeqCountsJob import SeqCountsJob
from functools import partial
import unittest
import pandas as pd
from pandas.testing import assert_frame_equal


class TestSeqCountsJob(unittest.TestCase):
    def setUp(self):
        package_root = "sequence_processing_pipeline"
        self.path = partial(join, package_root, "tests")
        # where 2caa8226-cf69-45a3-bd40-1e90ec3d18d0 is a random qiita job id.
        self.exp = self.path('data', 'tellseq_output', 'integrate_test.sbatch')

        # where 150629_SN1001_0511_AH5L7GBCXX is a run-directory that already
        # exists.
        self.run_dir = self.path('data', 'sample_run_directories',
                                 '150629_SN1001_0511_AH5L7GBCXX')

        self.output_path = self.path('2caa8226-cf69-45a3-bd40-1e90ec3d18d0')

        self.files_to_count_path = self.path("data", "files_to_count.txt")

        self.queue_name = "qiita"
        self.node_count = "1"
        self.wall_time_limit = "1440"
        self.jmem = "8"
        self.modules_to_load = []
        self.qiita_job_id = "2caa8226-cf69-45a3-bd40-1e90ec3d18d0"
        self.cores_per_task = "1"
        self.raw_fastq_dir = join(self.output_path, "TellReadJob", "Full")
        self.max_array_length = 100
        self.exp_sbatch_output = self.path("data", "seq_counts.sbatch")
        self.exp_results = self.path("data", "SeqCounts.csv")
        self.dummy_sample_sheet = self.path("data",
                                            "tellseq_metag_dummy_sample"
                                            "_sheet.csv")

    def test_creation(self):
        def compare_files(obs, exp):
            with open(obs, 'r') as f:
                obs_lines = f.readlines()
                obs_lines = [x.strip() for x in obs_lines]
                obs_lines = [x for x in obs_lines if x != '']

            with open(exp, 'r') as f:
                exp_lines = f.readlines()
                exp_lines = [x.strip() for x in exp_lines]
                exp_lines = [x for x in exp_lines if x != '']

            for obs_line, exp_line in zip(obs_lines, exp_lines):
                self.assertEqual(obs_line, exp_line)

        # test basic good-path
        job = SeqCountsJob(self.run_dir, self.output_path, self.queue_name,
                           self.node_count, self.wall_time_limit, self.jmem,
                           self.modules_to_load, self.qiita_job_id,
                           self.max_array_length, self.files_to_count_path,
                           self.cores_per_task)

        obs = job._generate_job_script()

        compare_files(obs, self.exp_sbatch_output)

        # hack log path so that it points to test data directory rather than
        # the output directory for a run we didn't run().
        job.log_path = self.path("data", "seq_counts_logs")

        obs = pd.read_csv(job._aggregate_counts(self.dummy_sample_sheet),
                          sep=',', dtype='str')
        exp = pd.read_csv(self.exp_results, sep=',', dtype='str')

        # assert_frame_equal will raise an AssertionError if the dfs are not
        # equal. The AssertionError message itself will be the have the
        # best description of the error so return it to the user.
        assert_frame_equal(obs, exp, check_like=True)


if __name__ == '__main__':
    unittest.main()
