from os.path import join
from sequence_processing_pipeline.TellReadJob import TellReadJob
from functools import partial
import unittest


class TestTellReadJob(unittest.TestCase):
    def setUp(self):
        package_root = "sequence_processing_pipeline"
        self.path = partial(join, package_root, "tests")
        # where 2caa8226-cf69-45a3-bd40-1e90ec3d18d0 is a random qiita job id.
        self.obs = self.path('2caa8226-cf69-45a3-bd40-1e90ec3d18d0',
                             'TellReadJob', 'tellread_test.sbatch')
        self.exp = self.path('data', 'tellseq_output', 'tellread_test.sbatch')

        # where 150629_SN1001_0511_AH5L7GBCXX is a run-directory that already
        # exists.
        self.run_dir = self.path('data', 'sample_run_directories',
                                 '150629_SN1001_0511_AH5L7GBCXX')

        self.output_path = self.path('2caa8226-cf69-45a3-bd40-1e90ec3d18d0')

        self.sample_sheet_path = self.path('data',
                                           'tellseq_metag_dummy_sample_'
                                           'sheet.csv')

        self.queue_name = "qiita"
        self.node_count = "1"
        self.wall_time_limit = "96:00:00"
        self.jmem = "16"
        self.modules_to_load = ["singularity_3.6.4"]
        self.qiita_job_id = "2caa8226-cf69-45a3-bd40-1e90ec3d18d0"
        self.label = "150629_SN1001_0511_AH5L7GBCXX-test"
        self.reference_base = ""
        self.reference_map = ""
        self.tmp1_path = join(self.output_path, "TellReadJob", "output",
                              "tmp1")
        # reflects location of script on host.
        self.sing_script_path = ("$HOME/qiita-spots/tellread-release-novaseqX/"
                                 "run_tellread_sing.sh")
        self.lane = "1"
        self.cores_per_task = "4"

    def test_creation(self):
        # test basic good-path
        job = TellReadJob(self.run_dir, self.output_path,
                          self.sample_sheet_path, self.queue_name,
                          self.node_count, self.wall_time_limit,
                          self.jmem, self.modules_to_load, self.qiita_job_id,
                          self.reference_base, self.reference_map,
                          self.sing_script_path, self.cores_per_task)

        job._generate_job_script()

        with open(self.obs, 'r') as f:
            obs_lines = f.readlines()

        with open(self.exp, 'r') as f:
            exp_lines = f.readlines()

        for obs_line, exp_line in zip(obs_lines, exp_lines):
            self.assertEqual(obs_line, exp_line)


if __name__ == '__main__':
    unittest.main()
