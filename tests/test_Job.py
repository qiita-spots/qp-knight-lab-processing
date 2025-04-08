import unittest
from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import PipelineError
from os.path import abspath, join, dirname, split, isdir, exists
from os import makedirs, chmod, remove
from functools import partial
from shutil import rmtree, copyfile
import re


class TestJob(unittest.TestCase):
    def setUp(self):
        self.remove_these = []

    def tearDown(self):
        for some_path in self.remove_these:
            if isdir(some_path):
                rmtree(some_path)
            else:
                remove(some_path)

    def test_system_call(self):
        package_root = abspath('./sequence_processing_pipeline')
        self.path = partial(join, package_root, 'tests', 'data')

        output_dir = self.path('my_output_dir')

        job = Job(self.path('211021_A00000_0000_SAMPLE'),
                  output_dir, '200nnn_xnnnnn_nnnn_xxxxxxxxxx',
                  ['ls'], 1000, None)

        self.assertTrue(job._which('ls') in ['/bin/ls', '/usr/bin/ls'])

        with self.assertRaises(PipelineError) as e:
            job._file_check('/does/not/exist')

        self.assertRegex(str(e.exception),
                         r"^file '/does/not/exist' does not exist.")

        obs = job._find_files(package_root)
        obs = [re.sub(r'^.*?/sequence_processing_pipeline',
                      r'sequence_processing_pipeline', x) for x in obs]
        self.assertIn('sequence_processing_pipeline/NuQCJob.py', obs)

        callback_results = []

        def my_callback(jid=None, status=None):
            callback_results.append((jid, status))

        obs = job._system_call('ls ' + join(package_root, 'tests', 'bin'),
                               callback=my_callback)

        exp = ['bcl2fastq\nbcl-convert\nfastqc\n',
               'bcl-convert\nbcl2fastq\nfastqc\n']

        self.assertIn(obs['stdout'], exp)
        self.assertEqual(obs['stderr'], '')
        self.assertEqual(obs['return_code'], 0)

        for item in callback_results:
            self.assertTrue(isinstance(item[0], int))
            self.assertIn(item[1], ['RUNNING', 'COMPLETED'])

        with self.assertRaises(PipelineError) as e:
            job._system_call('error')

        self.assertRegex(str(e.exception), (r'^Execute command-line statement'
                                            r' failure:\nCommand: '
                                            r'error\nreturn code:'
                                            r' 127'))

        self.remove_these.append(output_dir)

    def test_group_commands(self):
        package_root = abspath('./sequence_processing_pipeline')
        self.path = partial(join, package_root, 'tests', 'data')

        job = Job(self.path('211021_A00000_0000_SAMPLE'),
                  self.path('my_output_dir'), '200nnn_xnnnnn_nnnn_xxxxxxxxxx',
                  ['ls'], 2, None)

        cmds = list(range(1, 8))
        cmds = [str(x) for x in cmds]
        results = job._group_commands(cmds)
        self.assertEqual(results[0], '1;3;5;7')
        self.assertEqual(results[1], '2;4;6')
        self.assertEqual(len(results), 2)

    def test_extract_project_names_from_fastq_dir(self):
        package_root = abspath('./sequence_processing_pipeline')
        base_path = partial(join, package_root, 'tests', 'data')

        dummy_fastqs = [
            "ConvertJob/NPH_15288/359180337_S27_L007_R1_001.fastq.gz",
            "ConvertJob/NPH_15288/BLANK2_NPH_2_E_S13_L007_R1_001.fastq.gz",
            "ConvertJob/Undetermined_S0_L007_R2_001.fastq.gz",
            ("NuQCJob/only-adapter-filtered/NPH_15288/"
             "359180398_S9_L007_R1_001.fastq.gz"),
            ("NuQCJob/only-adapter-filtered/NPH_15288/"
             "BLANK2_NPH_9_E_S69_L007_R2_001.fastq.gz"),
            ("NuQCJob/NPH_15288/filtered_sequences/"
             "359180396_S7_L007_R2_001.trimmed.fastq.gz")]

        self.remove_these.append(base_path('7b9d7d9c-2cd4-4d54-94ac-'
                                           '40e07a713585'))

        for fastq in dummy_fastqs:
            full_path = base_path('7b9d7d9c-2cd4-4d54-94ac-40e07a713585',
                                  fastq)
            makedirs(dirname(full_path), exist_ok=True)
            with open(full_path, 'w') as f:
                f.write("This is a dummy file.")

        # fake a basic Job() object with enough metadata to test
        # extract_project_names_from_fastq_dir().
        job = Job(base_path('211021_A00000_0000_SAMPLE'),
                  base_path('7b9d7d9c-2cd4-4d54-94ac-40e07a713585'),
                  '200nnn_xnnnnn_nnnn_xxxxxxxxxx', ['ls'], 2, None)

        # results from ConvertJob and NuQCJob should be equal.
        tmp = base_path('7b9d7d9c-2cd4-4d54-94ac-40e07a713585', 'ConvertJob')
        obs = job.extract_project_names_from_fastq_dir(tmp)
        self.assertEqual(obs, ['NPH_15288'])

        tmp = base_path('7b9d7d9c-2cd4-4d54-94ac-40e07a713585', 'NuQCJob')
        obs = job.extract_project_names_from_fastq_dir(tmp)
        self.assertEqual(obs, ['NPH_15288'])

    def test_query_slurm(self):
        package_root = abspath('./sequence_processing_pipeline')
        base_path = partial(join, package_root, 'tests', 'data')

        # set up a fake job so that we can test the query_jobs() method.
        # it doesn't matter what the parameters are so long as the job
        # passes initialization.
        job = Job(base_path('211021_A00000_0000_SAMPLE'),
                  base_path('7b9d7d9c-2cd4-4d54-94ac-40e07a713585'),
                  '200nnn_xnnnnn_nnnn_xxxxxxxxxx', ['ls'], 2, None)

        # locate python binary path
        # we have a python script called fake_squeue.py that can simulate
        # repeated calls to squeue. It does this by generating a fake random
        # set of array job ids for each job id passed to it and records their
        # state in my_state.json. Each array job is set to change state from
        # RUNNING to either COMPLETED or FAILED between three to seven squeue
        # calls. The choice of which job-ids will succeed or fail, as is which
        # individual array-ids will succeed or fail is random.
        python_path = split(job._which('python'))[0]
        squeue_path = join(python_path, 'squeue')
        foo = join(package_root, 'scripts', 'fake_squeue.py')

        # place the fake squeue file in a place that's known to be in the
        # PATH. Make sure this file is removed after this test is complete.
        # Also make sure the saved state file is removed.
        copyfile(foo, squeue_path)
        chmod(squeue_path, 0o755)
        self.remove_these.append(squeue_path)
        self.remove_these.append(join(package_root, 'scripts',
                                      'my_state.json'))

        job_ids = ['1234567', '1234568', '1234569', '1234570']
        jobs = job._query_slurm(job_ids)

        # jobs is a dictionary of unique array_ids and/or job-ids for non-
        # array jobs. The faked squeue reports anywhere between five and
        # fifteen array jobs for a given job-id. After the first invocation
        # all processes should be in the 'RUNNING' state.
        # e.g.: "1234567_1": "RUNNING"

        for j in jobs:
            self.assertEqual(jobs[j], 'RUNNING')
            if '_' in j:
                jid, aid = j.split('_')
            else:
                jid = j
                aid = None

            # assert the job id component of the array-id is a valid job id.
            self.assertIn(jid, job_ids)

            if aid:
                # assert the array-id component of the array-id is between 0
                # and 15 as defined in the fake squeue script.
                aid = int(aid)
                self.assertLess(aid, 15)
                self.assertGreaterEqual(aid, 0)

    def test_query_slurm_single_job(self):
        # perform test_query_slurm() but with a single job only.
        package_root = abspath('./sequence_processing_pipeline')
        base_path = partial(join, package_root, 'tests', 'data')

        # set up a fake job so that we can test the query_jobs() method.
        # it doesn't matter what the parameters are so long as the job
        # passes initialization.
        job = Job(base_path('211021_A00000_0000_SAMPLE'),
                  base_path('7b9d7d9c-2cd4-4d54-94ac-40e07a713585'),
                  '200nnn_xnnnnn_nnnn_xxxxxxxxxx', ['ls'], 2, None)

        # locate python binary path
        # we have a python script called fake_squeue.py that can simulate
        # repeated calls to squeue. It does this by generating a fake random
        # set of array job ids for each job id passed to it and records their
        # state in my_state.json. Each array job is set to change state from
        # RUNNING to either COMPLETED or FAILED between three to seven squeue
        # calls. The choice of which job-ids will succeed or fail, as is which
        # individual array-ids will succeed or fail is random.
        python_path = split(job._which('python'))[0]
        squeue_path = join(python_path, 'squeue')
        foo = join(package_root, 'scripts', 'fake_squeue.py')

        # place the fake squeue file in a place that's known to be in the
        # PATH. Make sure this file is removed after this test is complete.
        # Also make sure the saved state file is removed.
        copyfile(foo, squeue_path)
        chmod(squeue_path, 0o755)
        self.remove_these.append(squeue_path)
        self.remove_these.append(join(package_root, 'scripts',
                                      'my_state.json'))

        job_ids = ['1234567']
        jobs = job._query_slurm(job_ids)

        # jobs is a dictionary of unique array_ids and/or job-ids for non-
        # array jobs. The faked squeue reports anywhere between five and
        # fifteen array jobs for a given job-id. After the first invocation
        # all processes should be in the 'RUNNING' state.
        # e.g.: "1234567_1": "RUNNING"

        for j in jobs:
            self.assertEqual(jobs[j], 'RUNNING')
            if '_' in j:
                jid, aid = j.split('_')
            else:
                jid = j
                aid = None

            # assert the job id component of the array-id is a valid job id.
            self.assertIn(jid, job_ids)

            if aid:
                # assert the array-id component of the array-id is between 0
                # and 15 as defined in the fake squeue script.
                aid = int(aid)
                self.assertLess(aid, 15)
                self.assertGreaterEqual(aid, 0)

    def test_wait_on_job_ids(self):
        package_root = abspath('./sequence_processing_pipeline')
        base_path = partial(join, package_root, 'tests', 'data')

        job = Job(base_path('211021_A00000_0000_SAMPLE'),
                  base_path('7b9d7d9c-2cd4-4d54-94ac-40e07a713585'),
                  '200nnn_xnnnnn_nnnn_xxxxxxxxxx', ['ls'], 2, None)

        python_path = split(job._which('python'))[0]
        squeue_path = join(python_path, 'squeue')
        foo = join(package_root, 'scripts', 'fake_squeue.py')
        copyfile(foo, squeue_path)
        chmod(squeue_path, 0o755)
        self.remove_these.append(squeue_path)
        self.remove_these.append(join(package_root, 'scripts',
                                      'my_state.json'))

        job_ids = ['1', '2', '3', '4']

        # to shorten the test time, set polling_interval_in_seconds to be
        # lower than one minute.
        Job.polling_interval_in_seconds = 10
        results = job.wait_on_job_ids(job_ids)

        # calling query_slurm one more time after wait_on_job_ids() is called
        # will technically advance the counter one more, which means that this
        # doesn't confirm that wait_on_job_ids() doesn't return before EVERY
        # single job is either COMPLETED or FAILED. However it does confirm
        # that wait_on_job_ids() doesn't return once the FIRST completed array
        # job is either COMPLETED or FAILED while others are still RUNNING.
        # This was previously an issue.
        obs = job._query_slurm(job_ids)

        for array_id in obs:
            state = obs[array_id]
            # w/out relying on states defined in Job, simply confirm all are
            # either COMPLETED or FAILED.
            self.assertIn(state, ['COMPLETED', 'FAILED'])

        # since wait_on_job_ids() now returns the same data structure as
        # query_slurm(), they should be equal.
        self.assertDictEqual(obs, results)

    def test_mark_completed_commands(self):
        package_root = abspath('./sequence_processing_pipeline')
        self.path = partial(join, package_root, 'tests', 'data')

        job = Job(self.path('211021_A00000_0000_SAMPLE'),
                  self.path('my_output_dir'), '200nnn_xnnnnn_nnnn_xxxxxxxxxx',
                  ['ls'], 2, None)

        fp = join(job.output_path, 'job_completed')

        self.assertFalse(exists(fp))

        job.mark_job_completed()

        self.assertTrue(exists(fp))


if __name__ == '__main__':
    unittest.main()
