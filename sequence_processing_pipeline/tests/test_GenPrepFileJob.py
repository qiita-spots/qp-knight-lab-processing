from sequence_processing_pipeline.GenPrepFileJob import GenPrepFileJob
from os.path import join
from os import makedirs
import unittest
import shutil
import re


class TestGenPrepFileJob(unittest.TestCase):
    def setUp(self):
        self.package_root = 'sequence_processing_pipeline'
        self.qiita_job_id = 'b197f317-1c06-4619-9af3-65721149c1e8'
        self.working_directory_root = join(self.package_root,
                                           self.qiita_job_id)
        try:
            shutil.rmtree(self.working_directory_root)
        except FileNotFoundError:
            # Clean up test directory just in case
            pass
        makedirs(self.working_directory_root)
        self.run_id = '210518_A00953_0305_TEST'
        self.run_dir = join(self.working_directory_root, self.run_id)
        self.convert_job_path = join(self.run_dir, 'ConvertJob')
        self.reports_path = join(self.convert_job_path, 'Reports')
        makedirs(self.reports_path, exist_ok=True)
        self.qc_job_path = join(self.run_dir, 'QCJob')
        self.project_list = ['Project1']
        makedirs(join(self.qc_job_path, self.project_list[0],
                      'filtered_sequences'))
        makedirs(join(self.qc_job_path, self.project_list[0],
                      'fastp_reports_dir', 'json'))

    def tearDown(self):
        shutil.rmtree(self.working_directory_root)

    def test_creation(self):
        sample_sheet_path = join(self.package_root, 'tests', 'data',
                                 'good-sample-sheet.csv')

        job = GenPrepFileJob(self.run_dir,
                             self.convert_job_path,
                             self.qc_job_path,
                             join(self.working_directory_root, 'OutputPath'),
                             sample_sheet_path,
                             'seqpro',
                             [],
                             'abcdabcdabcdabcdabcdabcdabcdabcd',
                             self.reports_path)

        results = job._system_call(f'find {self.run_dir}')
        lines = results['stdout'].split('\n')
        lines = [re.sub(r'^.*?sequence_processing_pipeline\/', '', x)
                 for x in lines]
        obs = [x for x in lines if x]
        obs.sort()

        exp = ['b197f317-1c06-4619-9af3-65721149c1e8/210518_A00953_0305_TEST',
               'b197f317-1c06-4619-9af3-65721149c1e8/210518_A00953_0305_TEST/'
               'ConvertJob',
               'b197f317-1c06-4619-9af3-65721149c1e8/210518_A00953_0305_TEST/'
               'ConvertJob/Reports',
               'b197f317-1c06-4619-9af3-65721149c1e8/210518_A00953_0305_TEST/'
               'QCJob',
               'b197f317-1c06-4619-9af3-65721149c1e8/210518_A00953_0305_TEST/'
               'QCJob/Project1',
               'b197f317-1c06-4619-9af3-65721149c1e8/210518_A00953_0305_TEST/'
               'QCJob/Project1/fastp_reports_dir',
               'b197f317-1c06-4619-9af3-65721149c1e8/210518_A00953_0305_TEST/'
               'QCJob/Project1/fastp_reports_dir/json',
               'b197f317-1c06-4619-9af3-65721149c1e8/210518_A00953_0305_TEST/'
               'QCJob/Project1/filtered_sequences']

        self.assertEquals(obs, exp)

    def test_get_prep_file_paths(self):
        self.maxDiff = None
        sample_sheet_path = join(self.package_root, 'tests', 'data',
                                 'good-sample-sheet.csv')

        # create a sample Job object
        job = GenPrepFileJob(self.run_dir,
                             self.convert_job_path,
                             self.qc_job_path,
                             join(self.working_directory_root, 'OutputPath'),
                             sample_sheet_path,
                             'seqpro',
                             [],
                             'abcdabcdabcdabcdabcdabcdabcdabcd',
                             self.reports_path)

        # We cannot run the object and test the output that is returned from
        # seqpro, but we can test the helper method against canned stdout and
        # verify the results.
        stdout = ('1111\tmetagenomics_pooling_notebook/metapool/tests/VFTEST/'
                  '200318_A00953_0082_AH5TWYDSXY.Project_1111.1.tsv\n'
                  '1111\tmetagenomics_pooling_notebook/metapool/tests/VFTEST/'
                  '200318_A00953_0082_AH5TWYDSXY.Project_1111.3.tsv\n'
                  '666\tmetagenomics_pooling_notebook/metapool/tests/VFTEST/'
                  '200318_A00953_0082_AH5TWYDSXY.Trojecp_666.3.tsv')

        obs = job._get_prep_file_paths(stdout)

        exp = {'1111': [('metagenomics_pooling_notebook/metapool/tests/VFTEST'
                         '/200318_A00953_0082_AH5TWYDSXY.Project_1111.1.tsv'),
                        ('metagenomics_pooling_notebook/metapool/tests/VFTEST'
                         '/200318_A00953_0082_AH5TWYDSXY.Project_1111.3.tsv')],
               '666': [('metagenomics_pooling_notebook/metapool/tests/VFTEST'
                        '/200318_A00953_0082_AH5TWYDSXY.Trojecp_666.3.tsv')]}

        self.assertDictEqual(obs, exp)


class TestReplication(unittest.TestCase):
    def setUp(self):
        self.package_root = 'sequence_processing_pipeline'
        self.qiita_job_id = 'b197f317-1c06-4619-9af3-65721149c1e8'
        self.working_directory_root = join(self.package_root,
                                           self.qiita_job_id)
        try:
            shutil.rmtree(self.working_directory_root)
        except FileNotFoundError:
            # Clean up test directory just in case
            pass
        makedirs(self.working_directory_root)
        self.run_id = '210518_A00953_0305_TEST'
        self.run_dir = join(self.working_directory_root, self.run_id)
        self.convert_job_path = join(self.run_dir, 'ConvertJob')
        self.reports_path = join(self.convert_job_path, 'Reports')
        makedirs(self.reports_path, exist_ok=True)
        self.qc_job_path = join(self.run_dir, 'QCJob')
        self.project_list = ['Project1']
        makedirs(join(self.qc_job_path, self.project_list[0],
                      'filtered_sequences'))
        makedirs(join(self.qc_job_path, self.project_list[0],
                      'fastp_reports_dir', 'json'))

    def tearDown(self):
        shutil.rmtree(self.working_directory_root)

    def test_sample_sheet_replicate_file_creation(self):
        # since the GenPrepFileJob tests do not run seqpro (this is done in
        # qp-klp tests), the main objective of these tests is to confirm that
        # given a known good input, the right methods from metapool are called
        # and the correct number of output files are created.
        sample_sheet_path = join(self.package_root, 'tests', 'data',
                                 'good_sheet_w_replicates.csv')

        job = GenPrepFileJob(self.run_dir,
                             self.convert_job_path,
                             self.qc_job_path,
                             join(self.working_directory_root, 'OutputPath'),
                             sample_sheet_path,
                             'seqpro',
                             [],
                             'abcdabcdabcdabcdabcdabcdabcdabcd',
                             self.reports_path)

        exp = [['seqpro', '--verbose',
                ('sequence_processing_pipeline/b197f317-1c06-4619-9af3-'
                 '65721149c1e8/OutputPath/GenPrepFileJob/210518_A00953_0305'
                 '_TEST'),
                ('"sequence_processing_pipeline/b197f317-1c06-4619-9af3-'
                 '65721149c1e8/OutputPath/GenPrepFileJob/'
                 'replicate_sheet_1.csv"'),
                ('sequence_processing_pipeline/b197f317-1c06-4619-9af3-'
                 '65721149c1e8/OutputPath/GenPrepFileJob/PrepFiles/1')],
               ['seqpro', '--verbose',
                ('sequence_processing_pipeline/b197f317-1c06-4619-9af3-'
                 '65721149c1e8/OutputPath/GenPrepFileJob/210518_A00953_0305'
                 '_TEST'),
                ('"sequence_processing_pipeline/b197f317-1c06-4619-9af3-'
                 '65721149c1e8/OutputPath/GenPrepFileJob/'
                 'replicate_sheet_2.csv"'),
                ('sequence_processing_pipeline/b197f317-1c06-4619-9af3-'
                 '65721149c1e8/OutputPath/GenPrepFileJob/PrepFiles/2')],
               ['seqpro', '--verbose',
                'sequence_processing_pipeline/b197f317-1c06-4619-9af3-'
                '65721149c1e8/OutputPath/GenPrepFileJob/210518_A00953_0305'
                '_TEST',
                '"sequence_processing_pipeline/b197f317-1c06-4619-9af3-'
                '65721149c1e8/OutputPath/GenPrepFileJob/'
                'replicate_sheet_3.csv"',
                'sequence_processing_pipeline/b197f317-1c06-4619-9af3-'
                '65721149c1e8/OutputPath/GenPrepFileJob/PrepFiles/3']]

        self.assertEqual(job.commands, exp)

    def test_pre_prep_replicate_file_creation(self):
        # since the GenPrepFileJob tests do not run seqpro (this is done in
        # qp-klp tests), the main objective of these tests is to confirm that
        # given a known good input, the right methods from metapool are called
        # and the correct number of output files are created.
        sample_sheet_path = join(self.package_root, 'tests', 'data',
                                 'pre_prep_w_replicates.csv')

        job = GenPrepFileJob(self.run_dir,
                             self.convert_job_path,
                             self.qc_job_path,
                             join(self.working_directory_root,
                                  'OutputPath'),
                             sample_sheet_path,
                             'seqpro',
                             [],
                             'abcdabcdabcdabcdabcdabcdabcdabcd',
                             self.reports_path,
                             is_amplicon=True)

        exp = [['seqpro', '--verbose',
                ('sequence_processing_pipeline/b197f317-1c06-4619-9af3-'
                 '65721149c1e8/OutputPath/GenPrepFileJob/210518_A00953_0305'
                 '_TEST'),
                ('"sequence_processing_pipeline/b197f317-1c06-4619-9af3-'
                 '65721149c1e8/OutputPath/GenPrepFileJob/'
                 'replicate_sheet_1.txt"'),
                ('sequence_processing_pipeline/b197f317-1c06-4619-9af3-'
                 '65721149c1e8/OutputPath/GenPrepFileJob/PrepFiles/1')],
               ['seqpro', '--verbose',
                ('sequence_processing_pipeline/b197f317-1c06-4619-9af3-'
                 '65721149c1e8/OutputPath/GenPrepFileJob/210518_A00953_0305'
                 '_TEST'),
                ('"sequence_processing_pipeline/b197f317-1c06-4619-9af3-'
                 '65721149c1e8/OutputPath/GenPrepFileJob/'
                 'replicate_sheet_2.txt"'),
                ('sequence_processing_pipeline/b197f317-1c06-4619-9af3-'
                 '65721149c1e8/OutputPath/GenPrepFileJob/PrepFiles/2')],
               ['seqpro', '--verbose',
                ('sequence_processing_pipeline/b197f317-1c06-4619-9af3-'
                 '65721149c1e8/OutputPath/GenPrepFileJob/210518_A00953_0305'
                 '_TEST'),
                ('"sequence_processing_pipeline/b197f317-1c06-4619-9af3-'
                 '65721149c1e8/OutputPath/GenPrepFileJob/'
                 'replicate_sheet_3.txt"'),
                ('sequence_processing_pipeline/b197f317-1c06-4619-9af3-'
                 '65721149c1e8/OutputPath/GenPrepFileJob/PrepFiles/3')]]

        self.assertEqual(job.commands, exp)


if __name__ == '__main__':
    unittest.main()
