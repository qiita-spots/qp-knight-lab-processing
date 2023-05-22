# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from unittest import TestCase
from qp_klp.Step import Step
from sequence_processing_pipeline.Pipeline import Pipeline
from os.path import join, abspath
from functools import partial
from os import makedirs, chmod
from shutil import rmtree
from os import environ, remove
from json import dumps


class BaseStepTests(TestCase):
    CONFIGURATION = {
        "configuration": {
            "pipeline": {
                "archive_path": ("sequence_processing_pipeline/tests/data/"
                                 "sequencing/knight_lab_completed_runs"),
                "search_paths": ["/tmp", "qp_klp/tests/data"],
                "amplicon_search_paths": ["/tmp", "qp_klp/tests/data"]
            },
            "bcl2fastq": {
                "nodes": 1,
                "nprocs": 16,
                "queue": "qiita",
                "wallclock_time_in_hours": 36,
                "modules_to_load": ["bcl2fastq_2.20.0.422"],
                "executable_path": "bcl2fastq",
                "per_process_memory_limit": "10gb"
            },
            "bcl-convert": {
                "nodes": 1,
                "nprocs": 16,
                "queue": "qiita",
                "wallclock_time_in_hours": 36,
                "modules_to_load": ["bclconvert_3.7.5"],
                "executable_path": "bcl-convert",
                "per_process_memory_limit": "10gb"
            },
            "qc": {
                "nodes": 1,
                "nprocs": 16,
                "queue": "qiita",
                "wallclock_time_in_hours": 1,
                "minimap_databases": ["/databases/minimap2/human-phix-db.mmi"],
                "kraken2_database": "/databases/minimap2/hp_kraken-db.mmi",
                "modules_to_load": ["fastp_0.20.1", "samtools_1.12",
                                    " minimap2_2.18"],
                "fastp_executable_path": "fastp",
                "minimap2_executable_path": "minimap2",
                "samtools_executable_path": "samtools",
                "job_total_memory_limit": "20gb",
                "job_pool_size": 30,
                "job_max_array_length": 1000
            },
            "seqpro": {
                "seqpro_path": "seqpro",
                "modules_to_load": []
            },
            "fastqc": {
                "nodes": 1,
                "nprocs": 16,
                "queue": "qiita",
                "nthreads": 16,
                "wallclock_time_in_hours": 1,
                "modules_to_load": ["fastqc_0.11.5"],
                "fastqc_executable_path": "fastqc",
                "multiqc_executable_path": "multiqc",
                "multiqc_config_file_path": ("sequence_processing_pipeline/"
                                             "multiqc-bclconvert-config.yaml"),
                "job_total_memory_limit": "20gb",
                "job_pool_size": 30,
                "job_max_array_length": 1000
            }
        }
    }

    def setUp(self):
        package_root = abspath('./qp_klp')
        cc_path = partial(join, package_root, 'tests', 'data')
        self.good_config_file = join(package_root, 'configuration.json')
        self.good_run_id = '211021_A00000_0000_SAMPLE'
        self.good_sample_sheet_path = cc_path('good-sample-sheet.csv')
        self.good_mapping_file_path = cc_path('good-mapping-file.txt')
        self.output_file_path = cc_path('output_dir')
        self.qiita_id = '077c4da8-74eb-4184-8860-0207f53623be'
        makedirs(self.output_file_path, exist_ok=True)

        self.pipeline = Pipeline(None, self.good_run_id,
                                 self.good_sample_sheet_path, None,
                                 self.output_file_path, self.qiita_id,
                                 'metagenomic', BaseStepTests.CONFIGURATION)

        self.config = BaseStepTests.CONFIGURATION['configuration']

        self.fake_bin_path = self._get_searchable_path()

        self.delete_these = []

    def _is_writable(self, a_path):
        try:
            tmp = join(a_path, 'qpklp_temp_file')
            with open(tmp, 'w') as f:
                f.write('this is a test\n')
            remove(tmp)
            return True
        except IOError:
            return False

    def _get_searchable_path(self):
        searchable_paths = []

        if 'CONDA_PREFIX' in environ:
            # create fake binaries in bin directory of Conda environment
            searchable_paths.append(environ['CONDA_PREFIX'] + '/bin')
        else:
            # if CONDA_PREFIX doesn't exist, select a path from a list of
            # searchable paths that contains 'env' and assume it's writable.
            tmp = environ['PATH']
            searchable_paths += tmp.split(':')

        for a_path in searchable_paths:
            if self._is_writable(a_path):
                return a_path

    def _create_fake_bin(self, name, content):
        tmp = join(self.fake_bin_path, name)
        with open(tmp, 'w') as f:
            f.write(f"#!/bin/sh\n{content}\n")
        chmod(tmp, 0o777)
        self.delete_these.append(tmp)
        return tmp

    def _create_test_input(self, stage):
        if stage >= 1:
            fake_path = join(self.output_file_path, 'ConvertJob', 'logs')
            makedirs(fake_path, exist_ok=True)

            self._create_fake_bin('sbatch', "echo 'Submitted "
                                            "batch job 9999999'")

            self._create_fake_bin('sacct', "echo '9999999|99999999-9999-9999"
                                           "-9999-999999999999.txt|COMPLETED"
                                           "|09:53:41|0:0'")

        if stage >= 2:
            fake_path = join(self.output_file_path, 'QCJob', 'logs')
            makedirs(fake_path, exist_ok=True)

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
        for fake_bin in self.delete_these:
            remove(fake_bin)

    def test_creation(self):
        # Test base-class creation method, even though base-class will never
        # be instantiated by itself in normal usage.
        self._delete_test_output()

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

        step = Step(self.pipeline, self.qiita_id, sn_tid_map_by_project, None)

        self.assertIsNotNone(step)

    def test_convert_bcl_to_fastq(self):
        self._delete_test_output()
        self._create_test_input(1)

        sn_tid_map_by_project = {}
        step = Step(self.pipeline, self.qiita_id, sn_tid_map_by_project, None)

        fake_path = join(self.output_file_path, 'ConvertJob', 'logs')
        makedirs(fake_path, exist_ok=True)

        step._convert_bcl_to_fastq(self.config['bcl-convert'],
                                   self.good_sample_sheet_path)

    def test_quality_control(self):
        self._delete_test_output()
        self._create_test_input(2)

        fake_path = join(self.output_file_path, 'QCJob', 'logs')
        makedirs(fake_path, exist_ok=True)

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

        sn_tid_map_by_project = {}
        step = Step(self.pipeline, self.qiita_id, sn_tid_map_by_project, None)
        step._quality_control(self.config['qc'], self.good_sample_sheet_path)

    def test_generate_pipeline(self):
        tmp = join('.', 'tmp.config')
        with open(tmp, 'w') as f:
            f.write(dumps(BaseStepTests.CONFIGURATION, indent=2))

        pipeline = Step.generate_pipeline('metagenomic',
                                          self.good_sample_sheet_path,
                                          1,
                                          tmp,
                                          self.good_run_id,
                                          self.output_file_path,
                                          self.qiita_id)

        self.assertIsNotNone(pipeline)

        pipeline = Step.generate_pipeline('amplicon',
                                          self.good_mapping_file_path,
                                          1,
                                          tmp,
                                          self.good_run_id,
                                          self.output_file_path,
                                          self.qiita_id)

        self.assertIsNotNone(pipeline)

        remove(tmp)

    def test_get_project_info(self):
        obs = self.pipeline.get_project_info()

        exp = [{'project_name': 'NYU_BMS_Melanoma_13059', 'qiita_id': '13059'},
               {'project_name': 'Feist_11661', 'qiita_id': '11661'},
               {'project_name': 'Gerwick_6123', 'qiita_id': '6123'}]

        self.assertEqual(obs, exp)

    def test_parse_prep_file(self):
        good_prep_file = join('qp_klp', 'tests', 'good-prep-file-small.txt')

        obs = Step.parse_prep_file(good_prep_file)

        # assert that prep-files that begin with sample-names of the form
        # '363192526', '1e-3', and '123.000' are parsed as strings instead of
        # numeric values.
        exp = {'363192526': {'experiment_design_description': 'sample project',
                             'library_construction_protocol': ('Knight Lab Kap'
                                                               'a HyperPlus'),
                             'platform': 'Illumina', 'run_center': 'KLM',
                             'run_date': '2022-04-18',
                             'run_prefix': '363192526_S9_L001',
                             'sequencing_meth': 'sequencing by synthesis',
                             'center_name': 'UCSD',
                             'center_project_name': 'Sample_Project',
                             'instrument_model': 'Illumina iSeq',
                             'runid': '20220101_FS10001776_07_ABC12345-4567',
                             'lane': '1', 'sample project': 'Sample_Project',
                             'well_description': ('Sample_Project_99999_1-'
                                                  '4.363192526.A3'),
                             'i5_index_id': 'iTru5_09_A',
                             'sample_plate': 'Sample_Project_99999_1-4',
                             'index2': 'TCTGAGAG', 'index': 'CATCTACG',
                             'sample_well': 'A3',
                             'i7_index_id': 'iTru7_114_05',
                             'raw_reads': '10749',
                             'quality_filtered_reads': '1',
                             'non_host_reads': '4'},
               '363192073': {'experiment_design_description': 'sample project',
                             'library_construction_protocol': ('Knight Lab Ka'
                                                               'pa HyperPlus'),
                             'platform': 'Illumina', 'run_center': 'KLM',
                             'run_date': '2022-04-18',
                             'run_prefix': '363192073_S195_L001',
                             'sequencing_meth': 'sequencing by synthesis',
                             'center_name': 'UCSD',
                             'center_project_name': 'Sample_Project',
                             'instrument_model': 'Illumina iSeq',
                             'runid': '20220101_FS10001776_07_ABC12345-4567',
                             'lane': '1', 'sample project': 'Sample_Project',
                             'well_description': ('Sample_Project_99999_1-'
                                                  '4.363192073.F1'),
                             'i5_index_id': 'iTru5_103_A',
                             'sample_plate': 'Sample_Project_99999_1-4',
                             'index2': 'TGGTCCTT', 'index': 'GCAATTCG',
                             'sample_well': 'F1',
                             'i7_index_id': 'iTru7_305_11',
                             'raw_reads': '16435',
                             'quality_filtered_reads': '2',
                             'non_host_reads': '5'},
               '363193755': {'experiment_design_description': 'sample project',
                             'library_construction_protocol': ('Knight Lab Ka'
                                                               'pa HyperPlus'),
                             'platform': 'Illumina', 'run_center': 'KLM',
                             'run_date': '2022-04-18',
                             'run_prefix': '363193755_S7_L001',
                             'sequencing_meth': 'sequencing by synthesis',
                             'center_name': 'UCSD',
                             'center_project_name': 'Sample_Project',
                             'instrument_model': 'Illumina iSeq',
                             'runid': '20220101_FS10001776_07_ABC12345-4567',
                             'lane': '1', 'sample project': 'Sample_Project',
                             'well_description': ('Sample_Project_99999_1-'
                                                  '4.363193755.M1'),
                             'i5_index_id': 'iTru5_07_A',
                             'sample_plate': 'Sample_Project_99999_1-4',
                             'index2': 'GGTGTCTT', 'index': 'GATTGCTC',
                             'sample_well': 'M1',
                             'i7_index_id': 'iTru7_114_03',
                             'raw_reads': '14303',
                             'quality_filtered_reads': '3',
                             'non_host_reads': '6'},
               '1e-3': {'experiment_design_description': 'sample project',
                        'library_construction_protocol': ('Knight Lab Kapa '
                                                          'HyperPlus'),
                        'platform': 'Illumina', 'run_center': 'KLM',
                        'run_date': '2022-04-18',
                        'run_prefix': '363192073_S195_L001',
                        'sequencing_meth': 'sequencing by synthesis',
                        'center_name': 'UCSD',
                        'center_project_name': 'Sample_Project',
                        'instrument_model': 'Illumina iSeq',
                        'runid': '20220101_FS10001776_07_ABC12345-4567',
                        'lane': '1', 'sample project': 'Sample_Project',
                        'well_description': ('Sample_Project_99999_1-'
                                             '4.363192073.F1'),
                        'i5_index_id': 'iTru5_103_A',
                        'sample_plate': 'Sample_Project_99999_1-4',
                        'index2': 'TGGTCCTT', 'index': 'GCAATTCG',
                        'sample_well': 'F1', 'i7_index_id': 'iTru7_305_11',
                        'raw_reads': '16435', 'quality_filtered_reads': '11',
                        'non_host_reads': '13'},
               '123.000': {'experiment_design_description': 'sample project',
                           'library_construction_protocol': ('Knight Lab Kapa'
                                                             ' HyperPlus'),
                           'platform': 'Illumina', 'run_center': 'KLM',
                           'run_date': '2022-04-18',
                           'run_prefix': '363193755_S7_L001',
                           'sequencing_meth': 'sequencing by synthesis',
                           'center_name': 'UCSD',
                           'center_project_name': 'Sample_Project',
                           'instrument_model': 'Illumina iSeq',
                           'runid': '20220101_FS10001776_07_ABC12345-4567',
                           'lane': '1', 'sample project': 'Sample_Project',
                           'well_description': ('Sample_Project_99999_1-'
                                                '4.363193755.M1'),
                           'i5_index_id': 'iTru5_07_A',
                           'sample_plate': 'Sample_Project_99999_1-4',
                           'index2': 'GGTGTCTT', 'index': 'GATTGCTC',
                           'sample_well': 'M1', 'i7_index_id': 'iTru7_114_03',
                           'raw_reads': '14303',
                           'quality_filtered_reads': '12',
                           'non_host_reads': '14'}}

        self.assertDictEqual(obs, exp)
