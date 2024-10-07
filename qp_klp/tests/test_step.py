# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from unittest import TestCase
from sequence_processing_pipeline.Pipeline import Pipeline, PipelineError
from os.path import join, abspath, exists, dirname
from functools import partial
from os import makedirs, chmod, access, W_OK
from shutil import rmtree, copy, which, copytree
from os import environ, remove, getcwd
import pandas as pd
from metapool import parse_prep
import re
from copy import deepcopy
from tempfile import TemporaryDirectory
from FailedSamplesRecord import FailedSamplesRecord


class FakeClient():
    def __init__(self):
        self.cwd = getcwd()
        self.base_path = join(self.cwd, 'qp_klp/tests/data/QDir')
        self.qdirs = {'Demultiplexed': 'Demultiplexed',
                      'beta_div_plots': 'analysis/beta_div_plots',
                      'rarefaction_curves': 'analysis/rarefaction_curves',
                      'taxa_summary': 'analysis/taxa_summary',
                      'q2_visualization': 'working_dir',
                      'distance_matrix': 'working_dir',
                      'ordination_results': 'working_dir',
                      'alpha_vector': 'working_dir',
                      'FASTQ': 'FASTQ',
                      'BIOM': 'BIOM',
                      'per_sample_FASTQ': 'per_sample_FASTQ',
                      'SFF': 'SFF',
                      'FASTA': 'FASTA',
                      'FASTA_Sanger': 'FASTA_Sanger',
                      'FeatureData': 'FeatureData',
                      'job-output-folder': 'job-output-folder',
                      'BAM': 'BAM',
                      'VCF': 'VCF',
                      'SampleData': 'SampleData',
                      'uploads': 'uploads'}

        self.samples_in_13059 = ['13059.SP331130A04', '13059.AP481403B02',
                                 '13059.LP127829A02', '13059.BLANK3.3B',
                                 '13059.EP529635B02', '13059.EP542578B04',
                                 '13059.EP446602B01', '13059.EP121011B01',
                                 '13059.EP636802A01', '13059.SP573843A04']

        # note these samples have known tids, but aren't in good-sample-sheet.
        self.samples_in_11661 = ['11661.1.24', '11661.1.57', '11661.1.86',
                                 '11661.10.17', '11661.10.41', '11661.10.64',
                                 '11661.11.18', '11661.11.43', '11661.11.64',
                                 '11661.12.15']

        self.samples_in_6123 = ['3A', '4A', '5B', '6A', 'BLANK.41.12G', '7A',
                                '8A', 'ISB', 'GFR', '6123']

        self.info_in_11661 = {'number-of-samples': 10,
                              'categories': ['sample_type', 'tube_id']}

        self.info_in_13059 = {'number-of-samples': 10,
                              'categories': ['anonymized_name',
                                             'collection_timestamp',
                                             'description',
                                             'dna_extracted',
                                             'elevation', 'empo_1',
                                             'empo_2', 'empo_3',
                                             'env_biome', 'env_feature',
                                             'env_material',
                                             'env_package',
                                             'geo_loc_name', 'host_age',
                                             'host_age_units',
                                             'host_body_habitat',
                                             'host_body_mass_index',
                                             'host_body_product',
                                             'host_body_site',
                                             'host_common_name',
                                             'host_height',
                                             'host_height_units',
                                             'host_life_stage',
                                             'host_scientific_name',
                                             'host_subject_id',
                                             'host_taxid', 'host_weight',
                                             'host_weight_units',
                                             'latitude', 'longitude',
                                             'nyuid',
                                             'physical_specimen_location',
                                             'physical_specimen_remaining',
                                             'predose_time',
                                             'sample_type',
                                             'scientific_name', 'sex',
                                             'subject_id', 'taxon_id',
                                             'title', 'tube_id']}

        # Study not in qiita-rc. Faking results.
        self.info_in_6123 = {'number-of-samples': 10,
                             'categories': ['sample_type', 'subject_id',
                                            'title']}

        self.tids_13059 = {"header": ["tube_id"],
                           "samples": {'13059.SP331130A04': ['SP331130A-4'],
                                       '13059.AP481403B02': ['AP481403B-2'],
                                       '13059.LP127829A02': ['LP127829A-2'],
                                       '13059.BLANK3.3B': ['BLANK3.3B'],
                                       '13059.EP529635B02': ['EP529635B-2'],
                                       '13059.EP542578B04': ['EP542578B-4'],
                                       '13059.EP446602B01': ['EP446602B-1'],
                                       '13059.EP121011B01': ['EP121011B-1'],
                                       '13059.EP636802A01': ['EP636802A-1'],
                                       '13059.SP573843A04': ['SP573843A-4']}}

        self.tids_11661 = {"header": ["tube_id"],
                           "samples": {"11661.1.24": ["1.24"],
                                       "11661.1.57": ["1.57"],
                                       "11661.1.86": ["1.86"],
                                       "11661.10.17": ["10.17"],
                                       "11661.10.41": ["10.41"],
                                       "11661.10.64": ["10.64"],
                                       "11661.11.18": ["11.18"],
                                       "11661.11.43": ["11.43"],
                                       "11661.11.64": ["11.64"],
                                       "11661.12.15": ["12.15"]}}

        for key in self.qdirs:
            self.qdirs[key] = join(self.base_path, self.qdirs[key])

        for qdir in self.qdirs:
            makedirs(self.qdirs[qdir], exist_ok=True)

        self.fake_id = 1000
        self._server_url = "some.server.url"
        self.saved_posts = {}

    def get(self, url):
        m = {'/api/v1/study/11661/samples': self.samples_in_11661,
             '/api/v1/study/11661/samples/categories=tube_id': self.tids_11661,
             '/api/v1/study/11661/samples/info': self.info_in_11661,
             '/api/v1/study/13059/samples': self.samples_in_13059,
             '/api/v1/study/13059/samples/categories=tube_id': self.tids_13059,
             '/api/v1/study/13059/samples/info': self.info_in_13059,
             '/api/v1/study/6123/samples': self.samples_in_6123,
             '/api/v1/study/6123/samples/info': self.info_in_6123,
             '/qiita_db/artifacts/types/': self.qdirs}

        if url in m:
            return m[url]

        return None

    def post(self, url, data=None):
        if '/qiita_db/prep_template/' == url:
            self.fake_id += 1
            return {'prep': self.fake_id}
        elif '/qiita_db/artifact/' == url:
            self.saved_posts[str(self.fake_id)] = data
            self.fake_id += 1
            return {'job_id': self.fake_id}
        else:
            raise ValueError("Unsupported URL")


class AnotherFakeClient():
    def __init__(self):
        self.cwd = getcwd()
        self.base_path = join(self.cwd, 'qp_klp/tests/data/QDir')
        self.qdirs = {'Demultiplexed': 'Demultiplexed',
                      'beta_div_plots': 'analysis/beta_div_plots',
                      'rarefaction_curves': 'analysis/rarefaction_curves',
                      'taxa_summary': 'analysis/taxa_summary',
                      'q2_visualization': 'working_dir',
                      'distance_matrix': 'working_dir',
                      'ordination_results': 'working_dir',
                      'alpha_vector': 'working_dir',
                      'FASTQ': 'FASTQ',
                      'BIOM': 'BIOM',
                      'per_sample_FASTQ': 'per_sample_FASTQ',
                      'SFF': 'SFF',
                      'FASTA': 'FASTA',
                      'FASTA_Sanger': 'FASTA_Sanger',
                      'FeatureData': 'FeatureData',
                      'job-output-folder': 'job-output-folder',
                      'BAM': 'BAM',
                      'VCF': 'VCF',
                      'SampleData': 'SampleData',
                      'uploads': 'uploads'}

        self.samples_in_99999 = ['99999.AAAAAAAAAAA',
                                 '99999.BBBBBBBBBBB',
                                 '99999.CCCCCCCCCCC',
                                 '99999.BLANK1.1BCD']

        self.info_in_99999 = {'number-of-samples': 10,
                              'categories': ['column1', 'column2', 'tube_id']}

        self.tids_99999 = {"header": ["tube_id"],
                           "samples": {'99999.AAAAAAAAAAA': ['1234567890a'],
                                       '99999.BBBBBBBBBBB': ['234567890ab'],
                                       '99999.CCCCCCCCCCC': ['34567890abc'],
                                       '99999.BLANK1.1BCD': ['BLANK1.1BCD']}}

        for key in self.qdirs:
            self.qdirs[key] = join(self.base_path, self.qdirs[key])

        for qdir in self.qdirs:
            makedirs(self.qdirs[qdir], exist_ok=True)

    def get(self, url):
        m = {'/api/v1/study/99999/samples': self.samples_in_99999,
             '/api/v1/study/99999/samples/categories=tube_id': self.tids_99999,
             '/api/v1/study/99999/samples/info': self.info_in_99999,
             '/qiita_db/artifacts/types/': self.qdirs}

        if url in m:
            return m[url]

        return None


class BaseStepTests(TestCase):
    '''
    BaseStepTests contains all the configuration information and helper
    functions used by every child StepTests class. This class does not
    include any tests. All tests defined in this class will be inherited by
    every child and will consequently be run multiple times. Hence, general
    functionality is instead tested by BasicStepSteps class.
    '''

    def setUp(self):
        package_root = abspath('./qp_klp')
        cc_path = partial(join, package_root, 'tests', 'data')
        self.good_config_file = join(package_root, 'configuration.json')
        self.good_run_id = '211021_A00000_0000_SAMPLE'
        self.good_sample_sheet_path = cc_path('good-sample-sheet.csv')
        self.another_good_sample_sheet_path = cc_path('another-good-sample-'
                                                      'sheet.csv')
        self.sheet_w_replicates_path = cc_path('good_sheet_w_replicates.csv')
        self.good_mapping_file_path = cc_path('good-mapping-file.txt')
        self.good_prep_info_file_path = cc_path('good-sample-prep.tsv')
        self.good_transcript_sheet_path = cc_path('good-sample-sheet-'
                                                  'transcriptomics.csv')
        self.output_file_path = cc_path('output_dir')
        self.process_shell_script = cc_path('process_all_fastq_files.sh')
        self.master_config_path = cc_path('configuration.json')
        self.dummy_fastq_file = cc_path('dummy.fastq.gz')
        self.mini_sheet_path = cc_path('mini-sample-sheet.csv')
        self.qiita_id = '077c4da8-74eb-4184-8860-0207f53623be'
        makedirs(self.output_file_path, exist_ok=True)

        self.pipeline = Pipeline(self.master_config_path,
                                 self.good_run_id,
                                 self.good_sample_sheet_path, None,
                                 self.output_file_path, self.qiita_id,
                                 Step.METAGENOMIC_TYPE)

        self.pipeline_mini = Pipeline(self.master_config_path,
                                      self.good_run_id,
                                      self.mini_sheet_path, None,
                                      self.output_file_path, self.qiita_id,
                                      Step.METAGENOMIC_TYPE)

        self.another_pipeline = Pipeline(self.master_config_path,
                                         self.good_run_id,
                                         self.another_good_sample_sheet_path,
                                         None, self.output_file_path,
                                         self.qiita_id, Step.METAGENOMIC_TYPE)

        self.pipeline_replicates = Pipeline(self.master_config_path,
                                            self.good_run_id,
                                            self.sheet_w_replicates_path, None,
                                            self.output_file_path,
                                            self.qiita_id,
                                            Step.METAGENOMIC_TYPE)

        self.amplicon_pipeline = Pipeline(self.master_config_path,
                                          self.good_run_id, None,
                                          self.good_mapping_file_path,
                                          self.output_file_path,
                                          self.qiita_id,
                                          Step.AMPLICON_TYPE)

        self.fake_bin_path = self._get_searchable_path()

        self.delete_these = []

    def tearDown(self):
        if exists(self.output_file_path):
            rmtree(self.output_file_path)
        for fake_bin in self.delete_these:
            if exists(fake_bin):
                remove(fake_bin)
        if exists('tmp.config'):
            remove('tmp.config')

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
            if access(a_path, W_OK):
                return a_path

    def _create_fake_bin(self, name, content):
        tmp = join(self.fake_bin_path, name)
        with open(tmp, 'w') as f:
            f.write(f"#!/bin/sh\n{content}\n")
        chmod(tmp, 0o777)
        self.delete_these.append(tmp)
        return tmp

    def _create_fake_file(self, path):
        with open(path, 'w') as f:
            f.write("This is a file.")

    def _create_test_input(self, stage):
        if stage >= 1:
            # create an empty ConvertJob directory to test initialization
            # with. Create fake binaries to test job submission.
            fake_path = join(self.output_file_path, 'ConvertJob', 'logs')
            makedirs(fake_path, exist_ok=True)
            fake_path = join(self.output_file_path, 'ConvertJob', 'Reports')
            makedirs(fake_path, exist_ok=True)

            self._create_fake_bin('sbatch', "echo 'Submitted "
                                            "batch job 9999999'")

            self._create_fake_bin('squeue', "echo 'ARRAY_JOB_ID,JOBID,STATE\n"
                                  "9999999,9999999,COMPLETED'")

        if stage >= 2:
            # generate dummy fastq files in ConvertJob and create an empty
            # NuQCJob directory to use for testing NuQCJob initialization.
            fake_path = join(self.output_file_path, 'NuQCJob', 'logs')
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
                        self._create_fake_file(file_path)

        if stage >= 3:
            # create a fake GenPrepFileJob directory.
            fake_path = join(self.output_file_path, 'GenPrepFileJob',
                             'PrepFiles')
            makedirs(fake_path, exist_ok=True)
            names = ['NYU_BMS_Melanoma_13059.1.tsv', 'Feist_11661.1.tsv',
                     'Gerwick_6123.1.tsv']

            for name in names:
                self._create_fake_file(join(fake_path, name))

            fake_paths = [join(self.output_file_path, 'NuQCJob',
                               'NYU_BMS_Melanoma_13059', 'fastp_reports_dir'),
                          join(self.output_file_path, 'NuQCJob',
                               'Feist_11661', 'fastp_reports_dir'),
                          join(self.output_file_path, 'NuQCJob',
                               'Gerwick_6123', 'fastp_reports_dir')
                          ]

            for fake_path in fake_paths:
                makedirs(fake_path, exist_ok=True)
                self._create_fake_file(join(fake_path, 'a_file'))

            names = ['NYU_BMS_Melanoma_13059', 'Feist_11661',
                     'Gerwick_6123']

            for project in names:
                file_name = f'{self.good_run_id}_{project}_blanks.tsv'
                fake_path = join(self.output_file_path, file_name)
                self._create_fake_file(fake_path)

            tarballs = ['logs-ConvertJob.tgz', 'logs-FastQCJob.tgz',
                        'logs-GenPrepFileJob.tgz', 'logs-QCJob.tgz',
                        'prep-files.tgz', 'reports-ConvertJob.tgz',
                        'reports-FastQCJob.tgz', 'reports-QCJob.tgz',
                        'sample-files.tgz']

            for file_name in tarballs:
                fake_path = join(self.output_file_path, file_name)
                self._create_fake_file(fake_path)

            suffixes = ['o1611416-26', 'e1611416-26']
            for file_name in suffixes:
                file_name = f'{self.good_run_id}_FastQCJob.{file_name}'
                fake_path = join(self.output_file_path, 'FastQCJob', 'logs')
                makedirs(fake_path, exist_ok=True)
                self._create_fake_file(join(fake_path, file_name))

            # we're just going to create a directory for FastQC results and
            # create a single file. We aren't going to replicate the entire
            # directory structure for now.
            fake_path = join(self.output_file_path, 'FastQCJob', 'fastqc')
            makedirs(fake_path, exist_ok=True)
            self._create_fake_file(join(fake_path, 'a_file.txt'))

            fake_path = join(self.output_file_path, 'GenPrepFileJob', 'logs')
            makedirs(fake_path, exist_ok=True)
            self._create_fake_file(join(fake_path, 'a_file.txt'))

            fake_path = join(self.output_file_path, 'failed_samples.html')
            self._create_fake_file(fake_path)

    def _create_alternate_test_input(self):
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
            trimmed_files_path = join(self.output_file_path, 'NuQCJob',
                                      project, 'filtered_sequences')
            empty_files_path = join(self.output_file_path, 'NuQCJob',
                                    project, 'zero_files')
            adapter_trimmed_files_path = join(self.output_file_path,
                                              'NuQCJob',
                                              'only-adapter-filtered',
                                              project)

            fake_paths = [trimmed_files_path, empty_files_path,
                          adapter_trimmed_files_path]

            for fake_path in fake_paths:
                makedirs(fake_path, exist_ok=True)

            empty_files = {
                    'Feist_11661': [
                        'CDPH-SAL_Salmonella_Typhi_MDL-150',
                        'CDPH-SAL_Salmonella_Typhi_MDL-151'
                        ],
                    'Gerwick_6123': ['8A', '9A', '10A'],
                    'NYU_BMS_Melanoma_13059': ['XX581451B02', 'XY256645B01',
                                               'XZ112567B02'
                                               ]}

            f_list = []
            for sample in exp[project]:
                f_list += [
                    join(trimmed_files_path,
                         f'{sample}_SXXX_L001_R1_001.trimmed.fastq.gz'),
                    join(trimmed_files_path,
                         f'{sample}_SXXX_L001_R2_001.trimmed.fastq.gz'),
                    join(trimmed_files_path,
                         f'{sample}_SXXX_L001_I1_001.trimmed.fastq.gz'),
                    join(trimmed_files_path,
                         f'{sample}_SXXX_L001_I2_001.trimmed.fastq.gz'),
                    join(adapter_trimmed_files_path,
                         f'{sample}_SXXX_L001_R1_001.fastq.gz'),
                    join(adapter_trimmed_files_path,
                         f'{sample}_SXXX_L001_R2_001.fastq.gz')
                ]

            for sample in empty_files[project]:
                f_list += [
                    join(trimmed_files_path,
                         f'{sample}_SXXX_L001_R1_001.trimmed.fastq.gz'),
                    join(trimmed_files_path,
                         f'{sample}_SXXX_L001_R2_001.trimmed.fastq.gz')
                ]

            for file_path in f_list:
                self._create_fake_file(file_path)


class BasicStepTests(BaseStepTests):
    def atest_creation(self):
        # Test base-class creation method, even though base-class will never
        # be instantiated by itself in normal usage.

        with self.assertRaisesRegex(ValueError, "A pipeline object is needed"
                                                " to initialize Step"):
            Step(None, self.qiita_id, None)

        with self.assertRaisesRegex(ValueError, "A Qiita job-id is needed to "
                                                "initialize Step"):
            Step(self.pipeline, None, None)

        step = Step(self.pipeline, self.qiita_id, None)

        self.assertIsNotNone(step)

    def atest_convert_bcl_to_fastq(self):
        self._create_test_input(1)

        step = Step(self.pipeline, self.qiita_id, None)

        fake_path = join(self.output_file_path, 'ConvertJob', 'logs')
        makedirs(fake_path, exist_ok=True)

        config = self.pipeline.config_profile['profile']['configuration']
        step._convert_bcl_to_fastq(config['bcl-convert'],
                                   self.good_sample_sheet_path)

    def atest_quality_control(self):
        self._create_test_input(2)

        fake_path = join(self.output_file_path, 'NuQCJob', 'logs')
        makedirs(fake_path, exist_ok=True)

        exp = {'Feist_11661': (['CDPH-SAL_Salmonella_Typhi_MDL-143',
                                'CDPH-SAL_Salmonella_Typhi_MDL-144',
                                'CDPH-SAL_Salmonella_Typhi_MDL-145',
                                'CDPH-SAL_Salmonella_Typhi_MDL-146',
                                'CDPH-SAL_Salmonella_Typhi_MDL-147'], False),
               'Gerwick_6123': (['3A', '4A', '5B', '6A', '7A'], True),
               'NYU_BMS_Melanoma_13059': (['AP581451B02', 'EP256645B01',
                                           'EP112567B02', 'EP337425B01',
                                           'LP127890A01'], False)}
        for project in exp:
            fake_path = join(self.output_file_path, 'ConvertJob', project)
            fake_path2 = join(self.output_file_path, 'NuQCJob', project)
            makedirs(fake_path, exist_ok=True)
            if exp[project][1]:
                makedirs(join(fake_path2, 'filtered_sequences'), exist_ok=True)
            else:
                makedirs(join(fake_path2, 'trimmed_sequences'), exist_ok=True)

            for sample in exp[project][0]:
                r1 = join(fake_path, f'{sample}_SXXX_L001_R1_001.fastq.gz')
                r2 = join(fake_path, f'{sample}_SXXX_L001_R2_001.fastq.gz')

                for file_path in [r1, r2]:
                    self._create_fake_file(file_path)

        step = Step(self.pipeline, self.qiita_id, None)
        config = self.pipeline.config_profile['profile']['configuration']
        step._quality_control(config['nu-qc'],
                              self.good_sample_sheet_path)

    def atest_generate_pipeline(self):
        pipeline = Step.generate_pipeline(Step.METAGENOMIC_TYPE,
                                          self.good_sample_sheet_path,
                                          1,
                                          self.master_config_path,
                                          self.good_run_id,
                                          self.output_file_path,
                                          self.qiita_id)

        self.assertIsNotNone(pipeline)

        pipeline = Step.generate_pipeline(Step.AMPLICON_TYPE,
                                          self.good_mapping_file_path,
                                          1,
                                          self.master_config_path,
                                          self.good_run_id,
                                          self.output_file_path,
                                          self.qiita_id)

        self.assertIsNotNone(pipeline)

        pipeline = Step.generate_pipeline(Step.METATRANSCRIPTOMIC_TYPE,
                                          self.good_transcript_sheet_path,
                                          1,
                                          self.master_config_path,
                                          self.good_run_id,
                                          self.output_file_path,
                                          self.qiita_id)

        self.assertIsNotNone(pipeline)

    def atest_get_project_info(self):
        obs = self.pipeline.get_project_info()

        exp = [{'project_name': 'NYU_BMS_Melanoma_13059', 'qiita_id': '13059',
                'contains_replicates': False},
               {'project_name': 'Feist_11661', 'qiita_id': '11661',
                'contains_replicates': False},
               {'project_name': 'Gerwick_6123', 'qiita_id': '6123',
                'contains_replicates': False}]

        self.assertEqual(obs, exp)

    def atest_parse_prep_file(self):
        good_prep_file = join('qp_klp', 'tests', 'good-prep-file-small.txt')

        obs = Step._parse_prep_file(good_prep_file)

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
                             'raw_reads_r1r2': '10749',
                             'quality_filtered_reads_r1r2': '1',
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
                             'raw_reads_r1r2': '16435',
                             'quality_filtered_reads_r1r2': '2',
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
                             'raw_reads_r1r2': '14303',
                             'quality_filtered_reads_r1r2': '3',
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
                        'raw_reads_r1r2': '16435',
                        'quality_filtered_reads_r1r2': '11',
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
                           'raw_reads_r1r2': '14303',
                           'quality_filtered_reads_r1r2': '12',
                           'non_host_reads': '14'}}

        self.assertDictEqual(obs, exp)

        # simply confirm that a DataFrame is returned when convert_to_dict is
        # False. We already know that the contents of obs will be correct.
        obs = Step._parse_prep_file(good_prep_file, convert_to_dict=False)
        self.assertIsInstance(obs, pd.DataFrame)

    def atest_generate_special_map(self):
        fake_client = FakeClient()
        step = Step(self.pipeline, self.qiita_id, None)
        step.generate_special_map(fake_client)
        obs = step.special_map

        exp = [('NYU_BMS_Melanoma_13059',
                join(fake_client.base_path, 'uploads/13059'), '13059'),
               ('Feist_11661',
                join(fake_client.base_path, 'uploads/11661'), '11661'),
               ('Gerwick_6123',
                join(fake_client.base_path, 'uploads/6123'), '6123')]

        self.assertEquals(obs, exp)

    def atest_get_samples_in_qiita(self):
        fake_client = FakeClient()
        step = Step(self.pipeline, self.qiita_id, None)
        obs_samples, obs_tids = step.get_samples_in_qiita(fake_client, '13059')

        exp_samples = {'EP121011B01', 'EP529635B02', 'EP542578B04',
                       'SP573843A04', 'SP331130A04', 'EP446602B01',
                       'BLANK3.3B', 'AP481403B02', 'LP127829A02',
                       'EP636802A01'}

        exp_tids = {'13059.SP331130A04': ['SP331130A-4'],
                    '13059.AP481403B02': ['AP481403B-2'],
                    '13059.LP127829A02': ['LP127829A-2'],
                    '13059.BLANK3.3B': ['BLANK3.3B'],
                    '13059.EP529635B02': ['EP529635B-2'],
                    '13059.EP542578B04': ['EP542578B-4'],
                    '13059.EP446602B01': ['EP446602B-1'],
                    '13059.EP121011B01': ['EP121011B-1'],
                    '13059.EP636802A01': ['EP636802A-1'],
                    '13059.SP573843A04': ['SP573843A-4']}

        self.assertEqual(obs_samples, exp_samples)
        self.assertDictEqual(obs_tids, exp_tids)

    def atest_get_tube_ids_from_qiita(self):
        fake_client = FakeClient()
        step = Step(self.pipeline, self.qiita_id, None)
        step._get_tube_ids_from_qiita(fake_client)
        obs = step.tube_id_map

        exp = {'13059': {'SP331130A04': 'SP331130A-4',
                         'AP481403B02': 'AP481403B-2',
                         'LP127829A02': 'LP127829A-2',
                         'BLANK3.3B': 'BLANK3.3B',
                         'EP529635B02': 'EP529635B-2',
                         'EP542578B04': 'EP542578B-4',
                         'EP446602B01': 'EP446602B-1',
                         'EP121011B01': 'EP121011B-1',
                         'EP636802A01': 'EP636802A-1',
                         'SP573843A04': 'SP573843A-4'},
               '11661': {'1.24': '1.24', '1.57': '1.57', '1.86': '1.86',
                         '10.17': '10.17', '10.41': '10.41', '10.64': '10.64',
                         '11.18': '11.18', '11.43': '11.43', '11.64': '11.64',
                         '12.15': '12.15'}}

        self.assertDictEqual(obs, exp)

    def atest_compare_samples_against_qiita(self):
        fake_client = FakeClient()
        step = Step(self.pipeline_mini, self.qiita_id, None)
        results = step._compare_samples_against_qiita(fake_client)

        # confirm projects in results match what's expected
        obs = [project['project_name'] for project in results]
        exp = ["NYU_BMS_Melanoma", "Gerwick"]

        self.assertEqual(obs, exp)

        # confirm projects using tube-ids match what's expected

        # results are a list of project dicts, rather than a dict of dicts.
        # however they are faked and can be expected to be returned in a
        # fixed order. Assert the order is as expected so the following tests
        # will be meaningful.

        # good_sample_sheet.csv will have some but not all sample-names
        # exchanged for tube-ids.
        self.assertCountEqual([proj['project_name'] for proj in results],
                              ['NYU_BMS_Melanoma', 'Gerwick'])

        # since Gerwick doesn't have tube-ids, it should always use sample-
        # names. NYU has tube-id in FakeQiita() so it's possible to test
        # tube-ids.
        self.assertCountEqual([proj['used_tids'] for proj in results],
                              [True, False])

        # 'NOTINQIITA1' is a sample-name from the sample-sheet and should not
        # be in fake-Qiita, as defined in FakeQiita() class. Therefore, it
        # should appear in the 'samples_not_in_qiita' list.
        self.assertIn('NOTINQIITA1', results[0]['samples_not_in_qiita'])

        # 'BLANK3.3B' is defined in the sample-sheet and also in FakeQiita,
        # both as a sample-name and as a tube-id (One of the few to be so
        # named). It shouldn't appear in 'samples_not_in_qiita' list.
        self.assertNotIn('BLANK3.3B', results[0]['samples_not_in_qiita'])

        # 'SP331130A-4' is a tube-id in qiita and should be present in the
        # 'examples_in_qiita' list
        # the tube-ids in 'examples_in_qiita' list should be a subset of all
        # the tube-ids in FakeQiita().
        exp = {'SP331130A-4', 'AP481403B-2', 'LP127829A-2', 'BLANK3.3B',
               'EP529635B-2', 'EP542578B-4', 'EP446602B-1', 'EP121011B-1',
               'EP636802A-1', 'SP573843A-4'}

        self.assertTrue(set(results[0]['examples_in_qiita']).issubset(exp))

        # Gerwick has a small number of samples in the sample-sheet, and all
        # of which are in FakeQiita().
        self.assertEqual(results[1]['samples_not_in_qiita'], set())

    def atest_generate_commands(self):
        self._create_test_input(3)

        fake_client = FakeClient()

        # self.pipeline represents a metagenomic pathway.
        step = Step(self.pipeline, self.qiita_id, None)

        # need to generate some metadata in order to generate commands.
        step.generate_special_map(fake_client)

        # test base _generate_commands() method; contains only commands used
        # across all pipeline types.
        step.generate_commands()

        exp = [
            (f'cd {self.output_file_path}; '
             'tar zcvf logs-ConvertJob.tgz ConvertJob/logs'),
            (f'cd {self.output_file_path}; '
             'tar zcvf reports-ConvertJob.tgz ConvertJob/Reports '
             'ConvertJob/logs'),
            (f'cd {self.output_file_path}; '
             'tar zcvf logs-NuQCJob.tgz NuQCJob/logs'),
            (f'cd {self.output_file_path}; '
             'tar zcvf logs-FastQCJob.tgz FastQCJob/logs'),
            (f'cd {self.output_file_path}; '
             'tar zcvf reports-FastQCJob.tgz FastQCJob/fastqc'),
            (f'cd {self.output_file_path}; '
             'tar zcvf logs-GenPrepFileJob.tgz GenPrepFileJob/logs'),
            (f'cd {self.output_file_path}; '
             'tar zcvf prep-files.tgz GenPrepFileJob/PrepFiles'),
            (f'cd {self.output_file_path}; '
             'mv failed_samples.html final_results'),
            (f'cd {self.output_file_path}; '
             'tar zcvf reports-NuQCJob.tgz NuQCJob/Feist_11661/'
             'fastp_reports_dir '
             'NuQCJob/Gerwick_6123/fastp_reports_dir '
             'NuQCJob/NYU_BMS_Melanoma_13059/fastp_reports_dir'),
            (f'cd {self.output_file_path}; '
             'tar zcvf sample-files.tgz 211021_A00000_0000_SAMPLE_Feist_11661'
             '_blanks.tsv 211021_A00000_0000_SAMPLE_Gerwick_6123_blanks.tsv '
             '211021_A00000_0000_SAMPLE_NYU_BMS_Melanoma_13059_blanks.tsv'),
            (f'cd {self.output_file_path}; (find *.tgz -maxdepth 1 -type f '
             '| xargs mv -t final_results) || true')]

        # replace unique string w/the base-directory path in the expected
        # output.
        for i in range(0, len(exp)):
            exp[i] = exp[i].replace('BASE_DIRECTORY', getcwd())

        self.assertEqual(step.cmds, exp)

    def atest_overwrite_prep_files(self):
        # use a prep-file specifically for modification by the
        # _overwrite_prep_files() method.
        fake_client = FakeClient()
        step = Step(self.pipeline, self.qiita_id, None)

        # copy the file so that we do not overwrite the original, which is
        # useful for other tests.

        sample_path = join(dirname(self.good_prep_info_file_path),
                           ('20230101_XX99999999_99_LOL99999-9999.'
                            'NYU_BMS_Melanoma_13059.1.tsv'))

        copy(self.good_prep_info_file_path, sample_path)

        # needed to prep for _overwrite_prep_files()
        step._get_tube_ids_from_qiita(fake_client)
        step._overwrite_prep_files([sample_path])

        # read in the changed prep-file and confirm that the sample_name
        # column contains sample-names instead of tube-ids and that the
        # tube-ids have been moved to a new column named 'old_sample_name'.
        df = pd.read_csv(sample_path, sep='\t', dtype=str, index_col=False)

        new_sample_names = set(df['sample_name'])

        # use the list of sample-names for the project stored in FakeClient()
        # as the expected set of metadata.
        exp = set([sample_name.replace('13059.', '') for sample_name in
                   fake_client.samples_in_13059])
        self.assertEqual(new_sample_names, exp)

        # confirm tids are where they're expected to be as well
        new_old_sample_names = set(df['old_sample_name'])
        tids = fake_client.tids_13059['samples']
        exp = set([tids[t][0] for t in tids])
        self.assertEqual(new_old_sample_names, exp)

    def atest_compare_samples_against_qiita_error_handling(self):
        fake_client = AnotherFakeClient()

        # In addition to the samples in AnotherFakeClient(),
        # another-good-sample-sheet.csv has an unregistered BLANK:
        # 'BLANK_UNREG'. _compare_samples_against_qiita() should detect
        # that it is a 'new BLANK', and not include it in the value for
        # 'samples_not_in_qiita.'

        # another-good-sample-sheet.csv contains one sample with a sample-name
        # corresponding to a registered tube-id, but with a leading zero:
        # '034567890abc'. We want to show that _compare_samples_against_qiita()
        # is able to resolve '034567890abc' to '34567890abc' and doesn't
        # return it in 'samples_not_in_qiita'.

        # another-good-sample-sheet.csv contains one sample that is truly
        # unregistered with AnotherFakeQiita(): '4567890abcd'. We want to
        # demonstrate that after checking for leading zeroes that the function
        # decides this is still an unregistered sample and returns it in
        # 'samples_not_in_qiita'.

        # self.pipeline represents a metagenomic pathway.
        step = Step(self.another_pipeline, self.qiita_id, None)

        obs = step._compare_samples_against_qiita(fake_client)

        exp = [{'samples_not_in_qiita': {'4567890abcd'},
                'examples_in_qiita': ['1234567890a', '234567890ab',
                                      '34567890abc', 'BLANK1.1BCD'],
                'project_name': 'TestProject', 'total_in_qiita': 4,
                'used_tids': True,
                'messages': [
                    "The total number of samples found in TestProject that"
                    " aren't BLANK is: 4",
                    "Number of values in sheet that aren't sample-names in"
                    " Qiita: 4",
                    "Number of values in sheet that aren't tube-ids in "
                    "Qiita: 1",
                    "More values in sheet matched tube-ids than sample-names"
                    " with TestProject"]}]

        self.assertEqual(obs, exp)

    def atest_precheck(self):
        fake_client = AnotherFakeClient()

        # test that Step.precheck() raises a PipelineError with the correct
        # message, given the configuration of Step() and AnotherFakeClient().

        step = Step(self.another_pipeline, self.qiita_id, None)

        msg = ("The total number of samples found in TestProject that aren't"
               " BLANK is: 4\nNumber of values in sheet that aren't sample-"
               "names in Qiita: 4\nNumber of values in sheet that aren't tube"
               "-ids in Qiita: 1\nMore values in sheet matched tube-ids than"
               " sample-names with TestProject")

        with self.assertRaisesRegex(PipelineError, msg):
            step.precheck(fake_client)

    def atest_project_metadata_check(self):
        fake_client = FakeClient()

        # self.pipeline represents a metagenomic pathway.
        step = Step(self.pipeline, self.qiita_id, None)

        # _project_metadata_check() should return w/out raising an Error if
        # step and fake_client is used.
        step._project_metadata_check(fake_client)

        fake_client.info_in_11661['categories'].append('well_id_384')
        fake_client.info_in_13059['categories'].append('well_id_384')

        msg = ("'well_id_384' exists in Qiita study 13059's sample metadata"
               "\n'well_id_384' exists in Qiita study 11661's sample metadata")
        with self.assertRaisesRegex(PipelineError, msg):
            step._project_metadata_check(fake_client)

    def atest_conditional_fastqc_finder(self):
        self._create_alternate_test_input()

        # For a metagenomic pipeline, we expect indexed files to be removed
        # from the results. We also expect only trimmed files from Feist_11661
        # retrieved, and none from other projects, adapter-trimmed-only files,
        # or zero-length files.
        step = Step(self.pipeline_replicates, self.qiita_id, None)
        results = step._get_postqc_fastq_files(self.output_file_path,
                                               'Feist_11661')

        exp = {
            "raw_forward_seqs": [
                "/NuQCJob/Feist_11661/filtered_sequences/CDPH-SAL_Salmonella"
                "_Typhi_MDL-143_SXXX_L001_R1_001.trimmed.fastq.gz",
                "/NuQCJob/Feist_11661/filtered_sequences/CDPH-SAL_Salmonella"
                "_Typhi_MDL-144_SXXX_L001_R1_001.trimmed.fastq.gz",
                "/NuQCJob/Feist_11661/filtered_sequences/CDPH-SAL_Salmonella"
                "_Typhi_MDL-145_SXXX_L001_R1_001.trimmed.fastq.gz",
                "/NuQCJob/Feist_11661/filtered_sequences/CDPH-SAL_Salmonella"
                "_Typhi_MDL-146_SXXX_L001_R1_001.trimmed.fastq.gz",
                "/NuQCJob/Feist_11661/filtered_sequences/CDPH-SAL_Salmonella"
                "_Typhi_MDL-147_SXXX_L001_R1_001.trimmed.fastq.gz",
                "/NuQCJob/Feist_11661/filtered_sequences/CDPH-SAL_Salmonella"
                "_Typhi_MDL-150_SXXX_L001_R1_001.trimmed.fastq.gz",
                "/NuQCJob/Feist_11661/filtered_sequences/CDPH-SAL_Salmonella"
                "_Typhi_MDL-151_SXXX_L001_R1_001.trimmed.fastq.gz"
            ],
            "raw_reverse_seqs": [
                "/NuQCJob/Feist_11661/filtered_sequences/CDPH-SAL_Salmonella"
                "_Typhi_MDL-143_SXXX_L001_R2_001.trimmed.fastq.gz",
                "/NuQCJob/Feist_11661/filtered_sequences/CDPH-SAL_Salmonella"
                "_Typhi_MDL-144_SXXX_L001_R2_001.trimmed.fastq.gz",
                "/NuQCJob/Feist_11661/filtered_sequences/CDPH-SAL_Salmonella"
                "_Typhi_MDL-145_SXXX_L001_R2_001.trimmed.fastq.gz",
                "/NuQCJob/Feist_11661/filtered_sequences/CDPH-SAL_Salmonella"
                "_Typhi_MDL-146_SXXX_L001_R2_001.trimmed.fastq.gz",
                "/NuQCJob/Feist_11661/filtered_sequences/CDPH-SAL_Salmonella"
                "_Typhi_MDL-147_SXXX_L001_R2_001.trimmed.fastq.gz",
                "/NuQCJob/Feist_11661/filtered_sequences/CDPH-SAL_Salmonella"
                "_Typhi_MDL-150_SXXX_L001_R2_001.trimmed.fastq.gz",
                "/NuQCJob/Feist_11661/filtered_sequences/CDPH-SAL_Salmonella"
                "_Typhi_MDL-151_SXXX_L001_R2_001.trimmed.fastq.gz"
            ]
        }

        # metagenomic runs shouldn't return a set of data like exp above.
        # It shouldn't include I1 and I2 files.
        self.assertEqual(set(results.keys()), {'raw_forward_seqs',
                                               'raw_reverse_seqs'})
        for key in results.keys():
            # remove base output_file_path from the results.
            obs = [x.replace(self.output_file_path, '')
                   for x in results[key]]
            self.assertEqual(set(obs), set(exp[key]))

        # Hack an amplicon pipeline. reuse project-names, sample-names and
        # qiita-ids. Expected results should be just as they are for
        # metagenomic pipelines, except the index files are included.
        step = Step(self.amplicon_pipeline, self.qiita_id, None)

        exp['raw_barcodes'] = [
            "/NuQCJob/Feist_11661/filtered_sequences/CDPH-SAL_Salmonella"
            "_Typhi_MDL-143_SXXX_L001_I1_001.trimmed.fastq.gz",
            "/NuQCJob/Feist_11661/filtered_sequences/CDPH-SAL_Salmonella"
            "_Typhi_MDL-143_SXXX_L001_I2_001.trimmed.fastq.gz",
            "/NuQCJob/Feist_11661/filtered_sequences/CDPH-SAL_Salmonella"
            "_Typhi_MDL-144_SXXX_L001_I1_001.trimmed.fastq.gz",
            "/NuQCJob/Feist_11661/filtered_sequences/CDPH-SAL_Salmonella"
            "_Typhi_MDL-144_SXXX_L001_I2_001.trimmed.fastq.gz",
            "/NuQCJob/Feist_11661/filtered_sequences/CDPH-SAL_Salmonella"
            "_Typhi_MDL-145_SXXX_L001_I1_001.trimmed.fastq.gz",
            "/NuQCJob/Feist_11661/filtered_sequences/CDPH-SAL_Salmonella"
            "_Typhi_MDL-145_SXXX_L001_I2_001.trimmed.fastq.gz",
            "/NuQCJob/Feist_11661/filtered_sequences/CDPH-SAL_Salmonella"
            "_Typhi_MDL-146_SXXX_L001_I1_001.trimmed.fastq.gz",
            "/NuQCJob/Feist_11661/filtered_sequences/CDPH-SAL_Salmonella"
            "_Typhi_MDL-146_SXXX_L001_I2_001.trimmed.fastq.gz",
            "/NuQCJob/Feist_11661/filtered_sequences/CDPH-SAL_Salmonella"
            "_Typhi_MDL-147_SXXX_L001_I1_001.trimmed.fastq.gz",
            "/NuQCJob/Feist_11661/filtered_sequences/CDPH-SAL_Salmonella"
            "_Typhi_MDL-147_SXXX_L001_I2_001.trimmed.fastq.gz"
        ]

        results = step._get_postqc_fastq_files(self.output_file_path,
                                               'Feist_11661')

        self.assertEqual(set(results.keys()), {'raw_barcodes',
                                               'raw_forward_seqs',
                                               'raw_reverse_seqs'})
        for key in results.keys():
            obs = [x.replace(self.output_file_path, '')
                   for x in results[key]]
            self.assertEqual(set(obs), set(exp[key]))


class ReplicateTests(BaseStepTests):
    def setUp(self):
        super().setUp()

        self._create_test_input(3)

        # Fake enough of a run so that GenPrepFileJob can generate
        # prep-info files based on real input.

        # seqpro path
        self.seqpro_path = which('seqpro')

        self.project_list = ['Feist_11661', 'NYU_BMS_Melanoma_13059']

        # create Job working directories as needed.
        data_dir = partial(join, 'qp_klp', 'tests', 'data')
        self.output_dir = partial(data_dir, 'output_dir')
        run_dir = partial(self.output_dir, 'GenPrepFileJob',
                          '211021_A00000_0000_SAMPLE')

        demultiplex_stats_path = data_dir('Demultiplex_Stats.csv')
        fastp_stats_path = data_dir('sample_fastp.json')

        convert_job_reports_dir = self.output_dir('ConvertJob', 'Reports')
        qc_job_reports_dir = self.output_dir('NuQCJob', 'Feist_11661',
                                             'fastp_reports_dir', 'json')

        run_dir_stats_dir = run_dir('Stats')
        json_dir_13059 = run_dir('NYU_BMS_Melanoma_13059', 'json')
        self.fastq_dir_11661 = run_dir('Feist_11661', 'filtered_sequences')
        self.fastq_dir_13059 = run_dir('NYU_BMS_Melanoma_13059',
                                       'filtered_sequences')

        self.qc_fastq_dir_11661 = self.output_dir('NuQCJob', 'Feist_11661',
                                                  'filtered_sequences')

        self.qc_fastq_dir_13059 = self.output_dir('NuQCJob',
                                                  'NYU_BMS_Melanoma_13059',
                                                  'filtered_sequences')

        create_these = [run_dir_stats_dir, convert_job_reports_dir,
                        qc_job_reports_dir, json_dir_13059,
                        self.fastq_dir_11661, self.fastq_dir_13059]

        for some_dir in create_these:
            makedirs(some_dir, exist_ok=True)

        # Copy pre-made files containing numbers of reads for each sample
        # into place.

        copy(demultiplex_stats_path,
             self.output_dir('ConvertJob', 'Reports', 'Demultiplex_Stats.csv'))

        samples_11661 = ['BLANK_43_12H_A4', 'JBI_KHP_HGL_022_A16',
                         'BLANK_43_12G_B2', 'JBI_KHP_HGL_023_B18',
                         'RMA_KHP_rpoS_Mage_Q97N_A9', 'JBI_KHP_HGL_023_A18',
                         'RMA_KHP_rpoS_Mage_Q97L_B8', 'JBI_KHP_HGL_022_B16',
                         'JBI_KHP_HGL_022_A15', 'RMA_KHP_rpoS_Mage_Q97E_A12',
                         'RMA_KHP_rpoS_Mage_Q97L_A8',
                         'RMA_KHP_rpoS_Mage_Q97N_A10', 'JBI_KHP_HGL_021_B14',
                         'BLANK_43_12G_A1', 'BLANK_43_12H_B4',
                         'RMA_KHP_rpoS_Mage_Q97N_B10', 'JBI_KHP_HGL_024_A19',
                         'RMA_KHP_rpoS_Mage_Q97D_B6',
                         'RMA_KHP_rpoS_Mage_Q97D_A6', 'JBI_KHP_HGL_023_A17',
                         'RMA_KHP_rpoS_Mage_Q97E_A11', 'BLANK_43_12G_A2',
                         'RMA_KHP_rpoS_Mage_Q97D_A5', 'JBI_KHP_HGL_024_A20',
                         'JBI_KHP_HGL_024_B20', 'RMA_KHP_rpoS_Mage_Q97L_A7',
                         'JBI_KHP_HGL_021_A14', 'RMA_KHP_rpoS_Mage_Q97E_B12',
                         'BLANK_43_12H_A3', 'JBI_KHP_HGL_021_A13']

        samples_13059 = ['AP581451B02_A21', 'EP256645B01_A23',
                         'EP112567B02_C1', 'EP337425B01_C3',
                         'LP127890A01_C5', 'EP159692B04_C7',
                         'EP987683A01_C9', 'AP959450A03_C11',
                         'SP464350A04_C13', 'EP121011B01_C15',
                         'AP581451B02_A22', 'EP256645B01_A24',
                         'EP112567B02_C2', 'EP337425B01_C4',
                         'LP127890A01_C6', 'EP159692B04_C8',
                         'EP987683A01_C10', 'AP959450A03_C12',
                         'SP464350A04_C14', 'EP121011B01_C16',
                         'AP581451B02_B22', 'EP256645B01_B24',
                         'EP112567B02_D2', 'EP337425B01_D4',
                         'LP127890A01_D6', 'EP159692B04_D8',
                         'EP987683A01_D10', 'AP959450A03_D12',
                         'SP464350A04_D14', 'EP121011B01_D16']

        for sample in samples_11661:
            copy(fastp_stats_path, join(qc_job_reports_dir,
                                        f'{sample}_S270_L001_R1_001.json'))

        for sample in samples_13059:
            copy(fastp_stats_path, join(json_dir_13059,
                                        f'{sample}_S270_L001_R1_001.json'))

        # create fake fastq files. Metapool checks to confirm that they are
        # gzipped, hence we copy a small dummy fastq.gz file w/countable
        # sequences here.

        for name in samples_11661:
            for n_file in [f"{name}_S270_L001_R1_001.trimmed.fastq.gz",
                           f"{name}_S270_L001_R2_001.trimmed.fastq.gz"]:
                copy(self.dummy_fastq_file, join(self.fastq_dir_11661, n_file))

        for name in samples_13059:
            for n_file in [f"{name}_S270_L001_R1_001.trimmed.fastq.gz",
                           f"{name}_S270_L001_R2_001.trimmed.fastq.gz"]:
                copy(self.dummy_fastq_file, join(self.fastq_dir_13059, n_file))

    def tearDown(self):
        if exists(self.output_file_path):
            rmtree(self.output_file_path)

    def atest_replicates(self):
        self.maxDiff = None

        # Create run_dir_stats_dir Step object and generate prep-files.
        step = Step(self.pipeline_replicates, self.qiita_id, None)

        # Fake an empty tube-id map so that _generate_prep_file() doesn't
        # abort early.
        step.tube_id_map = {}

        config = self.pipeline.config_profile['profile']['configuration']

        job = step._generate_prep_file(config['seqpro'],
                                       self.sheet_w_replicates_path,
                                       self.seqpro_path)

        # Metagenomic.generate_prep_file() and Amplicon.generate_prep_file()
        # both perform self.prep_file_paths = job.prep_file_paths after the
        # above completes. Recreate that behavior here.
        step.prep_file_paths = job.prep_file_paths

        prep_output_path = self.output_dir('GenPrepFileJob', 'PrepFiles')

        # Confirm that the generated prep-info files are found in the
        # correct location w/the correct names. For testing purposes, remove
        # the absolute path up to 'qp_klp' and test against relative paths
        # found in exp.

        obs = step.prep_file_paths

        for qiita_id in obs:
            obs[qiita_id] = [re.sub(r"^.*?\/qp_klp\/", "qp_klp/", x) for x in
                             obs[qiita_id]]

        # _generate_prep_file() should have created a prep-info file for each
        # of the three replicates and two projects defined in the original
        # sample-sheet, for a total of six files. When replicates are defined,
        # prep_output_path should contain numerically named directories, one
        # for each replicate, containing the prep-info files for that
        # replicate.

        exp = {
            "11661": [
                join(prep_output_path, "1",
                     "211021_A00000_0000_SAMPLE.Feist_11661.1.tsv"),
                join(prep_output_path, "2",
                     "211021_A00000_0000_SAMPLE.Feist_11661.1.tsv"),
                join(prep_output_path, "3",
                     "211021_A00000_0000_SAMPLE.Feist_11661.1.tsv")
            ],
            "13059": [
                join(prep_output_path, "1",
                     "211021_A00000_0000_SAMPLE.NYU_BMS_Melanoma_13059.1.tsv"),
                join(prep_output_path, "2",
                     "211021_A00000_0000_SAMPLE.NYU_BMS_Melanoma_13059.1.tsv"),
                join(prep_output_path, "3",
                     "211021_A00000_0000_SAMPLE.NYU_BMS_Melanoma_13059.1.tsv")
            ]
        }

        self.assertDictEqual(obs, exp)

        # verify each prep_info_file contains the expected number of rows and
        # more importantly are not empty.
        for qiita_id in obs:
            for prep_info_file in obs[qiita_id]:
                df = parse_prep(prep_info_file)
                self.assertEqual(10, df.shape[0])

        # Now, let's fake enough of the later steps to test the loading of
        # each prep-info file into Qiita and ensure that each post is unique
        # from a unique file.

        fake_client = FakeClient()
        step.update_prep_templates(fake_client)
        step.generate_special_map(fake_client)

        # Since load_preps_into_qiita() relies on the hierarchy of files in
        # NuQCJob()'s output to do it's work, copy the faked fastq files made
        # for GenPrepFileJob() into the right location in NuQCJob's output.
        # (We are skipping a number of steps, hence the need to do this.
        #  Ordinarily, the data would be created in NuQCJob() and migrated
        #  through the pipeline to GenPrepFileJob(), not the other way around.)
        copytree(self.fastq_dir_13059, self.qc_fastq_dir_13059)
        copytree(self.fastq_dir_11661, self.qc_fastq_dir_11661)
        results = step.load_preps_into_qiita(fake_client)

        obs = {
            'RMA_KHP_rpoS_Mage_Q97L_A7_S270_L001_R1_001.trimmed.fastq.gz': 0,
            'RMA_KHP_rpoS_Mage_Q97E_A12_S270_L001_R2_001.trimmed.fastq.gz': 0,
            'JBI_KHP_HGL_022_B16_S270_L001_R2_001.trimmed.fastq.gz': 0,
            'AP581451B02_A21_S270_L001_R2_001.trimmed.fastq.gz': 0,
            'EP159692B04_C8_S270_L001_R1_001.trimmed.fastq.gz': 0,
            'LP127890A01_D6_S270_L001_R1_001.trimmed.fastq.gz': 0}

        # take advantage of the posted-files listing being in string format
        # and test for the appearance of specific individual file names in
        # each post. Each fastq file should only appear in one and only one
        # post.
        for prep_id in fake_client.saved_posts:
            posted_files_listing = fake_client.saved_posts[prep_id]['files']
            for fastq in obs:
                if fastq in posted_files_listing:
                    obs[fastq] += 1

        exp = {
            'RMA_KHP_rpoS_Mage_Q97L_A7_S270_L001_R1_001.trimmed.fastq.gz': 1,
            'RMA_KHP_rpoS_Mage_Q97E_A12_S270_L001_R2_001.trimmed.fastq.gz': 1,
            'JBI_KHP_HGL_022_B16_S270_L001_R2_001.trimmed.fastq.gz': 1,
            'AP581451B02_A21_S270_L001_R2_001.trimmed.fastq.gz': 1,
            'EP159692B04_C8_S270_L001_R1_001.trimmed.fastq.gz': 1,
            'LP127890A01_D6_S270_L001_R1_001.trimmed.fastq.gz': 1}

        self.assertDictEqual(obs, exp)

        # confirm that six preps were created by confirming that
        # load_preps_into_qiita() returned six unique prep-ids, and six
        # unique job-ids.
        self.assertEqual(len(results['Linking JobID'].unique()), 6)
        self.assertEqual(len(results['Qiita Prep ID'].unique()), 6)

        # confirm that three of the preps belong to the Feist project, while
        # the other three belong to NYU.
        projects = ['Feist_11661', 'NYU_BMS_Melanoma_13059']
        for project in projects:
            self.assertEqual(results['Project'].value_counts()[project], 3)


class FailedSamplesRecordTests(TestCase):
    def setUp(self):
        class MockSample():
            def __init__(self, sample_id, project_name):
                self.Sample_ID = sample_id
                self.Sample_Project = project_name

        self.samples = [MockSample('A', 'ProjectOne'),
                        MockSample('B', 'ProjectTwo'),
                        MockSample('C', 'ProjectThree'),
                        MockSample('D', 'ProjectFour')]

        self.output = TemporaryDirectory()

    def atest_failed_samples_record(self):
        fsr = FailedSamplesRecord(self.output.name, self.samples)

        # assert that a state file doesn't already exist and attempt to load()
        # it. State should remain unchanged.
        exp = deepcopy(fsr.sample_state)
        self.assertFalse(exists(fsr.output_path))
        fsr.load()
        self.assertEqual(fsr.sample_state, exp)

        # confirm that dump() creates the appropriate file.
        self.assertFalse(exists(fsr.output_path))
        fsr.dump()
        self.assertTrue(exists(fsr.output_path))

        # load the dumped() file, and confirm that nothing changed since
        # the state wasn't update()d.
        fsr.load()
        self.assertEqual(fsr.sample_state, exp)

        # assert samples A and C failed at the ConvertJob stage, assert
        # state changed. dump() state and load() it. Confirm state on disk
        # reflects changes.
        fsr.update(['A', 'C'], 'ConvertJob')
        self.assertNotEqual(fsr.sample_state, exp)
        fsr.dump()
        fsr.load()
        exp = {'A': 'ConvertJob', 'B': None, 'C': 'ConvertJob', 'D': None}
        self.assertEqual(fsr.sample_state, exp)

        # B should be failed at FastQCJob but A and C should still be
        # failed at ConvertJob
        fsr.update(['A', 'B', 'C'], "FastQCJob")
        exp = {'A': 'ConvertJob', 'B': 'FastQCJob',
               'C': 'ConvertJob', 'D': None}
        self.assertEqual(fsr.sample_state, exp)

        # confirm html file exists. Assume pandas DataFrame.to_html() works
        # as intended.
        self.assertFalse(exists(fsr.report_path))
        fsr.generate_report()
        self.assertTrue(exists(fsr.report_path))
