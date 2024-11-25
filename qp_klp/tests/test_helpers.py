# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from unittest import TestCase
from os.path import join, abspath, exists
from os import makedirs
from shutil import rmtree
from os import remove, getcwd


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


class TestHelpers(TestCase):
    def setUp(self):
        self.fake_bin_path = ""
        self.delete_these_files = []
        self.delete_these_dirs = []
        # self.fake_bin_path = self.get_searchable_path()

        # self.output_dir represents a qiita working directory.
        package_root = abspath('./qp_klp/tests/data')
        self.output_dir = join(package_root,
                               "077c4da8-74eb-4184-8860-0207f53623be")
        self.delete_these_dirs = [self.output_dir]

        # We want a clean directory, nothing from leftover runs
        makedirs(self.output_dir, exist_ok=False)

        self.debug = False

    def tearDown(self):
        if not self.debug:
            rmtree(self.output_dir)

            for fp in self.delete_these_files:
                if exists(fp):
                    remove(fp)

    def test_foo(self):
        self.assertTrue(True)
