# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from unittest import TestCase
from qp_klp.klp_util import parse_prep_file
from os.path import join


class KLPUtilTests(TestCase):
    def test_parse_prep_file(self):
        good_prep_file = join('qp_klp', 'tests', 'good-prep-file-small.txt')

        obs = parse_prep_file(good_prep_file)

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
