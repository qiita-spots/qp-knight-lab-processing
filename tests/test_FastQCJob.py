import unittest
from os.path import join, exists, isfile
from functools import partial
from sequence_processing_pipeline.FastQCJob import FastQCJob
from sequence_processing_pipeline.PipelineError import (PipelineError,
                                                        JobFailedError)
from os import makedirs, listdir, mkdir
from shutil import rmtree, move
from json import load


class TestFastQCJob(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        package_root = ''
        self.path = partial(join, package_root, 'tests', 'data')
        self.qiita_job_id = 'abcdabcdabcdabcdabcdabcdabcdabcd'
        self.output_path = self.path('output_dir2')
        self.fastqc_log_path = join(self.output_path, 'logs')
        self.raw_fastq_files_path = ('tests/data'
                                     '/211021_A00000_0000_SAMPLE/Data/Fastq/p'
                                     'roject1')
        self.processed_fastq_files_path = ('tests'
                                           '/data/211021_A00000_0000_SAMPLE/sa'
                                           'mple-sequence-directory')
        self.config_yml = join(package_root, 'multiqc-bclconvert-config.yaml')
        self.qc_root_path = join(self.output_path, 'QCJob')
        makedirs(self.qc_root_path, exist_ok=True)

        self.StudyB_ids = [
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R08624',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_1_69',
            'JM-Metabolic__GN0_2148', 'JM-Metabolic__GN0_2183',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0400',
            'CDPH-SAL_Salmonella_Typhi_MDL-144',
            'P21_E_coli_ELI361', 'JBI_KHP_HGL_023',
            'Pputida_PALE__HGL_Pputida_158',
            'CDPH-SAL_Salmonella_Typhi_MDL-166',
            'CDPH-SAL_Salmonella_Typhi_MDL-147',
            'JM-Metabolic__GN03252',
            'stALE_E_coli_A11_F43_I1_R1',
            'P21_E_coli_ELI355',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_33_50',
            'stALE_E_coli_A14_F42_I1_R1',
            'Deoxyribose_PALE_ALE__MG1655_BOP27_10_51',
            'JM-Metabolic__GN0_2099', 'JM-Metabolic__GN02657',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_5_22',
            'stALE_E_coli_A1_F21_I1_R1',
            'JM-Metabolic__GN0_2005',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0421',
            'Pputida_TALE__HGL_Pputida_129',
            'Deoxyribose_PALE_ALE__MG1655_BOP27_10_28',
            'Pputida_TALE__HGL_Pputida_113',
            'JM-Metabolic__GN02487', 'P21_E_coli_ELI348',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_26_27',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0385',
            'stALE_E_coli_A10_F43_I1_R1',
            'stALE_E_coli_A11_F119_I1_R1',
            'JM-Metabolic__GN03132',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0381',
            'Pputida_TALE__HGL_Pputida_131',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11044',
            'Pputida_TALE__HGL_Pputida_122',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0419',
            'Pputida_TALE__HGL_Pputida_114', 'P21_E_coli_ELI366',
            'Pputida_TALE__HGL_Pputida_140', 'JBI_KHP_HGL_026',
            'Deoxyribose_PALE_ALE__MG1655_BOP27_6_21',
            'stALE_E_coli_A9_F21_I1_R1', 'P21_E_coli_ELI359',
            'JM-Metabolic__GN04306',
            'Pputida_JBEI__HGL_Pputida_108_BP7',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_20_71',
            'P21_E_coli_ELI345',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_20_43',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_18_59',
            'stALE_E_coli_A12_F136_I1_R1',
            'JM-Metabolic__GN04540',
            'Pputida_PALE__HGL_Pputida_153',
            'CDPH-SAL_Salmonella_Typhi_MDL-156',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_5_36',
            'RMA_KHP_rpoS_Mage_Q97L',
            'stALE_E_coli_A13_F42_I1_R1',
            'JM-Metabolic__GN02766',
            'stALE_E_coli_A7_F42_I1_R1',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_22_28',
            'Pputida_TALE__HGL_Pputida_144',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0520',
            'Pputida_TALE__HGL_Pputida_117',
            'stALE_E_coli_A5_F21_I1_R1',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0375',
            'stALE_E_coli_A4_F42_I1_R1',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_22_16',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_3_61',
            'Pputida_TALE__HGL_Pputida_116',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0353',
            'JM-Metabolic__GN02769', 'JM-Metabolic__GN04488',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0403',
            'Pputida_PALE__HGL_Pputida_149',
            'JM-Metabolic__GN0_2290', 'JM-Metabolic__GN02787',
            'JM-Metabolic__GN02449',
            'JBI_KHP_HGL_030_Amitesh_soxR_oxyR',
            'Pputida_JBEI__HGL_Pputida_110_M2',
            'JM-Metabolic__GN0_2169',
            'stALE_E_coli_A8_F42_I1_R1',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11103',
            'CDPH-SAL_Salmonella_Typhi_MDL-148',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0366',
            'Pputida_TALE__HGL_Pputida_125',
            'stALE_E_coli_A11_F21_I1_R1',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_19_25',
            'CDPH-SAL_Salmonella_Typhi_MDL-150',
            'stALE_E_coli_A15_F21_I1_R1',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_31_24',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_17_25',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0399',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0389',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0418',
            'Pputida_TALE__HGL_Pputida_130',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0486',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0373',
            'CDPH-SAL_Salmonella_Typhi_MDL-162',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_29_49',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_3_23',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_32_56',
            'CDPH-SAL_Salmonella_Typhi_MDL-160',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0474',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0354',
            'stALE_E_coli_A14_F133_I1_R1',
            'Pputida_PALE__HGL_Pputida_173',
            'JM-Metabolic__GN0_2175',
            'stALE_E_coli_A13_F20_I1_R1',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0395',
            'stALE_E_coli_A17_F118_I1_R1',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_17_64',
            'CDPH-SAL_Salmonella_Typhi_MDL-154',
            'Pputida_PALE__HGL_Pputida_156',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0326',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_32_6',
            'JM-Metabolic__GN0_2215',
            'Pputida_PALE__HGL_Pputida_162',
            'Pputida_TALE__HGL_Pputida_143',
            'stALE_E_coli_A12_F21_I1_R1',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0382',
            'JBI_KHP_HGL_022',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_17_58',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0524',
            'Deoxyribose_PALE_ALE__MG1655_BOP27_10_13',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0376',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_28_53',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0473',
            'JM-Metabolic__GN0_2380', 'RMA_KHP_rpoS_Mage_Q97D',
            'JBI_KHP_HGL_025', 'Pputida_TALE__HGL_Pputida_135',
            'stALE_E_coli_A12_F43_I1_R1',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_26_69',
            'JM-Metabolic__GN0_2172',
            'stALE_E_coli_A15_F42_I1_R1',
            'JM-Metabolic__GN0_2007',
            'Pputida_PALE__HGL_Pputida_154',
            'stALE_E_coli_A16_F134_I1_R1',
            'CDPH-SAL_Salmonella_Typhi_MDL-149',
            'P21_E_coli_ELI352',
            'CDPH-SAL_Salmonella_Typhi_MDL-164',
            'Pputida_PALE__HGL_Pputida_167',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_30_7',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0484',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0357',
            'stALE_E_coli_A4_F21_I1_R1', 'JM-Metabolic__GN02748',
            'JBI_KHP_HGL_021',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0372',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_7_51',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0352',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11102',
            'JM-Metabolic__GN05002',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0398',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0483',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0390',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_29_25',
            'JM-Metabolic__GN02590', 'stALE_E_coli_A8_F20_I1_R1',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0392',
            'stALE_E_coli_A6_F21_I1_R1',
            'JM-Metabolic__GN0_2317',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0364',
            'Pputida_TALE__HGL_Pputida_139',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0383',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_1_57',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0388',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_22_52',
            'P21_E_coli_ELI364', 'Pputida_PALE__HGL_Pputida_166',
            'Pputida_PALE__HGL_Pputida_175',
            'Pputida_PALE__HGL_Pputida_164',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0330',
            'CDPH-SAL_Salmonella_Typhi_MDL-143',
            'Pputida_PALE__HGL_Pputida_146',
            'Pputida_TALE__HGL_Pputida_121',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_21_46',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0397',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_20_16',
            'P21_E_coli_ELI350', 'Pputida_PALE__HGL_Pputida_161',
            'RMA_KHP_rpoS_Mage_Q97E',
            'stALE_E_coli_A5_F42_I1_R1',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R08704',
            'Pputida_PALE__HGL_Pputida_176',
            'Pputida_TALE__HGL_Pputida_137',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0409',
            'JM-Metabolic__GN04563',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_33_41',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0394',
            'Pputida_PALE__HGL_Pputida_172', 'P21_E_coli_ELI367',
            'JM-Metabolic__GN04014',
            'CDPH-SAL_Salmonella_Typhi_MDL-146',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_28_28',
            'Pputida_PALE__HGL_Pputida_160',
            'CDPH-SAL_Salmonella_Typhi_MDL-158',
            'P21_E_coli_ELI363',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0523',
            'JM-Metabolic__GN02567',
            'Pputida_JBEI__HGL_Pputida_111_M5',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0417',
            'CDPH-SAL_Salmonella_Typhi_MDL-155',
            'Pputida_TALE__HGL_Pputida_133',
            'Pputida_PALE__HGL_Pputida_145',
            'Pputida_JBEI__HGL_Pputida_107_BP6',
            'P21_E_coli_ELI368', 'JM-Metabolic__GN02531',
            'JM-Metabolic__GN0_2165',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_19_55',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0355',
            'JM-Metabolic__GN02446',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_24_24',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11135',
            'CDPH-SAL_Salmonella_Typhi_MDL-165',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_21_51',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_1_24',
            'Pputida_TALE__HGL_Pputida_132',
            'stALE_E_coli_A10_F131_I1_R1',
            'stALE_E_coli_A17_F21_I1_R1',
            'JM-Metabolic__GN02424',
            'JM-Metabolic__GN02529',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R10727',
            'Pputida_TALE__HGL_Pputida_119',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0378',
            'CDPH-SAL_Salmonella_Typhi_MDL-157',
            'JM-Metabolic__GN0_2354',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0407',
            'JBI_KHP_HGL_029_Amitesh_oxyR',
            'JM-Metabolic__GN04094',
            'Pputida_PALE__HGL_Pputida_170',
            'CDPH-SAL_Salmonella_Typhi_MDL-159',
            'Pputida_TALE__HGL_Pputida_118',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0393',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0408',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0522',
            'Deoxyribose_PALE_ALE__MG1655_BOP27_4_48',
            'JM-Metabolic__GN03409',
            'Pputida_TALE__HGL_Pputida_120',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0329',
            'P21_E_coli_ELI358', 'P21_E_coli_ELI344',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_5_46',
            'stALE_E_coli_A15_F117_I1_R1',
            'Pputida_PALE__HGL_Pputida_163',
            'Pputida_TALE__HGL_Pputida_112',
            'JM-Metabolic__GN05377',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_3_50',
            'stALE_E_coli_A7_F21_I1_R1',
            'stALE_E_coli_A9_F44_I1_R1',
            'JM-Metabolic__GN0_2009',
            'Pputida_PALE__HGL_Pputida_152',
            'CDPH-SAL_Salmonella_Typhi_MDL-152',
            'JM-Metabolic__GN04682',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_19_63',
            'JM-Metabolic__GN0_2277',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0525',
            'RMA_KHP_rpoS_Mage_Q97N',
            'Pputida_TALE__HGL_Pputida_128',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_7_23',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0377',
            'CDPH-SAL_Salmonella_Typhi_MDL-153',
            'P21_E_coli_ELI349', 'Pputida_PALE__HGL_Pputida_150',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0519',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_18_35',
            'Pputida_PALE__HGL_Pputida_151', 'JBI_KHP_HGL_027',
            'JBI_KHP_HGL_028_Amitesh_soxR',
            'stALE_E_coli_A3_F40_I1_R1',
            'CDPH-SAL_Salmonella_Typhi_MDL-168',
            'JM-Metabolic__GN0_2404',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11078',
            'stALE_E_coli_A2_F21_I1_R1',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0368',
            'Pputida_TALE__HGL_Pputida_134',
            'Pputida_TALE__HGL_Pputida_141',
            'stALE_E_coli_A16_F20_I1_R1',
            'CDPH-SAL_Salmonella_Typhi_MDL-145',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11101',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0401',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0391',
            'JM-Metabolic__GN05367',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0517',
            'JM-Metabolic__GN0_2375',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0420',
            'P21_E_coli_ELI365',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_7_41',
            'JM-Metabolic__GN0_2337',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_32_20',
            'Deoxyribose_PALE_ALE__MG1655_BOP27_4_14',
            'P21_E_coli_ELI354',
            'Deoxyribose_PALE_ALE__MG1655_BOP27_4_23',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_28_13',
            'Pputida_PALE__HGL_Pputida_171',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0396',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0380',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0518',
            'JM-Metabolic__GN02501',
            'Pputida_TALE__HGL_Pputida_138',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0387',
            'CDPH-SAL_Salmonella_Typhi_MDL-167',
            'CDPH-SAL_Salmonella_Typhi_MDL-161',
            'Pputida_PALE__HGL_Pputida_168',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0367',
            'stALE_E_coli_A3_F18_I1_R1',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_30_60',
            'Pputida_TALE__HGL_Pputida_124',
            'Pputida_TALE__HGL_Pputida_142',
            'JM-Metabolic__GN04255',
            'stALE_E_coli_A6_F43_I1_R1',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0406',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0328',
            'stALE_E_coli_A14_F20_I1_R1',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0370',
            'stALE_E_coli_A18_F130_I1_R1',
            'Pputida_PALE__HGL_Pputida_169',
            'P21_E_coli_ELI369',
            'P21_E_coli_ELI347', 'stALE_E_coli_A13_F121_I1_R1',
            'JM-Metabolic__GN0_2393',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_30_22',
            'P21_E_coli_ELI362',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0402',
            'JBI_KHP_HGL_024', 'stALE_E_coli_A18_F18_I1_R1',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0356',
            'JM-Metabolic__GN0_2094',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_33_21',
            'stALE_E_coli_A18_F39_I1_R1',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0404',
            'Pputida_PALE__HGL_Pputida_157', 'P21_E_coli_ELI351',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0386',
            'Pputida_PALE__HGL_Pputida_148',
            'Deoxyribose_PALE_ALE__MG1655_BOP27_6_35',
            'P21_E_coli_ELI353', 'Pputida_PALE__HGL_Pputida_147',
            'JM-Metabolic__GN05128',
            'Pputida_PALE__HGL_Pputida_174',
            'Pputida_TALE__HGL_Pputida_115',
            'Pputida_PALE__HGL_Pputida_159',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0485',
            'stALE_E_coli_A4_F21_I1_R2', 'JM-Metabolic__GN04665',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0371',
            'CDPH-SAL_Salmonella_Typhi_MDL-163',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0327',
            'stALE_E_coli_A16_F42_I1_R1',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0374',
            'JM-Metabolic__GN03218',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_31_62',
            'Pputida_PALE__HGL_Pputida_155',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11153',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_24_9',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0384',
            'JM-Metabolic__GN05109',
            'JBI_KHP_HGL_031_Amitesh_rpoS',
            'CDPH-SAL_Salmonella_Typhi_MDL-151',
            'JM-Metabolic__GN04428', 'JM-Metabolic__GN02514',
            'Pputida_JBEI__HGL_Pputida_109_BP8',
            'JM-Metabolic__GN0_2254',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_26_6',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_18_19',
            'Pputida_TALE__HGL_Pputida_123',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_31_42',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0405',
            'P21_E_coli_ELI357',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0369',
            'Pputida_TALE__HGL_Pputida_136',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0521',
            'Pputida_PALE__HGL_Pputida_165',
            'Pputida_TALE__HGL_Pputida_127',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_21_23',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_29_57',
            'stALE_E_coli_A10_F21_I1_R1',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11154',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0516',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_24_52',
            'Pputida_TALE__HGL_Pputida_126',
            'JM-Metabolic__GN04612']

        self.StudyC_ids = ['5B', '3A', '6A', '8A', '4A', '7A', 'GFR', 'ISB']

        self.nyu_ids = ['22_001_710_503_791_00', '22_001_801_552_503_00',
                        'AP006367B02', 'AP029018B01', 'AP032412B04',
                        'AP032413B04', 'AP046324B02', 'AP046327B02',
                        'AP062219B03', 'AP065292B01', 'AP103463B01',
                        'AP173299B04', 'AP173301B04', 'AP173305B04',
                        'AP223470B01', 'AP298002B02', 'AP309872B03',
                        'AP324642B04', 'AP470859B01', 'AP481403B02',
                        'AP531397B04', 'AP549678B01', 'AP549681B02',
                        'AP568785B04', 'AP568787B02', 'AP581451B02',
                        'AP616837B04', 'AP668628B04', 'AP668631B04',
                        'AP687591B04', 'AP696363B02', 'AP732307B04',
                        'AP744361A02', 'AP745799A04', 'AP771472A04',
                        'AP780167B02', 'AP787247B04', 'AP795068B04',
                        'AP891020A04', 'AP905750A02', 'AP911328B01',
                        'AP953594A02', 'AP959450A03', 'AP967057A04', 'C14',
                        'C18', 'C20', 'C3', 'C5', 'C6', 'C9', 'EP001624B01',
                        'EP001625B01', 'EP012991B03', 'EP023801B04',
                        'EP023808B02', 'EP032410B02', 'EP032412B02',
                        'EP043583B01', 'EP054632B01', 'EP061002B01',
                        'EP073160B01', 'EP073209B02', 'EP073216B01',
                        'EP087938B02', 'EP090129B04', 'EP112567B02',
                        'EP121011B01', 'EP121013B01', 'EP128904B02',
                        'EP128910B01', 'EP159692B04', 'EP159695B01',
                        'EP163771B01', 'EP182060B03', 'EP182065B04',
                        'EP182243B02', 'EP182346B04', 'EP184255B04',
                        'EP190307B01', 'EP202095B04', 'EP202452B01',
                        'EP207036B01', 'EP207041B01', 'EP207042B04',
                        'EP212214B01', 'EP216516B04', 'EP230245B01',
                        'EP238034B01', 'EP244360B01', 'EP244366B01',
                        'EP256644B01', 'EP256645B01', 'EP260543B04',
                        'EP260544B04', 'EP273332B04', 'EP282107B01',
                        'EP282108B01', 'EP282276B04', 'EP291979B04',
                        'EP291980B04', 'EP305735B04', 'EP316863B03',
                        'EP320438B01', 'EP333541B04', 'EP337325B04',
                        'EP337425B01', 'EP339053B02', 'EP339057B02',
                        'EP339059B02', 'EP339061B02', 'EP372981B04',
                        'EP379938B01', 'EP385379B01', 'EP385384B01',
                        'EP385387B01', 'EP393712B02', 'EP393714B01',
                        'EP393715B01', 'EP393717B01', 'EP393718B01',
                        'EP400447B04', 'EP400448B04', 'EP410041B01',
                        'EP410042B01', 'EP410046B01', 'EP422407B01',
                        'EP431562B04', 'EP431570B01', 'EP431575B01',
                        'EP446602B01', 'EP446604B03', 'EP446610B02',
                        'EP447926B04', 'EP447927B04', 'EP447928B04',
                        'EP447929B04', 'EP447940B04', 'EP447975B02',
                        'EP448041B04', 'EP451428B04', 'EP455757B04',
                        'EP455759B04', 'EP455763B04', 'EP479266B04',
                        'EP479270B03', 'EP479794B02', 'EP479894B04',
                        'EP483291B04', 'EP484973B04', 'EP487995B04',
                        'EP504030B04', 'EP529635B02', 'EP533388B01',
                        'EP533389B03', 'EP533426B03', 'EP533429B04',
                        'EP542577B04', 'EP542578B04', 'EP554501B04',
                        'EP554506B04', 'EP554513B02', 'EP554515B04',
                        'EP554518B04', 'EP573296B01', 'EP573310B01',
                        'EP573313B01', 'EP584756B04', 'EP587475B04',
                        'EP587476B04', 'EP587477B04', 'EP587478B04',
                        'EP606652B04', 'EP606656B03', 'EP606662B04',
                        'EP606663B04', 'EP617440B01', 'EP617441B01',
                        'EP617442B01', 'EP617443B01', 'EP636802A01',
                        'EP649418A02', 'EP649623A01', 'EP649653A04',
                        'EP649737A03', 'EP656055A04', 'EP657260A01',
                        'EP657385A04', 'EP657386A01', 'EP667743A04',
                        'EP675042B01', 'EP675044A01', 'EP675075A04',
                        'EP683835A01', 'EP685640B01', 'EP702221B04',
                        'EP718687A04', 'EP718688A01', 'EP721390A04',
                        'EP724905B01', 'EP727972A04', 'EP729433A02',
                        'EP729434A01', 'EP738468A01', 'EP738469A01',
                        'EP749735A07', 'EP759450A04', 'EP768164A02',
                        'EP768748A04', 'EP772143A02', 'EP772145A02',
                        'EP784608A01', 'EP786631A04', 'EP790019A01',
                        'EP790020A02', 'EP790021A04', 'EP790023A01',
                        'EP805337A01', 'EP808104A01', 'EP808105A01',
                        'EP808106A01', 'EP808109A01', 'EP808110A04',
                        'EP808111A03', 'EP808112A04', 'EP843906A04',
                        'EP846485A01', 'EP868682A01', 'EP872341A01',
                        'EP876243A04', 'EP882752A01', 'EP886422A01',
                        'EP890157A02', 'EP890158A02', 'EP899038A04',
                        'EP905975A04', 'EP915769A04', 'EP921593A04',
                        'EP921594A04', 'EP927458A04', 'EP927459A04',
                        'EP927461A04', 'EP927462A02', 'EP929277A02',
                        'EP940013A01', 'EP944059A02', 'EP970001A01',
                        'EP970005A01', 'EP980752B04', 'EP981129A02',
                        'EP987683A01', 'EP996831B04', 'LP127767A01',
                        'LP127890A01', 'LP128476A01', 'LP128479A01',
                        'LP128538A01', 'LP128539A01', 'LP128540A01',
                        'LP128541A01', 'LP128543A01', 'LP154981A01',
                        'LP154986A01', 'LP166715A01', 'LP169879A01',
                        'LP191039A01', 'LP196272A01', 'SP205732A02',
                        'SP205754A01', 'SP229387A04', 'SP230380A02',
                        'SP230381A01', 'SP230382A04', 'SP231628A02',
                        'SP231629A02', 'SP231630A02', 'SP231631A02',
                        'SP232077A04', 'SP232079A01', 'SP232114A04',
                        'SP232270A02', 'SP232309A01', 'SP232310A04',
                        'SP232311A04', 'SP235186A04', 'SP235189A01',
                        'SP246941A01', 'SP247340A04', 'SP257517A04',
                        'SP257519A04', 'SP280481A02', 'SP284095A03',
                        'SP284096A02', 'SP317293A02', 'SP317297A02',
                        'SP331134A04', 'SP335002A04', 'SP353893A02',
                        'SP365864A04', 'SP388683A02', 'SP399724A04',
                        'SP404403A02', 'SP404405A02', 'SP404409A02',
                        'SP404412A02', 'SP408629A01', 'SP410793A01',
                        'SP410796A02', 'SP415021A02', 'SP415023A02',
                        'SP415025A01', 'SP415030A01', 'SP416130A04',
                        'SP453872A01', 'SP464350A04', 'SP464352A03',
                        'SP471496A04', 'SP478193A02', 'SP490298A02',
                        'SP491897A02', 'SP491898A02', 'SP491900A02',
                        'SP491907A02', 'SP503615A02', 'SP506933A04',
                        'SP511289A02', 'SP511294A04', 'SP515443A04',
                        'SP515763A04', 'SP531696A04', 'SP561451A04',
                        'SP573823A04', 'SP573824A04', 'SP573843A04',
                        'SP573849A04', 'SP573859A04', 'SP573860A01',
                        'SP577399A02', 'SP584547A02', 'SP584551A08',
                        'SP612495A04', 'SP612496A01', 'SP631994A04',
                        'SP640978A02', 'SP641029A02', 'SP645141A03',
                        'SP681591A04', 'SP683466A02', 'SP704319A04',
                        'SP754514A04', 'ep256643b01', 'lp127896a01']

        self.sample_ids = self.StudyB_ids + self.StudyC_ids + self.nyu_ids

        base_path = join(self.output_path, 'FastQCJob')
        file_names = []
        makedirs(join(base_path, 'StudyB_11661', 'trimmed_sequences'),
                 exist_ok=True)
        makedirs(join(base_path, 'StudyC_6123', 'filtered_sequences'),
                 exist_ok=True)

        file_names.append(join(base_path,
                               'StudyB_11661',
                               'trimmed_sequences',
                               ('CDPH-SAL_Salmonella_Typhi_MDL-143_S1_L003_R1_'
                                '001.trimmed_fastqc.html')))

        file_names.append(join(base_path,
                               'StudyC_6123',
                               'filtered_sequences',
                               '3A_S169_L003_R2_001.trimmed_fastqc.html'))

        # write files out to disk
        for file_name in file_names:
            with open(file_name, 'w') as f2:
                f2.write("This is a file.")

        # set up dummy logs
        self.fastqc_log_path = join(self.output_path, "FastQCJob", "logs")
        makedirs(self.fastqc_log_path, exist_ok=True)

        log_files = {
            'slurm-9999999_35.out': ["---------------",
                                     "Run details:",
                                     ("hds-fe848a9e-c0e9-49d9-978d-"
                                      "27565a314e8b 1908305 b2-018"),
                                     "---------------",
                                     "+ this",
                                     "+ that",
                                     "+ blah",
                                     ("something error: Generic Standin Error"
                                      " (GSE).")],
            'slurm-9999999_17.out': ["---------------",
                                     "Run details:",
                                     ("hds-fe848a9e-c0e9-49d9-978d-"
                                      "27565a314e8b 1908305 b2-018"),
                                     "---------------",
                                     "+ this",
                                     "+ that",
                                     "+ blah",
                                     ("something error: Another Standin Error"
                                      " (ASE).")]
        }

        for log_file in log_files:
            fp = join(self.fastqc_log_path, log_file)

            with open(fp, 'w') as f:
                lines = log_files[log_file]
                for line in lines:
                    f.write(f"{line}\n")

    def tearDown(self):
        rmtree(self.output_path)

        zero_path = join(self.raw_fastq_files_path, 'zero_files')

        if exists(zero_path):
            for f in listdir(zero_path):
                source = join(zero_path, f)
                move(source, join(self.raw_fastq_files_path, f))

            rmtree(zero_path)

    def test_generate_job_scripts(self):
        job = FastQCJob(self.qc_root_path, self.output_path,
                        self.raw_fastq_files_path.replace('/project1', ''),
                        self.processed_fastq_files_path,
                        16, 16,
                        'tests/bin/fastqc',
                        ['my_module.1.1'], self.qiita_job_id, 'queue_name',
                        4, 23, '8g', 30, 1000, False)

        job_script_path = join(job.output_path, 'FastQCJob.sh')
        array_details_path = join(job.output_path, 'FastQCJob.array-details')
        self.assertEqual(exists(job_script_path), True)
        self.assertEqual(exists(array_details_path), True)

        exp = ["#!/bin/bash",
               "#SBATCH -J abcdabcdabcdabcdabcdabcdabcdabcd_FastQCJob",
               "#SBATCH -p queue_name", "#SBATCH -N 4", "#SBATCH -n 16",
               "#SBATCH --time 23", "#SBATCH --mem 8gG",
               "#SBATCH --array 1-4%30", "set -x", "set +e",
               "set -o pipefail", "date", "hostname",
               "echo ${SLURM_JOBID} ${SLURM_ARRAY_TASK_ID}",
               "cd tests/data/output_dir2/"
               "FastQCJob", "", "module load my_module.1.1", "",
               "step=${SLURM_ARRAY_TASK_ID}",
               "cmd0=$(head -n $step tests/data/"
               "output_dir2/FastQCJob/FastQCJob.array-details | tail -n 1)",
               "eval $cmd0", "echo \"Cmd Completed: $cmd0\" > logs/FastQCJob_"
               "$step.completed"]

        with open(job_script_path, 'r') as f:
            obs = f.readlines()
            obs = [x.strip() for x in obs]

        for a, b in zip(obs, exp):
            self.assertEqual(a, b)

        exp = ["fastqc --noextract -t 16 tests/da"
               "ta/211021_A00000_0000_SAMPLE/Data/Fastq/project1/sample1_R1_.f"
               "astq.gz tests/data/211021_A00000_"
               "0000_SAMPLE/Data/Fastq/project1/sample1_R2_.fastq.gz -o "
               "tests/data/output_dir2/FastQCJob/fastqc"
               "/project1/bclconvert",
               "fastqc --noextract -t 16 tests/da"
               "ta/211021_A00000_0000_SAMPLE/Data/Fastq/project1/sample2_R1_.f"
               "astq.gz tests/data/211021_A00000_"
               "0000_SAMPLE/Data/Fastq/project1/sample2_R2_.fastq.gz -o "
               "tests/data/output_dir2/FastQCJob/fastqc"
               "/project1/bclconvert",
               "fastqc --noextract -t 16 tests/da"
               "ta/211021_A00000_0000_SAMPLE/sample-sequence-directory/project"
               "1/filtered_sequences/sample1_R1_.trimmed.fastq.gz "
               "tests/data/211021_A00000_0000_SAMPLE/sample-s"
               "equence-directory/project1/filtered_sequences/sample1_R2_.trim"
               "med.fastq.gz -o tests/data/output"
               "_dir2/FastQCJob/fastqc/project1/filtered_sequences",
               "fastqc --noextract -t 16 tests/da"
               "ta/211021_A00000_0000_SAMPLE/sample-sequence-directory/project"
               "1/filtered_sequences/sample2_R1_.trimmed.fastq.gz "
               "tests/data/211021_A00000_0000_SAMPLE/sample-s"
               "equence-directory/project1/filtered_sequences/sample2_R2_.trim"
               "med.fastq.gz -o tests/data/output"
               "_dir2/FastQCJob/fastqc/project1/filtered_sequences"]

        with open(array_details_path, 'r') as f:
            obs = f.readlines()
            obs = [x.strip() for x in obs]

        for a, b in zip(obs, exp):
            self.assertEqual(a, b)

    def test_audit(self):
        job = FastQCJob(self.qc_root_path, self.output_path,
                        self.raw_fastq_files_path.replace('/project1', ''),
                        self.processed_fastq_files_path,
                        16, 16,
                        'tests/bin/fastqc', [],
                        self.qiita_job_id, 'queue_name', 4, 23, '8g', 30,
                        1000, False)

        obs = job.audit(self.sample_ids)

        exp = ['22_001_710_503_791_00', '22_001_801_552_503_00', '4A',
               '5B', '6A', '7A', '8A',
               'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_17_25',
               'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_17_58',
               'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_17_64',
               'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_19_25',
               'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_19_55',
               'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_19_63',
               'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_1_24',
               'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_1_57',
               'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_1_69',
               'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_21_23',
               'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_21_46',
               'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_21_51',
               'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_29_25',
               'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_29_49',
               'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_29_57',
               'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_31_24',
               'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_31_42',
               'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_31_62',
               'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_33_21',
               'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_33_41',
               'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_33_50',
               'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_3_23',
               'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_3_50',
               'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_3_61',
               'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_5_22',
               'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_5_36',
               'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_5_46',
               'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_7_23',
               'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_7_41',
               'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_7_51',
               'AP006367B02', 'AP029018B01', 'AP032412B04', 'AP032413B04',
               'AP046324B02', 'AP046327B02', 'AP062219B03', 'AP065292B01',
               'AP103463B01', 'AP173299B04', 'AP173301B04', 'AP173305B04',
               'AP223470B01', 'AP298002B02', 'AP309872B03', 'AP324642B04',
               'AP470859B01', 'AP481403B02', 'AP531397B04', 'AP549678B01',
               'AP549681B02', 'AP568785B04', 'AP568787B02', 'AP581451B02',
               'AP616837B04', 'AP668628B04', 'AP668631B04', 'AP687591B04',
               'AP696363B02', 'AP732307B04', 'AP744361A02', 'AP745799A04',
               'AP771472A04', 'AP780167B02', 'AP787247B04', 'AP795068B04',
               'AP891020A04', 'AP905750A02', 'AP911328B01', 'AP953594A02',
               'AP959450A03', 'AP967057A04', 'C14', 'C18', 'C20', 'C3', 'C5',
               'C6', 'C9',
               'CDPH-SAL_Salmonella_Typhi_MDL-144',
               'CDPH-SAL_Salmonella_Typhi_MDL-145',
               'CDPH-SAL_Salmonella_Typhi_MDL-146',
               'CDPH-SAL_Salmonella_Typhi_MDL-147',
               'CDPH-SAL_Salmonella_Typhi_MDL-148',
               'CDPH-SAL_Salmonella_Typhi_MDL-149',
               'CDPH-SAL_Salmonella_Typhi_MDL-150',
               'CDPH-SAL_Salmonella_Typhi_MDL-151',
               'CDPH-SAL_Salmonella_Typhi_MDL-152',
               'CDPH-SAL_Salmonella_Typhi_MDL-153',
               'CDPH-SAL_Salmonella_Typhi_MDL-154',
               'CDPH-SAL_Salmonella_Typhi_MDL-155',
               'CDPH-SAL_Salmonella_Typhi_MDL-156',
               'CDPH-SAL_Salmonella_Typhi_MDL-157',
               'CDPH-SAL_Salmonella_Typhi_MDL-158',
               'CDPH-SAL_Salmonella_Typhi_MDL-159',
               'CDPH-SAL_Salmonella_Typhi_MDL-160',
               'CDPH-SAL_Salmonella_Typhi_MDL-161',
               'CDPH-SAL_Salmonella_Typhi_MDL-162',
               'CDPH-SAL_Salmonella_Typhi_MDL-163',
               'CDPH-SAL_Salmonella_Typhi_MDL-164',
               'CDPH-SAL_Salmonella_Typhi_MDL-165',
               'CDPH-SAL_Salmonella_Typhi_MDL-166',
               'CDPH-SAL_Salmonella_Typhi_MDL-167',
               'CDPH-SAL_Salmonella_Typhi_MDL-168',
               'Deoxyribose_PALE_ALE__MG1655_BOP27_10_13',
               'Deoxyribose_PALE_ALE__MG1655_BOP27_10_28',
               'Deoxyribose_PALE_ALE__MG1655_BOP27_10_51',
               'Deoxyribose_PALE_ALE__MG1655_BOP27_4_14',
               'Deoxyribose_PALE_ALE__MG1655_BOP27_4_23',
               'Deoxyribose_PALE_ALE__MG1655_BOP27_4_48',
               'Deoxyribose_PALE_ALE__MG1655_BOP27_6_21',
               'Deoxyribose_PALE_ALE__MG1655_BOP27_6_35',
               'Deoxyribose_PALE_ALE__MG1655_Lib4_18_19',
               'Deoxyribose_PALE_ALE__MG1655_Lib4_18_35',
               'Deoxyribose_PALE_ALE__MG1655_Lib4_18_59',
               'Deoxyribose_PALE_ALE__MG1655_Lib4_20_16',
               'Deoxyribose_PALE_ALE__MG1655_Lib4_20_43',
               'Deoxyribose_PALE_ALE__MG1655_Lib4_20_71',
               'Deoxyribose_PALE_ALE__MG1655_Lib4_22_16',
               'Deoxyribose_PALE_ALE__MG1655_Lib4_22_28',
               'Deoxyribose_PALE_ALE__MG1655_Lib4_22_52',
               'Deoxyribose_PALE_ALE__MG1655_Lib4_24_24',
               'Deoxyribose_PALE_ALE__MG1655_Lib4_24_52',
               'Deoxyribose_PALE_ALE__MG1655_Lib4_24_9',
               'Deoxyribose_PALE_ALE__MG1655_Lib4_26_27',
               'Deoxyribose_PALE_ALE__MG1655_Lib4_26_6',
               'Deoxyribose_PALE_ALE__MG1655_Lib4_26_69',
               'Deoxyribose_PALE_ALE__MG1655_Lib4_28_13',
               'Deoxyribose_PALE_ALE__MG1655_Lib4_28_28',
               'Deoxyribose_PALE_ALE__MG1655_Lib4_28_53',
               'Deoxyribose_PALE_ALE__MG1655_Lib4_30_22',
               'Deoxyribose_PALE_ALE__MG1655_Lib4_30_60',
               'Deoxyribose_PALE_ALE__MG1655_Lib4_30_7',
               'Deoxyribose_PALE_ALE__MG1655_Lib4_32_20',
               'Deoxyribose_PALE_ALE__MG1655_Lib4_32_56',
               'Deoxyribose_PALE_ALE__MG1655_Lib4_32_6', 'EP001624B01',
               'EP001625B01', 'EP012991B03', 'EP023801B04', 'EP023808B02',
               'EP032410B02', 'EP032412B02', 'EP043583B01', 'EP054632B01',
               'EP061002B01', 'EP073160B01', 'EP073209B02', 'EP073216B01',
               'EP087938B02', 'EP090129B04', 'EP112567B02', 'EP121011B01',
               'EP121013B01', 'EP128904B02', 'EP128910B01', 'EP159692B04',
               'EP159695B01', 'EP163771B01', 'EP182060B03', 'EP182065B04',
               'EP182243B02', 'EP182346B04', 'EP184255B04', 'EP190307B01',
               'EP202095B04', 'EP202452B01', 'EP207036B01', 'EP207041B01',
               'EP207042B04', 'EP212214B01', 'EP216516B04', 'EP230245B01',
               'EP238034B01', 'EP244360B01', 'EP244366B01', 'EP256644B01',
               'EP256645B01', 'EP260543B04', 'EP260544B04', 'EP273332B04',
               'EP282107B01', 'EP282108B01', 'EP282276B04', 'EP291979B04',
               'EP291980B04', 'EP305735B04', 'EP316863B03', 'EP320438B01',
               'EP333541B04', 'EP337325B04', 'EP337425B01', 'EP339053B02',
               'EP339057B02', 'EP339059B02', 'EP339061B02', 'EP372981B04',
               'EP379938B01', 'EP385379B01', 'EP385384B01', 'EP385387B01',
               'EP393712B02', 'EP393714B01', 'EP393715B01', 'EP393717B01',
               'EP393718B01', 'EP400447B04', 'EP400448B04', 'EP410041B01',
               'EP410042B01', 'EP410046B01', 'EP422407B01', 'EP431562B04',
               'EP431570B01', 'EP431575B01', 'EP446602B01', 'EP446604B03',
               'EP446610B02', 'EP447926B04', 'EP447927B04', 'EP447928B04',
               'EP447929B04', 'EP447940B04', 'EP447975B02', 'EP448041B04',
               'EP451428B04', 'EP455757B04', 'EP455759B04', 'EP455763B04',
               'EP479266B04', 'EP479270B03', 'EP479794B02', 'EP479894B04',
               'EP483291B04', 'EP484973B04', 'EP487995B04', 'EP504030B04',
               'EP529635B02', 'EP533388B01', 'EP533389B03', 'EP533426B03',
               'EP533429B04', 'EP542577B04', 'EP542578B04', 'EP554501B04',
               'EP554506B04', 'EP554513B02', 'EP554515B04', 'EP554518B04',
               'EP573296B01', 'EP573310B01', 'EP573313B01', 'EP584756B04',
               'EP587475B04', 'EP587476B04', 'EP587477B04', 'EP587478B04',
               'EP606652B04', 'EP606656B03', 'EP606662B04', 'EP606663B04',
               'EP617440B01', 'EP617441B01', 'EP617442B01', 'EP617443B01',
               'EP636802A01', 'EP649418A02', 'EP649623A01', 'EP649653A04',
               'EP649737A03', 'EP656055A04', 'EP657260A01', 'EP657385A04',
               'EP657386A01', 'EP667743A04', 'EP675042B01', 'EP675044A01',
               'EP675075A04', 'EP683835A01', 'EP685640B01', 'EP702221B04',
               'EP718687A04', 'EP718688A01', 'EP721390A04', 'EP724905B01',
               'EP727972A04', 'EP729433A02', 'EP729434A01', 'EP738468A01',
               'EP738469A01', 'EP749735A07', 'EP759450A04', 'EP768164A02',
               'EP768748A04', 'EP772143A02', 'EP772145A02', 'EP784608A01',
               'EP786631A04', 'EP790019A01', 'EP790020A02', 'EP790021A04',
               'EP790023A01', 'EP805337A01', 'EP808104A01', 'EP808105A01',
               'EP808106A01', 'EP808109A01', 'EP808110A04', 'EP808111A03',
               'EP808112A04', 'EP843906A04', 'EP846485A01', 'EP868682A01',
               'EP872341A01', 'EP876243A04', 'EP882752A01', 'EP886422A01',
               'EP890157A02', 'EP890158A02', 'EP899038A04', 'EP905975A04',
               'EP915769A04', 'EP921593A04', 'EP921594A04', 'EP927458A04',
               'EP927459A04', 'EP927461A04', 'EP927462A02', 'EP929277A02',
               'EP940013A01', 'EP944059A02', 'EP970001A01', 'EP970005A01',
               'EP980752B04', 'EP981129A02', 'EP987683A01', 'EP996831B04',
               'GFR', 'ISB', 'JBI_KHP_HGL_021', 'JBI_KHP_HGL_022',
               'JBI_KHP_HGL_023', 'JBI_KHP_HGL_024', 'JBI_KHP_HGL_025',
               'JBI_KHP_HGL_026', 'JBI_KHP_HGL_027',
               'JBI_KHP_HGL_028_Amitesh_soxR', 'JBI_KHP_HGL_029_Amitesh_oxyR',
               'JBI_KHP_HGL_030_Amitesh_soxR_oxyR',
               'JBI_KHP_HGL_031_Amitesh_rpoS',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0326',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0327',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0328',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0329',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0330',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0352',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0353',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0354',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0355',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0356',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0357',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0364',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0366',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0367',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0368',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0369',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0370',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0371',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0372',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0373',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0374',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0375',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0376',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0377',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0378',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0380',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0381',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0382',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0383',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0384',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0385',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0386',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0387',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0388',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0389',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0390',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0391',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0392',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0393',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0394',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0395',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0396',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0397',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0398',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0399',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0400',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0401',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0402',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0403',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0404',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0405',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0406',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0407',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0408',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0409',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0417',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0418',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0419',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0420',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0421',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0473',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0474',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0483',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0484',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0485',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0486',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0516',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0517',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0518',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0519',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0520',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0521',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0522',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0523',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0524',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0525',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-R08624',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-R08704',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-R10727',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11044',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11078',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11101',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11102',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11103',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11135',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11153',
               'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11154',
               'JM-Metabolic__GN02424', 'JM-Metabolic__GN02446',
               'JM-Metabolic__GN02449', 'JM-Metabolic__GN02487',
               'JM-Metabolic__GN02501', 'JM-Metabolic__GN02514',
               'JM-Metabolic__GN02529', 'JM-Metabolic__GN02531',
               'JM-Metabolic__GN02567', 'JM-Metabolic__GN02590',
               'JM-Metabolic__GN02657', 'JM-Metabolic__GN02748',
               'JM-Metabolic__GN02766', 'JM-Metabolic__GN02769',
               'JM-Metabolic__GN02787', 'JM-Metabolic__GN03132',
               'JM-Metabolic__GN03218', 'JM-Metabolic__GN03252',
               'JM-Metabolic__GN03409', 'JM-Metabolic__GN04014',
               'JM-Metabolic__GN04094', 'JM-Metabolic__GN04255',
               'JM-Metabolic__GN04306', 'JM-Metabolic__GN04428',
               'JM-Metabolic__GN04488', 'JM-Metabolic__GN04540',
               'JM-Metabolic__GN04563', 'JM-Metabolic__GN04612',
               'JM-Metabolic__GN04665', 'JM-Metabolic__GN04682',
               'JM-Metabolic__GN05002', 'JM-Metabolic__GN05109',
               'JM-Metabolic__GN05128', 'JM-Metabolic__GN05367',
               'JM-Metabolic__GN05377', 'JM-Metabolic__GN0_2005',
               'JM-Metabolic__GN0_2007', 'JM-Metabolic__GN0_2009',
               'JM-Metabolic__GN0_2094', 'JM-Metabolic__GN0_2099',
               'JM-Metabolic__GN0_2148', 'JM-Metabolic__GN0_2165',
               'JM-Metabolic__GN0_2169', 'JM-Metabolic__GN0_2172',
               'JM-Metabolic__GN0_2175', 'JM-Metabolic__GN0_2183',
               'JM-Metabolic__GN0_2215', 'JM-Metabolic__GN0_2254',
               'JM-Metabolic__GN0_2277', 'JM-Metabolic__GN0_2290',
               'JM-Metabolic__GN0_2317', 'JM-Metabolic__GN0_2337',
               'JM-Metabolic__GN0_2354', 'JM-Metabolic__GN0_2375',
               'JM-Metabolic__GN0_2380', 'JM-Metabolic__GN0_2393',
               'JM-Metabolic__GN0_2404', 'LP127767A01', 'LP127890A01',
               'LP128476A01', 'LP128479A01', 'LP128538A01', 'LP128539A01',
               'LP128540A01', 'LP128541A01', 'LP128543A01', 'LP154981A01',
               'LP154986A01', 'LP166715A01', 'LP169879A01', 'LP191039A01',
               'LP196272A01', 'P21_E_coli_ELI344', 'P21_E_coli_ELI345',
               'P21_E_coli_ELI347', 'P21_E_coli_ELI348', 'P21_E_coli_ELI349',
               'P21_E_coli_ELI350', 'P21_E_coli_ELI351', 'P21_E_coli_ELI352',
               'P21_E_coli_ELI353', 'P21_E_coli_ELI354', 'P21_E_coli_ELI355',
               'P21_E_coli_ELI357', 'P21_E_coli_ELI358', 'P21_E_coli_ELI359',
               'P21_E_coli_ELI361', 'P21_E_coli_ELI362', 'P21_E_coli_ELI363',
               'P21_E_coli_ELI364', 'P21_E_coli_ELI365', 'P21_E_coli_ELI366',
               'P21_E_coli_ELI367', 'P21_E_coli_ELI368', 'P21_E_coli_ELI369',
               'Pputida_JBEI__HGL_Pputida_107_BP6',
               'Pputida_JBEI__HGL_Pputida_108_BP7',
               'Pputida_JBEI__HGL_Pputida_109_BP8',
               'Pputida_JBEI__HGL_Pputida_110_M2',
               'Pputida_JBEI__HGL_Pputida_111_M5',
               'Pputida_PALE__HGL_Pputida_145',
               'Pputida_PALE__HGL_Pputida_146',
               'Pputida_PALE__HGL_Pputida_147',
               'Pputida_PALE__HGL_Pputida_148',
               'Pputida_PALE__HGL_Pputida_149',
               'Pputida_PALE__HGL_Pputida_150',
               'Pputida_PALE__HGL_Pputida_151',
               'Pputida_PALE__HGL_Pputida_152',
               'Pputida_PALE__HGL_Pputida_153',
               'Pputida_PALE__HGL_Pputida_154',
               'Pputida_PALE__HGL_Pputida_155',
               'Pputida_PALE__HGL_Pputida_156',
               'Pputida_PALE__HGL_Pputida_157',
               'Pputida_PALE__HGL_Pputida_158',
               'Pputida_PALE__HGL_Pputida_159',
               'Pputida_PALE__HGL_Pputida_160',
               'Pputida_PALE__HGL_Pputida_161',
               'Pputida_PALE__HGL_Pputida_162',
               'Pputida_PALE__HGL_Pputida_163',
               'Pputida_PALE__HGL_Pputida_164',
               'Pputida_PALE__HGL_Pputida_165',
               'Pputida_PALE__HGL_Pputida_166',
               'Pputida_PALE__HGL_Pputida_167',
               'Pputida_PALE__HGL_Pputida_168',
               'Pputida_PALE__HGL_Pputida_169',
               'Pputida_PALE__HGL_Pputida_170',
               'Pputida_PALE__HGL_Pputida_171',
               'Pputida_PALE__HGL_Pputida_172',
               'Pputida_PALE__HGL_Pputida_173',
               'Pputida_PALE__HGL_Pputida_174',
               'Pputida_PALE__HGL_Pputida_175',
               'Pputida_PALE__HGL_Pputida_176',
               'Pputida_TALE__HGL_Pputida_112',
               'Pputida_TALE__HGL_Pputida_113',
               'Pputida_TALE__HGL_Pputida_114',
               'Pputida_TALE__HGL_Pputida_115',
               'Pputida_TALE__HGL_Pputida_116',
               'Pputida_TALE__HGL_Pputida_117',
               'Pputida_TALE__HGL_Pputida_118',
               'Pputida_TALE__HGL_Pputida_119',
               'Pputida_TALE__HGL_Pputida_120',
               'Pputida_TALE__HGL_Pputida_121',
               'Pputida_TALE__HGL_Pputida_122',
               'Pputida_TALE__HGL_Pputida_123',
               'Pputida_TALE__HGL_Pputida_124',
               'Pputida_TALE__HGL_Pputida_125',
               'Pputida_TALE__HGL_Pputida_126',
               'Pputida_TALE__HGL_Pputida_127',
               'Pputida_TALE__HGL_Pputida_128',
               'Pputida_TALE__HGL_Pputida_129',
               'Pputida_TALE__HGL_Pputida_130',
               'Pputida_TALE__HGL_Pputida_131',
               'Pputida_TALE__HGL_Pputida_132',
               'Pputida_TALE__HGL_Pputida_133',
               'Pputida_TALE__HGL_Pputida_134',
               'Pputida_TALE__HGL_Pputida_135',
               'Pputida_TALE__HGL_Pputida_136',
               'Pputida_TALE__HGL_Pputida_137',
               'Pputida_TALE__HGL_Pputida_138',
               'Pputida_TALE__HGL_Pputida_139',
               'Pputida_TALE__HGL_Pputida_140',
               'Pputida_TALE__HGL_Pputida_141',
               'Pputida_TALE__HGL_Pputida_142',
               'Pputida_TALE__HGL_Pputida_143',
               'Pputida_TALE__HGL_Pputida_144', 'RMA_KHP_rpoS_Mage_Q97D',
               'RMA_KHP_rpoS_Mage_Q97E', 'RMA_KHP_rpoS_Mage_Q97L',
               'RMA_KHP_rpoS_Mage_Q97N', 'SP205732A02', 'SP205754A01',
               'SP229387A04', 'SP230380A02', 'SP230381A01', 'SP230382A04',
               'SP231628A02', 'SP231629A02', 'SP231630A02', 'SP231631A02',
               'SP232077A04', 'SP232079A01', 'SP232114A04', 'SP232270A02',
               'SP232309A01', 'SP232310A04', 'SP232311A04', 'SP235186A04',
               'SP235189A01', 'SP246941A01', 'SP247340A04', 'SP257517A04',
               'SP257519A04', 'SP280481A02', 'SP284095A03', 'SP284096A02',
               'SP317293A02', 'SP317297A02', 'SP331134A04', 'SP335002A04',
               'SP353893A02', 'SP365864A04', 'SP388683A02', 'SP399724A04',
               'SP404403A02', 'SP404405A02', 'SP404409A02', 'SP404412A02',
               'SP408629A01', 'SP410793A01', 'SP410796A02', 'SP415021A02',
               'SP415023A02', 'SP415025A01', 'SP415030A01', 'SP416130A04',
               'SP453872A01', 'SP464350A04', 'SP464352A03', 'SP471496A04',
               'SP478193A02', 'SP490298A02', 'SP491897A02', 'SP491898A02',
               'SP491900A02', 'SP491907A02', 'SP503615A02', 'SP506933A04',
               'SP511289A02', 'SP511294A04', 'SP515443A04', 'SP515763A04',
               'SP531696A04', 'SP561451A04', 'SP573823A04', 'SP573824A04',
               'SP573843A04', 'SP573849A04', 'SP573859A04', 'SP573860A01',
               'SP577399A02', 'SP584547A02', 'SP584551A08', 'SP612495A04',
               'SP612496A01', 'SP631994A04', 'SP640978A02', 'SP641029A02',
               'SP645141A03', 'SP681591A04', 'SP683466A02', 'SP704319A04',
               'SP754514A04', 'ep256643b01', 'lp127896a01',
               'stALE_E_coli_A10_F131_I1_R1', 'stALE_E_coli_A10_F21_I1_R1',
               'stALE_E_coli_A10_F43_I1_R1', 'stALE_E_coli_A11_F119_I1_R1',
               'stALE_E_coli_A11_F21_I1_R1', 'stALE_E_coli_A11_F43_I1_R1',
               'stALE_E_coli_A12_F136_I1_R1', 'stALE_E_coli_A12_F21_I1_R1',
               'stALE_E_coli_A12_F43_I1_R1', 'stALE_E_coli_A13_F121_I1_R1',
               'stALE_E_coli_A13_F20_I1_R1', 'stALE_E_coli_A13_F42_I1_R1',
               'stALE_E_coli_A14_F133_I1_R1', 'stALE_E_coli_A14_F20_I1_R1',
               'stALE_E_coli_A14_F42_I1_R1', 'stALE_E_coli_A15_F117_I1_R1',
               'stALE_E_coli_A15_F21_I1_R1', 'stALE_E_coli_A15_F42_I1_R1',
               'stALE_E_coli_A16_F134_I1_R1', 'stALE_E_coli_A16_F20_I1_R1',
               'stALE_E_coli_A16_F42_I1_R1', 'stALE_E_coli_A17_F118_I1_R1',
               'stALE_E_coli_A17_F21_I1_R1', 'stALE_E_coli_A18_F130_I1_R1',
               'stALE_E_coli_A18_F18_I1_R1', 'stALE_E_coli_A18_F39_I1_R1',
               'stALE_E_coli_A1_F21_I1_R1', 'stALE_E_coli_A2_F21_I1_R1',
               'stALE_E_coli_A3_F18_I1_R1', 'stALE_E_coli_A3_F40_I1_R1',
               'stALE_E_coli_A4_F21_I1_R1', 'stALE_E_coli_A4_F21_I1_R2',
               'stALE_E_coli_A4_F42_I1_R1', 'stALE_E_coli_A5_F21_I1_R1',
               'stALE_E_coli_A5_F42_I1_R1', 'stALE_E_coli_A6_F21_I1_R1',
               'stALE_E_coli_A6_F43_I1_R1', 'stALE_E_coli_A7_F21_I1_R1',
               'stALE_E_coli_A7_F42_I1_R1', 'stALE_E_coli_A8_F20_I1_R1',
               'stALE_E_coli_A8_F42_I1_R1', 'stALE_E_coli_A9_F21_I1_R1',
               'stALE_E_coli_A9_F44_I1_R1']

        self.assertEqual(sorted(obs), sorted(exp))

        # these fake sample-ids should be returned by audit() as missing.
        obs = job.audit(self.sample_ids + ['not-a-sample', 'BLANK1'])
        obs = sorted(obs)
        exp = sorted(exp + ['BLANK1', 'not-a-sample'])

        self.assertListEqual(obs, exp)

    def test_all_zero_files(self):
        # create a 'zero_files' sub-folder and move all files there
        # temporarily.
        mkdir(join(self.raw_fastq_files_path, 'zero_files'))

        for f in listdir(self.raw_fastq_files_path):
            source = join(self.raw_fastq_files_path, f)
            if isfile(source):
                move(source, join(self.raw_fastq_files_path, 'zero_files', f))

        # attempt to create a FastQCJob when all input fastq files are in
        # zero_files. This will cause an Error if all fastq files for just
        # one project are 'zero-filed', even if there are other projects where
        # this is not the case.
        with self.assertRaises(PipelineError) as e:
            FastQCJob(self.qc_root_path, self.output_path,
                      self.raw_fastq_files_path.replace('/project1', ''),
                      self.processed_fastq_files_path,
                      16, 16,
                      'tests/bin/fastqc', [],
                      self.qiita_job_id, 'queue_name', 4, 23, '8g', 30,
                      1000, False)

        self.assertEqual(str(e.exception), "There are no fastq files for "
                                           "FastQCJob to process in "
                                           "tests/data/"
                                           "211021_A00000_0000_SAMPLE/Data/"
                                           "Fastq.")

    def test_completed_file_generation(self):
        job = FastQCJob(self.qc_root_path, self.output_path,
                        self.raw_fastq_files_path.replace('/project1', ''),
                        self.processed_fastq_files_path,
                        16, 16,
                        'tests/bin/fastqc', [],
                        self.qiita_job_id, 'queue_name', 4, 23, '8g', 30,
                        1000, False)

        my_path = join(self.output_path, 'FastQCJob', 'logs')

        # since .completed files are generated when jobs are submitted via
        # the run() method and completed successfully, we must manually
        # create the files _was_successful() expects to see.
        for i in range(1, 5):
            with open(join(my_path, f'project1_{i}.completed'), 'w') as f:
                f.write("This is a .completed file.")

        self.assertTrue(len(job._get_failed_indexes('1234.barnacle')) == 0)

    def test_completed_file_generation_some_failures(self):
        job = FastQCJob(self.qc_root_path, self.output_path,
                        self.raw_fastq_files_path.replace('/project1', ''),
                        self.processed_fastq_files_path,
                        16, 16,
                        'tests/bin/fastqc', [],
                        self.qiita_job_id, 'queue_name', 4, 23, '8g', 30,
                        1000, False)

        my_path = join(self.output_path, 'FastQCJob', 'logs')

        # simulate one job failing and not generating a .completed file by
        # setting range to (0,3) from (0,4). get_failed_indexes() should
        # return [3] as a result.
        for i in range(0, 3):
            with open(join(my_path, f'project1_{i}.completed'), 'w') as f:
                f.write("This is a .completed file.")

        failed_indexes = job._get_failed_indexes('2345.barnacle')
        self.assertTrue(failed_indexes, [3])

        # verify that when one or more array jobs have failed to complete,
        # a file is created that gives a job id and an array index.
        log_fp = join(self.output_path,
                      'FastQCJob',
                      'logs',
                      'failed_indexes_2345.barnacle.json')

        self.assertTrue(exists(log_fp))

        with open(log_fp, 'r') as f:
            obs = load(f)
            exp = {"job_id": "2345.barnacle",
                   "failed_indexes": [3, 4]}
            self.assertDictEqual(obs, exp)

    def test_error_msg_from_logs(self):
        job = FastQCJob(self.qc_root_path, self.output_path,
                        self.raw_fastq_files_path.replace('/project1', ''),
                        self.processed_fastq_files_path,
                        16, 16,
                        'tests/bin/fastqc', [],
                        self.qiita_job_id, 'queue_name', 4, 23, '8g', 30,
                        1000, False)

        self.assertFalse(job is None)

        # an internal method to force submit_job() to raise a JobFailedError
        # instead of submitting the job w/sbatch and waiting for a failed
        # job w/squeue.
        self.assertTrue(job._toggle_force_job_fail())

        try:
            job.run()
        except JobFailedError as e:
            # assert that the text of the original error message was
            # preserved, while including known strings from each of the
            # sample log-files.
            print(str(e))
            self.assertIn('This job died.', str(e))
            self.assertIn('something error: Generic Standin Error (GSE)',
                          str(e))
            self.assertIn('something error: Another Standin Error (ASE)',
                          str(e))


if __name__ == '__main__':
    unittest.main()
