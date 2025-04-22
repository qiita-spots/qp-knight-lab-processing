from os.path import join, abspath, exists, dirname
from os import makedirs
from sequence_processing_pipeline.ConvertJob import ConvertJob
from sequence_processing_pipeline.PipelineError import (PipelineError,
                                                        JobFailedError)
from functools import partial
import unittest
from shutil import rmtree


class TestConvertJob(unittest.TestCase):
    def setUp(self):
        # A list of sample-ids found within good-sample-sheet.csv.
        self.feist_ids = ['JM-Metabolic__GN04255',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_5_22',
                          'Pputida_PALE__HGL_Pputida_148',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0384',
                          'stALE_E_coli_A8_F42_I1_R1', 'JBI_KHP_HGL_022',
                          'Pputida_TALE__HGL_Pputida_126',
                          'Pputida_TALE__HGL_Pputida_115',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_28_28',
                          'CDPH-SAL_Salmonella_Typhi_MDL-144',
                          'Deoxyribose_PALE_ALE__MG1655_BOP27_4_48',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0357',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_30_22',
                          'CDPH-SAL_Salmonella_Typhi_MDL-159',
                          'JM-Metabolic__GN03409',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_32_20',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0364',
                          'JM-Metabolic__GN02531',
                          'stALE_E_coli_A10_F131_I1_R1',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0399',
                          'Pputida_PALE__HGL_Pputida_155',
                          'Pputida_PALE__HGL_Pputida_166',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0327',
                          'P21_E_coli_ELI350',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_17_64',
                          'JM-Metabolic__GN0_2393',
                          'JBI_KHP_HGL_029_Amitesh_oxyR',
                          'stALE_E_coli_A11_F21_I1_R1', 'P21_E_coli_ELI363',
                          'stALE_E_coli_A18_F39_I1_R1',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_7_51',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11044',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_32_6',
                          'JBI_KHP_HGL_023', 'stALE_E_coli_A12_F21_I1_R1',
                          'Pputida_PALE__HGL_Pputida_149',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0385',
                          'CDPH-SAL_Salmonella_Typhi_MDL-145',
                          'Pputida_TALE__HGL_Pputida_114',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0378',
                          'Pputida_TALE__HGL_Pputida_127', 'P21_E_coli_ELI362',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_30_60',
                          'P21_E_coli_ELI351',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0326',
                          'Pputida_JBEI__HGL_Pputida_108_BP7',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_21_23',
                          'CDPH-SAL_Salmonella_Typhi_MDL-158',
                          'JM-Metabolic__GN0_2007',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0356',
                          'Pputida_PALE__HGL_Pputida_167',
                          'JM-Metabolic__GN0_2165',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_24_9',
                          'stALE_E_coli_A13_F20_I1_R1',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0398',
                          'JM-Metabolic__GN04563',
                          'Pputida_PALE__HGL_Pputida_154', 'P21_E_coli_ELI344',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0419',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0370',
                          'Pputida_PALE__HGL_Pputida_172',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_21_46',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_5_46',
                          'JM-Metabolic__GN0_2354',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_22_52',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_1_24',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_18_59',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0404',
                          'P21_E_coli_ELI359', 'Pputida_TALE__HGL_Pputida_142',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0516',
                          'JM-Metabolic__GN0_2317',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_5_36',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0390',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0525',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R08624',
                          'stALE_E_coli_A10_F43_I1_R1',
                          'CDPH-SAL_Salmonella_Typhi_MDL-150',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0474',
                          'CDPH-SAL_Salmonella_Typhi_MDL-163',
                          'Pputida_TALE__HGL_Pputida_132',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_28_13',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_18_35',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0371',
                          'JM-Metabolic__GN02446',
                          'Pputida_PALE__HGL_Pputida_173',
                          'JM-Metabolic__GN04014',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0418',
                          'JM-Metabolic__GN02567', 'P21_E_coli_ELI345',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0524',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0391',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_31_62',
                          'Pputida_JBEI__HGL_Pputida_107_BP6',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0517',
                          'Pputida_TALE__HGL_Pputida_133',
                          'CDPH-SAL_Salmonella_Typhi_MDL-162',
                          'CDPH-SAL_Salmonella_Typhi_MDL-151',
                          'stALE_E_coli_A12_F136_I1_R1',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R08704',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11101',
                          'Pputida_TALE__HGL_Pputida_143', 'P21_E_coli_ELI358',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0405',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0368',
                          'CDPH-SAL_Salmonella_Typhi_MDL-155',
                          'CDPH-SAL_Salmonella_Typhi_MDL-166',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_20_16',
                          'Pputida_TALE__HGL_Pputida_137',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0520',
                          'Pputida_PALE__HGL_Pputida_159',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0395',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_3_23',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0401',
                          'stALE_E_coli_A10_F21_I1_R1',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_7_41',
                          'JM-Metabolic__GN05128',
                          'JBI_KHP_HGL_031_Amitesh_rpoS',
                          'JM-Metabolic__GN0_2175', 'JM-Metabolic__GN0_2290',
                          'JM-Metabolic__GN04540',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0388',
                          'CDPH-SAL_Salmonella_Typhi_MDL-148',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0375',
                          'Pputida_TALE__HGL_Pputida_119',
                          'Deoxyribose_PALE_ALE__MG1655_BOP27_10_13',
                          'JM-Metabolic__GN02769', 'stALE_E_coli_A1_F21_I1_R1',
                          'RMA_KHP_rpoS_Mage_Q97N',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0400',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_31_24',
                          'Pputida_TALE__HGL_Pputida_136',
                          'CDPH-SAL_Salmonella_Typhi_MDL-167',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_3_61',
                          'CDPH-SAL_Salmonella_Typhi_MDL-154',
                          'stALE_E_coli_A2_F21_I1_R1',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0369',
                          'Pputida_PALE__HGL_Pputida_158',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0394',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0521',
                          'JM-Metabolic__GN0_2169',
                          'stALE_E_coli_A4_F21_I1_R2',
                          'Deoxyribose_PALE_ALE__MG1655_BOP27_10_51',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0389',
                          'Pputida_PALE__HGL_Pputida_145',
                          'Pputida_PALE__HGL_Pputida_176',
                          'Pputida_TALE__HGL_Pputida_118',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0374',
                          'CDPH-SAL_Salmonella_Typhi_MDL-149',
                          'JM-Metabolic__GN04488',
                          'Pputida_JBEI__HGL_Pputida_109_BP8',
                          'P21_E_coli_ELI355',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0408',
                          'P21_E_coli_ELI366', 'stALE_E_coli_A9_F21_I1_R1',
                          'Pputida_PALE__HGL_Pputida_150',
                          'JM-Metabolic__GN04682',
                          'Pputida_JBEI__HGL_Pputida_111_M5',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0485',
                          'stALE_E_coli_A11_F43_I1_R1',
                          'Pputida_PALE__HGL_Pputida_163',
                          'stALE_E_coli_A17_F118_I1_R1',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0352',
                          'P21_E_coli_ELI348', 'JM-Metabolic__GN04094',
                          'JM-Metabolic__GN04612', 'stALE_E_coli_A8_F20_I1_R1',
                          'Pputida_TALE__HGL_Pputida_123',
                          'JM-Metabolic__GN02529',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0381',
                          'JBI_KHP_HGL_027', 'Pputida_PALE__HGL_Pputida_162',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0484',
                          'Pputida_PALE__HGL_Pputida_151',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0353',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_26_6',
                          'P21_E_coli_ELI367',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0409',
                          'stALE_E_coli_A13_F42_I1_R1', 'P21_E_coli_ELI354',
                          'Deoxyribose_PALE_ALE__MG1655_BOP27_10_28',
                          'Pputida_TALE__HGL_Pputida_122', 'JBI_KHP_HGL_026',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_24_52',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11153',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0380',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_29_57',
                          'P21_E_coli_ELI349', 'JM-Metabolic__GN0_2215',
                          'stALE_E_coli_A12_F43_I1_R1',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0383',
                          'JM-Metabolic__GN0_2337',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_33_41',
                          'JM-Metabolic__GN02449', 'JBI_KHP_HGL_025',
                          'JM-Metabolic__GN0_2183',
                          'stALE_E_coli_A15_F42_I1_R1',
                          'Pputida_TALE__HGL_Pputida_121',
                          'JM-Metabolic__GN02487',
                          'CDPH-SAL_Salmonella_Typhi_MDL-143',
                          'Pputida_TALE__HGL_Pputida_112',
                          'JM-Metabolic__GN03218',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0417',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_20_43',
                          'JM-Metabolic__GN05002', 'stALE_E_coli_A4_F42_I1_R1',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_1_69',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_29_49',
                          'Pputida_PALE__HGL_Pputida_152',
                          'JM-Metabolic__GN05377',
                          'Pputida_PALE__HGL_Pputida_161',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0518',
                          'P21_E_coli_ELI357', 'RMA_KHP_rpoS_Mage_Q97D',
                          'P21_E_coli_ELI364', 'JM-Metabolic__GN02424',
                          'JM-Metabolic__GN0_2375',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_29_25',
                          'stALE_E_coli_A7_F42_I1_R1', 'JM-Metabolic__GN03132',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_20_71',
                          'Pputida_JBEI__HGL_Pputida_110_M2',
                          'JBI_KHP_HGL_024',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_31_42',
                          'stALE_E_coli_A18_F18_I1_R1',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0382',
                          'JM-Metabolic__GN0_2254',
                          'Pputida_TALE__HGL_Pputida_113',
                          'JM-Metabolic__GN04428',
                          'Pputida_TALE__HGL_Pputida_120',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_32_56',
                          'P21_E_coli_ELI365', 'RMA_KHP_rpoS_Mage_Q97E',
                          'stALE_E_coli_A16_F42_I1_R1',
                          'stALE_E_coli_A6_F43_I1_R1',
                          'Deoxyribose_PALE_ALE__MG1655_BOP27_6_21',
                          'Deoxyribose_PALE_ALE__MG1655_BOP27_4_23',
                          'stALE_E_coli_A15_F117_I1_R1',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0486',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0519',
                          'Pputida_PALE__HGL_Pputida_160',
                          'Pputida_PALE__HGL_Pputida_153',
                          'stALE_E_coli_A4_F21_I1_R1',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11078',
                          'stALE_E_coli_A14_F20_I1_R1',
                          'JM-Metabolic__GN0_2380',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0377',
                          'Pputida_TALE__HGL_Pputida_128',
                          'Pputida_PALE__HGL_Pputida_175',
                          'JM-Metabolic__GN03252',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_26_69',
                          'Pputida_PALE__HGL_Pputida_146',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_19_25',
                          'JM-Metabolic__GN02657',
                          'stALE_E_coli_A15_F21_I1_R1',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0329',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0403',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_19_55',
                          'Pputida_PALE__HGL_Pputida_168',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0522',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0397',
                          'stALE_E_coli_A9_F44_I1_R1',
                          'CDPH-SAL_Salmonella_Typhi_MDL-157',
                          'Pputida_TALE__HGL_Pputida_135',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_22_16',
                          'CDPH-SAL_Salmonella_Typhi_MDL-164',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0473',
                          'Pputida_TALE__HGL_Pputida_129',
                          'Deoxyribose_PALE_ALE__MG1655_BOP27_6_35',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0376',
                          'Pputida_PALE__HGL_Pputida_147',
                          'Pputida_PALE__HGL_Pputida_174',
                          'stALE_E_coli_A3_F40_I1_R1',
                          'JM-Metabolic__GN0_2099',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0396',
                          'stALE_E_coli_A7_F21_I1_R1',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0523',
                          'Pputida_PALE__HGL_Pputida_169',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_3_50',
                          'JM-Metabolic__GN0_2009',
                          'CDPH-SAL_Salmonella_Typhi_MDL-165',
                          'Pputida_TALE__HGL_Pputida_134',
                          'CDPH-SAL_Salmonella_Typhi_MDL-156',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11135',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R10727',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0402',
                          'RMA_KHP_rpoS_Mage_Q97L',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0328',
                          'Pputida_TALE__HGL_Pputida_144',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_17_58',
                          'CDPH-SAL_Salmonella_Typhi_MDL-152',
                          'Pputida_TALE__HGL_Pputida_130',
                          'JM-Metabolic__GN0_2277',
                          'CDPH-SAL_Salmonella_Typhi_MDL-161',
                          'stALE_E_coli_A5_F42_I1_R1',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0392',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_33_50',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_19_63',
                          'P21_E_coli_ELI368', 'stALE_E_coli_A14_F133_I1_R1',
                          'Pputida_TALE__HGL_Pputida_140',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0406',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11102',
                          'Pputida_PALE__HGL_Pputida_170',
                          'JM-Metabolic__GN0_2172',
                          'stALE_E_coli_A14_F42_I1_R1',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0372',
                          'JM-Metabolic__GN02514',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_28_53',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0407',
                          'Pputida_TALE__HGL_Pputida_141', 'P21_E_coli_ELI369',
                          'stALE_E_coli_A18_F130_I1_R1',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_33_21',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11103',
                          'CDPH-SAL_Salmonella_Typhi_MDL-160',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_1_57',
                          'Pputida_TALE__HGL_Pputida_131',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_18_19',
                          'CDPH-SAL_Salmonella_Typhi_MDL-153',
                          'stALE_E_coli_A13_F121_I1_R1',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0393',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_30_7',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0330',
                          'P21_E_coli_ELI347', 'JM-Metabolic__GN05367',
                          'Pputida_PALE__HGL_Pputida_171',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0373',
                          'P21_E_coli_ELI352', 'P21_E_coli_ELI361',
                          'stALE_E_coli_A3_F18_I1_R1',
                          'Pputida_PALE__HGL_Pputida_157',
                          'Pputida_PALE__HGL_Pputida_164',
                          'JM-Metabolic__GN02787',
                          'CDPH-SAL_Salmonella_Typhi_MDL-168',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0355',
                          'Pputida_TALE__HGL_Pputida_139',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0366',
                          'Deoxyribose_PALE_ALE__MG1655_BOP27_4_14',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_17_25',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0421',
                          'JM-Metabolic__GN0_2094',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_24_24',
                          'JM-Metabolic__GN02590',
                          'JBI_KHP_HGL_030_Amitesh_soxR_oxyR',
                          'stALE_E_coli_A5_F21_I1_R1', 'JM-Metabolic__GN04306',
                          'Pputida_TALE__HGL_Pputida_124',
                          'CDPH-SAL_Salmonella_Typhi_MDL-146',
                          'Pputida_TALE__HGL_Pputida_117',
                          'JM-Metabolic__GN04665', 'JM-Metabolic__GN0_2148',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0386',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0483',
                          'Pputida_PALE__HGL_Pputida_165',
                          'Pputida_PALE__HGL_Pputida_156',
                          'JM-Metabolic__GN0_2404',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_21_51',
                          'JM-Metabolic__GN05109', 'JM-Metabolic__GN02748',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_22_28',
                          'stALE_E_coli_A16_F20_I1_R1',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0367',
                          'JM-Metabolic__GN02501',
                          'Pputida_TALE__HGL_Pputida_138',
                          'JBI_KHP_HGL_028_Amitesh_soxR',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0354',
                          'stALE_E_coli_A6_F21_I1_R1',
                          'JM-Metabolic__GN0_2005',
                          'stALE_E_coli_A16_F134_I1_R1', 'P21_E_coli_ELI353',
                          'Pputida_TALE__HGL_Pputida_116',
                          'stALE_E_coli_A17_F21_I1_R1',
                          'CDPH-SAL_Salmonella_Typhi_MDL-147',
                          'JM-Metabolic__GN02766',
                          'Pputida_TALE__HGL_Pputida_125', 'JBI_KHP_HGL_021',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_7_23',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0387',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11154',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0420',
                          'stALE_E_coli_A11_F119_I1_R1',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_26_27',
                          'P21_E_coli_ELI355',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0408',
                          'stALE_E_coli_A9_F21_I1_R1', 'P21_E_coli_ELI366',
                          'JM-Metabolic__GN04488',
                          'Pputida_JBEI__HGL_Pputida_109_BP8',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0352',
                          'stALE_E_coli_A17_F118_I1_R1',
                          'JM-Metabolic__GN04682',
                          'Pputida_PALE__HGL_Pputida_150',
                          'Pputida_JBEI__HGL_Pputida_111_M5',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0485',
                          'stALE_E_coli_A11_F43_I1_R1',
                          'Pputida_PALE__HGL_Pputida_163',
                          'JM-Metabolic__GN04094', 'P21_E_coli_ELI348',
                          'JM-Metabolic__GN04612', 'stALE_E_coli_A8_F20_I1_R1',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0381',
                          'JBI_KHP_HGL_027', 'Pputida_TALE__HGL_Pputida_123',
                          'JM-Metabolic__GN02529',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0353',
                          'Pputida_PALE__HGL_Pputida_162',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0484',
                          'Pputida_PALE__HGL_Pputida_151', 'P21_E_coli_ELI367',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0409',
                          'stALE_E_coli_A13_F42_I1_R1', 'P21_E_coli_ELI354',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_26_6',
                          'JBI_KHP_HGL_026',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_24_52',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0380',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11153',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_29_57',
                          'Deoxyribose_PALE_ALE__MG1655_BOP27_10_28',
                          'Pputida_TALE__HGL_Pputida_122',
                          'JM-Metabolic__GN0_2215', 'P21_E_coli_ELI349',
                          'stALE_E_coli_A12_F43_I1_R1',
                          'Pputida_PALE__HGL_Pputida_159',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0395',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0520',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0368',
                          'CDPH-SAL_Salmonella_Typhi_MDL-155',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_20_16',
                          'CDPH-SAL_Salmonella_Typhi_MDL-166',
                          'Pputida_TALE__HGL_Pputida_137',
                          'stALE_E_coli_A10_F21_I1_R1',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_7_41',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_3_23',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0401',
                          'CDPH-SAL_Salmonella_Typhi_MDL-148',
                          'Pputida_TALE__HGL_Pputida_119',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0375',
                          'Deoxyribose_PALE_ALE__MG1655_BOP27_10_13',
                          'JM-Metabolic__GN02769',
                          'JBI_KHP_HGL_031_Amitesh_rpoS',
                          'JM-Metabolic__GN05128', 'JM-Metabolic__GN0_2290',
                          'JM-Metabolic__GN04540', 'JM-Metabolic__GN0_2175',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0388',
                          'stALE_E_coli_A1_F21_I1_R1',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_31_24',
                          'RMA_KHP_rpoS_Mage_Q97N',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0400',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0521',
                          'Pputida_PALE__HGL_Pputida_158',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0394',
                          'JM-Metabolic__GN0_2169',
                          'Pputida_TALE__HGL_Pputida_136',
                          'CDPH-SAL_Salmonella_Typhi_MDL-167',
                          'CDPH-SAL_Salmonella_Typhi_MDL-154',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_3_61',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0369',
                          'stALE_E_coli_A2_F21_I1_R1',
                          'stALE_E_coli_A4_F21_I1_R2',
                          'Deoxyribose_PALE_ALE__MG1655_BOP27_10_51',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0374',
                          'Pputida_TALE__HGL_Pputida_118',
                          'CDPH-SAL_Salmonella_Typhi_MDL-149',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0389',
                          'Pputida_PALE__HGL_Pputida_145',
                          'Pputida_PALE__HGL_Pputida_176', 'P21_E_coli_ELI344',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0419',
                          'Pputida_PALE__HGL_Pputida_172',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_21_46',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0370',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_22_52',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0404',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_18_59',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_1_24',
                          'Pputida_TALE__HGL_Pputida_142', 'P21_E_coli_ELI359',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_5_46',
                          'JM-Metabolic__GN0_2354',
                          'CDPH-SAL_Salmonella_Typhi_MDL-150',
                          'stALE_E_coli_A10_F43_I1_R1',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R08624',
                          'CDPH-SAL_Salmonella_Typhi_MDL-163',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0474',
                          'Pputida_TALE__HGL_Pputida_132',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0516',
                          'JM-Metabolic__GN0_2317',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_5_36',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0525',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0390',
                          'JM-Metabolic__GN02446',
                          'Pputida_PALE__HGL_Pputida_173',
                          'JM-Metabolic__GN04014',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_18_35',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_28_13',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0371',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0418',
                          'P21_E_coli_ELI345', 'JM-Metabolic__GN02567',
                          'Pputida_TALE__HGL_Pputida_133',
                          'CDPH-SAL_Salmonella_Typhi_MDL-162',
                          'stALE_E_coli_A12_F136_I1_R1',
                          'CDPH-SAL_Salmonella_Typhi_MDL-151',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0391',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0524',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_31_62',
                          'Pputida_JBEI__HGL_Pputida_107_BP6',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0517',
                          'P21_E_coli_ELI358', 'Pputida_TALE__HGL_Pputida_143',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0405',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R08704',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11101',
                          'Pputida_TALE__HGL_Pputida_126',
                          'Pputida_TALE__HGL_Pputida_115',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_28_28',
                          'CDPH-SAL_Salmonella_Typhi_MDL-144',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_5_22',
                          'JM-Metabolic__GN04255',
                          'Pputida_PALE__HGL_Pputida_148',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0384',
                          'stALE_E_coli_A8_F42_I1_R1', 'JBI_KHP_HGL_022',
                          'Deoxyribose_PALE_ALE__MG1655_BOP27_4_48',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0399',
                          'stALE_E_coli_A10_F131_I1_R1',
                          'Pputida_PALE__HGL_Pputida_155',
                          'Pputida_PALE__HGL_Pputida_166',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0357',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_32_20',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_30_22',
                          'CDPH-SAL_Salmonella_Typhi_MDL-159',
                          'JM-Metabolic__GN03409', 'JM-Metabolic__GN02531',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0364',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0327',
                          'P21_E_coli_ELI350', 'JM-Metabolic__GN0_2393',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_17_64',
                          'JBI_KHP_HGL_029_Amitesh_oxyR', 'P21_E_coli_ELI363',
                          'stALE_E_coli_A11_F21_I1_R1',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11044',
                          'stALE_E_coli_A18_F39_I1_R1',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_7_51',
                          'CDPH-SAL_Salmonella_Typhi_MDL-145',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0378',
                          'Pputida_TALE__HGL_Pputida_114',
                          'Pputida_TALE__HGL_Pputida_127', 'JBI_KHP_HGL_023',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_32_6',
                          'stALE_E_coli_A12_F21_I1_R1',
                          'Pputida_PALE__HGL_Pputida_149',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0385',
                          'Pputida_JBEI__HGL_Pputida_108_BP7',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_21_23',
                          'P21_E_coli_ELI362',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_30_60',
                          'P21_E_coli_ELI351',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0326',
                          'JM-Metabolic__GN0_2165',
                          'Pputida_PALE__HGL_Pputida_167',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_24_9',
                          'stALE_E_coli_A13_F20_I1_R1',
                          'JM-Metabolic__GN04563',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0398',
                          'Pputida_PALE__HGL_Pputida_154',
                          'CDPH-SAL_Salmonella_Typhi_MDL-158',
                          'JM-Metabolic__GN0_2007',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0356',
                          'P21_E_coli_ELI352', 'P21_E_coli_ELI361',
                          'stALE_E_coli_A3_F18_I1_R1',
                          'CDPH-SAL_Salmonella_Typhi_MDL-168',
                          'Pputida_TALE__HGL_Pputida_139',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0355',
                          'Deoxyribose_PALE_ALE__MG1655_BOP27_4_14',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0366',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_17_25',
                          'Pputida_PALE__HGL_Pputida_157',
                          'Pputida_PALE__HGL_Pputida_164',
                          'JM-Metabolic__GN02787', 'JM-Metabolic__GN0_2094',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_24_24',
                          'JM-Metabolic__GN02590',
                          'JBI_KHP_HGL_030_Amitesh_soxR_oxyR',
                          'stALE_E_coli_A5_F21_I1_R1',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0421',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0386',
                          'JM-Metabolic__GN0_2148', 'JM-Metabolic__GN04306',
                          'Pputida_TALE__HGL_Pputida_124',
                          'CDPH-SAL_Salmonella_Typhi_MDL-146',
                          'JM-Metabolic__GN04665',
                          'Pputida_TALE__HGL_Pputida_117',
                          'JM-Metabolic__GN02748',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_22_28',
                          'stALE_E_coli_A16_F20_I1_R1',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0367',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0354',
                          'stALE_E_coli_A6_F21_I1_R1',
                          'Pputida_TALE__HGL_Pputida_138',
                          'JM-Metabolic__GN02501',
                          'JBI_KHP_HGL_028_Amitesh_soxR',
                          'JM-Metabolic__GN0_2005',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0483',
                          'Pputida_PALE__HGL_Pputida_165',
                          'JM-Metabolic__GN0_2404',
                          'Pputida_PALE__HGL_Pputida_156',
                          'JM-Metabolic__GN05109',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_21_51',
                          'P21_E_coli_ELI353', 'stALE_E_coli_A16_F134_I1_R1',
                          'JBI_KHP_HGL_021',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_7_23',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11154',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0387',
                          'stALE_E_coli_A17_F21_I1_R1',
                          'Pputida_TALE__HGL_Pputida_116',
                          'CDPH-SAL_Salmonella_Typhi_MDL-147',
                          'JM-Metabolic__GN02766',
                          'Pputida_TALE__HGL_Pputida_125',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_26_27',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0420',
                          'stALE_E_coli_A11_F119_I1_R1',
                          'stALE_E_coli_A5_F42_I1_R1',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0392',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_19_63',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_33_50',
                          'CDPH-SAL_Salmonella_Typhi_MDL-152',
                          'JM-Metabolic__GN0_2277',
                          'Pputida_TALE__HGL_Pputida_130',
                          'CDPH-SAL_Salmonella_Typhi_MDL-161',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11102',
                          'P21_E_coli_ELI368', 'stALE_E_coli_A14_F133_I1_R1',
                          'Pputida_TALE__HGL_Pputida_140',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0406',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0372',
                          'JM-Metabolic__GN02514', 'JM-Metabolic__GN0_2172',
                          'Pputida_PALE__HGL_Pputida_170',
                          'stALE_E_coli_A14_F42_I1_R1',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_28_53',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_33_21',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11103',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0407',
                          'Pputida_TALE__HGL_Pputida_141', 'P21_E_coli_ELI369',
                          'stALE_E_coli_A18_F130_I1_R1',
                          'stALE_E_coli_A13_F121_I1_R1',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0393',
                          'CDPH-SAL_Salmonella_Typhi_MDL-160',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_1_57',
                          'Pputida_TALE__HGL_Pputida_131',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_18_19',
                          'CDPH-SAL_Salmonella_Typhi_MDL-153',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0330',
                          'P21_E_coli_ELI347',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_30_7',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0373',
                          'JM-Metabolic__GN05367',
                          'Pputida_PALE__HGL_Pputida_171',
                          'stALE_E_coli_A4_F21_I1_R1',
                          'JM-Metabolic__GN0_2380',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11078',
                          'stALE_E_coli_A14_F20_I1_R1',
                          'Pputida_PALE__HGL_Pputida_175',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_26_69',
                          'JM-Metabolic__GN03252',
                          'Pputida_PALE__HGL_Pputida_146',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0377',
                          'Pputida_TALE__HGL_Pputida_128',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0329',
                          'stALE_E_coli_A15_F21_I1_R1',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0403',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_19_25',
                          'JM-Metabolic__GN02657',
                          'CDPH-SAL_Salmonella_Typhi_MDL-157',
                          'Pputida_TALE__HGL_Pputida_135',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0473',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_22_16',
                          'CDPH-SAL_Salmonella_Typhi_MDL-164',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_19_55',
                          'Pputida_PALE__HGL_Pputida_168',
                          'stALE_E_coli_A9_F44_I1_R1',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0397',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0522',
                          'Pputida_PALE__HGL_Pputida_147',
                          'Pputida_PALE__HGL_Pputida_174',
                          'Deoxyribose_PALE_ALE__MG1655_BOP27_6_35',
                          'Pputida_TALE__HGL_Pputida_129',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0376',
                          'JM-Metabolic__GN0_2099',
                          'stALE_E_coli_A3_F40_I1_R1',
                          'CDPH-SAL_Salmonella_Typhi_MDL-165',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_3_50',
                          'JM-Metabolic__GN0_2009',
                          'Pputida_TALE__HGL_Pputida_134',
                          'CDPH-SAL_Salmonella_Typhi_MDL-156',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0523',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0396',
                          'stALE_E_coli_A7_F21_I1_R1',
                          'Pputida_PALE__HGL_Pputida_169',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0402',
                          'Pputida_TALE__HGL_Pputida_144',
                          'RMA_KHP_rpoS_Mage_Q97L',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0328',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_17_58',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11135',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R10727',
                          'stALE_E_coli_A15_F42_I1_R1',
                          'JM-Metabolic__GN0_2183',
                          'Pputida_TALE__HGL_Pputida_121',
                          'JM-Metabolic__GN02487',
                          'CDPH-SAL_Salmonella_Typhi_MDL-143',
                          'Pputida_TALE__HGL_Pputida_112',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0383',
                          'JM-Metabolic__GN02449',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_33_41',
                          'JM-Metabolic__GN0_2337', 'JBI_KHP_HGL_025',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0417',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_20_43',
                          'JM-Metabolic__GN05002', 'JM-Metabolic__GN03218',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_29_49',
                          'JM-Metabolic__GN05377',
                          'Pputida_PALE__HGL_Pputida_152',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0518',
                          'Pputida_PALE__HGL_Pputida_161',
                          'stALE_E_coli_A4_F42_I1_R1',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_1_69',
                          'JM-Metabolic__GN02424', 'P21_E_coli_ELI357',
                          'RMA_KHP_rpoS_Mage_Q97D', 'P21_E_coli_ELI364',
                          'JM-Metabolic__GN03132',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_20_71',
                          'Pputida_JBEI__HGL_Pputida_110_M2',
                          'JM-Metabolic__GN0_2375',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_29_25',
                          'stALE_E_coli_A7_F42_I1_R1',
                          'Pputida_TALE__HGL_Pputida_113',
                          'JM-Metabolic__GN0_2254', 'JM-Metabolic__GN04428',
                          'Pputida_TALE__HGL_Pputida_120', 'JBI_KHP_HGL_024',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_31_42',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0382',
                          'stALE_E_coli_A18_F18_I1_R1',
                          'stALE_E_coli_A16_F42_I1_R1',
                          'stALE_E_coli_A6_F43_I1_R1',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_32_56',
                          'P21_E_coli_ELI365', 'RMA_KHP_rpoS_Mage_Q97E',
                          'stALE_E_coli_A15_F117_I1_R1',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0486',
                          'Pputida_PALE__HGL_Pputida_160',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0519',
                          'Pputida_PALE__HGL_Pputida_153',
                          'Deoxyribose_PALE_ALE__MG1655_BOP27_4_23',
                          'Deoxyribose_PALE_ALE__MG1655_BOP27_6_21']

        self.gerwick_ids = ['7A', 'ISB', '5B', '6A', '8A', '3A', 'GFR', '4A',
                            '5B', '6A', '8A', 'ISB', '7A', 'GFR', '3A', '4A']

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

        # the entire list of sample-ids found w/in good-sample-sheet.csv
        self.sample_ids = self.feist_ids + self.gerwick_ids + self.nyu_ids
        package_root = abspath(join(dirname(__file__), '..'))
        # this base path is used extensively throughout the tests.
        self.base_path = partial(join, package_root, 'tests', 'data')
        self.good_output_path = self.base_path('output_dir')
        self.sample_sheet_path = self.base_path('good-sample-sheet.csv')
        self.sample_sheet_w_context_path = self.base_path(
            'good-sample-sheet_w_context.csv')
        self.sheet_w_repl_path = self.base_path('good_sheet_w_replicates.csv')
        self.good_input_path = self.base_path('input_dir')

        # self.good_input_path doesn't need to have anything in it for the
        # purposes of testing, but it does need to exist or else the object
        # will raise an Error.
        makedirs(self.good_input_path, exist_ok=True)

        # because we can't run bcl2fastq/bcl-convert in a unit-test
        # environment, we need to simulate the file/directory structure of
        # the output in order to test the audit() method. ConvertJob assumes
        # the location of the output files to be under the supplied output
        # path/ConvertJob. Hence, our faked_output_path will follow the same
        # convention.
        faked_output_path = join(self.good_output_path, 'ConvertJob')

        # mimic the Data/Fastq/<Project> hierarchy.
        fastq_base = join(faked_output_path, 'Data', 'Fastq')
        makedirs(join(fastq_base, 'Feist_11661'), exist_ok=True)
        makedirs(join(fastq_base, 'Gerwick_6123'), exist_ok=True)
        makedirs(join(fastq_base, 'NYU_BMS_Melanoma_13059'), exist_ok=True)

        # generate filenames and paths for the dummy fastq files.
        file_names = [join(fastq_base,
                           'Feist_11661',
                           '%s_R1.fastq.gz' % x) for x in self.feist_ids]
        file_names += [join(fastq_base,
                            'Feist_11661',
                            '%s_R2.fastq.gz' % x) for x in self.feist_ids]
        file_names += [join(fastq_base,
                            'Gerwick_6123',
                            '%s_R1.fastq.gz' % x) for x in self.gerwick_ids]
        file_names += [join(fastq_base,
                            'Gerwick_6123',
                            '%s_R2.fastq.gz' % x) for x in self.gerwick_ids]
        file_names += [join(fastq_base,
                            'NYU_BMS_Melanoma_13059',
                            '%s_R1.fastq.gz' % x) for
                       x in self.nyu_ids]
        file_names += [join(fastq_base,
                            'NYU_BMS_Melanoma_13059',
                            '%s_R2.fastq.gz' % x) for
                       x in self.nyu_ids]

        # create the dummy fastq files.
        for line in file_names:
            # create a fake forward-read file.
            with open(line, 'w') as f2:
                f2.write("This is a file.")
            # create a fake reverse-read file.
            line = line.replace('_R1.fastq', '_R2.fastq')
            with open(line, 'w') as f2:
                f2.write("This is a file.")

        # convert to Log directory and file created by bcl-convert.
        # this is separate from the standard 'log' directory created
        # by all jobs.
        self.convert_log_path = join(faked_output_path, 'Logs')
        makedirs(self.convert_log_path, exist_ok=True)
        self.convert_log_path = join(self.convert_log_path, 'Errors.log')
        with open(self.convert_log_path, 'w') as f:
            f.write("2024-01-01T12:12:12Z thread 99999 ERROR: Sample Sheet "
                    "Error: in OverrideCycles: Read # 2 specified does not "
                    "add up to the 8 bases expected from RunInfo.xml\n")

    def tearDown(self):
        rmtree(self.good_input_path)
        rmtree(self.good_output_path)

    def test_creation(self):
        run_dir = self.base_path('211021_A00000_0000_SAMPLE')
        inv_input_directory = self.base_path('inv_input_directory')
        qiita_id = 'abcdabcdabcdabcdabcdabcdabcdabcd'

        # ConvertJob should assert due to invalid_input_directory.
        with self.assertRaises(PipelineError) as e:
            ConvertJob(inv_input_directory, self.good_output_path,
                       self.sample_sheet_path, 'qiita', 1, 16, 1440, '10gb',
                       'tests/bin/bcl-convert', [], qiita_id)

        self.assertEqual(str(e.exception),
                         "directory_path '%s' does not exist." %
                         self.base_path('inv_input_directory'))

        job = ConvertJob(run_dir, self.good_output_path,
                         self.sample_sheet_path, 'qiita', 1, 16, 1440, '10gb',
                         'tests/bin/bcl-convert', [], qiita_id)

        job._generate_job_script()

        with open(join(self.good_output_path, 'ConvertJob',
                       'ConvertJob.sh')) as f:
            obs = ''.join(f.readlines())

        # ssp should be just the value of the self.path() partial function by
        # itself. For readability, SCRIPT_EXP addresses the '/' separator.
        # Hence, the trailing '/' is redundant and should be removed here.
        self.assertEqual(obs,
                         SCRIPT_EXP.format(ssp=self.base_path('').rstrip('/'),
                                           gop=self.good_output_path,
                                           run_dir=run_dir))

    def test_error_msg_from_logs(self):
        run_dir = self.base_path('211021_A00000_0000_SAMPLE')
        qiita_id = 'abcdabcdabcdabcdabcdabcdabcdabcd'

        job = ConvertJob(run_dir, self.good_output_path,
                         self.sample_sheet_path, 'qiita', 1, 16, 1440, '10gb',
                         'tests/bin/bcl-convert', [], qiita_id)

        # an internal method to force submit_job() to raise a JobFailedError
        # instead of submitting the job w/sbatch and waiting for a failed
        # job w/squeue.
        self.assertTrue(job._toggle_force_job_fail())

        error_msg = ("This job died.\n2024-01-01T12:12:12Z thread 99999 ERROR:"
                     " Sample Sheet Error: in OverrideCycles: Read # 2 "
                     "specified does not add up to the 8 bases expected from"
                     " RunInfo.xml")

        with self.assertRaisesRegex(JobFailedError, error_msg):
            job.run()

    def test_audit(self):
        # the faked output should be in self.good_output_path/ConvertJob.
        # ConvertJob already takes into account 'ConvertJob' and so the
        # correct path to the faked output is self.good_output_path, rather
        # than faked_output_path.
        job = ConvertJob(self.good_input_path, self.good_output_path,
                         self.sample_sheet_path, 'qiita', 1, 16, 1440, '10gb',
                         'tests/bin/bcl-convert', [], 'some_qiita_id')

        obs = job.audit(self.sample_ids)
        # there shouldn't be any missing samples.
        self.assertEqual(obs, [])

        # these fake sample-ids should be returned by audit() as missing.
        obs = job.audit(self.sample_ids + ['not-a-sample', 'BLANK1'])
        self.assertListEqual(obs, ['BLANK1', 'not-a-sample'])

    def test_parse_sample_sheet(self):
        js_path = self.base_path('sample-convertjob.sh')
        obs = ConvertJob.parse_job_script(js_path)

        exp = {
            'run_directory': ('/sequencing/igm_runs/210820_A00953_0380_'
                              'BHJ53TDSX2'),
            'sample_sheet_path': ('/qmounts/qiita_data/working_dir/'
                                  '3f6b1fe3-1f5d-4fec-af47-31a2e35fef91/'
                                  '2024-02-08_U19_Wisconsin_15445_reruns'
                                  '_NovaSeq_nonNA.csv')
        }

        self.assertDictEqual(obs, exp)

    def test_copy_sequences_bad_parameters(self):
        run_dir = self.base_path('211021_A00000_0000_SAMPLE')
        qiita_id = 'abcdabcdabcdabcdabcdabcdabcdabcd'

        job = ConvertJob(run_dir, self.good_output_path,
                         self.sample_sheet_path, 'qiita', 1, 16, 1440, '10gb',
                         'tests/bin/bcl-convert', [], qiita_id)

        # instead of calling run() and faking an entire ConvertJob run,
        # manually call _get_sample_sheet_info(), which is typically called
        # once a job has completed, to gather the metadata needed to properly
        # run copy_sequences() method.

        job._get_sample_sheet_info()

        sample_name = 'CDPH-SAL.Salmonella.Typhi.MDL-154'
        source_project = 'Feist_11661'
        other_projects = ['NYU_BMS_Melanoma_13059', 'Gerwick_6123']
        dest_project = other_projects[0]
        not_source_project = other_projects[1]
        not_a_sample_name = 'NOT_A_SAMPLE_NAME'
        not_a_project = 'NOT_A_PROJECT'

        err_msg = ("'NOT_A_SAMPLE_NAME' did not match any 'sample_name' values"
                   " in project 'Feist_11661'.")
        with self.assertRaisesRegex(ValueError, err_msg):
            job.copy_sequences(not_a_sample_name,
                               source_project,
                               dest_project)

        err_msg = ("'CDPH-SAL.Salmonella.Typhi.MDL-154' did not match any "
                   "'sample_name' values in project 'Gerwick_6123'.")
        with self.assertRaisesRegex(ValueError, err_msg):
            job.copy_sequences(sample_name,
                               not_source_project,
                               dest_project)

        with self.assertRaisesRegex(ValueError, "'NOT_A_PROJECT' is not "
                                                "defined in sample sheet"):
            job.copy_sequences(sample_name,
                               not_a_project,
                               dest_project)

        with self.assertRaisesRegex(ValueError, "'NOT_A_PROJECT' is not "
                                                "defined in sample sheet"):
            job.copy_sequences(sample_name,
                               source_project,
                               not_a_project)

        with self.assertRaisesRegex(ValueError, "source 'Feist_11661' and "
                                                "destination 'Feist_11661' "
                                                "projects are the same"):
            job.copy_sequences(sample_name,
                               source_project,
                               source_project)

    def test_copy_sequences_success(self):
        run_dir = self.base_path('211021_A00000_0000_SAMPLE')
        qiita_id = 'abcdabcdabcdabcdabcdabcdabcdabcd'

        job = ConvertJob(run_dir, self.good_output_path,
                         self.sample_sheet_path, 'qiita', 1, 16, 1440, '10gb',
                         'tests/bin/bcl-convert', [], qiita_id)

        sample_name = 'CDPH-SAL.Salmonella.Typhi.MDL-154'
        source_project = 'Feist_11661'
        dest_project = 'NYU_BMS_Melanoma_13059'
        projects = ['NYU_BMS_Melanoma_13059', 'Gerwick_6123', 'Feist_11661']

        # since we can't perform a real run, let's manually create a fake
        # fastq file and project directories in the 'output_dir' directory and
        # manually call job._get_sample_sheet_info() to obtain all of the
        # metadata needed to copy a sequence from one project into another.

        for some_project in projects:
            # fake the fastq file directories in ConvertJob, one for each
            # project defined in the sample-sheet.
            makedirs(join(self.good_output_path, 'ConvertJob', some_project))

        # fake a fastq file in the 'Feist_11661' directory for the purposes of
        # copying it into the 'NYU_BMS_Melanoma_13059' project.
        with open(join(self.good_output_path, 'ConvertJob', source_project,
                       'CDPH-SAL_Salmonella_Typhi_MDL-154_S1_L001_R1_001.'
                       'fastq.gz'), 'w') as f:
            f.write("Hello World!\n")

        # manually call the functionality that reads the sample-sheet and
        # attempts to associate samples with the fastq files generated by
        # bcl-convert.
        job._get_sample_sheet_info()

        # copy all fastq files associated w/
        # 'CDPH-SAL.Salmonella.Typhi.MDL-154' from 'Feist_11661' to
        # 'NYU_BMS_Melanoma_13059'. 'Gerwick_6123' should remain empty; the
        # code shouldn't copy anything into that project.
        job.copy_sequences(sample_name, source_project, dest_project)

        sample_info = job.info[source_project]['samples'][sample_name]

        # get the path for the source fastq file we created above and swap out
        # the project-level directory names to confirm and deny the existence
        # of the fastq file in other locations.
        source_file = sample_info['matching_files'][0]

        # file should have been copied here.
        dst_file = source_file.replace('Feist_11661', 'NYU_BMS_Melanoma_13059')
        self.assertTrue(exists(dst_file))

        # file should not have been copied here.
        dst_file = source_file.replace('Feist_11661', 'Gerwick_6123')
        self.assertFalse(exists(dst_file))

    def test_copy_sequences_success_w_replicates(self):
        # perform a similar test to the one above, but w/replicate samples.

        run_dir = self.base_path('211021_A00000_0000_SAMPLE')
        qiita_id = 'abcdabcdabcdabcdabcdabcdabcdabcd'

        job = ConvertJob(run_dir, self.good_output_path,
                         self.sheet_w_repl_path, 'qiita', 1, 16, 1440, '10gb',
                         'tests/bin/bcl-convert', [], qiita_id)

        sample_name = 'RMA.KHP.rpoS.Mage.Q97D'
        source_project = 'Feist_11661'
        dest_project = 'NYU_BMS_Melanoma_13059'
        projects = ['NYU_BMS_Melanoma_13059', 'Feist_11661']

        for some_project in projects:
            makedirs(
                join(self.good_output_path, 'ConvertJob', some_project))

        # fake a fastq file in the 'Feist_11661' directory for the purposes of
        # copying it into the 'NYU_BMS_Melanoma_13059' project.

        # instead of faking just a single fastq file, fake an R1, R2 and I1
        # fastq file for all three replicates of 'RMA.KHP.rpoS.Mage.Q97D'.

        fastqs = ['RMA_KHP_rpoS_Mage_Q97D_A5_S1_L001_R1_001.fastq.gz',
                  'RMA_KHP_rpoS_Mage_Q97D_A5_S1_L001_I1_001.fastq.gz',
                  'RMA_KHP_rpoS_Mage_Q97D_A5_S1_L001_R2_001.fastq.gz',

                  'RMA_KHP_rpoS_Mage_Q97D_A6_S1_L001_R1_001.fastq.gz',
                  'RMA_KHP_rpoS_Mage_Q97D_A6_S1_L001_I1_001.fastq.gz',
                  'RMA_KHP_rpoS_Mage_Q97D_A6_S1_L001_R2_001.fastq.gz',

                  'RMA_KHP_rpoS_Mage_Q97D_B6_S1_L001_R1_001.fastq.gz',
                  'RMA_KHP_rpoS_Mage_Q97D_B6_S1_L001_I1_001.fastq.gz',
                  'RMA_KHP_rpoS_Mage_Q97D_B6_S1_L001_R2_001.fastq.gz']

        for fastq in fastqs:
            with open(join(self.good_output_path, 'ConvertJob', source_project,
                           fastq), 'w') as f:
                f.write("Hello World!\n")

        job._get_sample_sheet_info()

        job.copy_sequences(sample_name, source_project, dest_project)

        files_to_match = []

        for smpl in job.info[source_project]['samples']:
            smpl_info = job.info[source_project]['samples'][smpl]
            if smpl_info['orig_name'] == 'RMA.KHP.rpoS.Mage.Q97D':
                files_to_match += smpl_info['matching_files']

        files_to_match = [fp.replace('Feist_11661',
                                     'NYU_BMS_Melanoma_13059')
                          for fp in files_to_match]

        for dst_file in files_to_match:
            self.assertTrue(exists(dst_file))

    def test_copy_controls_between_projects(self):
        # There is a lot of heavy-weight set-up for this test, so I am not
        # testing copying more than one control.  However, this is itself a
        # good check: there are a lot of samples named "BLANK.<etc>" in the
        # sample sheet, but only ONE (which is NOT named "BLANK.<etc>") is
        # listed in the SampleContext section. The fact that the code
        # copies just that one shows it is reading what it should.

        run_dir = self.base_path('211021_A00000_0000_SAMPLE')
        qiita_id = 'abcdabcdabcdabcdabcdabcdabcdabcd'

        job = ConvertJob(
            run_dir, self.good_output_path, self.sample_sheet_w_context_path,
            'qiita', 1, 16, 1440, '10gb', 'tests/bin/bcl-convert', [],
            qiita_id)

        sample_name = 'CDPH-SAL.Salmonella.Typhi.MDL-154'
        source_project = 'Feist_11661'
        projects = ['NYU_BMS_Melanoma_13059', 'Gerwick_6123', 'Feist_11661']

        # since we can't perform a real run, let's manually create a fake
        # fastq file and project directories in the 'output_dir' directory and
        # manually call job._get_sample_sheet_info() to obtain all of the
        # metadata needed to copy a sequence from one project into another.

        for some_project in projects:
            # fake the fastq file directories in ConvertJob, one for each
            # project defined in the sample-sheet.
            makedirs(join(self.good_output_path, 'ConvertJob', some_project),
                     exist_ok=True)

        # fake a fastq file in the 'Feist_11661' directory for the purposes of
        # copying it into the 'NYU_BMS_Melanoma_13059' project.
        with open(join(self.good_output_path, 'ConvertJob', source_project,
                       'CDPH-SAL_Salmonella_Typhi_MDL-154_S1_L001_R1_001.'
                       'fastq.gz'), 'w') as f:
            f.write("Hello World!\n")

        job.copy_controls_between_projects()

        sample_info = job.info[source_project]['samples'][sample_name]

        # get the path for the source fastq file we created above and swap out
        # the project-level directory names to confirm and deny the existence
        # of the fastq file in other locations.
        source_file = sample_info['matching_files'][0]

        # file should have been copied here.
        dst_file = source_file.replace('Feist_11661', 'NYU_BMS_Melanoma_13059')
        self.assertTrue(exists(dst_file))

        # file should not have been copied here.
        dst_file = source_file.replace('Feist_11661', 'Gerwick_6123')
        self.assertFalse(exists(dst_file))


SCRIPT_EXP = ''.join([
    '#!/bin/bash\n',
    '#SBATCH --job-name abcdabcdabcdabcdabcdabcdabcdabcd_ConvertJob\n',
    '#SBATCH -p qiita\n',
    '#SBATCH -N 1\n',
    '#SBATCH -n 16\n',
    '#SBATCH --time 1440\n',
    '#SBATCH --mail-type=ALL\n',
    '#SBATCH --mail-user qiita.help@gmail.com\n',
    '#SBATCH --mem-per-cpu 10gb\n',
    'set -x\n',
    'date\n',
    'hostname\n',
    'cd {run_dir}\n',
    'tests/bin/bcl-convert --sample-sheet "{ssp}/good-sample-sheet.csv" '
    '--output-directory {gop}/ConvertJob '
    '--bcl-input-directory . '
    '--bcl-num-decompression-threads 16 --bcl-num-conversion-threads 16 '
    '--bcl-num-compression-threads 16 --bcl-num-parallel-tiles 16 '
    '--bcl-sampleproject-subdirectories true --force\n'])


if __name__ == '__main__':
    unittest.main()
