#!/bin/bash -l
#SBATCH -J tellread
#SBATCH -p qiita
#SBATCH -N 1
#SBATCH -c 999
#SBATCH --mem 10G
#SBATCH --time 216

#SBATCH --output tellread_%x-%A.out
#SBATCH --error tellread_%x-%A.err

set -x

export TMPDIR=/my/tmp1/path
mkdir -p ${TMPDIR}
export TMPDIR=$(mktemp -d)

mkdir -p /Users/ccowart/SAMPLE_SHEET_MADNESSS/one/qp_klp/tests/data/077c4da8-74eb-4184-8860-0207f53623be/TellReadJob/output

module load foo_module
/my/sing/script/path \
    -i qp_klp/tests/data/211021_A00000_0000_SAMPLE \
    -o /Users/ccowart/SAMPLE_SHEET_MADNESSS/one/qp_klp/tests/data/077c4da8-74eb-4184-8860-0207f53623be/TellReadJob/output \
    -s $(echo LS_8_10_2013_SRE,LS_12_17_2014_SRE,LS_4_4_2015_SRE,LS_2_23_2015_SRE,LS_9_28_2014_SRE,LS_12_14_2013_SRE,LS_4_7_2013_SRE,LS_7_14_2013_SRE,LS_10_27_2013_SRE,LS_1_19_2014_SRE,LS_9_3_2013_SRE,LS_2_25_2013_SRE,LS_7_26_2015_SRE,LS_2_17_2014_SRE,LS_6_29_2015_SRE,LS_3_24_2015_SRE,LS_1_6_2015_SRE,T_LS_7_15_15B_SRE,LS_6_9_2013_SRE,Person_A_SRE,LS_8_22_2014_R2_SRE,Person_B_SRE,LS_8_22_2014_R1_SRE,Person_C_SRE,LS_12_28_2011_SRE,Person_D_SRE,LS_5_4_2014_SRE,45208_1_1,LS_11_6_2012_SRE,45248_2_2,LS_4_3_2012_SRE,45261_2_1,45272_11_2,T_LS_7_12_15A,45316_8_1,T_LS_7_8_15A,45327_7_2,LS_8_10_2013,45272_1_swab_2,LS_6_29_2015,45326_1_swab_2,LS_3_8_2015,T_LS_7_19_15A,LS_4_29_2013,T_LS_7_15_15B,LS_11_16_2014,T_LS_7_19_15B,LS_1_19_2014,LS_3_24_2015,LS_2_8_2013,LS_11_10_2013,Marine_Sediment_0_2cm_R1,LS_3_23_2014,Marine_Sediment_5_7cm_R1,LS_1_14_2015,Marine_Sediment_10_12cm_R2,LS_8_25_2014,Marine_Sediment_15_17cm_R1,LS_1_26_2013,Marine_Sediment_20_22cm_R1,LS_6_16_2014,Marine_Sediment_25_27cm_R2,LS_7_27_2014,Marine_Sediment_30_32cm_R3,Person_A_R3,Soil_SynCom_T4_2_Tube5,Person_B_R2,A21,Person_C_R4,A23,Person_D_R2,A27,Soil_SynCom_T1_2_Tube1,A30,Soil_SynCom_T2_2_Tube2,A31,Soil_SynCom_T3_2_Tube3,S1_T1_A,Soil_SynCom_T4_1_Tube4,S2_T1_B_A,S2_T1_01BH1_Y_A,S1_T1_1CIM_A,S2_MT1_1HBI_Y_A,S1_M1_B_1CIM_A,S1_T1_B_LBM_A,BLANK_K15_cancer_patient,S2_MT1_LBM_A,BLANK_M15_cancer_patient,S2_T1_A,BLANK_O15_cancer_patient,1CIM_M_CNTL_A,BLANK_A17_cancer_patient,1CIM_G_CNTL_A,BLANK_C17_cancer_patient,GC_1HCOM_A,BLANK_E17_cancer_patient | tr -d '"') \
    -g $(echo NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE | tr -d '"') \
    -j ${SLURM_JOB_CPUS_PER_NODE}  \
    -l s_4

if [[ -d /Users/ccowart/SAMPLE_SHEET_MADNESSS/one/qp_klp/tests/data/077c4da8-74eb-4184-8860-0207f53623be/TellReadJob/output/Full ]]; then
    echo "Run appears successful"
elif [[ -d /Users/ccowart/SAMPLE_SHEET_MADNESSS/one/qp_klp/tests/data/077c4da8-74eb-4184-8860-0207f53623be/TellReadJob/output/1_demult/Full ]]; then
    echo "Run appears unsuccessful but has output"
    exit 1
else
    echo "Run appears unsuccessful"
    exit 1
fi