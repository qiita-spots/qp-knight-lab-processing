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

export TMPDIR=/Users/ccowart/SAMPLE_SHEET_MADNESSS/one/qp_klp/tests/data/077c4da8-74eb-4184-8860-0207f53623be/TellReadJob/output/tmp1
mkdir -p ${TMPDIR}
export TMPDIR=$(mktemp -d)

mkdir -p /Users/ccowart/SAMPLE_SHEET_MADNESSS/one/qp_klp/tests/data/077c4da8-74eb-4184-8860-0207f53623be/TellReadJob/output

module load foo_module
/my/sing/script/path \
    -i qp_klp/tests/data/211021_A00000_0000_SAMPLE \
    -o /Users/ccowart/SAMPLE_SHEET_MADNESSS/one/qp_klp/tests/data/077c4da8-74eb-4184-8860-0207f53623be/TellReadJob/output \
    -s $(echo C591,C550,C503,C566,C556,C578,C592 | tr -d '"') \
    -g $(echo NONE,NONE,NONE,NONE,NONE,NONE,NONE | tr -d '"') \
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