#!/bin/bash -l
#SBATCH -J seq_counts
#SBATCH --time 1440
#SBATCH --mem 8G
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -p qiita
#SBATCH --array=1-8

#SBATCH --output tests/2caa8226-cf69-45a3-bd40-1e90ec3d18d0/SeqCountsJob/logs/%x_%A_%a.out
#SBATCH --error tests/2caa8226-cf69-45a3-bd40-1e90ec3d18d0/SeqCountsJob/logs/%x_%A_%a.err

set -x
set -e
set -o pipefail
mkdir -p tests/2caa8226-cf69-45a3-bd40-1e90ec3d18d0/SeqCountsJob/logs

files=($(cat tests/data/files_to_count.txt))
my_file=${files[$((${SLURM_ARRAY_TASK_ID} - 1))]}

echo "${my_file}"

conda activate qp-knight-lab-processing-2022.03

seqtk size ${my_file}
