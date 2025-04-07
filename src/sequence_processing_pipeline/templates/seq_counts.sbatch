#!/bin/bash -l
#SBATCH -J {{job_name}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -N {{node_count}}
#SBATCH -c {{cores_per_task}}
#SBATCH -p {{queue_name}}
#SBATCH --array=1-{{file_count}}

#SBATCH --output {{output_path}}/logs/%x_%A_%a.out
#SBATCH --error {{output_path}}/logs/%x_%A_%a.err

set -x
set -e
set -o pipefail

mkdir -p {{output_path}}/logs

files=($(cat {{files_to_count_path}}))
my_file=${files[$((${SLURM_ARRAY_TASK_ID} - 1))]}

echo "${my_file}"

conda activate qp-knight-lab-processing-2022.03

seqtk size ${my_file}
