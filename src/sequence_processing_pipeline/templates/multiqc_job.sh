#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -p {{queue_name}}
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
# fastqc/multiqc use the same mem value in the templates (2gb); however multiqc requires ~12.
#SBATCH --mem 1{{mem_in_gb}}G
#SBATCH --array {{array_params}}
set -x
set +e
set -o pipefail
date
hostname
echo ${SLURM_JOBID} ${SLURM_ARRAY_TASK_ID}
cd {{output_path}}
{% if modules_to_load is defined %}
    module load {{modules_to_load}}
{% endif %}
step=${SLURM_ARRAY_TASK_ID}
cmd0=$(head -n $step {{array_details}} | tail -n 1)
eval $cmd0
if [ $? -eq 1 ]; then
    echo "multiqc failed."
    exit 1
fi
echo "Cmd Completed: $cmd0" > logs/MultiQCJob_$step.completed
