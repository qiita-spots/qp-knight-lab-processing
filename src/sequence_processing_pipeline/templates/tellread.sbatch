#!/bin/bash -l
#SBATCH -J {{job_name}}
#SBATCH -p {{queue_name}}
#SBATCH -N {{node_count}}
#SBATCH -c {{cores_per_task}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH --time {{wall_time_limit}}

#SBATCH --output {{output}}/logs/tellread_%x-%A.out
#SBATCH --error {{output}}/logs/tellread_%x-%A.err

set -x

module load {{modules_to_load}}
{{sing_script_path}} \
    -i {{rundir_path}} \
    -o {{output}} \
    -s $(echo {{samples}} | tr -d '"') \
    -g $(echo {{refs}} | tr -d '"') \
    -j ${SLURM_JOB_CPUS_PER_NODE} {{extra}} \
    -l {{lane}}

# get the timestamp for the most recently changed file in directory '.'

# hard-limit for wait time set to ~ 8 hours.
# (4 checks per hour, for 8 hours equals 32 iterations)
for i in $(seq 1 32);
do
    before="$(find {{output}}/Full -type f -printf '%T@\n' | sort -n | tail -1)"
    # assume TellReadJob is finished if ctime hasn't changed in 15 minutes
    # for any fastq file in the directory.
    sleep 900
    after="$(find {{output}}/Full -type f -printf '%T@\n' | sort -n | tail -1)"

    echo "$before   $after"

    if [[ "$before" == "$after" ]]; then
        echo "DONE"
        exit 0
    else
        echo "NOT DONE"
    fi
done

# if we've reached this point then we've exceeded our hard-limit for waiting.
# return w/an error.
exit 1

