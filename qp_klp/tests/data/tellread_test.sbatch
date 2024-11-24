#!/bin/bash -l
#SBATCH -J tellread
#SBATCH -p qiita
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem 16G
#SBATCH --time 1440

#SBATCH --output /home/runner/qp-knight-lab-processing/qp_klp/tests/data/077c4da8-74eb-4184-8860-0207f53623be/TellReadJob/logs/tellread_%x-%A.out
#SBATCH --error /home/runner/qp-knight-lab-processing/qp_klp/tests/data/077c4da8-74eb-4184-8860-0207f53623be/TellReadJob/logs/tellread_%x-%A.err

set -x

module load singularity_3.6.4
/home/qiita_test/qiita-spots/ADDITIONAL_REQUIREMENTS/tellread-release-novaseqX/run_tellread_sing.sh \
    -i qp_klp/tests/data/211021_A00000_0000_SAMPLE \
    -o /home/runner/qp-knight-lab-processing/qp_klp/tests/data/077c4da8-74eb-4184-8860-0207f53623be/TellReadJob \
    -s $(echo C591,C550,C503,C566,C556,C578,C592 | tr -d '"') \
    -g $(echo NONE,NONE,NONE,NONE,NONE,NONE,NONE | tr -d '"') \
    -j ${SLURM_JOB_CPUS_PER_NODE}  \
    -l s_4

# get the timestamp for the most recently changed file in directory '.'

# hard-limit for wait time set to ~ 8 hours.
# (4 checks per hour, for 8 hours equals 32 iterations)
for i in $(seq 1 32);
do
    before="$(find /home/runner/qp-knight-lab-processing/qp_klp/tests/data/077c4da8-74eb-4184-8860-0207f53623be/TellReadJob/Full -type f -printf '%T@\n' | sort -n | tail -1)"
    # assume TellReadJob is finished if ctime hasn't changed in 15 minutes
    # for any fastq file in the directory.
    sleep 900
    after="$(find /home/runner/qp-knight-lab-processing/qp_klp/tests/data/077c4da8-74eb-4184-8860-0207f53623be/TellReadJob/Full -type f -printf '%T@\n' | sort -n | tail -1)"

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
