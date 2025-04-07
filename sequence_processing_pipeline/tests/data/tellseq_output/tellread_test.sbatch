#!/bin/bash -l
#SBATCH -J tellread
#SBATCH -p qiita
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem 16G
#SBATCH --time 96:00:00

#SBATCH --output sequence_processing_pipeline/tests/2caa8226-cf69-45a3-bd40-1e90ec3d18d0/TellReadJob/logs/tellread_%x-%A.out
#SBATCH --error sequence_processing_pipeline/tests/2caa8226-cf69-45a3-bd40-1e90ec3d18d0/TellReadJob/logs/tellread_%x-%A.err

set -x

module load singularity_3.6.4
$HOME/qiita-spots/tellread-release-novaseqX/run_tellread_sing.sh \
    -i sequence_processing_pipeline/tests/data/sample_run_directories/150629_SN1001_0511_AH5L7GBCXX \
    -o sequence_processing_pipeline/tests/2caa8226-cf69-45a3-bd40-1e90ec3d18d0/TellReadJob \
    -s $(echo C501,C509,C502,C510,C503,C511,C504,C512,C505,C513,C506,C514,C507,C515,C508,C516,C517,C525,C518,C526,C519,C527,C520,C528,C521,C529,C522,C530,C523,C531,C524,C532,C533,C541,C534,C542,C535,C543,C536,C544,C537,C545,C538,C546,C539,C547,C540,C548,C549,C557,C550,C558,C551,C559,C552,C560,C553,C561,C554,C562,C555,C563,C556,C564,C565,C573,C566,C574,C567,C575,C568,C576,C569,C577,C570,C578,C571,C579,C572,C580,C581,C589,C582,C590,C583,C591,C584,C592,C585,C593,C586,C594,C587,C595,C588,C596 | tr -d '"') \
    -g $(echo NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE | tr -d '"') \
    -j ${SLURM_JOB_CPUS_PER_NODE}  \
    -l s_4

# get the timestamp for the most recently changed file in directory '.'

# hard-limit for wait time set to ~ 8 hours.
# (4 checks per hour, for 8 hours equals 32 iterations)
for i in $(seq 1 32);
do
    before="$(find sequence_processing_pipeline/tests/2caa8226-cf69-45a3-bd40-1e90ec3d18d0/TellReadJob/Full -type f -printf '%T@\n' | sort -n | tail -1)"
    # assume TellReadJob is finished if ctime hasn't changed in 15 minutes
    # for any fastq file in the directory.
    sleep 900
    after="$(find sequence_processing_pipeline/tests/2caa8226-cf69-45a3-bd40-1e90ec3d18d0/TellReadJob/Full -type f -printf '%T@\n' | sort -n | tail -1)"

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
