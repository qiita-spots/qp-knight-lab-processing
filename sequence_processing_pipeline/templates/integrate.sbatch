#!/bin/bash -l
#SBATCH -J {{job_name}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -N {{node_count}}
#SBATCH -c {{cores_per_task}}
#SBATCH -p {{queue_name}}
#SBATCH --array=1-{{barcode_id_count}}

#SBATCH --output {{output_dir}}/logs/integrate_%x_%A_%a.out
#SBATCH --error {{output_dir}}/logs/integrate_%x_%A_%a.err

set -x
set -e

samples=($(cat {{output_dir}}/sample_index_list.txt | cut -f 2))
sample=${samples[$((${SLURM_ARRAY_TASK_ID} - 1))]} 

export TMPDIR={{tmp_dir}}

# get list of samples and determine which sample this array instance will work
# on.
samples=($(cat {{output_dir}}/sample_index_list.txt | cut -f 2))
sample=${samples[$((${SLURM_ARRAY_TASK_ID} - 1))]}

echo "Processing sample ${sample}..."

# make temp directory
export TMPDIR={{tmp_dir}}
mkdir -p $TMPDIR


# TODO: All three input files must be non-zero in length.
# If possible, do this check as part of normal FSR operation.
# Previously this was done right here BEFORE integrating, rather
# than after.

# NB: non-zero file-length check removed for now. This should be performed
# by FSR after processing is done.
# TODO: Make sure raw_fastq_dir is TellReadJob/Full
r1_in={{raw_fastq_dir}}/TellReadJob_R1_${sample}.fastq.gz.corrected.err_barcode_removed.fastq
r2_in={{raw_fastq_dir}}/TellReadJob_R2_${sample}.fastq.gz.corrected.err_barcode_removed.fastq
i1_in={{raw_fastq_dir}}/TellReadJob_I1_${sample}.fastq.gz.corrected.err_barcode_removed.fastq

# create output directory
mkdir -p {{output_dir}}/integrated

# generate output file names
r1_out={{output_dir}}/integrated/${sample}.R1.fastq.gz
r2_out={{output_dir}}/integrated/${sample}.R2.fastq.gz
i1_out={{output_dir}}/integrated/${sample}.I1.fastq.gz

# generate 'integrated' I1 fastq.gz file. We do this as part of each array so
# they're done in parallel.
gzip -c ${i1_in} > ${i1_out}

# generate integrated R1 and R2 fastq.gz files.
conda activate qp-knight-lab-processing-2022.03

python {{integrate_script_path}} integrate \
--no-sort \
--r1-in ${r1_in} \
--r2-in ${r2_in} \
--i1-in ${i1_in} \
--r1-out ${r1_out} \
--r2-out ${r2_out} \
--threads {{cores_per_task}}
