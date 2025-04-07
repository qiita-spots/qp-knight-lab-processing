#!/bin/bash -l
#SBATCH -J integrate
#SBATCH --time 96:00:00
#SBATCH --mem 16G
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -p qiita
#SBATCH --array=1-96

#SBATCH --output sequence_processing_pipeline/tests/2caa8226-cf69-45a3-bd40-1e90ec3d18d0/TRIntegrateJob/logs/integrate_%x_%A_%a.out
#SBATCH --error sequence_processing_pipeline/tests/2caa8226-cf69-45a3-bd40-1e90ec3d18d0/TRIntegrateJob/logs/integrate_%x_%A_%a.err

set -x
set -e

samples=($(cat sequence_processing_pipeline/tests/2caa8226-cf69-45a3-bd40-1e90ec3d18d0/TRIntegrateJob/sample_index_list.txt | cut -f 2))
sample=${samples[$((${SLURM_ARRAY_TASK_ID} - 1))]} 

export TMPDIR=sequence_processing_pipeline/tests/2caa8226-cf69-45a3-bd40-1e90ec3d18d0/TRIntegrateJob/tmp

# get list of samples and determine which sample this array instance will work
# on.
samples=($(cat sequence_processing_pipeline/tests/2caa8226-cf69-45a3-bd40-1e90ec3d18d0/TRIntegrateJob/sample_index_list.txt | cut -f 2))
sample=${samples[$((${SLURM_ARRAY_TASK_ID} - 1))]}

echo "Processing sample ${sample}..."

# make temp directory
export TMPDIR=sequence_processing_pipeline/tests/2caa8226-cf69-45a3-bd40-1e90ec3d18d0/TRIntegrateJob/tmp
mkdir -p $TMPDIR


# TODO: All three input files must be non-zero in length.
# If possible, do this check as part of normal FSR operation.
# Previously this was done right here BEFORE integrating, rather
# than after.

# NB: non-zero file-length check removed for now. This should be performed
# by FSR after processing is done.
# TODO: Make sure raw_fastq_dir is TellReadJob/Full
r1_in=sequence_processing_pipeline/tests/2caa8226-cf69-45a3-bd40-1e90ec3d18d0/TellReadJob/Full/TellReadJob_R1_${sample}.fastq.gz.corrected.err_barcode_removed.fastq
r2_in=sequence_processing_pipeline/tests/2caa8226-cf69-45a3-bd40-1e90ec3d18d0/TellReadJob/Full/TellReadJob_R2_${sample}.fastq.gz.corrected.err_barcode_removed.fastq
i1_in=sequence_processing_pipeline/tests/2caa8226-cf69-45a3-bd40-1e90ec3d18d0/TellReadJob/Full/TellReadJob_I1_${sample}.fastq.gz.corrected.err_barcode_removed.fastq

# create output directory
mkdir -p sequence_processing_pipeline/tests/2caa8226-cf69-45a3-bd40-1e90ec3d18d0/TRIntegrateJob/integrated

# generate output file names
r1_out=sequence_processing_pipeline/tests/2caa8226-cf69-45a3-bd40-1e90ec3d18d0/TRIntegrateJob/integrated/${sample}.R1.fastq.gz
r2_out=sequence_processing_pipeline/tests/2caa8226-cf69-45a3-bd40-1e90ec3d18d0/TRIntegrateJob/integrated/${sample}.R2.fastq.gz
i1_out=sequence_processing_pipeline/tests/2caa8226-cf69-45a3-bd40-1e90ec3d18d0/TRIntegrateJob/integrated/${sample}.I1.fastq.gz

# generate 'integrated' I1 fastq.gz file. We do this as part of each array so
# they're done in parallel.
gzip -c ${i1_in} > ${i1_out}

# generate integrated R1 and R2 fastq.gz files.
conda activate qp-knight-lab-processing-2022.03

python sequence_processing_pipeline/contrib/integrate-indices-np.py integrate \
--no-sort \
--r1-in ${r1_in} \
--r2-in ${r2_in} \
--i1-in ${i1_in} \
--r1-out ${r1_out} \
--r2-out ${r2_out} \
--threads 4