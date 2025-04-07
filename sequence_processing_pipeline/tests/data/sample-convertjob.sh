#!/bin/bash
#SBATCH --job-name 3f6b1fe3-1f5d-4fec-af47-31a2e35fef91_ConvertJob
#SBATCH -p qiita
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time 180
#SBATCH --mail-type=ALL
#SBATCH --mail-user qiita.help@gmail.com
#SBATCH --mem-per-cpu 12gb
set -x
date
hostname
cd /sequencing/igm_runs/210820_A00953_0380_BHJ53TDSX2
module load bclconvert_3.7.5
bcl-convert --sample-sheet "/qmounts/qiita_data/working_dir/3f6b1fe3-1f5d-4fec-af47-31a2e35fef91/2024-02-08_U19_Wisconsin_15445_reruns_NovaSeq_nonNA.csv" --output-directory /qmounts/qiita_data/working_dir/3f6b1fe3-1f5d-4fec-af47-31a2e35fef91/ConvertJob --bcl-input-directory . --bcl-num-decompression-threads 16 --bcl-num-conversion-threads 16 --bcl-num-compression-threads 16 --bcl-num-parallel-tiles 16 --bcl-sampleproject-subdirectories true --force
