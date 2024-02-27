#!/bin/bash -l
#SBATCH -J 077c4da8-74eb-4184-8860-0207f53623be_NuQCJob
#SBATCH -p qiita
### wall-time-limit in minutes
#SBATCH --time 2028
#SBATCH --mem 20G
#SBATCH -N 2
### Note cores_per_task maps to fastp & minimap2 thread counts
### as well as sbatch -c. demux threads remains fixed at 1.
### Note -c set to 4 and thread counts set to 7 during testing.
#SBATCH -c 2

echo "---------------"
echo "Run details:"
echo "$SLURM_JOB_NAME $SLURM_JOB_ID $SLURMD_NODENAME $SLURM_ARRAY_TASK_ID"
echo "---------------"

if [[ -z "${SLURM_ARRAY_TASK_ID}" ]]; then
    echo "Not operating within an array"
    exit 1
fi
if [[ -z ${MMI} ]]; then
    echo "MMI is not set"
    exit 1
fi
if [[ -z ${PREFIX} ]]; then
    echo "PREFIX is not set"
    exit 1
fi
if [[ -z ${OUTPUT} ]]; then
    echo "OUTPUT is not set"
    exit 1
fi

echo "MMI is ${MMI}"

conda activate qp-knight-lab-processing-2022.03
module load fastp_0.20.1 samtools_1.12 minimap2_2.18

set -x
set -e
set -o pipefail

export FILES=$(printf "%s-%d" ${PREFIX} ${SLURM_ARRAY_TASK_ID})
if [[ ! -f ${FILES} ]]; then
    logger ${FILES} not found
    exit 1
fi
# set a temp directory, make a new unique one under it and
# make sure we clean up as we're dumping to shm
# DO NOT do this casually. Only do a clean up like this if
# you know for sure TMPDIR is what you want.

TMPDIR=/dev/shm
export TMPDIR=${TMPDIR}
export TMPDIR=$(mktemp -d)
echo $TMPDIR

mkdir -p REMOVED/qp-knight-lab-processing/qp_klp/tests/data/output_dir/NuQCJob/fastp_reports_dir/html
mkdir -p REMOVED/qp-knight-lab-processing/qp_klp/tests/data/output_dir/NuQCJob/fastp_reports_dir/json

export ADAPTER_ONLY_OUTPUT=${OUTPUT}/only-adapter-filtered
mkdir -p ${ADAPTER_ONLY_OUTPUT}

function cleanup {
  echo "Removing $TMPDIR"
  rm -fr $TMPDIR
  unset TMPDIR
}
trap cleanup EXIT

export delimiter=::MUX::
export r1_tag=/1
export r2_tag=/2
function mux-runner () {
    n=$(wc -l ${FILES} | cut -f 1 -d" ")

    jobd=${TMPDIR}
    id_map=${jobd}/id_map
    seqs_r1=${jobd}/seqs.r1.fastq
    seqs_r2=${jobd}/seqs.r2.fastq
    r1_filt=${jobd}/seqs.r1.adapter-removed.fastq
    r2_filt=${jobd}/seqs.r2.adapter-removed.fastq

    for i in $(seq 1 ${n})
    do
        line=$(head -n ${i} ${FILES} | tail -n 1)
        r1=$(echo ${line} | cut -f 1 -d" ")
        r2=$(echo ${line} | cut -f 2 -d" ")
        base=$(echo ${line} | cut -f 3 -d" ")
        r1_name=$(basename ${r1} .fastq.gz)
        r2_name=$(basename ${r2} .fastq.gz)
        r1_adapter_only=${ADAPTER_ONLY_OUTPUT}/${r1_name}.fastq.gz
        r2_adapter_only=${ADAPTER_ONLY_OUTPUT}/${r2_name}.fastq.gz

        s_name=$(basename "${r1}" | sed -r 's/\.fastq\.gz//')
        html_name=$(echo "$s_name.html")
        json_name=$(echo "$s_name.json")

        echo "${i}	${r1_name}	${r2_name}	${base}" >> ${id_map}

        fastp \
            -l 100 \
            -i ${r1} \
            -I ${r2} \
            -w 2 \
            --adapter_fasta fastp_known_adapters_formatted.fna \
            --html REMOVED/qp-knight-lab-processing/qp_klp/tests/data/output_dir/NuQCJob/fastp_reports_dir/html/${html_name} \
            --json REMOVED/qp-knight-lab-processing/qp_klp/tests/data/output_dir/NuQCJob/fastp_reports_dir/json/${json_name} \
            --out1 ${r1_filt} \
            --out2 ${r2_filt}

        # multiplex and write adapter filtered data all at once
        cat ${r1_filt} | \
            sed -r "1~4s/^@(.*)/@${i}${delimiter}\1/" \
            >> ${seqs_r1} &
        cat ${r2_filt} | \
            sed -r "1~4s/^@(.*)/@${i}${delimiter}\1/" \
            >> ${seqs_r2} &
        cat ${r1_filt} | \
            gzip -c > ${r1_adapter_only} &
        cat ${r2_filt} | \
            gzip -c > ${r2_adapter_only} &
        wait

        rm ${r1_filt} &
        rm ${r2_filt} &
        wait
    done
}
export -f mux-runner

function minimap2-runner () {
    if [[ -z ${fcurrent} ]]; then
        echo "fcurrent not set"
        exit 1
    fi

    jobd=${TMPDIR}

    id_map=${jobd}/id_map
    if [[ ! -f ${id_map} ]]; then
        echo "No samples..."
        return
    fi

    seqs_r1=${jobd}/seqs.r1.fastq
    seqs_r2=${jobd}/seqs.r2.fastq
    seqs_mmpe_r1=${jobd}/seqs.mmpe.r1.fastq
    seqs_mmpe_r2=${jobd}/seqs.mmpe.r2.fastq
    seqs_mmpese=${jobd}/seqs.mmpese.fastq
    seqs_mmpese_r1=${jobd}/seqs.mmpese.r1.fastq
    seqs_mmpese_r2=${jobd}/seqs.mmpese.r2.fastq

    # PE operation
    minimap2 -2 -ax sr -t 12 \
        ${fcurrent} \
        ${seqs_r1} \
        ${seqs_r2} | \
            samtools fastq \
                -@ 4 \
                -f 12 \
                -F 256 \
                -N \
                -1 ${seqs_mmpe_r1} \
                -2 ${seqs_mmpe_r2}

    rm ${seqs_r1} &
    rm ${seqs_r2} &
    # no need to block

    # SE operation
    # run r1/r2 serially to avoid minimap2 picking up on interleve
    minimap2 -2 -ax sr -t 12 \
        ${fcurrent} \
        <(cat ${seqs_mmpe_r1} ${seqs_mmpe_r2}) | \
            samtools fastq \
                -@ 4 \
                -f 4 \
                -F 256 \
                -0 ${seqs_mmpese}

    rm ${seqs_mmpe_r1} &
    rm ${seqs_mmpe_r2} &
    # no need to block

    REMOVED/sequence_processing_pipeline/scripts/splitter ${seqs_mmpese} ${seqs_mmpese_r1} ${delimiter} ${r1_tag} &
    REMOVED/sequence_processing_pipeline/scripts/splitter ${seqs_mmpese} ${seqs_mmpese_r2} ${delimiter} ${r2_tag} &
    wait

    rm ${seqs_mmpese} &

    fastq_pair -t 50000000 ${seqs_mmpese_r1} ${seqs_mmpese_r2}
    rm ${seqs_mmpese_r1}.single.fq &
    rm ${seqs_mmpese_r2}.single.fq &
    rm ${seqs_mmpese_r1} &
    rm ${seqs_mmpese_r2} &
    wait

    mv ${seqs_mmpese_r1}.paired.fq ${seqs_r1}
    mv ${seqs_mmpese_r2}.paired.fq ${seqs_r2}
}
export -f minimap2-runner

function demux-runner () {
    n_demux_jobs=${SLURM_CPUS_PER_TASK}
    jobd=${TMPDIR}
    id_map=${jobd}/id_map
    seqs_r1=${jobd}/seqs.r1.fastq
    seqs_r2=${jobd}/seqs.r2.fastq

    id_map=${jobd}/id_map
    if [[ ! -f ${id_map} ]]; then
        echo "No samples..."
        return
    fi

    for idx in $(seq 0 ${n_demux_jobs})
    do
        REMOVED/demux \
            --id-map ${id_map} \
            --infile <(cat ${seqs_r1} ${seqs_r2}) \
            --output ${OUTPUT} \
            --task ${idx} \
            --maxtask ${n_demux_jobs} &
    done
    wait
}
export -f demux-runner

mmi_files=${TMPDIR}/mmi-files
if [[ -d ${MMI} ]]; then
    /bin/ls -1 ${MMI}/*.mmi > ${mmi_files}
else
    echo ${MMI} > ${mmi_files}
fi
n_files=$(wc -l ${mmi_files} | cut -f 1 -d" ")

mux-runner

# process with minimap2
for idx in $(seq 1 ${n_files})
do
    fcurrent=$(head -n ${idx} ${mmi_files} | tail -n 1)
    export fcurrent
    echo "$(date) :: $(basename ${fcurrent})"
    minimap2-runner
done

mkdir -p ${OUTPUT}

echo "$(date) :: demux start"
demux-runner
echo "$(date) :: demux stop"

touch ${OUTPUT}/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.completed
