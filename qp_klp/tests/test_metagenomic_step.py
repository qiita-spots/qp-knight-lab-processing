# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from qp_klp.tests.test_step import BaseStepTests
from qp_klp.Metagenomic import Metagenomic
from sequence_processing_pipeline.Pipeline import Pipeline
from os import makedirs
from os.path import join, split, exists, basename, dirname
from random import randrange, randint, choices
from glob import glob
from shutil import copyfile
import sys
import re


class MetagenomicTests(BaseStepTests):
    def setUp(self):
        super().setUp()

    def test_metagenomic_creation(self):
        # Test base-class creation method, even though base-class will never
        # be instantiated by itself in normal usage.

        with self.assertRaisesRegex(ValueError, "A pipeline object is needed"
                                                " to initialize Step"):
            Metagenomic(None, self.qiita_id, None)

        with self.assertRaisesRegex(ValueError, "A Qiita job-id is needed to "
                                                "initialize Step"):
            Metagenomic(self.pipeline, None, None)

        step = Metagenomic(self.pipeline, self.qiita_id, None)

        self.assertIsNotNone(step)

        makedirs(self.output_file_path, exist_ok=True)

        trans_pipeline = Pipeline(None,
                                  self.good_run_id,
                                  self.good_transcript_sheet_path,
                                  None,
                                  self.output_file_path,
                                  self.qiita_id,
                                  Metagenomic.METATRANSCRIPTOMIC_TYPE,
                                  BaseStepTests.CONFIGURATION)

        step = Metagenomic(trans_pipeline, self.qiita_id, None)

        self.assertIsNotNone(step)

    def test_metagenomic_convert_bcl_to_fastq(self):
        self._create_test_input(1)

        step = Metagenomic(self.pipeline, self.qiita_id, None)
        step.convert_bcl_to_fastq()

    def _generate_fake_fastq_files(self, output_path, file_count,
                                   max_fastq_size_in_mb):
        if not exists(output_path):
            raise ValueError("%s doesn't exist" % output_path)

        makedirs(output_path, exist_ok=True)

        # generate a random subset of file_count files to be 'small-ones' that
        # fail zero-length threshold. Select no more than 1/3rd of the files
        # requested and return the filenames to the user for testing.
        small_ones = choices(range(1, file_count), k=randint(1, file_count//3))
        return_these = []

        for i in range(1, file_count):
            forward_read_fp = join(output_path,
                                   "SAMPLE%d_L001_R1_001.fastq.gz" % i)
            reverse_read_fp = join(output_path,
                                   "SAMPLE%d_L001_R2_001.fastq.gz" % i)

            if i in small_ones:
                # NuQC's minimum threshold for moving files to zero-files
                # directory is under 3100 bytes.
                file_size = randrange(1024, 3099)
            else:
                # where 1048576 equals 1MB in bytes
                file_size = randrange(1048576, max_fastq_size_in_mb * 1048576)

            if file_size < 3100:
                return_these.append(forward_read_fp)
                return_these.append(reverse_read_fp)

            with open(forward_read_fp, 'w') as f:
                # assume 'A' is one byte in size on disk.
                f.write("A" * file_size)

            # for convenience, let reverse read be equal in size to forward
            # read.
            with open(reverse_read_fp, 'w') as f:
                f.write("A" * file_size)

        return return_these

    def test_metagenomic_quality_control(self):
        self._create_test_input(2)

        metadata = {'NYU_BMS_Melanoma_13059': {'needs_filtering': False,
                                               'samples': []},
                    'Feist_11661': {'needs_filtering': False, 'samples': []},
                    'Gerwick_6123': {'needs_filtering': True, 'samples': []}}

        # since Slurm isn't available for tests, and process_all_fastq_files.sh
        # can't actually be executed, generate a fake set of raw fastq files
        # in a 'ConvertJob' directory. In place of actual trimmed/filtered
        # files, copy the faked files into their expected location under
        # 'NuQCJob'.
        # This allows us to still confirm that the run() method and functions
        # like zero-length file filtering still work as intended. We can also
        # still confirm that process_all_fastq_files.sh is created and looks
        # as expected.

        small_ones = []
        for project_name in metadata:
            small_ones += self._generate_fake_fastq_files(
                join(self.output_file_path, 'ConvertJob', project_name), 25,
                10)

        step = Metagenomic(self.pipeline, self.qiita_id, None)

        for project_name in metadata:
            # after initialization of step object but before run() is called,
            # copy the raw files and rename them into QC'ed files.
            src = join(self.output_file_path, 'ConvertJob', project_name)
            sub_dir = 'filtered_sequences' if metadata[project_name][
                'needs_filtering'] else 'trimmed_sequences'
            dst = join(self.output_file_path, 'NuQCJob', project_name, sub_dir)

            makedirs(dst, exist_ok=True)

            faked_files = glob(src + '/SAMPLE*.fastq.gz')
            for fp in faked_files:
                _, file_name = split(fp)
                file_path = join(dst, file_name.replace('.fastq.gz',
                                                        '.trimmed.fastq.gz'))
                copyfile(fp, file_path)

        # Since the 'sbatch' and 'sacct' commands don't exist in the testing
        # environment, fake them by creating fakes that will output to stdout
        # the metadata needed to keep run() working as intended. These faked
        # binary files overwrite the versions created by test_step.py and are
        # automatically deleted after testing. test_step.py sets chmod for us.
        with open(join(dirname(sys.executable), 'sbatch'), 'w') as f:
            # code will extract 777 for use as the fake job-id in slurm.
            f.write("echo Hello 777")

        with open(join(dirname(sys.executable), 'sacct'), 'w') as f:
            # fake sacct will return job-id 777 completed successfully.
            # faked output files created in test method() will generate
            # faked results.
            f.write("echo \"777|cc_fake_job.sh|COMPLETED|00:10:00|0:0\"")

        # execute the quality_control() method, which will in turn call NuQC's
        # run() method.
        step.quality_control()

        # after step.quality_control() executes, process_all_fastq_files.sh
        # should be created. Confirm the output of this file is as expected.

        exp = ["#!/bin/bash -l",
               "#SBATCH -J 077c4da8-74eb-4184-8860-0207f53623be_NuQCJob",
               "#SBATCH -p qiita", "### wall-time-limit in minutes",
               "#SBATCH --time 60", "#SBATCH --mem 20G",
               "#SBATCH -N 1", "#SBATCH -c 4", "",
               "if [[ -z \"${SLURM_ARRAY_TASK_ID}\" ]]; then",
               "    echo \"Not operating within an array\"", "    exit 1",
               "fi", "", "if [[ \"${SLURM_ARRAY_TASK_MIN}\" -ne 1 ]]; then",
               "    echo \"Min array ID is not 1\"", "    exit 1", "fi", "",
               "if [[ -z ${MMI} ]]; then", "    echo \"MMI is not set\"",
               "    exit 1", "fi", "", "if [[ -z ${PREFIX} ]]; then",
               "    echo \"PREFIX is not set\"", "    exit 1", "fi", "",
               "if [[ -z ${OUTPUT} ]]; then", "    echo \"OUTPUT is not set\"",
               "    exit 1", "fi", "", "if [[ -z ${TMPDIR} ]]; then",
               "    echo \"TMPDIR is not set\"", "    exit 1", "fi", "",
               "echo \"MMI is ${MMI}\"", "",
               "conda activate human-depletion", "", "set -x", "set -e", "",
               "date", "hostname",
               "echo ${SLURM_JOBID} ${SLURM_ARRAY_TASK_ID}",
               "### output_path = output_path passed to Job objects + "
               "'NuQCJob'",
               "### e.g.: working-directory/ConvertJob, working-directory/"
               "QCJob...",
               "cd REMOVED/qp-knight-lab-processing/qp_klp/tests"
               "/data/output_dir/NuQCJob/tmp", "",
               "### set a temp directory, make a new unique one under it and",
               "### make sure we clean up as we're dumping to shm",
               "### DO NOT do this casually. Only do a clean up like this if",
               "### you know for sure TMPDIR is what you want.", "",
               "mkdir -p ${TMPDIR}", "export TMPDIR=${TMPDIR}",
               "export TMPDIR=$(mktemp -d)", "echo $TMPDIR", "",
               ("mkdir -p REMOVED/qp-knight-lab-processing/qp_klp/tests/data"
                "/output_dir/NuQCJob/fastp_reports_dir/html"),
               ("mkdir -p REMOVED/qp-knight-lab-processing/qp_klp/tests/data"
                "/output_dir/NuQCJob/fastp_reports_dir/json"),
               "",
               "function cleanup {", "  echo \"Removing $TMPDIR\"",
               "  rm -fr $TMPDIR", "  unset TMPDIR", "}", "trap cleanup EXIT",
               "",
               "export FILES=$(pwd)/$(printf \"%s-%d\" ${PREFIX} ${SLURM_ARRAY"
               "_TASK_ID})",
               "if [[ ! -f ${FILES} ]]; then", "    logger ${FILES} not found",
               "    exit 1", "fi", "", "delimiter=::MUX::",
               "n=$(wc -l ${FILES} | cut -f 1 -d\" \")", "",
               "for i in $(seq 1 ${n})", "do", "    echo \"Beginning loop "
                                               "iteration ${i}\"",
               "    line=$(head -n ${i} ${FILES} | tail -n 1)",
               "    r1=$(echo ${line} | cut -f 1 -d\" \")",
               "    r2=$(echo ${line} | cut -f 2 -d\" \")",
               "    base=$(echo ${line} | cut -f 3 -d\" \")",
               "    r1_name=$(basename ${r1} .fastq.gz)",
               "    r2_name=$(basename ${r2} .fastq.gz)", "",
               "    # for now just make sure each file is saved and we can "
               "read the data inside",
               "    # to sort them out later.",
               "    html_name=$(echo \"$r1_name.html\")",
               "    json_name=$(echo \"$r1_name.json\")",
               "",
               "    echo \"${i}	${r1_name}	${r2_name}	${base}\" >> ${TMPDIR}"
               "/id_map",
               "", "    fastp \\", "        -l 45 \\", "        -i ${r1} \\",
               "        -I ${r2} \\", "        -w 7 \\",
               "        --adapter_fasta fastp_known_adapters_formatted.fna \\",
               ("        --html REMOVED/qp-knight-lab-processing/qp_klp/tests"
                "/data/output_dir/NuQCJob/fastp_reports_dir/html/"
                "${html_name} \\"),
               ("        --json REMOVED/qp-knight-lab-processing/qp_klp/tests"
                "/data/output_dir/NuQCJob/fastp_reports_dir/json/"
                "${json_name} \\"),
               "        --stdout | \\",
               "            sed -r \"1~4s/^@(.*)/@${i}${delimiter}\\1/\"",
               "done > ${TMPDIR}/seqs.fastq", "",
               "function minimap2_runner () {",
               "    mmi=$1", "", "    echo \"$(date) :: $(basename ${mmi})\"",
               "    minimap2 -2 -ax sr -t 7 ${mmi} ${TMPDIR}/seqs.fastq | \\",
               "        samtools fastq -@ 1 -f 12 -F 256 > ${TMPDIR}/seqs_new."
               "fastq",
               "    mv ${TMPDIR}/seqs_new.fastq ${TMPDIR}/seqs.fastq", "}", "",
               "function runner () {", "    i=${1}", "",
               "    # with the current proposed resources, we likely do not "
               "get",
               "    # benefit for parallel gzip write",
               "    REMOVED/demux \\",
               "        --id-map ${TMPDIR}/id_map \\",
               "        --infile ${TMPDIR}/seqs.fastq \\",
               "        --output ${OUTPUT} \\",
               "        --encoded-id ${i} \\",
               "        --threads 1", "}",
               "export -f runner", "", "if [[ -f ${MMI} ]]; then",
               "    minimap2_runner ${MMI}", "else",
               "    for mmi in ${MMI}/*.mmi", "    do",
               "        minimap2_runner ${mmi}", "    done", "fi", "",
               "mkdir -p ${OUTPUT}", "", "jobs=${SLURM_JOB_CPUS_PER_NODE}", "",
               "echo \"$(date) :: demux start\"", "",
               "seq 1 ${n} | parallel -j ${jobs} runner", "",
               "echo \"$(date) :: demux stop\"",
               "",
               "touch ${OUTPUT}/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.completed",
               ""]

        with open(join(self.output_file_path, 'NuQCJob',
                       'process_all_fastq_files.sh')) as f:
            obs = f.readlines()
            obs = [line.rstrip() for line in obs]

            # remove part of the absolute path so that comparison test is
            # valid across multiple installations.
            p = re.compile(r"^cd (.*)/qp-knight-lab-processing/qp_klp/tests"
                           r"/data/output_dir/NuQCJob/tmp$")
            m = p.match(obs[51])
            obs[51] = obs[51].replace(m[1], 'REMOVED')

            p = re.compile(r"mkdir -p (.*)/qp-knight-lab-processing/qp_klp/"
                           r"tests/data/output_dir/NuQCJob/fastp_reports_dir/"
                           r"html")
            m = p.match(obs[63])
            obs[63] = obs[63].replace(m[1], 'REMOVED')

            p = re.compile(r"mkdir -p (.*)/qp-knight-lab-processing/qp_klp/"
                           r"tests/data/output_dir/NuQCJob/fastp_reports_dir/"
                           r"json")
            m = p.match(obs[64])
            obs[64] = obs[64].replace(m[1], 'REMOVED')

            p = re.compile(r"^\s+\-\-html (.*)/qp-knight-lab-processing/"
                           r"qp_klp/tests/data/output_dir/NuQCJob/fastp_"
                           r"reports_dir/html/\${html_name} \\$")
            m = p.match(obs[105])
            obs[105] = obs[105].replace(m[1], 'REMOVED')

            p = re.compile(r"^\s+\-\-json (.*)/qp-knight-lab-processing/"
                           r"qp_klp/tests/data/output_dir/NuQCJob/fastp_"
                           r"reports_dir/json/\${json_name} \\$")
            m = p.match(obs[106])
            obs[106] = obs[106].replace(m[1], 'REMOVED')

            p = re.compile(r"^\s+(.*/bin)/demux \\")
            m = p.match(obs[125])
            obs[125] = obs[125].replace(m[1], 'REMOVED')

            for obs_line, exp_line in zip(obs, exp):
                self.assertEqual(obs_line, exp_line)

        # Lastly, for the random subset of zero-length files 'small_ones',
        # generate a list of files that run() failed and moved into the
        # zero_files sub-folder and compare them for equality.
        zf_files = []
        for project_name in metadata:
            zf_path = join(self.output_file_path, 'NuQCJob', project_name,
                           'zero_files')
            zf_files += glob(join(zf_path, 'SAMPLE*.fastq.gz'))

        # small_ones lists the original raw files when they were created.
        # However, zf_files contains files that have been QC'ed. For
        # comparison purposes, remove the .trimmed or .filtered substring
        # from the names.
        small_ones = [basename(fp) for fp in small_ones]
        zf_files = [
            basename(fp).replace('.trimmed.', '.').replace('.filtered.', '.')
            for fp in zf_files]

        self.assertEqual(set(small_ones), set(zf_files))

        # if end of this function is reached without an Error raised, then
        # run() ran correctly.
