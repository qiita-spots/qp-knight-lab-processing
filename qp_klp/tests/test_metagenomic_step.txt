# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from qp_klp.tests.test_step import BaseStepTests
from qp_klp.Assays import Metagenomic
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

    def atest_metagenomic_creation(self):
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

        trans_pipeline = Pipeline('qp_klp/tests/data/configuration.json',
                                  self.good_run_id,
                                  self.good_transcript_sheet_path,
                                  None,
                                  self.output_file_path,
                                  self.qiita_id,
                                  Metagenomic.METATRANSCRIPTOMIC_TYPE)

        step = Metagenomic(trans_pipeline, self.qiita_id, None)

        self.assertIsNotNone(step)

    def atest_metagenomic_convert_bcl_to_fastq(self):
        self._create_test_input(1)

        step = Metagenomic(self.pipeline, self.qiita_id, None)
        step.convert_raw_to_fastq()

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

    def atest_metagenomic_quality_control(self):
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

        # Since the 'sbatch' and 'squeue' commands don't exist in the testing
        # environment, fake them by creating fakes that will output to stdout
        # the metadata needed to keep run() working as intended. These faked
        # binary files overwrite the versions created by test_step.py and are
        # automatically deleted after testing. test_step.py sets chmod for us.
        with open(join(dirname(sys.executable), 'sbatch'), 'w') as f:
            # code will extract 777 for use as the fake job-id in slurm.
            f.write("echo Hello 777")

        with open(join(dirname(sys.executable), 'squeue'), 'w') as f:
            # fake squeue will return job-id 777 completed successfully.
            # faked output files created in test method() will generate
            # faked results.
            f.write("echo \"ARRAY_JOB_ID,JOBID,STATE\n777,777,COMPLETED\"")

        # execute the quality_control() method, which will in turn call NuQC's
        # run() method.
        step.quality_control()

        # after step.quality_control() executes, process_all_fastq_files.sh
        # should be created. Confirm the output of this file is as expected.
        with open(self.process_shell_script, 'r') as f:
            exp = f.readlines()
            exp = [line.rstrip() for line in exp]

        with open(join(self.output_file_path, 'NuQCJob',
                       'process_all_fastq_files.sh')) as f:
            obs = f.readlines()
            obs = [line.rstrip() for line in obs]
            with open('tmpfoobar', 'w') as foo:
                for line in obs:
                    foo.write("%s\n" % line)

            # remove part of the absolute path so that comparison test is
            # valid across multiple installations.
            patterns = [
                re.compile(r"^\s+\-\-html (.*)/qp-knight-lab-processing/"
                           r"qp_klp/tests/data/output_dir/NuQCJob/fastp_"
                           r"reports_dir/html/\${html_name} \\$"),

                re.compile(r"^\s+\-\-json (.*)/qp-knight-lab-processing/"
                           r"qp_klp/tests/data/output_dir/NuQCJob/fastp_"
                           r"reports_dir/json/\${json_name} \\$"),

                re.compile(r"^\s+python (.*/bin)/demux \\"),

                re.compile(r"    (.*)/sequence_processing_pipeline/scripts"
                           r"/splitter "),

                re.compile(r"^\s+mv \$\{jobd}/seqs.movi.txt.gz (.*)/qp-knight"
                           r"-lab-processing")
            ]

            for pattern in patterns:
                for i in range(0, len(obs)):
                    m = pattern.match((obs[i]))
                    if m:
                        obs[i] = obs[i].replace(m[1], 'REMOVED')

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
