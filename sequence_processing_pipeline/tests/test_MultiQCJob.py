import unittest
from os.path import join, exists
from functools import partial
from sequence_processing_pipeline.MultiQCJob import MultiQCJob
from sequence_processing_pipeline.PipelineError import JobFailedError
from os import makedirs, listdir
from shutil import rmtree, move


class TestMultiQCJob(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        package_root = 'sequence_processing_pipeline'
        self.path = partial(join, package_root, 'tests', 'data')
        self.qiita_job_id = 'abcdabcdabcdabcdabcdabcdabcdabcd'
        self.output_path = self.path('output_dir2')
        self.multiqc_log_path = join(self.output_path, 'logs')
        self.raw_fastq_files_path = ('sequence_processing_pipeline/tests/data'
                                     '/211021_A00000_0000_SAMPLE/Data/Fastq/p'
                                     'roject1')
        self.processed_fastq_files_path = ('sequence_processing_pipeline/tests'
                                           '/data/211021_A00000_0000_SAMPLE/sa'
                                           'mple-sequence-directory')
        self.config_yml = join(package_root, 'multiqc-bclconvert-config.yaml')
        self.qc_root_path = join(self.output_path, 'QCJob')
        self.fastqc_root_path = join(self.output_path, 'FastQCJob')
        makedirs(self.qc_root_path, exist_ok=True)

        base_path = join(self.output_path, 'MultiQCJob')
        file_names = []
        makedirs(join(base_path, 'Feist_11661', 'trimmed_sequences'),
                 exist_ok=True)
        makedirs(join(base_path, 'Gerwick_6123', 'filtered_sequences'),
                 exist_ok=True)

        file_names.append(join(base_path,
                               'Feist_11661',
                               'trimmed_sequences',
                               ('CDPH-SAL_Salmonella_Typhi_MDL-143_S1_L003_R1_'
                                '001.trimmed_fastqc.html')))

        file_names.append(join(base_path,
                               'Gerwick_6123',
                               'filtered_sequences',
                               '3A_S169_L003_R2_001.trimmed_fastqc.html'))

        # write files out to disk
        for file_name in file_names:
            with open(file_name, 'w') as f2:
                f2.write("This is a file.")

        # set up dummy logs
        self.multiqc_log_path = join(self.output_path, "MultiQCJob", "logs")
        makedirs(self.multiqc_log_path, exist_ok=True)

        log_files = {
            'slurm-9999999_35.out': ["---------------",
                                     "Run details:",
                                     ("hds-fe848a9e-c0e9-49d9-978d-"
                                      "27565a314e8b 1908305 b2-018"),
                                     "---------------",
                                     "+ this",
                                     "+ that",
                                     "+ blah",
                                     ("something error: Generic Standin Error"
                                      " (GSE).")],
            'slurm-9999999_17.out': ["---------------",
                                     "Run details:",
                                     ("hds-fe848a9e-c0e9-49d9-978d-"
                                      "27565a314e8b 1908305 b2-018"),
                                     "---------------",
                                     "+ this",
                                     "+ that",
                                     "+ blah",
                                     ("something error: Another Standin Error"
                                      " (ASE).")]
        }

        for log_file in log_files:
            fp = join(self.multiqc_log_path, log_file)

            with open(fp, 'w') as f:
                lines = log_files[log_file]
                for line in lines:
                    f.write(f"{line}\n")

    def tearDown(self):
        rmtree(self.output_path)

        zero_path = join(self.raw_fastq_files_path, 'zero_files')

        if exists(zero_path):
            for f in listdir(zero_path):
                source = join(zero_path, f)
                move(source, join(self.raw_fastq_files_path, f))

            rmtree(zero_path)

    def test_generate_job_scripts(self):
        job = MultiQCJob(self.qc_root_path, self.output_path,
                         self.raw_fastq_files_path.replace('/project1', ''),
                         self.processed_fastq_files_path,
                         16, 16,
                         'sequence_processing_pipeline/tests/bin/multiqc',
                         ['multiqc.2.0'], self.qiita_job_id, 'queue_name', 4,
                         23, '8g', 30, self.fastqc_root_path, 1000,
                         "sequence_processing_pipeline/"
                         "multiqc-bclconvert-config.yaml", False)

        job_script_path = join(job.output_path, 'MultiQCJob.sh')
        array_details_path = join(job.output_path, 'MultiQCJob.array-details')
        self.assertEqual(exists(job_script_path), True)
        self.assertEqual(exists(array_details_path), True)

        exp = ["#!/bin/bash",
               "#SBATCH -J abcdabcdabcdabcdabcdabcdabcdabcd_MultiQCJob",
               "#SBATCH -p queue_name", "#SBATCH -N 4", "#SBATCH -n 16",
               "#SBATCH --time 23", "#SBATCH --mem 8gG",
               "#SBATCH --array 1-1%30", "set -x", "set +e",
               "set -o pipefail", "date", "hostname",
               "echo ${SLURM_JOBID} ${SLURM_ARRAY_TASK_ID}",
               "cd sequence_processing_pipeline/tests/data/output_dir2/"
               "MultiQCJob", "", "module load multiqc.2.0", "",
               "step=${SLURM_ARRAY_TASK_ID}",
               "cmd0=$(head -n $step sequence_processing_pipeline/tests/data/"
               "output_dir2/MultiQCJob/MultiQCJob.array-details | tail -n 1)",
               "eval $cmd0", "echo \"Cmd Completed: $cmd0\" > logs/MultiQCJob_"
               "$step.completed"]

        with open(job_script_path, 'r') as f:
            obs = f.readlines()
            obs = [x.strip() for x in obs]

        for a, b in zip(obs, exp):
            self.assertEqual(a, b)

        exp = ["multiqc -c sequence_processing_pipeline/multiqc-bclconvert-con"
               "fig.yaml --fullnames --force -o sequence_processing_pipeline/t"
               "ests/data/output_dir2/MultiQCJob/multiqc/project1"]

        with open(array_details_path, 'r') as f:
            obs = f.readlines()
            obs = [x.strip() for x in obs]

        for a, b in zip(obs, exp):
            self.assertEqual(a, b)

    def test_error_msg_from_logs(self):
        job = MultiQCJob(self.qc_root_path, self.output_path,
                         self.raw_fastq_files_path.replace('/project1', ''),
                         self.processed_fastq_files_path,
                         16, 16,
                         'sequence_processing_pipeline/tests/bin/multiqc',
                         ['multiqc.2.0'], self.qiita_job_id, 'queue_name', 4,
                         23, '8g', 30, self.fastqc_root_path, 1000,
                         "sequence_processing_pipeline/"
                         "multiqc-bclconvert-config.yaml", False)

        self.assertFalse(job is None)

        # an internal method to force submit_job() to raise a JobFailedError
        # instead of submitting the job w/sbatch and waiting for a failed
        # job w/squeue.
        self.assertTrue(job._toggle_force_job_fail())

        try:
            job.run()
        except JobFailedError as e:
            # assert that the text of the original error message was
            # preserved, while including known strings from each of the
            # sample log-files.
            print(str(e))
            self.assertIn('This job died.', str(e))
            self.assertIn('something error: Generic Standin Error (GSE)',
                          str(e))
            self.assertIn('something error: Another Standin Error (ASE)',
                          str(e))


if __name__ == '__main__':
    unittest.main()
