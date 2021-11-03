# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from unittest import main
from os import remove, makedirs
from shutil import rmtree
from json import dumps
from tempfile import mkdtemp
from os.path import exists, isdir, join, realpath, dirname
from qiita_client.testing import PluginTestCase
from qiita_client import ArtifactInfo
from qp_klp import __version__, plugin
from qp_klp.klp import sequence_processing_pipeline
from time import sleep
from os import environ


class KLPTests(PluginTestCase):
    def _make_temp_dir(self):
        new_dir = mkdtemp()
        self._clean_up_files.append(new_dir)
        return new_dir

    def setUp(self):
        # this will allow us to see the full errors
        self.maxDiff = None

        plugin('https://localhost:8383', 'register', 'ignored')
        sleep(2)

        self._clean_up_files = []

        self.basedir = dirname(realpath(__file__))

        self.multiqc_config_filepath = join(self._make_temp_dir(),
                                            'multiqc-config.yaml')
        self.multiqc_config_data = [
            "title: 'Sequence processing summaries'\n",
            "output_fn_name: 'index.html'\n",
            "show_analysis_paths: False\n",
            "show_analysis_time: False\n",
            "run_modules:\n",
            "    - bclconvert\n",
            "    - fastqc\n",
            "    - fastp\n",
            "sample_names_ignore:\n",
            "    - 'blank*'\n",
            "    - 'BLANK*'\n",
            "no_version_check: False\n",
            "plots_force_interactive: True\n",
            "num_datasets_plot_limit: 1\n",
            "max_table_rows: 10000\n",
            "table_columns_visible:\n",
            "    'Sequence Quality (bclconvert raw)':\n",
            "        percent_fails: False\n",
            "        percent_duplicates: False\n",
            "        percent_gc: False\n",
            "        avg_sequence_length: True\n",
            "        total_sequences: True\n",
            "    'Trimming':\n",
            "        input_format: False\n",
            "        avg_sequence_length: False\n",
            "        total_record_count: False\n",
            "        mean_sequence_length: False\n",
            "        fraction_bp_trimmed: False\n",
            "        fraction_records_with_adapters: False\n",
            "    'Sequence Quality (trimmed)':\n",
            "        percent_fails: False\n",
            "        percent_duplicates: False\n",
            "        percent_gc: False\n",
            "        avg_sequence_length: True\n",
            "        total_sequences: True\n",
            "    'Human Filtering':\n",
            "        overall_alignment_rate: True\n",
            "    'Sequence Quality (filtered)':\n",
            "        percent_fails: False\n",
            "        percent_duplicates: False\n",
            "        percent_gc: False\n",
            "        avg_sequence_length: True\n",
            "        total_sequences: True\n",
            "module_order:\n",
            "    - bclconvert:\n",
            "         name: 'Base Calling'\n",
            "         info: 'Conversion from BCL files to FASTQ files.'\n",
            "    - fastqc:\n",
            "        name: 'Sequence Quality (raw)'\n",
            "        info: 'Sequence quality and summary statistics for raw s",
            "equences.'\n",
            "        path_filters:\n",
            "            - '*fastqc.zip'\n",
            "            - '*.csv'\n",
            "    - fastp:\n",
            "        name: 'Sequence Quality (adapter trimmed)'\n",
            "        info: 'Summary statistics from adapter trimming and qual",
            "ity control with fastp.'\n",
            "        fn: '*.json'\n",
            "    - fastqc:\n",
            "        name: 'Sequence Quality (trimmed)'\n",
            "        info: 'Sequence quality and summary statistics after qua",
            "lity-control and adapter trimming.'\n",
            "        path_filters:\n",
            "          - '*trimmed_fastqc.zip'\n",
            "          - '*fastp_fastqc.zip'\n"
        ]

        self.sample_csv_data = [
            "[Header],,,,,,,,,,\n",
            "IEMFileVersion,4,,,,,,,,,\n",
            "Investigator Name,Knight,,,,,,,,,\n",
            "Experiment Name,RKL0042,,,,,,,,,\n",
            "Date,2/26/20,,,,,,,,,\n",
            "Workflow,GenerateFASTQ,,,,,,,,,\n",
            "Application,FASTQ Only,,,,,,,,,\n",
            "Assay,Metagenomics,,,,,,,,,\n",
            "Description,,,,,,,,,,\n",
            "Chemistry,Default,,,,,,,,,\n",
            ",,,,,,,,,,\n",
            "[Reads],,,,,,,,,,\n",
            "150,,,,,,,,,,\n",
            "150,,,,,,,,,,\n",
            ",,,,,,,,,,\n",
            "[Settings],,,,,,,,,,\n",
            "ReverseComplement,0,,,,,,,,,\n",
            ",,,,,,,,,,\n",
            "[Data],,,,,,,,,,\n",
            "Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,"
            "index,I5_Index_ID,index2,Sample_Project,Well_description\n",
            "1,CDPH-SAL_Salmonella_Typhi_MDL-143,CDPH-SAL_Salmonella_Typhi_MD"
            "L-143,Feist_11661_P40,A1,iTru7_107_07,CCGACTAT,iTru5_01_A,ACCGAC"
            "AA,Feist_11661,CDPH-SAL_Salmonella Typhi_MDL-143\n",
            "1,CDPH-SAL_Salmonella_Typhi_MDL-144,CDPH-SAL_Salmonella_Typhi_MD"
            "L-144,Feist_11661_P40,C1,iTru7_107_08,CCGACTAT,iTru5_02_A,CTTCGC"
            "AA,Feist_11661,CDPH-SAL_Salmonella Typhi_MDL-144\n",
            ",,,,,,,,,,\n",
            "[Bioinformatics],,,,,,,,,,\n",
            "Sample_Project,QiitaID,BarcodesAreRC,ForwardAdapter,ReverseAdapt"
            "er,HumanFiltering,library_construction_protocol,experiment_desig"
            "n_description,,,\n",
            "Feist_11661,11661,FALSE,AACC,GGTT,FALSE,Knight Lab Kapa HP,Eqiip"
            "eriment,,,\n",
            ",,,,,,,,,,\n",
            "[Contact],,,,,,,,,,\n",
            "Email,Sample_Project,,,,,,,,,\n",
            "test@lol.com,Feist_11661,,,,,,,,,\n",
            ",,,,,,,,,,\n",
        ]

        self.out_dir = self._make_temp_dir()
        self.search_dir = self._make_temp_dir()

        spp_config = {
              "configuration": {
                "pipeline": {
                  "younger_than": 90,
                  "older_than": 24,
                  "archive_path": '/not/out_dir',
                  "search_paths": [self.search_dir]
                },
                "bcl2fastq": {
                  "nodes": 1,
                  "nprocs": 16,
                  "queue": "qiita",
                  "wallclock_time_in_hours": 36,
                  "modules_to_load": ["bcl2fastq_2.20.0.422"],
                  "executable_path": "bcl2fastq",
                  "per_process_memory_limit": "10gb"
                },
                "bcl-convert": {
                  "nodes": 1,
                  "nprocs": 16,
                  "queue": "qiita",
                  "wallclock_time_in_hours": 36,
                  "modules_to_load": ["bclconvert_3.7.5"],
                  "executable_path": "bcl-convert",
                  "per_process_memory_limit": "10gb"
                },
                "qc": {
                  "nodes": 1,
                  "nprocs": 16,
                  "queue": "qiita",
                  "wallclock_time_in_hours": 1,
                  "mmi_db": "/databases/minimap2/human-phix-db.mmi",
                  "modules_to_load": ["fastp_0.20.1", "samtools_1.12",
                                      "minimap2_2.18"],
                  "fastp_executable_path": "fastp",
                  "minimap2_executable_path": "minimap2",
                  "samtools_executable_path": "samtools",
                  "job_total_memory_limit": "20gb",
                  "job_pool_size": 30
                },
                "seqpro": {
                  "seqpro_path": "seqpro",
                  "modules_to_load": []
                },
                "fastqc": {
                  "nodes": 1,
                  "nprocs": 16,
                  "queue": "qiita",
                  "nthreads": 16,
                  "wallclock_time_in_hours": 1,
                  "modules_to_load": ["fastqc_0.11.5"],
                  "fastqc_executable_path": "fastqc",
                  "multiqc_executable_path": "multiqc",
                  "multiqc_config_file_path": self.multiqc_config_filepath,
                  "job_total_memory_limit": "20gb",
                  "job_pool_size": 30
                }
              }
            }

        # create the file and write the configuration out to disk
        # for use by sequence_processing_pipeline().
        config_filepath = environ['QP_KLP_CONFIG_FP']

        with open(config_filepath, 'w') as f:
            f.write(dumps(spp_config, indent=2))

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def test_sequence_processing_pipeline(self):
        # not a valid run_identifier folder and sample_sheet
        params = {"run_identifier": "NOT_A_RUN_IDENTIFIER",
                  "sample_sheet": "NA"}

        data = {
            "user": "demo@microbio.me",
            "command": dumps(["qp-klp", __version__,
                              "Sequence Processing Pipeline"]),
            "status": "running",
            "parameters": dumps(params),
        }

        job_id = self.qclient.post("/apitest/processing_job/",
                                   data=data)["job"]

        success, ainfo, msg = sequence_processing_pipeline(
            self.qclient, job_id, params, self.out_dir
        )
        self.assertFalse(success)
        self.assertEqual(msg, "This doesn't appear to be a valid sample sheet"
                              "; please review.")

        test_dir = join(self.search_dir, "200318_A00953_0082_AH5TWYDSXY")
        makedirs(test_dir)

        # create the sentinel files ConvertJob will check for.
        with open(join(test_dir, 'RTAComplete.txt'), 'w') as f:
            f.write("Hello World\n")

        with open(join(test_dir, 'RunInfo.xml'), 'w') as f:
            f.write("Hello World\n")

        # create the project directory QCJob will expect to find fastq files
        # in.
        fastq_dir = join(test_dir, 'Data', 'Fastq', 'Feist_11661')
        makedirs(fastq_dir)
        file_list = ["CDPH-SAL_Salmonella_Typhi_MDL-143_R1_.fastq.gz",
                     "CDPH-SAL_Salmonella_Typhi_MDL-143_R2_.fastq.gz",
                     "CDPH-SAL_Salmonella_Typhi_MDL-144_R1_.fastq.gz",
                     "CDPH-SAL_Salmonella_Typhi_MDL-144_R2_.fastq.gz"]

        for fastq_file in file_list:
            with open(join(fastq_dir, fastq_file), 'w') as f:
                f.write("Hello World\n")

        # write multi-qc config file to a known location
        with open(self.multiqc_config_filepath, 'w') as f:
            for line in self.multiqc_config_data:
                f.write(f"{line}\n")

        # create the Reports directory in the location GenPrepFileJob
        # expects.
        reports_dir = join(self.out_dir, 'ConvertJob', 'Reports')
        makedirs(reports_dir, exist_ok=True)

        # create QCJobs output directory for use by GenPrepFileJob
        qcj_output_fp = join(self.out_dir, 'QCJob', 'Feist_11661')
        makedirs(join(qcj_output_fp, 'filtered_sequences'))
        makedirs(join(qcj_output_fp, 'fastp_reports_dir', 'json'))

        # valid run_identifier folder but not sample_sheet
        # NOTE: we are not creating a new job for this test, which is fine
        params = {"run_identifier": "200318_A00953_0082_AH5TWYDSXY",
                  "sample_sheet": "NA"}

        success, ainfo, msg = sequence_processing_pipeline(
            self.qclient, job_id, params, self.out_dir
        )
        self.assertFalse(success)

        self.assertEqual(msg, "This doesn't appear to be a valid sample sheet"
                              "; please review.")

        # test success
        # both valid run_identifier and sample_sheet
        # NOTE: we are no creating a new job for this test, which is fine

        params = {
            "run_identifier": "200318_A00953_0082_AH5TWYDSXY",
            "sample_sheet": {
                "body": ''.join(self.sample_csv_data),
                "content_type": "text/plain",
                "filename": "prep_16S.txt",
            }
        }

        success, ainfo, msg = sequence_processing_pipeline(
            self.qclient, job_id, params, self.out_dir
        )

        if msg:
            # if success is True, msg should be None.
            print("Message returned: %s" % msg)

        self.assertTrue(success)

        exp = [ArtifactInfo("output",
                            "job-output-folder",
                            [(f"{self.out_dir}/final_results/", "directory")]
                            )
               ]

        self.assertEqual(ainfo, exp)


if __name__ == "__main__":
    main()
