# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from unittest import main
from os import remove, makedirs, walk
from shutil import rmtree
from json import dumps
from tempfile import mkdtemp
from os.path import exists, isdir, join, realpath, dirname
from qiita_client.testing import PluginTestCase
from qiita_client import ArtifactInfo
from qp_klp import __version__, plugin
from qp_klp.klp_util import FailedSamplesRecord, map_sample_names_to_tube_ids
from qp_klp.klp import sequence_processing_pipeline
from time import sleep
from os import environ
import logging
import re
from metapool import KLSampleSheet
from shutil import copy


logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(message)s',
                    datefmt='%d-%b-%y %H:%M:%S')


class KLPTests(PluginTestCase):
    def _make_temp_dir(self):
        new_dir = mkdtemp()
        self._clean_up_files.append(new_dir)
        return new_dir

    def setUp(self):
        # this will allow us to see the full errors
        self.maxDiff = None
        self.logger = logging.getLogger(__name__)

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
            "1,CDPH-SAL_Salmonella_Typhi_MDL-143,SKD7.640191"
            ",Feist_1_P40,A1,iTru7_107_07,CCGACTAT,iTru5_01_A,ACCGAC"
            "AA,Feist_1,CDPH-SAL_Salmonella Typhi_MDL-143\n",
            "1,CDPH-SAL_Salmonella_Typhi_MDL-144,SKB8.640193XX"
            ",Feist_1_P40,C1,iTru7_107_08,CCGACTAT,iTru5_02_A,CTTCGC"
            "AA,Feist_1,CDPH-SAL_Salmonella Typhi_MDL-144\n",
            ",,,,,,,,,,\n",
            "[Bioinformatics],,,,,,,,,,\n",
            "Sample_Project,QiitaID,BarcodesAreRC,ForwardAdapter,ReverseAdapt"
            "er,HumanFiltering,library_construction_protocol,experiment_desig"
            "n_description,,,\n",
            "Feist_1,11661,FALSE,AACC,GGTT,FALSE,Knight Lab Kapa HP,Equip"
            "eriment,,,\n",
            ",,,,,,,,,,\n",
            "[Contact],,,,,,,,,,\n",
            "Email,Sample_Project,,,,,,,,,\n",
            "test@lol.com,Feist_1,,,,,,,,,\n",
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
                    "search_paths": [self.search_dir],
                    "amplicon_search_paths": [self.search_dir]
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
                    "minimap_databases": [("/databases/minimap2/human-phix-db."
                                           "mmi")],
                    "kraken2_database": "/databases/minimap2/hp_kraken-db.mmi",
                    "modules_to_load": ["fastp_0.20.1", "samtools_1.12",
                                        "minimap2_2.18"],
                    "fastp_executable_path": "fastp",
                    "minimap2_executable_path": "minimap2",
                    "samtools_executable_path": "samtools",
                    "job_total_memory_limit": "20gb",
                    "job_pool_size": 30,
                    "job_max_array_length": 1000
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
                    "job_pool_size": 30,
                    "job_max_array_length": 1000
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
                  "sample_sheet": "NA",
                  "lane_number": 1}

        data = {
            "user": "demo@microbio.me",
            "command": dumps(["qp-klp", __version__,
                              "Sequence Processing Pipeline"]),
            "status": "running",
            "parameters": dumps(params),
        }

        job_id = self.qclient.post("/apitest/processing_job/",
                                   data=data)["job"]

        success, _, msg = sequence_processing_pipeline(
            self.qclient, job_id, params, self.out_dir
        )
        self.assertFalse(success)
        self.assertEqual(msg, "This doesn't appear to be a valid sample sheet"
                              " or mapping file; please review.")

        test_dir = join(self.search_dir, "200318_A00953_0082_AH5TWYDSXY")
        makedirs(test_dir)

        # create the sentinel files ConvertJob will check for.
        with open(join(test_dir, 'RTAComplete.txt'), 'w') as f:
            f.write("Hello World\n")

        with open(join(test_dir, 'RunInfo.xml'), 'w') as f:
            f.write("Hello World\n")

        # create the project directory sequence_processing_pipeline() will
        # expect to find fastq files in.
        fastq_dir = join(self.out_dir, 'ConvertJob', 'Feist_1')
        makedirs(fastq_dir)

        file_list = ["CDPH-SAL_Salmonella_Typhi_MDL-143_R1_.fastq.gz",
                     "CDPH-SAL_Salmonella_Typhi_MDL-143_R2_.fastq.gz",
                     "CDPH-SAL_Salmonella_Typhi_MDL-144_R1_.fastq.gz",
                     "CDPH-SAL_Salmonella_Typhi_MDL-144_R2_.fastq.gz"]

        for fastq_file in file_list:
            fp = join(fastq_dir, fastq_file)
            self.logger.debug("FASTQ FILE PATH: %s" % fp)
            with open(fp, 'w') as f:
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
        qcj_output_fp = join(self.out_dir, 'QCJob', 'Feist_1')
        qcj_filtered_sequences = join(qcj_output_fp, 'filtered_sequences')
        makedirs(qcj_filtered_sequences)

        # create GenPrepFileJob output directory for use by packaging code and
        # downstream-testing.
        gp_root_fp = join(self.out_dir, 'GenPrepFileJob')
        prep_files_root_fp = join(gp_root_fp, 'PrepFiles')
        makedirs(prep_files_root_fp, exist_ok=True)
        # TODO: Fix name
        prep_file_name = ('230224_M05314_0347_000000000-KVMH3.'
                          'ABTX_20230227_11052.1.tsv')

        # copy sample prep-info file into position.
        copy(join(self.basedir, prep_file_name),
             join(prep_files_root_fp, prep_file_name))

        # valid run_identifier folder but not sample_sheet
        # NOTE: we are not creating a new job for this test, which is fine
        params = {"run_identifier": "200318_A00953_0082_AH5TWYDSXY",
                  "sample_sheet": "NA",
                  "lane_number": 1}

        success, _, msg = sequence_processing_pipeline(
            self.qclient, job_id, params, self.out_dir
        )
        self.assertFalse(success)

        self.assertEqual(msg, "This doesn't appear to be a valid sample sheet"
                              " or mapping file; please review.")

        # test error due to missing sample_names

        params = {
            "run_identifier": "200318_A00953_0082_AH5TWYDSXY",
            "sample_sheet": {
                "body": ''.join(self.sample_csv_data),
                "content_type": "text/plain",
                # verify sequence_processing_pipeline() will convert spaces
                # to underscores ('_').
                "filename": "A sample sheet.csv",
            },
            "lane_number": 2
        }

        success, _, msg = sequence_processing_pipeline(
            self.qclient, job_id, params, self.out_dir)

        self.assertFalse(success)
        self.assertTrue(
            msg.startswith("Feist_1 has 1 missing samples (i.e. "
                           "SKB8.640193XX). Some samples from Qiita:"))
        self.assertTrue(
           msg.endswith(". No tube_id column in Qiita."))

        # test success
        # both valid run_identifier and sample_sheet
        # NOTE: we are not creating a new job for this test, which is fine

        # fix self.sample_csv_data
        self.sample_csv_data[21] = self.sample_csv_data[21].replace(
            'SKB8.640193XX', 'SKB8.640193')
        params = {
            "run_identifier": "200318_A00953_0082_AH5TWYDSXY",
            "sample_sheet": {
                "body": ''.join(self.sample_csv_data),
                "content_type": "text/plain",
                # verify sequence_processing_pipeline() will convert spaces
                # to underscores ('_').
                "filename": "A sample sheet.csv",
            },
            "lane_number": 2
        }

        # add files to QCJob/filtered...
        new_files = ['file1_R1_file1.trimmed.fastq.gz',
                     'file1_R2_file1.trimmed.fastq.gz',
                     'file1_I1_file1.trimmed.fastq.gz',
                     'file1_I2_file1.trimmed.fastq.gz',
                     'file2_R1_file2.trimmed.fastq.gz',
                     'file2_R2_file2.trimmed.fastq.gz',
                     'file2_I1_file2.trimmed.fastq.gz',
                     'file2_I2_file2.trimmed.fastq.gz']

        for new_file in new_files:
            nfp = join(qcj_filtered_sequences, new_file)
            with open(nfp, 'w') as nf:
                nf.write("Hello World!\n")

        success, ainfo, msg = sequence_processing_pipeline(
            self.qclient, job_id, params, self.out_dir
        )
        self.assertEquals(msg, 'Main Pipeline Finished, processing results')
        self.assertTrue(success)

        exp = [ArtifactInfo("output",
                            "job-output-folder",
                            [(f"{self.out_dir}/final_results/", "directory")]
                            )
               ]

        self.assertEqual(ainfo, exp)

        # confirm that 'cmds.log' exists.
        cmd_log_fp = join(self.out_dir, 'cmds.log')
        self.assertTrue(exists(cmd_log_fp))

        # confirm that fastq files were copied to uploads directory.
        uploads_fp = ('/home/runner/work/qp-knight-lab-processing/'
                      'qp-knight-lab-processing/qiita-dev/qiita_db/'
                      'support_files/test_data/uploads/11661')

        for some_file in file_list:
            some_path = join(uploads_fp, some_file)
            print("checking '%s' exists..." % some_path)
            self.assertTrue(exists(some_path))

        # confirm that an output directory named 'final_results' was created
        # by the pipeline and that 'prep_files.tgz' is one of the products
        # inside.
        self.assertTrue(exists(join(self.out_dir, 'final_results',
                                    'prep-files.tgz')))

        # confirm touched_studies.html was generated.
        ts_fp = join(self.out_dir, 'final_results', 'touched_studies.html')
        self.assertTrue(exists(ts_fp))

        # verify sequence_processing_pipeline() will convert spaces
        # to underscores ('_').
        self.assertTrue(exists(join(f"{self.out_dir}", 'A_sample_sheet.csv')))

        # verify cmd.log
        exp = ['cd OUT_DIR; tar zcvf logs-ConvertJob.tgz ConvertJob/logs',
               ('cd OUT_DIR; tar zcvf reports-ConvertJob.tgz ConvertJob/Repor'
                'ts ConvertJob/Logs'),
               'cd OUT_DIR; tar zcvf logs-QCJob.tgz QCJob/logs',
               'cd OUT_DIR; tar zcvf logs-FastQCJob.tgz FastQCJob/logs',
               'cd OUT_DIR; tar zcvf reports-FastQCJob.tgz FastQCJob/fastqc',
               ('cd OUT_DIR; tar zcvf logs-GenPrepFileJob.tgz GenPrepFileJob/'
                'logs'),
               'cd OUT_DIR; tar zcvf prep-files.tgz GenPrepFileJob/PrepFiles',
               ('cd OUT_DIR; tar zcvf reports-QCJob.tgz QCJob/Feist_1/fas'
                'tp_reports_dir'),
               'cd PREFIX/support_files/test_data/uploads/11661',
               'cd OUT_DIR; mv *.tgz final_results',
               'cd OUT_DIR; mv FastQCJob/multiqc final_results',
               'cd OUT_DIR; mv touched_studies.html final_results']

        cmdslog_fp = join(self.out_dir, 'cmds.log')
        with open(cmdslog_fp, 'r') as f:
            # read all lines into a list
            cmds = f.readlines()
            # remove newlines
            cmds = [x.strip() for x in cmds]
            # replace randomly-generated tmp directory with fixed text.
            cmds = [re.sub(r'^cd .*?;', r'cd OUT_DIR;', x) for x in cmds]

            cmds = [re.sub(r' .*\/support_files\/test_data\/uploads\/11661$', # noqa
                           r' PREFIX/support_files/test_data/uploads/11661',
                           x) for x in cmds]

            self.assertEqual(exp, cmds)

        # Note that because we are using self.sample_csv_data instead of
        # good-sample-sheet.csv as our sample-sheet, touched_studies.html
        # will include only the one project Feist_1, instead of all
        # three studies found in good-sample-sheet.csv.
        with open(ts_fp, 'r') as f:
            obs = f.readlines()
            obs = [x.strip() for x in obs]
            obs = ''.join(obs)
            exp = ('<table border="2" class="dataframe"><thead><tr style="text'
                   '-align: left;"><th>Project</th><th>Qiita Study ID</th><th>'
                   'Qiita URL</th></tr></thead><tbody><tr><td>Feist_1</td>'
                   '<td>11661</td><td><a href="https://localhost:21174/study/'
                   'description/11661" target="_blank">https://localhost:21174'
                   '/study/description/11661</a></td></tr></tbody></table>')
            self.assertEqual(obs, exp)

    def test_failed_samples_recorder(self):
        # since unittests can't run third-party code like bcl2fastq and
        # skip_exec will bypass a Job's run() and audit() commands, we will
        # instead test the FailedSamplesRecord class to confirm that it works
        # as expected.

        # self.basedir = .../qp-knight-lab-processing/qp-knight-lab-processing
        # /qp_klp/tests/
        sheet = KLSampleSheet(f'{self.basedir}/good-sample-sheet.csv')
        fsr = FailedSamplesRecord(self.basedir, sheet.samples)

        # we want to include samples from all projects in the sample-sheet.
        # order of projects listed is Feist_1, NYU_BMS_Melanoma_13059, and
        # Gerwick_6123.
        fail_set1 = ['Pputida_TALE__HGL_Pputida_121', 'EP073160B01', '5B']
        fail_set2 = ['Deoxyribose_PALE_ALE__MG1655_Lib4_20_16', 'EP202095B04',
                     '4A']
        fail_set3 = ['JM-MEC__Staphylococcus_aureusstrain_BERTI-R10727',
                     'EP159695B01', '6A']

        # simulate three write calls out to file. Each successive call should
        # append to the information already written out to file.
        fsr.write(fail_set1, 'ConvertJob')

        with open(f'{self.basedir}/failed_samples.html', 'r') as f:
            obs1 = f.readlines()
            obs1 = [x.strip() for x in obs1]
            obs1 = ''.join(obs1)
            exp1 = ('<table border="2" class="dataframe"><thead><tr style="tex'
                    't-align: left;"><th>Project</th><th>Sample ID</th><th>Fai'
                    'led at</th></tr></thead><tbody><tr><td>Feist_11661</td><t'
                    'd>Pputida_TALE__HGL_Pputida_121</td><td>ConvertJob</td></'
                    'tr><tr><td>Gerwick_6123</td><td>5B</td><td>ConvertJob</td'
                    '></tr><tr><td>NYU_BMS_Melanoma_13059</td><td>EP073160B01<'
                    '/td><td>ConvertJob</td></tr></tbody></table>')
            self.assertEqual(obs1, exp1)

        fsr.write(fail_set2, 'QCJob')
        with open(f'{self.basedir}/failed_samples.html', 'r') as f:
            obs2 = f.readlines()
            obs2 = [x.strip() for x in obs2]
            obs2 = ''.join(obs2)
            exp2 = ('<table border="2" class="dataframe"><thead><tr style="tex'
                    't-align: left;"><th>Project</th><th>Sample ID</th><th>Fai'
                    'led at</th></tr></thead><tbody><tr><td>Gerwick_6123</td><'
                    'td>4A</td><td>QCJob</td></tr><tr><td>Feist_11661</td><td>'
                    'Pputida_TALE__HGL_Pputida_121</td><td>ConvertJob</td></tr'
                    '><tr><td>Gerwick_6123</td><td>5B</td><td>ConvertJob</td><'
                    '/tr><tr><td>Feist_11661</td><td>Deoxyribose_PALE_ALE__MG1'
                    '655_Lib4_20_16</td><td>QCJob</td></tr><tr><td>NYU_BMS_Mel'
                    'anoma_13059</td><td>EP202095B04</td><td>QCJob</td></tr><t'
                    'r><td>NYU_BMS_Melanoma_13059</td><td>EP073160B01</td><td>'
                    'ConvertJob</td></tr></tbody></table>')
            self.assertEqual(obs2, exp2)

        fsr.write(fail_set3, 'FastQCJob')
        with open(f'{self.basedir}/failed_samples.html', 'r') as f:
            obs3 = f.readlines()
            obs3 = [x.strip() for x in obs3]
            obs3 = ''.join(obs3)
            exp3 = ('<table border="2" class="dataframe"><thead><tr style="tex'
                    't-align: left;"><th>Project</th><th>Sample ID</th><th>Fai'
                    'led at</th></tr></thead><tbody><tr><td>Gerwick_6123</td><'
                    'td>4A</td><td>QCJob</td></tr><tr><td>Feist_11661</td><td>'
                    'Pputida_TALE__HGL_Pputida_121</td><td>ConvertJob</td></tr'
                    '><tr><td>Gerwick_6123</td><td>5B</td><td>ConvertJob</td><'
                    '/tr><tr><td>Gerwick_6123</td><td>6A</td><td>FastQCJob</td'
                    '></tr><tr><td>Feist_11661</td><td>Deoxyribose_PALE_ALE__M'
                    'G1655_Lib4_20_16</td><td>QCJob</td></tr><tr><td>Feist_116'
                    '61</td><td>JM-MEC__Staphylococcus_aureusstrain_BERTI-R107'
                    '27</td><td>FastQCJob</td></tr><tr><td>NYU_BMS_Melanoma_13'
                    '059</td><td>EP159695B01</td><td>FastQCJob</td></tr><tr><t'
                    'd>NYU_BMS_Melanoma_13059</td><td>EP202095B04</td><td>QCJo'
                    'b</td></tr><tr><td>NYU_BMS_Melanoma_13059</td><td>EP07316'
                    '0B01</td><td>ConvertJob</td></tr></tbody></table>')
            self.assertEqual(obs3, exp3)

    def test_map_sample_names_to_tube_ids(self):
        sn_tid_map_by_proj = {'Sample_Project': {}}
        sn_tid_map_by_proj['Sample_Project']['363192526'] = 'tube_id01'
        sn_tid_map_by_proj['Sample_Project']['363192073'] = 'tube_id02'
        sn_tid_map_by_proj['Sample_Project']['363193755'] = 'tube_id03'
        sn_tid_map_by_proj['Sample_Project']['363192568'] = 'tube_id04'
        sn_tid_map_by_proj['Sample_Project']['363193764'] = 'tube_id05'
        sn_tid_map_by_proj['Sample_Project']['363197837'] = 'tube_id06'
        sn_tid_map_by_proj['Sample_Project']['363193058'] = 'tube_id07'
        sn_tid_map_by_proj['Sample_Project']['363192067'] = 'tube_id08'
        sn_tid_map_by_proj['Sample_Project']['363192066'] = 'tube_id09'
        sn_tid_map_by_proj['Sample_Project']['363193001'] = 'tube_id10'
        sn_tid_map_by_proj['Sample_Project']['363197803'] = 'tube_id11'
        sn_tid_map_by_proj['Sample_Project']['363192078'] = 'tube_id12'
        sn_tid_map_by_proj['Sample_Project']['363192050'] = 'tube_id13'
        sn_tid_map_by_proj['Sample_Project']['363197086'] = 'tube_id14'
        sn_tid_map_by_proj['Sample_Project']['363197089'] = 'tube_id15'
        sn_tid_map_by_proj['Sample_Project']['363193005'] = 'tube_id16'
        sn_tid_map_by_proj['Sample_Project']['363193007'] = 'tube_id17'
        sn_tid_map_by_proj['Sample_Project']['363197110'] = 'tube_id18'
        sn_tid_map_by_proj['Sample_Project']['363192124'] = 'tube_id19'
        sn_tid_map_by_proj['Sample_Project']['363193040'] = 'tube_id20'
        sn_tid_map_by_proj['Sample_Project']['363192553'] = 'tube_id21'
        sn_tid_map_by_proj['Sample_Project']['363197824'] = 'tube_id22'
        sn_tid_map_by_proj['Sample_Project']['363192082'] = 'tube_id23'
        sn_tid_map_by_proj['Sample_Project']['363192558'] = 'tube_id24'
        sn_tid_map_by_proj['Sample_Project']['363192598'] = 'tube_id25'
        sn_tid_map_by_proj['Sample_Project']['363192112'] = 'tube_id26'
        sn_tid_map_by_proj['Sample_Project']['363197075'] = 'tube_id27'
        sn_tid_map_by_proj['Sample_Project']['363192059'] = 'tube_id28'
        sn_tid_map_by_proj['Sample_Project']['363192578'] = 'tube_id29'
        sn_tid_map_by_proj['Sample_Project']['363192096'] = 'tube_id30'
        sn_tid_map_by_proj['Sample_Project']['363192105'] = 'tube_id31'
        sn_tid_map_by_proj['Sample_Project']['363193000'] = 'tube_id32'
        sn_tid_map_by_proj['Sample_Project']['363193014'] = 'tube_id33'
        sn_tid_map_by_proj['Sample_Project']['363193028'] = 'tube_id34'
        sn_tid_map_by_proj['Sample_Project']['363192559'] = 'tube_id35'
        sn_tid_map_by_proj['Sample_Project']['363193016'] = 'tube_id36'
        sn_tid_map_by_proj['Sample_Project']['363193032'] = 'tube_id37'
        sn_tid_map_by_proj['Sample_Project']['363191841'] = 'tube_id38'
        sn_tid_map_by_proj['Sample_Project']['363192566'] = 'tube_id39'
        sn_tid_map_by_proj['Sample_Project']['363193762'] = 'tube_id40'
        sn_tid_map_by_proj['Sample_Project']['363192796'] = 'tube_id41'
        sn_tid_map_by_proj['Sample_Project']['363192074'] = 'tube_id42'
        sn_tid_map_by_proj['Sample_Project']['363193071'] = 'tube_id43'
        sn_tid_map_by_proj['Sample_Project']['363192094'] = 'tube_id44'
        sn_tid_map_by_proj['Sample_Project']['363193024'] = 'tube_id45'
        sn_tid_map_by_proj['Sample_Project']['363197031'] = 'tube_id46'
        sn_tid_map_by_proj['Sample_Project']['363193054'] = 'tube_id47'
        sn_tid_map_by_proj['Sample_Project']['363192054'] = 'tube_id48'
        sn_tid_map_by_proj['Sample_Project']['363192547'] = 'tube_id49'

        output_dir = 'qp_klp/tests'

        copy(join(output_dir, 'good-prep-file.txt'),
             join(output_dir, 'Sample_Project.tsv'))

        results = map_sample_names_to_tube_ids(sn_tid_map_by_proj, output_dir)

        for project in results:
            for prep_file in results[project]:
                df = results[project][prep_file]

                # confirm .tsv file loaded everything as strings
                types = dict(df.dtypes)
                unique_types = set([types[k] for k in types])
                bad_types_present = False
                for type in unique_types:
                    if type != object:
                        bad_types_present = True
                        break

                self.assertFalse(bad_types_present)
                res = df['sample_name'].tolist()
                res = [x for x in res if not x.startswith('tube_id')]

                self.assertEquals(len(res), 0)

                # write modified results back out to file
                df.to_csv(prep_file, index=False, sep="\t")


class KLPAmpliconTests(PluginTestCase):
    def _make_temp_dir(self):
        new_dir = mkdtemp()
        self._clean_up_files.append(new_dir)
        return new_dir

    def setUp(self):
        # this will allow us to see the full errors
        self.maxDiff = None
        self.logger = logging.getLogger(__name__)

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

        self.out_dir = self._make_temp_dir()
        self.search_dir = self._make_temp_dir()

        spp_config = {
            "configuration": {
                "pipeline": {
                    "younger_than": 90,
                    "older_than": 24,
                    "archive_path": '/not/out_dir',
                    "search_paths": [self.search_dir],
                    "amplicon_search_paths": [self.search_dir]
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
                    "minimap_databases": [("/databases/minimap2/human-phix-db."
                                           "mmi")],
                    "kraken2_database": "/databases/minimap2/hp_kraken-db.mmi",
                    "modules_to_load": ["fastp_0.20.1", "samtools_1.12",
                                        "minimap2_2.18"],
                    "fastp_executable_path": "fastp",
                    "minimap2_executable_path": "minimap2",
                    "samtools_executable_path": "samtools",
                    "job_total_memory_limit": "20gb",
                    "job_pool_size": 30,
                    "job_max_array_length": 1000
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
                    "job_pool_size": 30,
                    "job_max_array_length": 1000
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
        test_dir = join(self.search_dir, "230224_M05314_0347_000000000-KVMH3")
        makedirs(test_dir)

        # create the sentinel files ConvertJob will check for.
        with open(join(test_dir, 'RTAComplete.txt'), 'w') as f:
            f.write("Hello World\n")

        # copy example RunInfo.xml into its proper location for testing.
        copy(f'{self.basedir}/RunInfo.xml', join(test_dir, 'RunInfo.xml'))

        # create the project directory sequence_processing_pipeline() will
        # expect to find fastq files in. Note amplicon pipeline expects
        # fastq files in ConvertJob folder, rather than ConvertJob/{Project}
        # folders.
        fastq_dir = join(self.out_dir, 'ConvertJob')
        makedirs(fastq_dir)

        file_list = [("230224_M05314_0347_000000000-KVMH3_SMPL1_S1_"
                      "L001_I1_001.fastq.gz"),
                     ("230224_M05314_0347_000000000-KVMH3_SMPL1_S1_"
                      "L001_R1_001.fastq.gz"),
                     ("230224_M05314_0347_000000000-KVMH3_SMPL1_S1_"
                      "L001_R2_001.fastq.gz")]

        for fastq_file in file_list:
            fp = join(fastq_dir, fastq_file)
            with open(fp, 'w') as f:
                f.write("Hello World\n")

        # write multi-qc config file to a known location
        with open(self.multiqc_config_filepath, 'w') as f:
            for line in self.multiqc_config_data:
                f.write(f"{line}\n")

        # create the Reports directory in the location GenPrepFileJob
        # expects.
        reports_dir = join(self.out_dir, 'ConvertJob', 'Reports')
        makedirs(reports_dir, exist_ok=True)

        # create GenPrepFileJob output directory for use by packaging code and
        # downstream-testing.
        gp_root_fp = join(self.out_dir, 'GenPrepFileJob')
        prep_files_root_fp = join(gp_root_fp, 'PrepFiles')
        makedirs(prep_files_root_fp, exist_ok=True)
        prep_file_name = ('230224_M05314_0347_000000000-KVMH3.'
                          'ABTX_20230227_11052.1.tsv')

        # copy sample prep-info file into position.
        copy(join(self.basedir, prep_file_name),
             join(prep_files_root_fp, prep_file_name))

        # the only difference between this test and test_spp_no_qiita_id_error
        # is project_names are missing qiita_id in bad_mapping_file.txt.
        with open(f'{self.basedir}/good_mapping_file.txt', 'r') as f:
            mapping_file = f.readlines()
            mapping_file = ''.join(mapping_file)

        params = {"run_identifier": "230224_M05314_0347_000000000-KVMH3",
                  "sample_sheet": {
                      "body": mapping_file,
                      "content_type": "text/plain",
                      # verify sequence_processing_pipeline() will convert
                      # spaces to underscores ('_').
                      "filename": "A sample sheet.csv",
                  },
                  "lane_number": 1}

        data = {
            "user": "demo@microbio.me",
            "command": dumps(["qp-klp", __version__,
                              "Sequence Processing Pipeline"]),
            "status": "running",
            "parameters": dumps(params),
        }

        job_id = self.qclient.post("/apitest/processing_job/",
                                   data=data)["job"]

        success, _, msg = sequence_processing_pipeline(
            self.qclient, job_id, params, self.out_dir
        )

        # on error, it is beneficial to report the msg.
        self.assertEqual(msg, 'Main Pipeline Finished, processing results')
        self.assertTrue(success)

        # confirm that 'cmds.log' exists.
        cmd_log_fp = join(self.out_dir, 'cmds.log')
        self.assertTrue(exists(cmd_log_fp))

        # confirm that fastq files were copied to uploads directory.
        uploads_fp = ('/home/runner/work/qp-knight-lab-processing/'
                      'qp-knight-lab-processing/qiita-dev/qiita_db/'
                      'support_files/test_data/uploads/1')

        for some_file in file_list:
            some_path = join(uploads_fp, some_file)
            self.assertTrue(exists(some_path))

        # confirm that an output directory named 'final_results' was created
        # by the pipeline and that 'prep_files.tgz' is one of the products
        # inside.
        self.assertTrue(exists(join(self.out_dir, 'final_results',
                                    'prep-files.tgz')))

        # confirm touched_studies.html was generated.
        self.assertTrue(exists(join(self.out_dir, 'final_results',
                                    'touched_studies.html')))

    def test_spp_no_qiita_id_error(self):
        test_dir = join(self.search_dir, "230224_M05314_0347_000000000-KVMH3")
        makedirs(test_dir)

        # create the sentinel files ConvertJob will check for.
        with open(join(test_dir, 'RTAComplete.txt'), 'w') as f:
            f.write("Hello World\n")

        # copy example RunInfo.xml into its proper location for testing.
        copy(f'{self.basedir}/RunInfo.xml', join(test_dir, 'RunInfo.xml'))

        # create the project directory sequence_processing_pipeline() will
        # expect to find fastq files in. Note amplicon pipeline expects
        # fastq files in ConvertJob folder, rather than ConvertJob/{Project}
        # folders.
        fastq_dir = join(self.out_dir, 'ConvertJob')
        makedirs(fastq_dir)

        file_list = ["CDPH-SAL_Salmonella_Typhi_MDL-143_R1_.fastq.gz",
                     "CDPH-SAL_Salmonella_Typhi_MDL-143_R2_.fastq.gz",
                     "CDPH-SAL_Salmonella_Typhi_MDL-144_R1_.fastq.gz",
                     "CDPH-SAL_Salmonella_Typhi_MDL-144_R2_.fastq.gz"]

        for fastq_file in file_list:
            fp = join(fastq_dir, fastq_file)
            with open(fp, 'w') as f:
                f.write("Hello World\n")

        # write multi-qc config file to a known location
        with open(self.multiqc_config_filepath, 'w') as f:
            for line in self.multiqc_config_data:
                f.write(f"{line}\n")

        # create the Reports directory in the location GenPrepFileJob
        # expects.
        reports_dir = join(self.out_dir, 'ConvertJob', 'Reports')
        makedirs(reports_dir, exist_ok=True)

        with open(f'{self.basedir}/bad_mapping_file.txt', 'r') as f:
            mapping_file = f.readlines()
            mapping_file = ''.join(mapping_file)

        params = {"run_identifier": "230224_M05314_0347_000000000-KVMH3",
                  "sample_sheet": {
                      "body": mapping_file,
                      "content_type": "text/plain",
                      # verify sequence_processing_pipeline() will convert
                      # spaces to underscores ('_').
                      "filename": "A sample sheet.csv",
                  },
                  "lane_number": 1}

        data = {
            "user": "demo@microbio.me",
            "command": dumps(["qp-klp", __version__,
                              "Sequence Processing Pipeline"]),
            "status": "running",
            "parameters": dumps(params),
        }

        job_id = self.qclient.post("/apitest/processing_job/",
                                   data=data)["job"]

        success, _, msg = sequence_processing_pipeline(
            self.qclient, job_id, params, self.out_dir
        )

        self.assertEqual(msg, "Values in the project_name column must be "
                              "appended with a Qiita ID.")
        self.assertFalse(success)


if __name__ == "__main__":
    main()
