# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from unittest import main
from os import remove
from shutil import rmtree
from json import dumps
from tempfile import mkdtemp
from os.path import exists, isdir, join, realpath, dirname
from qiita_client.testing import PluginTestCase
from qiita_client import ArtifactInfo
from qp_klp import __version__, plugin
from qp_klp.klp import sequence_processing_pipeline
from time import sleep


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

        spp_config = {
              "configuration": {
                "pipeline": {
                  "younger_than": 90,
                  "older_than": 24,
                  "archive_path": '/not/out_dir',
                  "search_paths": [self.out_dir]
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
                  "multiqc_config_file_path": None,
                  "job_total_memory_limit": "20gb",
                  "job_pool_size": 30
                }
              }
            }

        # use out_dir to store the configuration.json file as well.
        # create the file and write the configuration out to disk
        # for use by sequence_processing_pipeline().
        self.config_filepath = join(self.out_dir, 'configuration.json')

        with open(self.config_filepath, 'w') as f:
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
                  "config_filepath": self.config_filepath}

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

        # valid run_identifier folder but not sample_sheet
        # NOTE: we are not creating a new job for this test, which is fine
        params = {"run_identifier": "200318_A00953_0082_AH5TWYDSXY",
                  "sample_sheet": "NA",
                  "config_filepath": self.config_filepath}

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
            },
            "config_filepath": self.config_filepath
        }

        success, ainfo, msg = sequence_processing_pipeline(
            self.qclient, job_id, params, self.out_dir
        )

        if msg:
            # if success is True, msg should be None.
            print("Message returned: %s" % msg)

        self.assertTrue(success)
        exp = [
            ArtifactInfo(
                "output",
                "job-output-folder",
                [(f"{self.out_dir}/", "directory")]
            )
        ]

        self.assertEqual(ainfo, exp)


if __name__ == "__main__":
    main()
