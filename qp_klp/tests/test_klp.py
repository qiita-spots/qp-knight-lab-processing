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
from os.path import exists, isdir, realpath, dirname

from qiita_client.testing import PluginTestCase
from qiita_client import ArtifactInfo

from qp_klp import __version__, plugin

from qp_klp.klp import sequence_processing_pipeline


class KLPTests(PluginTestCase):
    def setUp(self):
        # this will allow us to see the full errors
        self.maxDiff = None

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

        plugin("https://localhost:8383", "register", "ignored")
        self._clean_up_files = []

        self.basedir = dirname(realpath(__file__))

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def test_sequence_processing_pipeline(self):
        # not a valid run_identifier folder and sample_sheet
        params = {"run_identifier": "/this/path/doesnt/exist",
                  "sample_sheet": "NA"}

        data = {
            "user": "demo@microbio.me",
            "command": dumps(["qp-klp", __version__,
                              "Sequence Processing Pipeline"]),
            "status": "running",
            "parameters": dumps(params),
        }
        jid = self.qclient.post("/apitest/processing_job/", data=data)["job"]
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        success, ainfo, msg = sequence_processing_pipeline(
            self.qclient, jid, params, out_dir
        )
        self.assertFalse(success)
        self.assertEqual(msg, "The path doesn't exist or is not a folder")

        # valid run_identifier folder but not sample_sheet
        # NOTE: we are no creating a new job for this test, which is fine
        params = {"run_identifier": out_dir, "sample_sheet": "NA"}

        success, ainfo, msg = sequence_processing_pipeline(
            self.qclient, jid, params, out_dir
        )
        self.assertFalse(success)
        self.assertEqual(msg, "Doesn't look like a valid uploaded file; "
                         "please review.")

        # test success
        # both valid run_identifier and sample_sheet
        # NOTE: we are no creating a new job for this test, which is fine

        params = {
            "run_identifier": out_dir,
            "sample_sheet": {
                "body": self.sample_csv_data,
                "content_type": "text/plain",
                "filename": "prep_16S.txt",
            },
        }

        success, ainfo, msg = sequence_processing_pipeline(
            self.qclient, jid, params, out_dir
        )
        self.assertTrue(success)
        exp = [
            ArtifactInfo(
                "output",
                "job-output-folder",
                [(f"{out_dir}/", "directory")]
            )
        ]

        self.assertEqual(ainfo, exp)


if __name__ == "__main__":
    main()
