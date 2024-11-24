# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from qp_klp.WorkflowFactory import WorkflowFactory
from unittest import TestCase
from os import makedirs
from qp_klp.SequencingTech import SEQTECH_NAME_ILLUMINA
from qp_klp.Assays import ASSAY_NAME_METAGENOMIC
from shutil import rmtree

'''

class MetagenomicWorkflowTests(TestCase):
    def setUp(self):
        self.output_dir = "qp_klp/tests/211021_A00000_0000_SAMPLE"
        self.remove_these = []

        kwargs = {"uif_path": "qp_klp/tests/data/sample-sheets/metagenomic/"
                  "illumina/good_sheet1.csv",
                  "qclient": None,
                  "lane_number": "1",
                  "config_fp": "qp_klp/tests/data/configuration.json",
                  "run_identifier": "211021_A00000_0000_SAMPLE",
                  "output_dir": self.output_dir,
                  "job_id": "78901",
                  "is_restart": False
                  }

        self._create_directory(kwargs['output_dir'])
        self.wf = WorkflowFactory.generate_workflow(**kwargs)

    def tearDown(self):
        for fp in self.remove_these:
            rmtree(fp)

    def _create_directory(self, fp):
        makedirs(fp, exist_ok=True)
        self.remove_these.append(fp)

    def test_execute_pipeline(self):
        print(self.wf.what_am_i())
        print(dir(self.wf))

        # write a test that tests confirm_mandatory_attributes by creating with multiple subsets of missing attributes?
        # write a test that tests determine_steps_to_skip by setting up run-directories in different states of completion.
        # write a test to execute pipeline all the way through.
        # confirm failed samples record is included in the end to end test.

        self.assertTrue(False)


    # do we need a separate FSR test?
'''
