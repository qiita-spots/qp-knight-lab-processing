# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from os import makedirs
from os.path import abspath, dirname, exists, join
from pathlib import Path
from shutil import copyfile, rmtree
from unittest import main

import pandas as pd
from metapool.sample_sheet import PROTOCOL_NAME_ILLUMINA, PROTOCOL_NAME_TELLSEQ
from qiita_client.testing import PluginTestCase

from qp_klp.Assays import (
    ASSAY_NAME_AMPLICON,
    ASSAY_NAME_METAGENOMIC,
    ASSAY_NAME_METATRANSCRIPTOMIC,
)
from qp_klp.WorkflowFactory import WorkflowFactory


class WorkflowFactoryTests(PluginTestCase):
    def setUp(self):
        self.base_dir = dirname(abspath(__file__))
        self.output_dir = join(self.base_dir, "test_output")
        self.remove_these = [self.output_dir]

        self.gz_source = f"{self.base_dir}/data/dummy.fastq.gz"

    def tearDown(self):
        for fp in self.remove_these:
            if exists(fp):
                rmtree(fp)

    def _create_directory(self, fp):
        # just making sure that we always start with a clean folder:
        if exists(fp):
            rmtree(fp)
        makedirs(fp)
        self.remove_these.append(fp)

    def test_no_parameters(self):
        # confirm factory doesn't create workflows when given no parameters.
        msg = "kwargs must not be None and must define 'uif_path'"
        with self.assertRaisesRegex(ValueError, msg):
            WorkflowFactory.generate_workflow()

        with self.assertRaisesRegex(ValueError, msg):
            WorkflowFactory.generate_workflow(**{})

        with self.assertRaisesRegex(ValueError, msg):
            WorkflowFactory.generate_workflow(**{"foo": "bar"})

    def test_workflow_performs_parameter_checking(self):
        # confirm that once factory is given enough information to select a
        # Workflow() class, the class itself confirms that it has all the
        # parameters it needs.
        kwargs = {
            "uif_path": "tests/data/sample-sheets/metagenomic/illumina/good_sheet1.csv"
        }

        msg = (
            "The following values must also be defined in kwargs for "
            "StandardMetagenomicWorkflow workflows: qclient, lane_number,"
            " config_fp, run_identifier, output_dir, job_id, lane_number,"
            " is_restart"
        )

        with self.assertRaisesRegex(ValueError, msg):
            WorkflowFactory.generate_workflow(**kwargs)

        kwargs = {"uif_path": "tests/data/pre-preps/good_pre_prep1.txt"}

        msg = (
            "The following values must also be defined in kwargs for "
            "StandardAmpliconWorkflow workflows: qclient, config_fp, "
            "run_identifier, output_dir, job_id, is_restart"
        )

        with self.assertRaisesRegex(ValueError, msg):
            WorkflowFactory.generate_workflow(**kwargs)

    def test_invalid_sample_sheets(self):
        # confirm an otherwise good-sample-sheet w/a bad SheetType is going
        # fail because it doesn't pass validation. SheetType directly
        # determines the Instrument Mixin to be used.
        kwargs = {
            "uif_path": "tests/data/sample-sheets/metagenomic/illumina/bad_sheet1.csv"
        }

        msg = "'not_a_metag' is an unrecognized SheetType"

        with self.assertRaisesRegex(ValueError, msg):
            WorkflowFactory.generate_workflow(**kwargs)

        # confirm an otherwise good-sample-sheet w/a bad Assay value is going
        # to fail because it doesn't pass validation.
        kwargs = {
            "uif_path": "tests/data/sample-sheets/metagenomic/illumina/bad_sheet2.csv"
        }

        msg = "'NotMetagenomic' is an unrecognized Assay type"

        with self.assertRaisesRegex(ValueError, msg):
            WorkflowFactory.generate_workflow(**kwargs)

        # confirm a file that is obviously not a mapping-file or sample-sheet
        # file will not follow through.
        kwargs = {"uif_path": "tests/data/Demultiplex_Stats.csv"}

        msg = (
            "Your uploaded file doesn't appear to be a sample-sheet or a mapping-file."
        )
        with self.assertRaisesRegex(ValueError, msg):
            WorkflowFactory.generate_workflow(**kwargs)

    def _inject_data(self, wf):
        """This is a helper method for testing that all the steps are run
        when testing wf.execute_pipeline()
        """
        samples = wf.pipeline.sample_sheet.samples
        tellseq = True
        if "index" in dict(samples[0]).keys():
            tellseq = False

        # inject Convert/NuQC/FastQC/MultiQC/GenPrepFileJob files so we can
        # move down the pipeline; first let's create the base folders
        if tellseq:
            reports_dir = f"{self.output_dir}/TellReadJob/"
            convert_dir = f"{self.output_dir}/TRIntegrateJob/integrated"
            makedirs(convert_dir, exist_ok=True)
        else:
            convert_dir = f"{self.output_dir}/ConvertJob"
            reports_dir = f"{self.output_dir}/ConvertJob/Reports"
            makedirs(convert_dir, exist_ok=True)
            Path(f"{convert_dir}/job_completed").touch()
        tellread_dir = f"{self.output_dir}/TellReadJob"
        nuqc_dir = f"{self.output_dir}/NuQCJob"
        fastqc_dir = f"{self.output_dir}/FastQCJob/logs/"
        multiqc_dir = f"{self.output_dir}/MultiQCJob/logs/"
        genprep_dir = f"{self.output_dir}/GenPrepFileJob/211021_A00000_0000_SAMPLE/"
        makedirs(reports_dir, exist_ok=True)
        makedirs(nuqc_dir, exist_ok=True)
        makedirs(fastqc_dir, exist_ok=True)
        makedirs(multiqc_dir, exist_ok=True)
        makedirs(genprep_dir, exist_ok=True)
        # now let's create the required project folders
        for project in wf.pipeline.get_project_info():
            sp = project["project_name"]
            makedirs(f"{convert_dir}/{sp}", exist_ok=True)
            makedirs(f"{nuqc_dir}/filtered_sequences/{sp}", exist_ok=True)
            makedirs(f"{genprep_dir}/{sp}/filtered_sequences/", exist_ok=True)

        # then loop over samples and stage all fastq.gz files
        dstats = []
        for i, sample in enumerate(samples):
            rp = sample["Sample_ID"]
            sp = sample["Sample_Project"]

            # ConvertJob
            if tellseq:
                bid = sample["barcode_id"]
                dstats.append(
                    {
                        "Lane": sample["Lane"],
                        "SampleID": rp,
                        "Sample_Project": sp,
                        "Barcode": bid,
                        "# Reads": 2,
                    }
                )
                dname = f"{convert_dir}/{sp}"
                copyfile(self.gz_source, f"{dname}/{bid}.R1.fastq.gz")
                copyfile(self.gz_source, f"{dname}/{bid}.R2.fastq.gz")
            else:
                dstats.append(
                    {
                        "Lane": sample["Lane"],
                        "SampleID": rp,
                        "Sample_Project": sp,
                        "Index": sample["index"],
                        "# Reads": 2,
                    }
                )
                dname = f"{convert_dir}/{sp}"
                Path(f"{dname}/{rp}_L001_R1_001.fastq.gz").touch()
                Path(f"{dname}/{rp}_L001_R2_001.fastq.gz").touch()

            # NuQCJob
            dname = f"{nuqc_dir}/filtered_sequences/{sp}"
            copyfile(self.gz_source, f"{dname}/{rp}_L001_R1_001.fastq.gz")
            copyfile(self.gz_source, f"{dname}/{rp}_L001_R2_001.fastq.gz")

            # GenPrepFileJob
            gprep_base = f"{genprep_dir}/{sp}/filtered_sequences/{rp}"
            Path(f"{gprep_base}_L001_R1_001.fastq.gz").touch()
            Path(f"{gprep_base}_L001_R2_001.fastq.gz").touch()

        # this is required by the Convert step
        if tellseq:
            pd.DataFrame(dstats).set_index("SampleID").to_csv(
                f"{tellread_dir}/sample_index_list_TellReadJob.txt"
            )
        else:
            pd.DataFrame(dstats).set_index("SampleID").to_csv(
                f"{reports_dir}/Demultiplex_Stats.csv"
            )

        # generating the "*.completed" files
        for i in range(len(samples) * 3):
            Path(f"{fastqc_dir}/FastQCJob_{i}.completed").touch()
            Path(f"{multiqc_dir}/MultiQCJob_{i}.completed").touch()

    def test_metagenomic_workflow_creation(self):
        kwargs = {
            "uif_path": "tests/data/sample-sheets/metagenomic/illumina/good_sheet1.csv",
            "qclient": self.qclient,
            "lane_number": "1",
            "config_fp": "tests/configuration.json",
            "run_identifier": "211021_A00000_0000_SAMPLE",
            "output_dir": self.output_dir,
            "job_id": "78901",
            "is_restart": False,
        }

        self._create_directory(kwargs["output_dir"])

        wf = WorkflowFactory.generate_workflow(**kwargs)

        # confirm that the proper type of workflow was generated.
        self.assertEqual(wf.protocol_type, PROTOCOL_NAME_ILLUMINA)
        self.assertEqual(wf.assay_type, ASSAY_NAME_METAGENOMIC)

        self._inject_data(wf)

        # Metagenomic is a valid data type in the default qiita test
        # database but job-id: 78901 doesn't exist; however, if we get
        # to here, it means that all steps have run to completion
        # and the system is trying to create the preps.
        # Note: Qiita job_id's are UUID in the database and this tests
        # uses 78901 as the job_id so the db complains about the format
        with self.assertRaisesRegex(
            RuntimeError, 'invalid input syntax for type uuid: "78901"'
        ):
            wf.execute_pipeline()

    def test_metatranscriptomic_workflow_creation(self):
        kwargs = {
            "uif_path": "tests/data/sample-sheets/"
            "metatranscriptomic/illumina/good_sheet1.csv",
            "qclient": self.qclient,
            "lane_number": "1",
            "config_fp": "tests/configuration.json",
            "run_identifier": "211021_A00000_0000_SAMPLE",
            "output_dir": self.output_dir,
            "job_id": "78901",
            "is_restart": False,
        }

        self._create_directory(kwargs["output_dir"])

        wf = WorkflowFactory.generate_workflow(**kwargs)

        # confirm that the proper type of workflow was generated.
        self.assertEqual(wf.protocol_type, PROTOCOL_NAME_ILLUMINA)
        self.assertEqual(wf.assay_type, ASSAY_NAME_METATRANSCRIPTOMIC)

        self._inject_data(wf)
        # Metatranscriptomic is not a valid data type in the default qiita test
        # database; however, if we get to here, it means that all steps have
        # ran to completion and the system is trying to create the preps.
        with self.assertRaisesRegex(
            RuntimeError, "Metatranscriptomic not valid for table data_type"
        ):
            wf.execute_pipeline()

    def test_amplicon_workflow_creation(self):
        runid = "211021_A00000_0000_SAMPLE"
        kwargs = {
            "uif_path": "tests/data/pre-preps/good_pre_prep1.txt",
            "qclient": self.qclient,
            "config_fp": "tests/configuration.json",
            "run_identifier": runid,
            "output_dir": self.output_dir,
            "job_id": "78901",
            "is_restart": False,
        }

        self._create_directory(kwargs["output_dir"])

        wf = WorkflowFactory.generate_workflow(**kwargs)

        # confirm that the proper type of workflow was generated.
        self.assertEqual(wf.protocol_type, PROTOCOL_NAME_ILLUMINA)
        self.assertEqual(wf.assay_type, ASSAY_NAME_AMPLICON)

        convert_dir = f"{self.output_dir}/ConvertJob/"
        convertr_dir = f"{self.output_dir}/ConvertJob/{runid}_SMPL1/"
        makedirs(convertr_dir, exist_ok=True)
        copyfile(
            self.gz_source, f"{convertr_dir}/{runid}_SMPL1_S00_L001_R1_001.fastq.gz"
        )
        copyfile(
            self.gz_source, f"{convertr_dir}/{runid}_SMPL1_S00_L001_R2_001.fastq.gz"
        )
        Path(f"{convert_dir}/job_completed").touch()

        nuqc_dir = f"{self.output_dir}/NuQCJob/filtered_sequences/{runid}_SMPL1"
        makedirs(nuqc_dir, exist_ok=True)
        copyfile(self.gz_source, f"{nuqc_dir}/{runid}_SMPL1_S00_L001_R1_001.fastq.gz")
        copyfile(self.gz_source, f"{nuqc_dir}/{runid}_SMPL1_S00_L001_R2_001.fastq.gz")

        for project in ["TestProj_1", "TestProj_2"]:
            genprep_dir = (
                f"{self.output_dir}/GenPrepFileJob/"
                f"211021_A00000_0000_SAMPLE/{project}/amplicon/"
            )
            makedirs(genprep_dir, exist_ok=True)
            copyfile(
                self.gz_source, f"{genprep_dir}/{runid}_SMPL1_S00_L001_R1_001.fastq.gz"
            )
            copyfile(
                self.gz_source, f"{genprep_dir}/{runid}_SMPL1_S00_L001_R2_001.fastq.gz"
            )

        fastqc_dir = f"{self.output_dir}/FastQCJob/logs/"
        multiqc_dir = f"{self.output_dir}/MultiQCJob/logs/"
        makedirs(fastqc_dir, exist_ok=True)
        makedirs(multiqc_dir, exist_ok=True)
        # generating the "*.completed" files
        for i in range(3):
            Path(f"{fastqc_dir}/FastQCJob_{i}.completed").touch()
            Path(f"{multiqc_dir}/MultiQCJob_{i}.completed").touch()

        # Amplicon is a valid data type in the default qiita test
        # database but job-id: 78901 doesn't exist; however, if we get
        # to here, it means that all steps have run to completion
        # and the system is trying to create the preps.
        # Note: Qiita job_id's are UUID in the database and this tests
        # uses 78901 as the job_id so the db complains about the format
        with self.assertRaisesRegex(
            RuntimeError, 'invalid input syntax for type uuid: "78901"'
        ):
            wf.execute_pipeline()

    def test_tellseq_workflow_creation(self):
        kwargs = {
            "uif_path": "tests/data/sample-sheets/metagenomic/tellseq/good_sheet1.csv",
            "qclient": self.qclient,
            "config_fp": "tests/configuration.json",
            "run_identifier": "211021_A00000_0000_SAMPLE",
            "output_dir": self.output_dir,
            "job_id": "78901",
            "lane_number": "1",
            "is_restart": False,
        }

        self._create_directory(kwargs["output_dir"])

        wf = WorkflowFactory.generate_workflow(**kwargs)

        # confirm that the proper type of workflow was generated.
        self.assertEqual(wf.protocol_type, PROTOCOL_NAME_TELLSEQ)
        self.assertEqual(wf.assay_type, ASSAY_NAME_METAGENOMIC)

        self._inject_data(wf)

        # Metagenomic is a valid data type in the default qiita test
        # database but job-id: 78901 doesn't exist; however, if we get
        # to here, it means that all steps have run to completion
        # and the system is trying to create the preps.
        # Note: Qiita job_id's are UUID in the database and this tests
        # uses 78901 as the job_id so the db complains about the format
        with self.assertRaisesRegex(
            RuntimeError, 'invalid input syntax for type uuid: "78901"'
        ):
            wf.execute_pipeline()


if __name__ == "__main__":
    main()
