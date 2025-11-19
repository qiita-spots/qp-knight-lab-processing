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
from unittest import TestCase, main

import metapool
import pandas as pd
from click.testing import CliRunner
from metapool.sample_sheet import PROTOCOL_NAME_PACBIO_SMRT
from qiita_client.testing import PluginTestCase

from qp_klp.Assays import ASSAY_NAME_METAGENOMIC
from qp_klp.scripts.pacbio_commands import generate_bam2fastq_commands
from qp_klp.WorkflowFactory import WorkflowFactory


class WorkflowBaseFactoryTests(TestCase):
    def setUp(self):
        self.base_dir = dirname(abspath(__file__))
        self.output_dir = join(self.base_dir, "test_output")
        self.remove_these = [self.output_dir]
        self._create_directory(self.output_dir)

        self.gz_source = f"{self.base_dir}/data/dummy.fastq.gz"

        metapool_base_dir = dirname(abspath(metapool.__file__))
        self.sample_sheets_base_dir = join(metapool_base_dir, "tests", "data")

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


class WorkflowFactoryTests(WorkflowBaseFactoryTests, PluginTestCase):
    def _inject_data(self, wf):
        """This is a helper method for testing that all the steps are run
        when testing wf.execute_pipeline()
        """
        samples = wf.pipeline.sample_sheet.samples

        convert_dir = f"{self.output_dir}/ConvertJob"
        reports_dir = f"{self.output_dir}/ConvertJob/Reports"
        makedirs(convert_dir, exist_ok=True)
        makedirs(reports_dir, exist_ok=True)
        Path(f"{convert_dir}/job_completed").touch()
        nuqc_dir = f"{self.output_dir}/NuQCJob"
        fastqc_dir = f"{self.output_dir}/FastQCJob/logs/"
        multiqc_dir = f"{self.output_dir}/MultiQCJob/logs/"
        genprep_dir = f"{self.output_dir}/GenPrepFileJob/r11111_20250101_111111/"
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
        for sample in samples:
            sn = sample["Sample_ID"]
            sp = sample["Sample_Project"]
            dstats.append({"SampleID": sn, "# Reads": 2})
            dname = f"{convert_dir}/{sp}"
            copyfile(self.gz_source, f"{dname}/{sn}_S000_L001_R1_001.fastq.gz")
            with open(f"{dname}/{sn}_S000_L001_R1_001.counts.txt", "w") as f:
                f.write("2")

            # NuQCJob
            dname = f"{nuqc_dir}/filtered_sequences/{sp}"
            copyfile(self.gz_source, f"{dname}/{sn}_S000_L001_R1_001.trimmed.fastq.gz")

            # GenPrepFileJob
            gprep_base = f"{genprep_dir}/{sp}/filtered_sequences/{sn}"
            copyfile(self.gz_source, f"{gprep_base}_S000_L001_R1_001.trimmed.fastq.gz")

        pd.DataFrame(dstats).set_index("SampleID").to_csv(
            f"{reports_dir}/Demultiplex_Stats.csv"
        )

        # generating the "*.completed" files
        for i in range(len(samples) * 3):
            Path(f"{fastqc_dir}/FastQCJob_{i}.completed").touch()
            Path(f"{multiqc_dir}/MultiQCJob_{i}.completed").touch()

    def test_pacbio_metagenomic_workflow_creation(self):
        kwargs = {
            "uif_path": f"{self.sample_sheets_base_dir}/good_pacbio_metagv10.csv",
            "qclient": self.qclient,
            "lane_number": "1",
            "config_fp": "tests/configuration.json",
            "run_identifier": "r11111_20250101_111111",
            "output_dir": self.output_dir,
            "job_id": "78901",
            "is_restart": False,
        }

        wf = WorkflowFactory.generate_workflow(**kwargs)
        # confirm that the proper type of workflow was generated.
        self.assertEqual(wf.protocol_type, PROTOCOL_NAME_PACBIO_SMRT)
        self.assertEqual(wf.assay_type, ASSAY_NAME_METAGENOMIC)

        self._inject_data(wf)

        # we can add some checks/tests of the initial scripts (mainly Convert)
        # but not doing now as is not required

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

    def test_pacbio_metagenomic_workflow_creation_twist_adaptors(self):
        kwargs = {
            "uif_path": f"{self.sample_sheets_base_dir}/good_pacbio_metagv11.csv",
            "qclient": self.qclient,
            "lane_number": "1",
            "config_fp": "tests/configuration.json",
            "run_identifier": "r11111_20250101_111111",
            "output_dir": self.output_dir,
            "job_id": "78901",
            "is_restart": False,
        }

        self._create_directory(kwargs["output_dir"])

        WorkflowFactory.generate_workflow(**kwargs)

        with open(f"{self.output_dir}/sample_list.tsv", "r") as fp:
            obs = fp.readlines()

        self.assertCountEqual(SAMPLE_LIST_11, obs)


class WorkflowFactoryTestsNoPlugin(WorkflowBaseFactoryTests):
    def _inject_data(self, add_files=False, remove_twisted=True):
        self.run_dir = f"{self.output_dir}/r11111_20250101_111111"
        self._create_directory(self.run_dir)

        if remove_twisted:
            _data = [line.rsplit("\t", 1)[0] + "\n" for line in SAMPLE_LIST_11]
        else:
            # note that we are actually adding both twisted and not twisted here
            _data = [line for line in SAMPLE_LIST_11]
            _data[-1] = SAMPLE_LIST_11[-1].rsplit("\t", 1)[0] + "\n"

        # making sure the header is intact
        _data[0] = SAMPLE_LIST_11[0]

        self.sample_list_fp = f"{self.output_dir}/sample_list.tsv"
        with open(self.sample_list_fp, "w") as fp:
            fp.write("".join(_data))

        if add_files:
            bd = f"{self.run_dir}/1_A01/hifi_reads/"
            self._create_directory(bd)
            if remove_twisted:
                for bc in ["bc3011", "bc0112", "bc9992"]:
                    Path(f"{bd}/mXXXX.hifi_reads.{bc}.bam").touch()
            else:
                # creating non twisted file
                Path(f"{bd}/mXXXX.hifi_reads.bc9992.bam").touch()
                # writting twisted
                bd = f"{self.run_dir}/1_A01/--MyProject_99999--/call-lima/execution/"
                makedirs(bd)
                for tb in [
                    "16_UDI_1_A01_F--16_UDI_1_A01_R",
                    "16_UDI_2_B01_F--16_UDI_2_B01_R",
                ]:
                    self._create_directory(f"{bd}/{tb}")
                    Path(f"{bd}/{tb}/mXXXX.hifi_reads.{tb}.bam").touch()

    def _execute_command(self):
        runner = CliRunner()
        result = runner.invoke(
            generate_bam2fastq_commands,
            [
                self.sample_list_fp,
                self.run_dir,
                self.output_dir,
            ],
        )
        return result

    def test_generate_bam2fastq_commands(self):
        # test missing barcodes
        self._inject_data()
        result = self._execute_command()
        self.assertRegex(
            str(result.exception),
            "r11111_20250101_111111 is missing barcodes: \['bc3011', 'bc0112', 'bc9992'\]",
        )

        # testing pacbio barcodes
        self._inject_data(add_files=True)
        result = self._execute_command()
        exp = [
            f"bam2fastq -j 1 -o {self.output_dir}/MyProject_99999/sample_1_S000_L001_R1_001 -c 9 "
            f"{self.run_dir}/1_A01/hifi_reads/mXXXX.hifi_reads.bc3011.bam; fqtools count "
            f"{self.output_dir}/MyProject_99999/sample_1_S000_L001_R1_001.fastq.gz > "
            f"{self.output_dir}/MyProject_99999/sample_1_S000_L001_R1_001.counts.txt",
            f"bam2fastq -j 1 -o {self.output_dir}/MyProject_99999/sample_2_S000_L001_R1_001 -c 9 "
            f"{self.run_dir}/1_A01/hifi_reads/mXXXX.hifi_reads.bc0112.bam; fqtools count "
            f"{self.output_dir}/MyProject_99999/sample_2_S000_L001_R1_001.fastq.gz > "
            f"{self.output_dir}/MyProject_99999/sample_2_S000_L001_R1_001.counts.txt",
            f"bam2fastq -j 1 -o {self.output_dir}/MyProject_99999/sample_3_S000_L001_R1_001 -c 9 "
            f"{self.run_dir}/1_A01/hifi_reads/mXXXX.hifi_reads.bc9992.bam; fqtools count "
            f"{self.output_dir}/MyProject_99999/sample_3_S000_L001_R1_001.fastq.gz > "
            f"{self.output_dir}/MyProject_99999/sample_3_S000_L001_R1_001.counts.txt",
            "",
        ]
        self.assertCountEqual(result.output.split("\n"), exp)

        # testing pacbio and twist adaptors combined
        self._inject_data(add_files=True, remove_twisted=False)
        result = self._execute_command()
        exp = [
            f"bam2fastq -j 1 -o {self.output_dir}/MyProject_99999/sample_1_S000_L001_R1_001 -c 9 "
            f"{self.run_dir}/1_A01/--MyProject_99999--/call-lima/execution/16_UDI_1_A01_F--16_UDI_1_A01_R/mXXXX.hifi_reads.16_UDI_1_A01_F--16_UDI_1_A01_R.bam; "
            f"fqtools count {self.output_dir}/MyProject_99999/sample_1_S000_L001_R1_001.fastq.gz > "
            f"{self.output_dir}/MyProject_99999/sample_1_S000_L001_R1_001.counts.txt",
            f"bam2fastq -j 1 -o {self.output_dir}/MyProject_99999/sample_2_S000_L001_R1_001 -c 9 "
            f"{self.run_dir}/1_A01/--MyProject_99999--/call-lima/execution/16_UDI_2_B01_F--16_UDI_2_B01_R/mXXXX.hifi_reads.16_UDI_2_B01_F--16_UDI_2_B01_R.bam; "
            f"fqtools count {self.output_dir}/MyProject_99999/sample_2_S000_L001_R1_001.fastq.gz > "
            f"{self.output_dir}/MyProject_99999/sample_2_S000_L001_R1_001.counts.txt",
            f"bam2fastq -j 1 -o {self.output_dir}/MyProject_99999/sample_3_S000_L001_R1_001 -c 9 "
            f"{self.run_dir}/1_A01/hifi_reads/mXXXX.hifi_reads.bc9992.bam; fqtools count "
            f"{self.output_dir}/MyProject_99999/sample_3_S000_L001_R1_001.fastq.gz > "
            f"{self.output_dir}/MyProject_99999/sample_3_S000_L001_R1_001.counts.txt",
            "",
        ]
        self.assertCountEqual(result.output.split("\n"), exp)


SAMPLE_LIST_11 = [
    "barcode\tsample_name\tproject_name\tlane\ttwist_adaptor_id\n",
    "bc3011\tsample_1\tMyProject_99999\t1\t16_UDI_1_A01_F--16_UDI_1_A01_R\n",
    "bc0112\tsample_2\tMyProject_99999\t1\t16_UDI_2_B01_F--16_UDI_2_B01_R\n",
    "bc9992\tsample_3\tMyProject_99999\t1\t16_UDI_5_E01_F--16_UDI_5_E01_R\n",
]

if __name__ == "__main__":
    main()
