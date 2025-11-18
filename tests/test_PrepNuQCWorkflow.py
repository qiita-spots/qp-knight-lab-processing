from functools import partial
from json import dumps
from os import environ, makedirs, remove
from os.path import basename, exists, isdir, join
from pathlib import Path
from shutil import copyfile, rmtree
from tempfile import mkdtemp
from unittest import main

from qiita_client.testing import PluginTestCase

from qp_klp import __version__ as version_qp_klp
from qp_klp import plugin
from qp_klp.klp import PrepNuQCJob


class TestPrepNuQCWorkflow(PluginTestCase):
    def setUp(self):
        plugin("https://localhost:21174", "register", "ignored")

        out_dir = mkdtemp()
        self.maxDiff = None
        self.out_dir = out_dir
        self._clean_up_files = []
        self._clean_up_files.append(out_dir)
        self.source_gz = "tests/data/dummy.fastq.gz"

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def _setup_test(self):
        prep_info_dict = {
            "SKB8.640193": {
                "run_prefix": "S22205_S104",
                "index": "GTCACTGT",
                "index2": "GCAAGATC",
            },
            "SKD8.640184": {
                "run_prefix": "S22282_S102",
                "index": "GCCATAAC",
                "index2": "AGTCTCAC",
            },
        }
        data = {
            "prep_info": dumps(prep_info_dict),
            # magic #1 = testing study
            "study": 1,
            "data_type": "Metagenomic",
        }
        pid = self.qclient.post("/apitest/prep_template/", data=data)["prep"]

        # inserting artifacts
        in_dir = mkdtemp()
        self._clean_up_files.append(in_dir)

        fps = [
            join(in_dir, "S22205_S104_L001_R1_001.fastq.gz"),
            join(in_dir, "S22205_S104_L001_R2_001.fastq.gz"),
            join(in_dir, "S22282_S102_L001_R1_001.fastq.gz"),
            join(in_dir, "S22282_S102_L001_R2_001.fastq.gz"),
        ]
        for fp in fps:
            copyfile(self.source_gz, fp)
        fp_summary = join(in_dir, "summary.html")
        copyfile("tests/data/summary.html", fp_summary)

        data = {
            "filepaths": dumps(
                [
                    (fps[0], "raw_forward_seqs"),
                    (fps[1], "raw_reverse_seqs"),
                    (fps[2], "raw_forward_seqs"),
                    (fps[3], "raw_reverse_seqs"),
                    (fp_summary, "html_summary"),
                ]
            ),
            "type": "per_sample_FASTQ",
            "name": "Test Woltka artifact",
            "prep": pid,
        }
        self.qclient.post("/apitest/artifact/", data=data)["artifact"]

        data = {
            "user": "demo@microbio.me",
            "command": dumps(
                ["qp-klp", version_qp_klp, "Human Filter & QC existing Prep"]
            ),
            "status": "running",
            "parameters": dumps({"prep_id": pid}),
        }
        job_id = self.qclient.post("/apitest/processing_job/", data=data)["job"]

        return pid, job_id, fps

    def test_prep_nuqcjob(self):
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # test failure
        jid = "bcc7ebcd-39c1-43e4-af2d-822e3589f14d"
        success, ainfo, msg = PrepNuQCJob(self.qclient, jid, {"prep_id": "1"}, out_dir)
        self.assertEqual(msg, "Prep 1 has a not valid data type: 18S")
        self.assertFalse(success)

        # let's test without setting the env PrepNuQCJob_TEST; the workflow
        # will fail due to FastQCJob
        pid, job_id, _ = self._setup_test()
        success, ainfo, msg = PrepNuQCJob(
            self.qclient, job_id, {"prep_id": pid}, out_dir
        )
        self.assertRegex(msg, r"There are no fastq files for FastQCJob")
        self.assertFalse(success)

        # now, let's test all the way to the end by using PrepNuQCJob_TEST
        pid, job_id, fps = self._setup_test()
        environ["PrepNuQCJob_TEST"] = "true"

        # inserting halfway files
        out_path = partial(join, out_dir)
        project_names = ["qiita-3-10_1", "qiita-4-11_1"]
        run_identifier = "250225_LH00444_0301_B22N7T2LT4"
        for project_name in project_names:
            nuqc_dir = out_path("NuQCJob", project_name, "filtered_sequences")
            makedirs(nuqc_dir, exist_ok=True)
            fastqc_dir = out_path("FastQCJob")
            makedirs(f"{fastqc_dir}/logs", exist_ok=True)
            genprep_dir = out_path(
                "GenPrepFileJob", run_identifier, project_name, "filtered_sequences"
            )
            makedirs(genprep_dir)

            for f in fps:
                bn = basename(f).replace(".fastq.gz", ".trimmed.fastq.gz")
                copyfile(self.source_gz, f"{nuqc_dir}/{bn}")
                copyfile(self.source_gz, f"{genprep_dir}/{bn}")
        for i in range(1, 7):
            Path(f"{fastqc_dir}/logs/FastQCJob_{i}.completed").touch()
        Path(f"{fastqc_dir}/job_completed").touch()

        success, ainfo, msg = PrepNuQCJob(
            self.qclient, job_id, {"prep_id": pid}, out_dir
        )
        environ["PrepNuQCJob_TEST"] = ""
        # we'll consider this error a success as it is failing during the
        # new prep insert as that user doesn't exist
        self.assertRegex(msg, r"ID 'qiita.help@gmail.com' does not exists")
        self.assertFalse(success)


if __name__ == "__main__":
    main()
