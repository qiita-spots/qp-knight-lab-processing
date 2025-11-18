# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from copy import deepcopy
from os import getcwd, makedirs, remove
from os.path import abspath, exists, join
from shutil import rmtree
from tempfile import TemporaryDirectory
from unittest import TestCase, main

from qp_klp.FailedSamplesRecord import FailedSamplesRecord
from qp_klp.WorkflowFactory import WorkflowFactory
from qp_klp.Workflows import WorkflowError


class FakeClient:
    def __init__(self):
        self.cwd = getcwd()
        self.base_path = join(self.cwd, "tests/data/QDir")
        self.qdirs = {
            "Demultiplexed": "Demultiplexed",
            "beta_div_plots": "analysis/beta_div_plots",
            "rarefaction_curves": "analysis/rarefaction_curves",
            "taxa_summary": "analysis/taxa_summary",
            "q2_visualization": "working_dir",
            "distance_matrix": "working_dir",
            "ordination_results": "working_dir",
            "alpha_vector": "working_dir",
            "FASTQ": "FASTQ",
            "BIOM": "BIOM",
            "per_sample_FASTQ": "per_sample_FASTQ",
            "SFF": "SFF",
            "FASTA": "FASTA",
            "FASTA_Sanger": "FASTA_Sanger",
            "FeatureData": "FeatureData",
            "job-output-folder": "job-output-folder",
            "BAM": "BAM",
            "VCF": "VCF",
            "SampleData": "SampleData",
            "uploads": "uploads",
        }

        self.samples_in_13059 = [
            "13059.SP331130A04",
            "13059.AP481403B02",
            "13059.LP127829A02",
            "13059.BLANK3.3B",
            "13059.EP529635B02",
            "13059.EP542578B04",
            "13059.EP446602B01",
            "13059.EP121011B01",
            "13059.EP636802A01",
            "13059.SP573843A04",
        ]

        # note these samples have known tids, but aren't in good-sample-sheet.
        self.samples_in_11661 = [
            "11661.1.24",
            "11661.1.57",
            "11661.1.86",
            "11661.10.17",
            "11661.10.41",
            "11661.10.64",
            "11661.11.18",
            "11661.11.43",
            "11661.11.64",
            "11661.12.15",
        ]

        self.samples_in_6123 = [
            "3A",
            "4A",
            "5B",
            "6A",
            "BLANK.41.12G",
            "7A",
            "8A",
            "ISB",
            "GFR",
            "6123",
        ]

        self.info_in_11661 = {
            "number-of-samples": 10,
            "categories": ["sample_type", "tube_id"],
        }

        self.info_in_13059 = {
            "number-of-samples": 10,
            "categories": [
                "anonymized_name",
                "collection_timestamp",
                "description",
                "dna_extracted",
                "elevation",
                "empo_1",
                "empo_2",
                "empo_3",
                "env_biome",
                "env_feature",
                "env_material",
                "env_package",
                "geo_loc_name",
                "host_age",
                "host_age_units",
                "host_body_habitat",
                "host_body_mass_index",
                "host_body_product",
                "host_body_site",
                "host_common_name",
                "host_height",
                "host_height_units",
                "host_life_stage",
                "host_scientific_name",
                "host_subject_id",
                "host_taxid",
                "host_weight",
                "host_weight_units",
                "latitude",
                "longitude",
                "nyuid",
                "physical_specimen_location",
                "physical_specimen_remaining",
                "predose_time",
                "sample_type",
                "scientific_name",
                "sex",
                "subject_id",
                "taxon_id",
                "title",
                "tube_id",
            ],
        }

        # Study not in qiita-rc. Faking results.
        self.info_in_6123 = {
            "number-of-samples": 10,
            "categories": ["sample_type", "subject_id", "title"],
        }

        self.tids_13059 = {
            "header": ["tube_id"],
            "samples": {
                "13059.SP331130A04": ["SP331130A-4"],
                "13059.AP481403B02": ["AP481403B-2"],
                "13059.LP127829A02": ["LP127829A-2"],
                "13059.BLANK3.3B": ["BLANK3.3B"],
                "13059.EP529635B02": ["EP529635B-2"],
                "13059.EP542578B04": ["EP542578B.4"],
                "13059.EP446602B01": ["EP446602B.1"],
                "13059.EP121011B01": ["EP121011B.1"],
                "13059.EP636802A01": ["EP636802A.1"],
                "13059.SP573843A04": ["SP573843A.4"],
            },
        }

        self.tids_11661 = {
            "header": ["tube_id"],
            "samples": {
                "11661.1.24": ["1.24"],
                "11661.1.57": ["1.57"],
                "11661.1.86": ["1.86"],
                "11661.10.17": ["10.17"],
                "11661.10.41": ["10.41"],
                "11661.10.64": ["10.64"],
                "11661.11.18": ["11.18"],
                "11661.11.43": ["11.43"],
                "11661.11.64": ["11.64"],
                "11661.12.15": ["12.15"],
            },
        }

        for key in self.qdirs:
            self.qdirs[key] = join(self.base_path, self.qdirs[key])

        for qdir in self.qdirs:
            makedirs(self.qdirs[qdir], exist_ok=True)

        self.fake_id = 1000
        self._server_url = "some.server.url"
        self.saved_posts = {}

    def get(self, url):
        m = {
            "/api/v1/study/11661/samples": self.samples_in_11661,
            "/api/v1/study/11661/samples/categories=tube_id": self.tids_11661,
            "/api/v1/study/11661/samples/info": self.info_in_11661,
            "/api/v1/study/13059/samples": self.samples_in_13059,
            "/api/v1/study/13059/samples/categories=tube_id": self.tids_13059,
            "/api/v1/study/13059/samples/info": self.info_in_13059,
            "/api/v1/study/6123/samples": self.samples_in_6123,
            "/api/v1/study/6123/samples/info": self.info_in_6123,
            "/qiita_db/artifacts/types/": self.qdirs,
        }

        if url in m:
            return m[url]

        return None

    def post(self, url, data=None):
        if "/qiita_db/prep_template/" == url:
            self.fake_id += 1
            return {"prep": self.fake_id}
        elif "/qiita_db/artifact/" == url:
            self.saved_posts[str(self.fake_id)] = data
            self.fake_id += 1
            return {"job_id": self.fake_id}
        else:
            raise ValueError("Unsupported URL")


class AnotherFakeClient:
    def __init__(self):
        self.cwd = getcwd()
        self.base_path = join(self.cwd, "tests/data/QDir")
        self.qdirs = {
            "Demultiplexed": "Demultiplexed",
            "beta_div_plots": "analysis/beta_div_plots",
            "rarefaction_curves": "analysis/rarefaction_curves",
            "taxa_summary": "analysis/taxa_summary",
            "q2_visualization": "working_dir",
            "distance_matrix": "working_dir",
            "ordination_results": "working_dir",
            "alpha_vector": "working_dir",
            "FASTQ": "FASTQ",
            "BIOM": "BIOM",
            "per_sample_FASTQ": "per_sample_FASTQ",
            "SFF": "SFF",
            "FASTA": "FASTA",
            "FASTA_Sanger": "FASTA_Sanger",
            "FeatureData": "FeatureData",
            "job-output-folder": "job-output-folder",
            "BAM": "BAM",
            "VCF": "VCF",
            "SampleData": "SampleData",
            "uploads": "uploads",
        }

        self.samples_in_99999 = [
            "99999.AAAAAAAAAAA",
            "99999.BBBBBBBBBBB",
            "99999.CCCCCCCCCCC",
            "99999.BLANK1.1BCD",
        ]

        self.info_in_99999 = {
            "number-of-samples": 10,
            "categories": ["column1", "column2", "tube_id"],
        }

        self.tids_99999 = {
            "header": ["tube_id"],
            "samples": {
                "99999.AAAAAAAAAAA": ["1234567890a"],
                "99999.BBBBBBBBBBB": ["234567890ab"],
                "99999.CCCCCCCCCCC": ["34567890abc"],
                "99999.BLANK1.1BCD": ["BLANK1.1BCD"],
            },
        }

        for key in self.qdirs:
            self.qdirs[key] = join(self.base_path, self.qdirs[key])

        for qdir in self.qdirs:
            makedirs(self.qdirs[qdir], exist_ok=True)

    def get(self, url):
        m = {
            "/api/v1/study/99999/samples": self.samples_in_99999,
            "/api/v1/study/99999/samples/categories=tube_id": self.tids_99999,
            "/api/v1/study/99999/samples/info": self.info_in_99999,
            "/qiita_db/artifacts/types/": self.qdirs,
        }

        if url in m:
            return m[url]

        return None


class TestHelpers(TestCase):
    def setUp(self):
        self.fake_bin_path = ""
        self.delete_these_files = []
        self.delete_these_dirs = []
        # self.fake_bin_path = self.get_searchable_path()

        # self.output_dir represents a qiita working directory.
        package_root = abspath("../tests/data")
        self.output_dir = join(package_root, "077c4da8-74eb-4184-8860-0207f53623be")
        self.delete_these_dirs = [self.output_dir]

        # We want a clean directory, nothing from leftover runs
        makedirs(self.output_dir, exist_ok=False)

        self.debug = False
        self.fake_client = FakeClient()
        self.kwargs = {
            "uif_path": "tests/data/sample-sheets/metagenomic/illumina/good_sheet1.csv",
            "qclient": self.fake_client,
            "lane_number": "1",
            "config_fp": "tests/configuration.json",
            "run_identifier": "211021_A00000_0000_SAMPLE",
            "output_dir": self.output_dir,
            "job_id": "077c4da8-74eb-4184-8860-0207f53623be",
            "is_restart": False,
        }

    def tearDown(self):
        if not self.debug:
            rmtree(self.output_dir)

            for fp in self.delete_these_files:
                if exists(fp):
                    remove(fp)

    def test_generate_special_map(self):
        wf = WorkflowFactory.generate_workflow(**self.kwargs)
        wf.generate_special_map()

        obs = wf.special_map
        exp = [
            (
                "StudyA_13059",
                join(self.fake_client.base_path, "uploads/13059"),
                "13059",
            ),
            (
                "StudyB_11661",
                join(self.fake_client.base_path, "uploads/11661"),
                "11661",
            ),
            ("StudyC_6123", join(self.fake_client.base_path, "uploads/6123"), "6123"),
        ]

        self.assertEqual(obs, exp)

    def test_get_project_info(self):
        wf = WorkflowFactory.generate_workflow(**self.kwargs)
        wf.generate_special_map()
        obs = wf.pipeline.get_project_info()

        exp = [
            {
                "project_name": "StudyA_13059",
                "qiita_id": "13059",
                "contains_replicates": False,
            },
            {
                "project_name": "StudyB_11661",
                "qiita_id": "11661",
                "contains_replicates": False,
            },
            {
                "project_name": "StudyC_6123",
                "qiita_id": "6123",
                "contains_replicates": False,
            },
        ]

        self.assertEqual(obs, exp)

    def test_get_samples_in_qiita(self):
        wf = WorkflowFactory.generate_workflow(**self.kwargs)
        obs_samples, obs_tids = wf.get_samples_in_qiita(self.fake_client, "13059")

        exp_samples = {
            "EP121011B01",
            "EP529635B02",
            "EP542578B04",
            "SP573843A04",
            "SP331130A04",
            "EP446602B01",
            "BLANK3.3B",
            "AP481403B02",
            "LP127829A02",
            "EP636802A01",
        }

        exp_tids = {
            "13059.SP331130A04": ["SP331130A-4"],
            "13059.AP481403B02": ["AP481403B-2"],
            "13059.LP127829A02": ["LP127829A-2"],
            "13059.BLANK3.3B": ["BLANK3.3B"],
            "13059.EP529635B02": ["EP529635B-2"],
            "13059.EP542578B04": ["EP542578B.4"],
            "13059.EP446602B01": ["EP446602B.1"],
            "13059.EP121011B01": ["EP121011B.1"],
            "13059.EP636802A01": ["EP636802A.1"],
            "13059.SP573843A04": ["SP573843A.4"],
        }

        self.assertEqual(obs_samples, exp_samples)
        self.assertDictEqual(obs_tids, exp_tids)

    def test_get_tube_ids_from_qiita(self):
        wf = WorkflowFactory.generate_workflow(**self.kwargs)
        wf._get_tube_ids_from_qiita()

        obs = wf.tube_id_map
        exp = {
            "13059": {
                "SP331130A04": "SP331130A-4",
                "AP481403B02": "AP481403B-2",
                "LP127829A02": "LP127829A-2",
                "BLANK3.3B": "BLANK3.3B",
                "EP529635B02": "EP529635B-2",
                "EP542578B04": "EP542578B.4",
                "EP446602B01": "EP446602B.1",
                "EP121011B01": "EP121011B.1",
                "EP636802A01": "EP636802A.1",
                "SP573843A04": "SP573843A.4",
            },
            "11661": {
                "1.24": "1.24",
                "1.57": "1.57",
                "1.86": "1.86",
                "10.17": "10.17",
                "10.41": "10.41",
                "10.64": "10.64",
                "11.18": "11.18",
                "11.43": "11.43",
                "11.64": "11.64",
                "12.15": "12.15",
            },
        }

        self.assertDictEqual(obs, exp)

    def test_project_metadata_check(self):
        wf = WorkflowFactory.generate_workflow(**self.kwargs)

        # _project_metadata_check() should return w/out raising an Error if
        # step and fake_client is used.
        wf._project_metadata_check()

        self.fake_client.info_in_11661["categories"].append("sample_well")
        self.fake_client.info_in_13059["categories"].append("sample_well")

        msg = (
            "'sample_well' exists in Qiita study 13059's sample metadata"
            "\n'sample_well' exists in Qiita study 11661's sample metadata"
        )
        with self.assertRaisesRegex(WorkflowError, msg):
            wf._project_metadata_check()


class FailedSamplesRecordTests(TestCase):
    def setUp(self):
        class MockSample:
            def __init__(self, sample_id, project_name):
                self.Sample_ID = sample_id
                self.Sample_Project = project_name

        self.samples = [
            MockSample("A", "ProjectOne"),
            MockSample("B", "ProjectTwo"),
            MockSample("C", "ProjectThree"),
            MockSample("D", "ProjectFour"),
        ]

        self.output = TemporaryDirectory()

    def test_failed_samples_record(self):
        fsr = FailedSamplesRecord(self.output.name, self.samples)

        # assert that a state file doesn't already exist and attempt to load()
        # it. State should remain unchanged.
        exp = deepcopy(fsr.sample_state)
        self.assertFalse(exists(fsr.output_path))
        fsr.load()
        self.assertEqual(fsr.sample_state, exp)

        # confirm that dump() creates the appropriate file.
        self.assertFalse(exists(fsr.output_path))
        fsr.dump()
        self.assertTrue(exists(fsr.output_path))

        # load the dumped() file, and confirm that nothing changed since
        # the state wasn't update()d.
        fsr.load()
        self.assertEqual(fsr.sample_state, exp)

        # assert samples A and C failed at the ConvertJob stage, assert
        # state changed. dump() state and load() it. Confirm state on disk
        # reflects changes.
        fsr.update(["A", "C"], "ConvertJob")
        self.assertNotEqual(fsr.sample_state, exp)
        fsr.dump()
        fsr.load()
        exp = {"A": "ConvertJob", "B": None, "C": "ConvertJob", "D": None}
        self.assertEqual(fsr.sample_state, exp)

        # B should be failed at FastQCJob but A and C should still be
        # failed at ConvertJob
        fsr.update(["A", "B", "C"], "FastQCJob")
        exp = {"A": "ConvertJob", "B": "FastQCJob", "C": "ConvertJob", "D": None}
        self.assertEqual(fsr.sample_state, exp)

        # confirm html file exists. Assume pandas DataFrame.to_html() works
        # as intended.
        self.assertFalse(exists(fsr.report_path))
        fsr.generate_report()
        self.assertTrue(exists(fsr.report_path))


if __name__ == "__main__":
    main()
