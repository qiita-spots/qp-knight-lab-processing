# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from unittest import TestCase
from sequence_processing_pipeline.Pipeline import Pipeline, PipelineError
from os.path import join, abspath, exists, dirname, split
from functools import partial
from os import makedirs, chmod, access, W_OK, walk
from shutil import rmtree, copy, which, copytree
from os import environ, remove, getcwd
import pandas as pd
from metapool import parse_prep
import re
from copy import deepcopy
from tempfile import TemporaryDirectory
from qp_klp.Workflows import StandardMetagenomicWorkflow, WorkflowFactory
from metapool import load_sample_sheet
from collections import defaultdict
from hashlib import md5
from random import randint


class FakeClient():
    def __init__(self):
        self.cwd = getcwd()
        self.base_path = join(self.cwd, 'qp_klp/tests/data/QDir')
        self.qdirs = {'Demultiplexed': 'Demultiplexed',
                      'beta_div_plots': 'analysis/beta_div_plots',
                      'rarefaction_curves': 'analysis/rarefaction_curves',
                      'taxa_summary': 'analysis/taxa_summary',
                      'q2_visualization': 'working_dir',
                      'distance_matrix': 'working_dir',
                      'ordination_results': 'working_dir',
                      'alpha_vector': 'working_dir',
                      'FASTQ': 'FASTQ',
                      'BIOM': 'BIOM',
                      'per_sample_FASTQ': 'per_sample_FASTQ',
                      'SFF': 'SFF',
                      'FASTA': 'FASTA',
                      'FASTA_Sanger': 'FASTA_Sanger',
                      'FeatureData': 'FeatureData',
                      'job-output-folder': 'job-output-folder',
                      'BAM': 'BAM',
                      'VCF': 'VCF',
                      'SampleData': 'SampleData',
                      'uploads': 'uploads'}

        self.samples_in_13059 = ['13059.SP331130A04', '13059.AP481403B02',
                                 '13059.LP127829A02', '13059.BLANK3.3B',
                                 '13059.EP529635B02', '13059.EP542578B04',
                                 '13059.EP446602B01', '13059.EP121011B01',
                                 '13059.EP636802A01', '13059.SP573843A04']

        # note these samples have known tids, but aren't in good-sample-sheet.
        self.samples_in_11661 = ['11661.1.24', '11661.1.57', '11661.1.86',
                                 '11661.10.17', '11661.10.41', '11661.10.64',
                                 '11661.11.18', '11661.11.43', '11661.11.64',
                                 '11661.12.15']

        self.samples_in_6123 = ['3A', '4A', '5B', '6A', 'BLANK.41.12G', '7A',
                                '8A', 'ISB', 'GFR', '6123']

        self.info_in_11661 = {'number-of-samples': 10,
                              'categories': ['sample_type', 'tube_id']}

        self.info_in_13059 = {'number-of-samples': 10,
                              'categories': ['anonymized_name',
                                             'collection_timestamp',
                                             'description',
                                             'dna_extracted',
                                             'elevation', 'empo_1',
                                             'empo_2', 'empo_3',
                                             'env_biome', 'env_feature',
                                             'env_material',
                                             'env_package',
                                             'geo_loc_name', 'host_age',
                                             'host_age_units',
                                             'host_body_habitat',
                                             'host_body_mass_index',
                                             'host_body_product',
                                             'host_body_site',
                                             'host_common_name',
                                             'host_height',
                                             'host_height_units',
                                             'host_life_stage',
                                             'host_scientific_name',
                                             'host_subject_id',
                                             'host_taxid', 'host_weight',
                                             'host_weight_units',
                                             'latitude', 'longitude',
                                             'nyuid',
                                             'physical_specimen_location',
                                             'physical_specimen_remaining',
                                             'predose_time',
                                             'sample_type',
                                             'scientific_name', 'sex',
                                             'subject_id', 'taxon_id',
                                             'title', 'tube_id']}

        # Study not in qiita-rc. Faking results.
        self.info_in_6123 = {'number-of-samples': 10,
                             'categories': ['sample_type', 'subject_id',
                                            'title']}

        self.tids_13059 = {"header": ["tube_id"],
                           "samples": {'13059.SP331130A04': ['SP331130A-4'],
                                       '13059.AP481403B02': ['AP481403B-2'],
                                       '13059.LP127829A02': ['LP127829A-2'],
                                       '13059.BLANK3.3B': ['BLANK3.3B'],
                                       '13059.EP529635B02': ['EP529635B-2'],
                                       '13059.EP542578B04': ['EP542578B-4'],
                                       '13059.EP446602B01': ['EP446602B-1'],
                                       '13059.EP121011B01': ['EP121011B-1'],
                                       '13059.EP636802A01': ['EP636802A-1'],
                                       '13059.SP573843A04': ['SP573843A-4']}}

        self.tids_11661 = {"header": ["tube_id"],
                           "samples": {"11661.1.24": ["1.24"],
                                       "11661.1.57": ["1.57"],
                                       "11661.1.86": ["1.86"],
                                       "11661.10.17": ["10.17"],
                                       "11661.10.41": ["10.41"],
                                       "11661.10.64": ["10.64"],
                                       "11661.11.18": ["11.18"],
                                       "11661.11.43": ["11.43"],
                                       "11661.11.64": ["11.64"],
                                       "11661.12.15": ["12.15"]}}

        for key in self.qdirs:
            self.qdirs[key] = join(self.base_path, self.qdirs[key])

        for qdir in self.qdirs:
            makedirs(self.qdirs[qdir], exist_ok=True)

        self.fake_id = 1000
        self._server_url = "some.server.url"
        self.saved_posts = {}

    def get(self, url):
        m = {'/api/v1/study/11661/samples': self.samples_in_11661,
             '/api/v1/study/11661/samples/categories=tube_id': self.tids_11661,
             '/api/v1/study/11661/samples/info': self.info_in_11661,
             '/api/v1/study/13059/samples': self.samples_in_13059,
             '/api/v1/study/13059/samples/categories=tube_id': self.tids_13059,
             '/api/v1/study/13059/samples/info': self.info_in_13059,
             '/api/v1/study/6123/samples': self.samples_in_6123,
             '/api/v1/study/6123/samples/info': self.info_in_6123,
             '/qiita_db/artifacts/types/': self.qdirs}

        if url in m:
            return m[url]

        return None

    def post(self, url, data=None):
        if '/qiita_db/prep_template/' == url:
            self.fake_id += 1
            return {'prep': self.fake_id}
        elif '/qiita_db/artifact/' == url:
            self.saved_posts[str(self.fake_id)] = data
            self.fake_id += 1
            return {'job_id': self.fake_id}
        else:
            raise ValueError("Unsupported URL")


class AnotherFakeClient():
    def __init__(self):
        self.cwd = getcwd()
        self.base_path = join(self.cwd, 'qp_klp/tests/data/QDir')
        self.qdirs = {'Demultiplexed': 'Demultiplexed',
                      'beta_div_plots': 'analysis/beta_div_plots',
                      'rarefaction_curves': 'analysis/rarefaction_curves',
                      'taxa_summary': 'analysis/taxa_summary',
                      'q2_visualization': 'working_dir',
                      'distance_matrix': 'working_dir',
                      'ordination_results': 'working_dir',
                      'alpha_vector': 'working_dir',
                      'FASTQ': 'FASTQ',
                      'BIOM': 'BIOM',
                      'per_sample_FASTQ': 'per_sample_FASTQ',
                      'SFF': 'SFF',
                      'FASTA': 'FASTA',
                      'FASTA_Sanger': 'FASTA_Sanger',
                      'FeatureData': 'FeatureData',
                      'job-output-folder': 'job-output-folder',
                      'BAM': 'BAM',
                      'VCF': 'VCF',
                      'SampleData': 'SampleData',
                      'uploads': 'uploads'}

        self.samples_in_99999 = ['99999.AAAAAAAAAAA',
                                 '99999.BBBBBBBBBBB',
                                 '99999.CCCCCCCCCCC',
                                 '99999.BLANK1.1BCD']

        self.info_in_99999 = {'number-of-samples': 10,
                              'categories': ['column1', 'column2', 'tube_id']}

        self.tids_99999 = {"header": ["tube_id"],
                           "samples": {'99999.AAAAAAAAAAA': ['1234567890a'],
                                       '99999.BBBBBBBBBBB': ['234567890ab'],
                                       '99999.CCCCCCCCCCC': ['34567890abc'],
                                       '99999.BLANK1.1BCD': ['BLANK1.1BCD']}}

        for key in self.qdirs:
            self.qdirs[key] = join(self.base_path, self.qdirs[key])

        for qdir in self.qdirs:
            makedirs(self.qdirs[qdir], exist_ok=True)

    def get(self, url):
        m = {'/api/v1/study/99999/samples': self.samples_in_99999,
             '/api/v1/study/99999/samples/categories=tube_id': self.tids_99999,
             '/api/v1/study/99999/samples/info': self.info_in_99999,
             '/qiita_db/artifacts/types/': self.qdirs}

        if url in m:
            return m[url]

        return None


class Basics(TestCase):
    def setUp(self):
        self.fake_bin_path = ""
        self.delete_these_files = []
        self.delete_these_dirs = []
        self.fake_bin_path = self.get_searchable_path()

        # self.output_dir represents a qiita working directory.
        package_root = abspath('./qp_klp/tests/data')
        self.output_dir = join(package_root, "077c4da8-74eb-4184-8860-0207f53623be")
        self.delete_these_dirs = [self.output_dir]

        # We want a clean directory, nothing from leftover runs
        makedirs(self.output_dir, exist_ok=False)

    def tearDown(self):
        for fp in self.delete_these_files:
            if exists(fp):
                #remove(fp)
                pass

    def create_fake_bin(self, name, content):
        tmp = join(self.fake_bin_path, name)
        with open(tmp, 'w') as f:
            f.write(f"#!/bin/sh\n{content}\n")
        chmod(tmp, 0o777)
        self.delete_these_files.append(tmp)
        return tmp

    def create_fake_file(self, fp):
        with open(fp, 'w') as f:
            f.write("This is a file.")

    def get_searchable_path(self):
        searchable_paths = []

        if 'CONDA_PREFIX' in environ:
            # create fake binaries in bin directory of Conda environment
            searchable_paths.append(environ['CONDA_PREFIX'] + '/bin')
        else:
            # if CONDA_PREFIX doesn't exist, select a path from a list of
            # searchable paths that contains 'env' and assume it's writable.
            tmp = environ['PATH']
            searchable_paths += tmp.split(':')

        for a_path in searchable_paths:
            if access(a_path, W_OK):
                return a_path

    def test_partial_pipeline(self):
        # Tests convert_raw_to_fastq() and quality_control() steps of
        # StandardMetagenomicWorkflow(), which in turn exercises
        # FailedSamplesRecord, Metagenomic and Illumina mixins, and the
        # base Workflow class.

        # create a shell-script mimicking what the real sbatch prints to
        # stdout. The code in Job() will extract the job-id (9999999) which
        # it will poll for using our fake squeue script.
        #
        # after printing to stdout, our fake sbatch will simulate the results
        # of bcl-convert generating fastqs by creating a directory full of
        # fastq files for each project.

        # the line returning the fake Slurm job-id.
        cmds = ["echo 'Submitted batch job 9999999'"]

        # the lines to recreate the directories a standard Job() object creates.
        cmds.append("mkdir -p %s" % join(self.output_dir, 'ConvertJob', 'logs'))
        cmds.append("mkdir -p %s" % join(self.output_dir, 'ConvertJob', 'Reports'))

        # use the list of sample_ids found in the sample-sheet to generate
        # fake fastq files for convert_raw_to_fastq() to find and manipulate.
        sheet = load_sample_sheet("qp_klp/tests/data/sample-sheets/metagenomic"
                                  "/illumina/good_sheet1.csv")
        exp = defaultdict(list)
        for sample in sheet.samples:
            sample = sample.to_json()
            exp[sample['Sample_Project']].append(sample['Sample_ID'])

        # in order to test post-process auditing, we will manually remove a
        # single sample_id from each project in order to simulate failed
        # conversions.
        simulated_failed_samples = []

        for project in exp:
            # sort the list so that files are created in a predictable order.
            exp[project].sort()

            simulated_failed_samples.append(exp[project].pop())

            fake_path = join(self.output_dir, 'ConvertJob', project)
            cmds.append(f"mkdir -p {fake_path}")

            for sample in exp[project]:
                r1 = join(fake_path, f'{sample}_S123_L001_R1_001.fastq.gz')
                r2 = join(fake_path, f'{sample}_S123_L001_R2_001.fastq.gz')

                # let r1 and r2 be the same size.
                size = randint(1, 10)
                for file_path in [r1, r2]:
                    cmds.append(f"dd if=/dev/zero of={file_path} bs={size}MB count=1 2>/dev/null")

        # write all the statements out into a bash-script named 'sbatch' and
        # place it somewhere in the PATH. (It will be removed on tearDown()).
        self.create_fake_bin('sbatch', "\n".join(cmds))

        # create fake squeue binary that writes to stdout what Job() needs to
        # see to think that bcl-convert has completed successfully.
        self.create_fake_bin('squeue', "echo 'ARRAY_JOB_ID,JOBID,STATE\n"
                              "9999999,9999999,COMPLETED'")

        # Create a Workflow object using WorkflowFactory(). No need to confirm
        # it is the correct one; that is confirmed in other tests. Run
        # convert_raw_to_fastq() and confirm that the directory structure and
        # the audit results match what is expected.
        kwargs = {"uif_path": "qp_klp/tests/data/sample-sheets/metagenomic/"
                  "illumina/good_sheet1.csv",
                  "qclient": FakeClient(),
                  "lane_number": "1",
                  "config_fp": "qp_klp/tests/data/configuration.json",
                  "run_identifier": '211021_A00000_0000_SAMPLE',
                  "output_dir": self.output_dir,
                  "job_id": "077c4da8-74eb-4184-8860-0207f53623be",
                  "is_restart": False
                  }

        wf = WorkflowFactory.generate_workflow(**kwargs)

        # Illumina.convert_raw_to_fastq() calls Job.audit() after bcl-convert
        # exists and will identify the samples that failed to process. Confirm
        # the values are correct here.
        audit_results = sorted(wf.convert_raw_to_fastq())

        self.assertEqual(audit_results, sorted(simulated_failed_samples))

        # NB: bcl-convert's presence in ConvertJob.sh confirms it is getting
        # path and binary name from the correct configuration file. Note it
        # doesn't check to see that the binary exists because where Job() runs
        # is not the same location as where bcl-convert will run (compute-
        # node.)

        # confirm ConvertJob.sh Slurm job script looks exactly as intended by
        # confirming its digest.

        # TODO: Redo this test since full path is now present
        #exp = "8b3354068d30142ff4158cd049808afe"
        #script_path = join(self.output_dir, 'ConvertJob', 'ConvertJob.sh')
        #obs = md5(open(script_path, 'rb').read()).hexdigest()
        #self.assertEqual(obs, exp)

        # ConvertJob successful.
        cmds = []

        # the lines to recreate the directories a standard Job() object creates.
        cmds.append("mkdir -p %s" % join(self.output_dir, 'NuQCJob', 'logs'))
        cmds.append("mkdir -p %s" % join(self.output_dir, 'NuQCJob', 'only-adapter-filtered'))
        cmds.append("mkdir -p %s" % join(self.output_dir, 'NuQCJob', 'fastp_reports_dir', 'html'))
        cmds.append("mkdir -p %s" % join(self.output_dir, 'NuQCJob', 'fastp_reports_dir', 'json'))
        cmds.append("mkdir -p %s" % join(self.output_dir, 'NuQCJob', 'tmp'))
        cmds.append("mkdir -p %s" % join(self.output_dir, 'NuQCJob', 'tmp.564341'))

        #fake_path = join(self.output_dir, 'ConvertJob', project)

        # simulate host-filtering scripts by scanning contents of ConvertJob
        # directory and copying the files into the expected location for post-
        # processing.
        for root, dirs, files in walk(join(self.output_dir, 'ConvertJob')):
            for _file in files:
                # don't process anything from ConvertJob directory that isn't
                # a simulated fastq file. Since these are fake files we can
                # assume 'Undetermined' fastq files are not present.
                if not _file.endswith('.fastq.gz'):
                    continue

                raw_file = join(root, _file)
                _, project_name = split(root)

                new_name = _file.replace('.fastq.gz', '.interleave.fastq.gz')
                file_path = join(self.output_dir, 'NuQCJob', 'only-adapter-filtered', new_name)
                cmds.append(f"cp {raw_file} {file_path}")

                cmds.append("mkdir -p %s" % join(self.output_dir, 'NuQCJob', project_name, 'filtered_sequences'))

                file_path = join(self.output_dir, 'NuQCJob', project_name, 'filtered_sequences', _file)
                cmds.append(f"cp {raw_file} {file_path}")

                new_name = _file.replace('.fastq.gz', '.html')
                file_path = join(self.output_dir, 'NuQCJob', 'fastp_reports_dir', 'html', new_name)
                cmds.append(f"echo 'This is an html file.' > {file_path}")

                new_name = _file.replace('.fastq.gz', '.json')
                file_path = join(self.output_dir, 'NuQCJob', 'fastp_reports_dir', 'json', new_name)
                cmds.append(f"echo 'This is a json file.' > {file_path}")

        cmds.append("echo 'Submitted batch job 9999999'")
        # write all the statements out into a bash-script named 'sbatch' and
        # place it somewhere in the PATH. (It will be removed on tearDown()).
        self.create_fake_bin('sbatch', "\n".join(cmds))
        audit_results = sorted(wf.quality_control())

        # add tests to test audit results, modify test to introduce some 'zero-length' files.
        # add some assertions to show that the post-processing step is munging the correct
        # directory structure needed for subsequent steps.

        self.assertTrue(False)
