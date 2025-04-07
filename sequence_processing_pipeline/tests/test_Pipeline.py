import json
import os
from sequence_processing_pipeline.PipelineError import PipelineError
from sequence_processing_pipeline.Pipeline import Pipeline, InstrumentUtils
import unittest
from os import makedirs, walk
from os.path import abspath, basename, join, exists
from functools import partial
import re
from shutil import copy
import pandas as pd
from tempfile import NamedTemporaryFile


class TestPipeline(unittest.TestCase):
    def setUp(self):
        package_root = abspath('./sequence_processing_pipeline')
        self.path = partial(join, package_root, 'tests', 'data')
        self.good_config_file = join(package_root, 'configuration.json')
        self.bad_config_file = self.path('bad_configuration.json')
        self.invalid_config_file = 'does/not/exist/configuration.json'
        self.good_run_id = '211021_A00000_0000_SAMPLE'
        # qiita_id is randomly-generated and does not match any known
        # existing qiita job_id.
        self.qiita_id = '077c4da8-74eb-4184-8860-0207f53623be'
        self.invalid_run_id = 'not-sample-sequence-directory'
        self.output_file_path = self.path('output_dir')
        makedirs(self.output_file_path, exist_ok=True)
        self.maxDiff = None
        self.good_sample_sheet_path = self.path('good-sample-sheet.csv')
        self.good_legacy_sheet_path = self.path('mgv90_test_sheet.csv')
        self.mp_sheet_path = self.path('multi-project-sheet.csv')
        self.bad_sample_sheet_path = self.path('duplicate_sample-sample-sheet'
                                               '.csv')
        self.bad_assay_type_path = self.path('bad-sample-sheet-metagenomics'
                                             '.csv')
        self.good_run_dir = self.path(self.good_run_id)
        self.runinfo_file = self.path(self.good_run_id, 'RunInfo.xml')
        self.rtacomplete_file = self.path(self.good_run_id, 'RTAComplete.txt')
        self.good_sheet_w_replicates = self.path('good_sheet_w_replicates.csv')

        # most of the tests here were written with the assumption that these
        # files already exist.
        self.create_runinfo_file()
        self.create_rtacomplete_file()

        # read good configuration file at initialization to avoid putting
        # most of the test code within 'with open()' expressions.
        with open(self.good_config_file, 'r') as f:
            self.good_config = json.load(f)

        self.delete_these = []

    def tearDown(self):
        # Pipeline is now the only class aware of these files, hence they
        # can be deleted at the end of testing.
        self.delete_runinfo_file()
        self.delete_rtacomplete_file()
        self.delete_more_files()

    def make_runinfo_file_unreadable(self):
        os.chmod(self.runinfo_file, 0o000)

    def make_runinfo_file_readable(self):
        os.chmod(self.runinfo_file, 0o777)

    def create_runinfo_file(self, four_reads=False):
        # since good sample RunInfo.xml files already exist to support
        # other tests, reuse them here.

        f_name = 'RunInfo_Good2.xml' if four_reads else 'RunInfo_Good1.xml'
        copy(join('sequence_processing_pipeline/tests/data/', f_name),
             self.runinfo_file)

    def delete_runinfo_file(self):
        try:
            os.remove(self.runinfo_file)
        except FileNotFoundError:
            # make method idempotent
            pass

    def create_rtacomplete_file(self):
        with open(self.rtacomplete_file, 'w') as f:
            f.write("")

    def delete_rtacomplete_file(self):
        try:
            os.remove(self.rtacomplete_file)
        except FileNotFoundError:
            # make method idempotent
            pass

    def delete_more_files(self):
        for file_path in self.delete_these:
            if exists(file_path):
                # if file no longer exists, that's okay.
                os.remove(file_path)

    def _make_mapping_file(self, output_file_path):
        cols = ('sample_name', 'barcode', 'library_construction_protocol',
                'mastermix_lot', 'sample_plate', 'center_project_name',
                'instrument_model', 'tm1000_8_tool', 'well_id_96',
                'tm50_8_tool', 'tm10_8_tool', 'well_id_384',
                'well_description', 'run_prefix', 'run_date', 'center_name',
                'tm300_8_tool', 'extraction_robot', 'qiita_prep_id',
                'experiment_design_description', 'platform', 'water_lot',
                'project_name', 'pcr_primers', 'sequencing_meth', 'plating',
                'orig_name', 'linker', 'runid', 'target_subfragment', 'primer',
                'primer_plate', 'run_center', 'primer_date', 'target_gene',
                'processing_robot', 'extractionkit_lot')
        rows = [['1.0', ], ['1e-3', ]]
        rows[0].extend(['foo'] * (len(cols) - 1))
        rows[1].extend(['foo'] * (len(cols) - 1))
        df = pd.DataFrame(rows, columns=cols)
        df.to_csv(output_file_path, sep='\t', index=False, header=True)

    def test_make_sif_fname(self):
        exp = '211021_A00000_0000_SAMPLE_NYU_BMS_Melanoma_13059_blanks.tsv'
        obs = Pipeline.make_sif_fname('211021_A00000_0000_SAMPLE',
                                      'NYU_BMS_Melanoma_13059')
        self.assertEqual(exp, obs)

    def test_is_sif_fp(self):
        obs1 = Pipeline.is_sif_fp("/path/to/sifs/211021_A00000_0000_SAMPLE_"
                                  "NYU_BMS_Melanoma_13059_blanks.tsv")
        self.assertTrue(obs1)

        obs2 = Pipeline.is_sif_fp("/path/to/sifs/211021_A00000_0000_SAMPLE_"
                                  "NYU_BMS_Melanoma_13059_lord_of_the.sif")
        self.assertFalse(obs2)

    def test_get_qiita_id_from_sif_fp(self):
        exp = "13059"
        obs = Pipeline.get_qiita_id_from_sif_fp(
            "/path/to/sifs/211021_A00000_0000_SAMPLE_"
            "NYU_BMS_Melanoma_13059_blanks.tsv")
        self.assertEqual(exp, obs)

    def test_validate_mapping_file_numeric_ids(self):
        with NamedTemporaryFile() as tmp:
            self._make_mapping_file(tmp.name)
            exp = ['1.0', '1e-3']
            pipeline = Pipeline(self.good_config_file, self.good_run_id,
                                tmp.name, self.output_file_path, self.qiita_id,
                                Pipeline.AMPLICON_PTYPE)

            # explicitly call _validate_mapping_file() for testing purposes.
            # Note that currently this method was already called once during
            # the Pipeline object's initialization, implying that tmp.name was
            # validated successfully.
            obs_df = pipeline._validate_mapping_file(tmp.name)
            self.assertEqual(list(obs_df['sample_name']), exp)

    def test_get_sample_names_from_sample_sheet(self):
        pipeline = Pipeline(self.good_config_file, self.good_run_id,
                            self.mp_sheet_path,
                            self.output_file_path, self.qiita_id,
                            Pipeline.METAGENOMIC_PTYPE)

        # get all names from all projects in the sample-sheet.
        # get all names in project 'Wisconsin_U19_15445'
        # get all names in project 'Wisconsin_U19_NA_15446'
        # attempt to get names from a project not in the sheet.
        # what happens when the fully-qualified project name is used
        # (includes qiita_id)?
        # get all names in project 'Wisconsin_U19_NA_15446'

        params = [None, 'Wisconsin_U19', 'Wisconsin_U19_NA', 'NotAProject',
                  'Wisconsin_U19_15445', 'Wisconsin_U19_NA_15446']

        exps = [{'3A', '4A', '5B', '6A', '7A', '8A'}, {'3A', '4A', '5B'},
                {'6A', '8A', '7A'}, set(), {'3A', '4A', '5B'},
                {'6A', '8A', '7A'}]

        for param, exp in zip(params, exps):
            obs = set(pipeline._get_sample_names_from_sample_sheet(param))
            self.assertEqual(obs, exp)

    def test_get_orig_names_from_sheet_with_replicates(self):
        pipeline = Pipeline(self.good_config_file, self.good_run_id,
                            self.good_sheet_w_replicates,
                            self.output_file_path, self.qiita_id,
                            Pipeline.METAGENOMIC_PTYPE)

        obs = pipeline.get_orig_names_from_sheet('Feist_11661')
        exp = {'BLANK.43.12G', 'BLANK.43.12H', 'JBI.KHP.HGL.021',
               'JBI.KHP.HGL.022', 'JBI.KHP.HGL.023', 'JBI.KHP.HGL.024',
               'RMA.KHP.rpoS.Mage.Q97D', 'RMA.KHP.rpoS.Mage.Q97E',
               'RMA.KHP.rpoS.Mage.Q97L', 'RMA.KHP.rpoS.Mage.Q97N'}

        self.assertEqual(set(obs), exp)

    def test_required_file_checks(self):
        # begin this test by deleting the RunInfo.txt file and verifying that
        # Pipeline object will raise an Error.
        self.delete_runinfo_file()

        with self.assertRaisesRegex(PipelineError, "required file 'RunInfo.xml"
                                                   "' is not present."):
            Pipeline(self.good_config_file, self.good_run_id,
                     self.good_sample_sheet_path,
                     self.output_file_path,
                     self.qiita_id, Pipeline.METAGENOMIC_PTYPE)

        # delete RTAComplete.txt and recreate RunInfo.txt file to verify that
        # an Error is raised when only RTAComplete.txt is missing.
        self.delete_rtacomplete_file()
        self.create_runinfo_file()

        with self.assertRaisesRegex(PipelineError, "required file 'RTAComplete"
                                                   ".txt' is not present."):
            Pipeline(self.good_config_file, self.good_run_id,
                     self.good_sample_sheet_path,
                     self.output_file_path,
                     self.qiita_id, Pipeline.METAGENOMIC_PTYPE)

        # make RunInfo.xml file unreadable and verify that Pipeline object
        # raises the expected Error.
        self.create_rtacomplete_file()
        self.make_runinfo_file_unreadable()

        with self.assertRaisesRegex(PipelineError, "RunInfo.xml is present, bu"
                                                   "t not readable"):
            Pipeline(self.good_config_file, self.good_run_id,
                     self.good_sample_sheet_path,
                     self.output_file_path,
                     self.qiita_id, Pipeline.METAGENOMIC_PTYPE)
        self.make_runinfo_file_readable()

    def test_creation(self):
        # Pipeline should assert due to config_file
        with self.assertRaises(PipelineError) as e:
            Pipeline(self.bad_config_file,
                     self.good_run_id,
                     self.good_sample_sheet_path,
                     self.output_file_path,
                     self.qiita_id, Pipeline.METAGENOMIC_PTYPE)

        msg = re.sub(r'not a key in .*?/sequence_processing_pipeline',
                     r'not a key in sequence_processing_pipeline',
                     str(e.exception))
        self.assertEqual(msg, "'search_paths' is not a key in "
                              "sequence_processing_pipeline/tests"
                              "/data/bad_configuration.json")

        # Pipeline should assert due to Assay having a bad value.
        with self.assertRaisesRegex(ValueError, "bad-sample-sheet-metagenomics"
                                                ".csv' does not appear to be a"
                                                " valid sample-sheet."):
            Pipeline(self.good_config_file,
                     self.good_run_id,
                     self.bad_assay_type_path,
                     self.output_file_path,
                     self.qiita_id, Pipeline.METAGENOMIC_PTYPE)

        # Pipeline should assert due to an invalid config file path.
        with self.assertRaises(PipelineError) as e:
            Pipeline(self.invalid_config_file,
                     self.good_run_id,
                     self.good_sample_sheet_path,
                     self.output_file_path,
                     self.qiita_id, Pipeline.METAGENOMIC_PTYPE)

        self.assertEqual(str(e.exception), 'does/not/exist/configuration.json '
                                           'does not exist.')

        # Pipeline should assert on config_file = None
        with self.assertRaises(PipelineError) as e:
            Pipeline(None,
                     self.good_run_id,
                     self.good_sample_sheet_path,
                     self.output_file_path,
                     self.qiita_id, Pipeline.METAGENOMIC_PTYPE)

        self.assertEqual(str(e.exception), 'configuration_file_path cannot be '
                                           'None')

        # Pipeline should assert due to invalid_run_id
        with self.assertRaises(PipelineError) as e:
            Pipeline(self.good_config_file,
                     self.invalid_run_id,
                     self.good_sample_sheet_path,
                     self.output_file_path,
                     self.qiita_id, Pipeline.METAGENOMIC_PTYPE)

        self.assertEqual(str(e.exception), "A run-dir for 'not-sample-sequence"
                                           "-directory' could not be found")

        # Pipeline should assert on run_id = None
        with self.assertRaises(PipelineError) as e:
            Pipeline(self.good_config_file,
                     None,
                     self.good_sample_sheet_path,
                     self.output_file_path,
                     self.qiita_id, Pipeline.METAGENOMIC_PTYPE)

        # Pipeline should assert if config file is not valid JSON.
        # good_sample_sheet_path is obviously not a valid JSON file.
        with self.assertRaisesRegex(PipelineError, "good-sample-sheet.csv is "
                                                   "not a valid json file"):
            Pipeline(self.good_sample_sheet_path,
                     self.good_run_id,
                     self.good_sample_sheet_path,
                     self.output_file_path,
                     self.qiita_id, Pipeline.METAGENOMIC_PTYPE)

        bad_json_file = self.path('configuration_profiles', 'bad.json')
        self.delete_these.append(bad_json_file)

        # test Error returned when root attribute 'profile' does not exist.
        with open(bad_json_file, 'w') as f:
            f.write('{ "not_profile": { "instrument_type": "default", '
                    '"configuration": { "bcl2fastq": { "nodes": 1, "nprocs": '
                    '16, "queue": "qiita", "wallclock_time_in_minutes": 216, '
                    '"modules_to_load": [ "bcl2fastq_2.20.0.422" ], '
                    '"executable_path": "bcl2fastq", '
                    '"per_process_memory_limit": "10gb" } } } }')

        with self.assertRaisesRegex(ValueError, "'profile' is not an attribute"
                                                " in 'sequence_processing_"
                                                "pipeline/tests/data/"
                                                "configuration_profiles/"
                                                "bad.json'"):
            Pipeline(self.good_config_file,
                     self.good_run_id,
                     self.good_sample_sheet_path,
                     self.output_file_path,
                     self.qiita_id, Pipeline.METAGENOMIC_PTYPE)

        # test Error returned when 'instrument_type' does not exist.
        bad_json_file = self.path('configuration_profiles', 'bad.json')

        with open(bad_json_file, 'w') as f:
            f.write('{ "profile": { "not_instrument_type": "default", '
                    '"configuration": { "bcl2fastq": { "nodes": 1, "nprocs": '
                    '16, "queue": "qiita", "wallclock_time_in_minutes": 216, '
                    '"modules_to_load": [ "bcl2fastq_2.20.0.422" ], '
                    '"executable_path": "bcl2fastq", '
                    '"per_process_memory_limit": "10gb" } } } }')

        with self.assertRaisesRegex(ValueError, "'instrument_type' is not an "
                                                "attribute in 'sequence_"
                                                "processing_pipeline/tests/"
                                                "data/configuration_profiles/"
                                                "bad.json'"):
            Pipeline(self.good_config_file,
                     self.good_run_id,
                     self.good_sample_sheet_path,
                     self.output_file_path,
                     self.qiita_id, Pipeline.METAGENOMIC_PTYPE)

        # test Error returned when a non-default profile is missing assay_type
        bad_json_file = self.path('configuration_profiles', 'bad.json')
        self.delete_these.append(bad_json_file)

        with open(bad_json_file, 'w') as f:
            f.write('{ "profile": { "instrument_type": "MiSeq", '
                    '"configuration": { "bcl2fastq": { "nodes": 1, "nprocs": '
                    '16, "queue": "qiita", "wallclock_time_in_minutes": 216, '
                    '"modules_to_load": [ "bcl2fastq_2.20.0.422" ], '
                    '"executable_path": "bcl2fastq", '
                    '"per_process_memory_limit": "10gb" } } } }')

        with self.assertRaisesRegex(ValueError, "'assay_type' is not an "
                                                "attribute in 'sequence_"
                                                "processing_pipeline/tests/"
                                                "data/configuration_profiles/"
                                                "bad.json'"):
            Pipeline(self.good_config_file,
                     self.good_run_id,
                     self.good_sample_sheet_path,
                     self.output_file_path,
                     self.qiita_id, Pipeline.METAGENOMIC_PTYPE)

    def test_sample_sheet_validation(self):
        # test successful validation of a good sample-sheet.
        # if self.good_sample_sheet_path points to a bad sample-sheet, then
        # Pipeline would raise a PipelineError w/warnings and error messages
        # contained w/in its 'message' member.
        try:
            Pipeline(self.good_config_file, self.good_run_id,
                     self.good_sample_sheet_path,
                     self.output_file_path,
                     self.qiita_id, Pipeline.METAGENOMIC_PTYPE)
        except PipelineError as e:
            self.fail(("test_filter_directories_for_time failed w/PipelineEr"
                       f"ror: {e.message}"))

        # test unsuccessful validation of a bad sample-sheet
        with self.assertRaises(PipelineError) as e:
            Pipeline(self.good_config_file, self.good_run_id,
                     self.bad_sample_sheet_path,
                     self.output_file_path,
                     self.qiita_id, Pipeline.METAGENOMIC_PTYPE)
        self.assertEqual(str(e.exception), ('Sample-sheet contains errors:\n'
                                            'A sample already exists with '
                                            'lane 1 and sample-id '
                                            'EP479894B04'))

    def test_generate_sample_information_files(self):
        # test sample-information-file generation.
        pipeline = Pipeline(self.good_config_file, self.good_run_id,
                            self.good_sample_sheet_path,
                            self.output_file_path, self.qiita_id,
                            Pipeline.METAGENOMIC_PTYPE)

        paths = pipeline.generate_sample_info_files()

        # confirm files exist in the expected location and with the expected
        # filenames.
        obs = [x.split('sequence_processing_pipeline/')[1] for x in paths]
        exp = [(f'tests/data/output_dir/{self.good_run_id}'
                '_NYU_BMS_Melanoma_13059_blanks.tsv'),
               (f'tests/data/output_dir/{self.good_run_id}'
                '_Feist_11661_blanks.tsv'),
               (f'tests/data/output_dir/{self.good_run_id}'
                '_Gerwick_6123_blanks.tsv')]

        # sort the lists to ensure both are in a fixed order.
        obs.sort()
        exp.sort()

        self.assertEqual(obs, exp)

        # confirm files contain the expected number of lines.
        # This is going to be based on the number of samples named 'BLANK*'
        # in good-sample-sheet.csv.
        exp_lines = {f'{self.good_run_id}_NYU_BMS_Melanoma_13059_blanks.tsv':
                     33,
                     f'{self.good_run_id}_Feist_11661_blanks.tsv': 8,
                     f'{self.good_run_id}_Gerwick_6123_blanks.tsv': 2}

        exp_first_lines = {
            f'{self.good_run_id}_NYU_BMS_Melanoma_13059_blanks.tsv':
            'BLANK1.1A\t2021-10-21\t193\t'
            'Control\tNegative\tSterile w'
            'ater blank\tSterile water blank\turban biome\tres'
            'earch facility\tsterile wate'
            'r\tmisc environment\tUSA:CA:'
            'San Diego\tBLANK1.1A\t32.5\t'
            '-117.25\tcontrol blank\tmeta'
            'genome\t256318\tBLANK1.1A\tN'
            'YU_BMS_Melanoma\tTRUE\t'
            'UCSD\tFALSE',
            f'{self.good_run_id}_Feist_11661_blanks.tsv':
            'BLANK.40.12G\t2021-10-21\t193\tControl'
            '\tNegative\tSterile water blank\tSterile water blank\turban '
            'biome\tresearch facility\tsterile water'
            '\tmisc environment\tUSA:CA:San Diego\tB'
            'LANK.40.12G\t32.5\t-117.25\tcontrol bla'
            'nk\tmetagenome\t256318\tBLANK.40.12G\t'
            'Feist\tTRUE\tUCSD\tFALSE',
            f'{self.good_run_id}_Gerwick_6123_blanks.tsv':
            'BLANK.41.12G\t2021-10-21\t193\tControl'
            '\tNegative\tSterile water blank\tSterile water blank\turban'
            ' biome\tresearch facility\tsterile wat'
            'er\tmisc environment\tUSA:CA:San Diego'
            '\tBLANK.41.12G\t32.5\t-117.25\tcontrol'
            ' blank\tmetagenome\t256318\tBLANK.41.1'
            '2G\tGerwick\tTRUE\tUCSD\tFALSE'
        }

        exp_last_lines = {
            f'{self.good_run_id}_NYU_BMS_Melanoma_13059_blanks.tsv':
            'BLANK4.4H\t2021-10-21\t193\t'
            'Control\tNegative\tSterile w'
            'ater blank\tSterile water blank\turban biome\tres'
            'earch facility\tsterile wate'
            'r\tmisc environment\tUSA:CA:'
            'San Diego\tBLANK4.4H\t32.5\t'
            '-117.25\tcontrol blank\tmeta'
            'genome\t256318\tBLANK4.4H\tN'
            'YU_BMS_Melanoma\tTRUE\t'
            'UCSD\tFALSE',
            f'{self.good_run_id}_Feist_11661_blanks.tsv':
            'BLANK.43.12H\t2021-10-21\t193\tControl'
            '\tNegative\tSterile water blank\tSterile water blank\turban'
            ' biome\tresearch facility\tsterile wat'
            'er\tmisc environment\tUSA:CA:San Diego'
            '\tBLANK.43.12H\t32.5\t-117.25\tcontrol'
            ' blank\tmetagenome\t256318\tBLANK.43.1'
            '2H\tFeist\tTRUE\tUCSD\tFALSE',
            f'{self.good_run_id}_Gerwick_6123_blanks.tsv':
            'BLANK.41.12G\t2021-10-21\t193\tContro'
            'l\tNegative\tSterile water blank\tSterile water blank\turb'
            'an biome\tresearch facility\tsterile '
            'water\tmisc environment\tUSA:CA:San D'
            'iego\tBLANK.41.12G\t32.5\t-117.25\tco'
            'ntrol blank\tmetagenome\t256318\tBLAN'
            'K.41.12G\tGerwick\tTRUE\tUCSD\t'
            'FALSE'
        }

        for some_path in paths:
            some_name = basename(some_path)
            with open(some_path, 'r') as f:
                obs_lines = f.readlines()
                self.assertEqual(len(obs_lines), exp_lines[some_name])
                # confirm that each file contains the expected header.
                header = obs_lines[0].strip()
                self.assertEqual(header, '\t'.join(Pipeline.sif_header))
                # confirm that the first line of each file is as expected.
                obs = obs_lines[1].strip()
                exp = exp_first_lines[some_name]

                self.assertEqual(obs, exp)

                # confirm that the last line of each file is as expected.
                obs = obs_lines[-1].strip()
                exp = exp_last_lines[some_name]
                self.assertEqual(obs, exp)

    def test_generate_sample_information_files_with_additional_meta(self):
        # TODO: With recent changes is this needed?
        # test sample-information-file generation.
        pipeline = Pipeline(self.good_config_file, self.good_run_id,
                            self.good_sample_sheet_path,
                            self.output_file_path, self.qiita_id,
                            Pipeline.METAGENOMIC_PTYPE)

        # create a dataframe with duplicate information to pass to
        # generate_sample_information_files(). Confirm that the duplicates
        # are dropped. Confirm 'NOTBLANK_999A' is also filtered out.
        df = pd.DataFrame(data=[('BLANK999_999A', 'NYU_BMS_Melanoma_13059'),
                                ('BLANK999_999A', 'NYU_BMS_Melanoma_13059'),
                                ('NOTBLANK_999A', 'NYU_BMS_Melanoma_13059')],
                          columns=['sample_name', 'project_name'])

        sif_path = pipeline.generate_sample_info_files(addl_info=df)

        # get the path for the NYU_BMS_Melanoma dataset.
        sif_path = [x for x in sif_path if 'NYU_BMS_Melanoma' in x][0]

        exp_first_line = ("BLANK1.1A\t2021-10-21\t193\t"
                          "Control\tNegative\tSterile water blank\t"
                          "Sterile water blank\turban biome\t"
                          "research facility\tsterile water\t"
                          "misc environment\tUSA:CA:San Diego\t"
                          "BLANK1.1A\t32.5\t-117.25\tcontrol blank\t"
                          "metagenome\t256318\tBLANK1.1A\t"
                          "NYU_BMS_Melanoma\tTRUE\tUCSD\tFALSE")

        exp_last_line = ("BLANK4.4H\t2021-10-21\t193\tControl\tNegative\t"
                         "Sterile water blank\tSterile water blank\t"
                         "urban biome\tresearch facility\tsterile water\t"
                         "misc environment\tUSA:CA:San Diego\tBLANK4.4H\t"
                         "32.5\t-117.25\tcontrol blank\tmetagenome\t256318\t"
                         "BLANK4.4H\tNYU_BMS_Melanoma\tTRUE\tUCSD\tFALSE")

        with open(sif_path, 'r') as f:
            obs_lines = f.readlines()

            # confirm that each file contains the expected header.
            header = obs_lines[0].strip()
            self.assertEqual(header, '\t'.join(Pipeline.sif_header))

            # confirm that the first line of each file is as expected.
            obs = obs_lines[1].strip()
            exp = exp_first_line

            self.assertEqual(obs, exp)

            # confirm that the last line of each file is as expected.
            obs = obs_lines[-1].strip()
            exp = exp_last_line

            self.assertEqual(obs, exp)

    def test_get_sample_ids(self):
        exp_sample_ids = ['CDPH-SAL__Salmonella__Typhi__MDL-143',
                          'CDPH-SAL_Salmonella_Typhi_MDL-144',
                          'CDPH-SAL_Salmonella_Typhi_MDL-145',
                          'CDPH-SAL_Salmonella_Typhi_MDL-146',
                          'CDPH-SAL_Salmonella_Typhi_MDL-147',
                          'CDPH-SAL_Salmonella_Typhi_MDL-148',
                          'CDPH-SAL_Salmonella_Typhi_MDL-149',
                          'CDPH-SAL_Salmonella_Typhi_MDL-150',
                          'CDPH-SAL_Salmonella_Typhi_MDL-151',
                          'CDPH-SAL_Salmonella_Typhi_MDL-152',
                          'CDPH-SAL_Salmonella_Typhi_MDL-153',
                          'CDPH-SAL_Salmonella_Typhi_MDL-154',
                          'CDPH-SAL_Salmonella_Typhi_MDL-155',
                          'CDPH-SAL_Salmonella_Typhi_MDL-156',
                          'CDPH-SAL_Salmonella_Typhi_MDL-157',
                          'CDPH-SAL_Salmonella_Typhi_MDL-158',
                          'CDPH-SAL_Salmonella_Typhi_MDL-159',
                          'CDPH-SAL_Salmonella_Typhi_MDL-160',
                          'CDPH-SAL_Salmonella_Typhi_MDL-161',
                          'CDPH-SAL_Salmonella_Typhi_MDL-162',
                          'CDPH-SAL_Salmonella_Typhi_MDL-163',
                          'CDPH-SAL_Salmonella_Typhi_MDL-164',
                          'CDPH-SAL_Salmonella_Typhi_MDL-165',
                          'CDPH-SAL_Salmonella_Typhi_MDL-166',
                          'CDPH-SAL_Salmonella_Typhi_MDL-167',
                          'CDPH-SAL_Salmonella_Typhi_MDL-168',
                          'P21_E_coli_ELI344', 'P21_E_coli_ELI345',
                          'P21_E_coli_ELI347', 'P21_E_coli_ELI348',
                          'P21_E_coli_ELI349', 'P21_E_coli_ELI350',
                          'P21_E_coli_ELI351', 'P21_E_coli_ELI352',
                          'P21_E_coli_ELI353', 'P21_E_coli_ELI354',
                          'P21_E_coli_ELI355', 'P21_E_coli_ELI357',
                          'P21_E_coli_ELI358', 'P21_E_coli_ELI359',
                          'P21_E_coli_ELI361', 'P21_E_coli_ELI362',
                          'P21_E_coli_ELI363', 'P21_E_coli_ELI364',
                          'P21_E_coli_ELI365', 'P21_E_coli_ELI366',
                          'P21_E_coli_ELI367', 'P21_E_coli_ELI368',
                          'P21_E_coli_ELI369', 'stALE_E_coli_A1_F21_I1_R1',
                          'stALE_E_coli_A2_F21_I1_R1',
                          'stALE_E_coli_A3_F18_I1_R1',
                          'stALE_E_coli_A3_F40_I1_R1',
                          'stALE_E_coli_A4_F21_I1_R1',
                          'stALE_E_coli_A4_F21_I1_R2',
                          'stALE_E_coli_A4_F42_I1_R1',
                          'stALE_E_coli_A5_F21_I1_R1',
                          'stALE_E_coli_A5_F42_I1_R1',
                          'stALE_E_coli_A6_F21_I1_R1',
                          'stALE_E_coli_A6_F43_I1_R1',
                          'stALE_E_coli_A7_F21_I1_R1',
                          'stALE_E_coli_A7_F42_I1_R1',
                          'stALE_E_coli_A8_F20_I1_R1',
                          'stALE_E_coli_A8_F42_I1_R1',
                          'stALE_E_coli_A9_F21_I1_R1',
                          'stALE_E_coli_A9_F44_I1_R1',
                          'stALE_E_coli_A10_F21_I1_R1',
                          'stALE_E_coli_A10_F43_I1_R1',
                          'stALE_E_coli_A10_F131_I1_R1',
                          'stALE_E_coli_A11_F21_I1_R1',
                          'stALE_E_coli_A11_F43_I1_R1',
                          'stALE_E_coli_A11_F119_I1_R1',
                          'stALE_E_coli_A12_F21_I1_R1',
                          'stALE_E_coli_A12_F43_I1_R1',
                          'stALE_E_coli_A12_F136_I1_R1',
                          'stALE_E_coli_A13_F20_I1_R1',
                          'stALE_E_coli_A13_F42_I1_R1',
                          'stALE_E_coli_A13_F121_I1_R1',
                          'stALE_E_coli_A14_F20_I1_R1',
                          'stALE_E_coli_A14_F42_I1_R1',
                          'stALE_E_coli_A14_F133_I1_R1',
                          'stALE_E_coli_A15_F21_I1_R1',
                          'stALE_E_coli_A15_F42_I1_R1',
                          'stALE_E_coli_A15_F117_I1_R1',
                          'stALE_E_coli_A16_F20_I1_R1',
                          'stALE_E_coli_A16_F42_I1_R1',
                          'stALE_E_coli_A16_F134_I1_R1',
                          'stALE_E_coli_A17_F21_I1_R1',
                          'stALE_E_coli_A17_F118_I1_R1',
                          'stALE_E_coli_A18_F18_I1_R1',
                          'stALE_E_coli_A18_F39_I1_R1',
                          'stALE_E_coli_A18_F130_I1_R1', '3A', '4A',
                          'BLANK_40_12G', 'BLANK_40_12H',
                          'Pputida_JBEI__HGL_Pputida_107_BP6',
                          'Pputida_JBEI__HGL_Pputida_108_BP7',
                          'Pputida_JBEI__HGL_Pputida_109_BP8',
                          'Pputida_JBEI__HGL_Pputida_110_M2',
                          'Pputida_JBEI__HGL_Pputida_111_M5',
                          'Pputida_TALE__HGL_Pputida_112',
                          'Pputida_TALE__HGL_Pputida_113',
                          'Pputida_TALE__HGL_Pputida_114',
                          'Pputida_TALE__HGL_Pputida_115',
                          'Pputida_TALE__HGL_Pputida_116',
                          'Pputida_TALE__HGL_Pputida_117',
                          'Pputida_TALE__HGL_Pputida_118',
                          'Pputida_TALE__HGL_Pputida_119',
                          'Pputida_TALE__HGL_Pputida_120',
                          'Pputida_TALE__HGL_Pputida_121',
                          'Pputida_TALE__HGL_Pputida_122',
                          'Pputida_TALE__HGL_Pputida_123',
                          'Pputida_TALE__HGL_Pputida_124',
                          'Pputida_TALE__HGL_Pputida_125',
                          'Pputida_TALE__HGL_Pputida_126',
                          'Pputida_TALE__HGL_Pputida_127',
                          'Pputida_TALE__HGL_Pputida_128',
                          'Pputida_TALE__HGL_Pputida_129',
                          'Pputida_TALE__HGL_Pputida_130',
                          'Pputida_TALE__HGL_Pputida_131',
                          'Pputida_TALE__HGL_Pputida_132',
                          'Pputida_TALE__HGL_Pputida_133',
                          'Pputida_TALE__HGL_Pputida_134',
                          'Pputida_TALE__HGL_Pputida_135',
                          'Pputida_TALE__HGL_Pputida_136',
                          'Pputida_TALE__HGL_Pputida_137',
                          'Pputida_TALE__HGL_Pputida_138',
                          'Pputida_TALE__HGL_Pputida_139',
                          'Pputida_TALE__HGL_Pputida_140',
                          'Pputida_TALE__HGL_Pputida_141',
                          'Pputida_TALE__HGL_Pputida_142',
                          'Pputida_TALE__HGL_Pputida_143',
                          'Pputida_TALE__HGL_Pputida_144',
                          'Pputida_PALE__HGL_Pputida_145',
                          'Pputida_PALE__HGL_Pputida_146',
                          'Pputida_PALE__HGL_Pputida_147',
                          'Pputida_PALE__HGL_Pputida_148',
                          'Pputida_PALE__HGL_Pputida_149',
                          'Pputida_PALE__HGL_Pputida_150',
                          'Pputida_PALE__HGL_Pputida_151',
                          'Pputida_PALE__HGL_Pputida_152',
                          'Pputida_PALE__HGL_Pputida_153',
                          'Pputida_PALE__HGL_Pputida_154',
                          'Pputida_PALE__HGL_Pputida_155',
                          'Pputida_PALE__HGL_Pputida_156',
                          'Pputida_PALE__HGL_Pputida_157',
                          'Pputida_PALE__HGL_Pputida_158',
                          'Pputida_PALE__HGL_Pputida_159',
                          'Pputida_PALE__HGL_Pputida_160',
                          'Pputida_PALE__HGL_Pputida_161',
                          'Pputida_PALE__HGL_Pputida_162',
                          'Pputida_PALE__HGL_Pputida_163',
                          'Pputida_PALE__HGL_Pputida_164',
                          'Pputida_PALE__HGL_Pputida_165',
                          'Pputida_PALE__HGL_Pputida_166',
                          'Pputida_PALE__HGL_Pputida_167',
                          'Pputida_PALE__HGL_Pputida_168',
                          'Pputida_PALE__HGL_Pputida_169',
                          'Pputida_PALE__HGL_Pputida_170',
                          'Pputida_PALE__HGL_Pputida_171',
                          'Pputida_PALE__HGL_Pputida_172',
                          'Pputida_PALE__HGL_Pputida_173',
                          'Pputida_PALE__HGL_Pputida_174',
                          'Pputida_PALE__HGL_Pputida_175',
                          'Pputida_PALE__HGL_Pputida_176',
                          'JM-Metabolic__GN0_2005', 'JM-Metabolic__GN0_2007',
                          'JM-Metabolic__GN0_2009', 'JM-Metabolic__GN0_2094',
                          'JM-Metabolic__GN0_2099', 'JM-Metabolic__GN0_2148',
                          'JM-Metabolic__GN0_2165', 'JM-Metabolic__GN0_2169',
                          'JM-Metabolic__GN0_2172', 'JM-Metabolic__GN0_2175',
                          'JM-Metabolic__GN0_2183', 'JM-Metabolic__GN0_2215',
                          'JM-Metabolic__GN0_2254', 'JM-Metabolic__GN0_2277',
                          'JM-Metabolic__GN0_2290', 'JM-Metabolic__GN0_2337',
                          'JM-Metabolic__GN0_2317', 'JM-Metabolic__GN0_2354',
                          'JM-Metabolic__GN0_2375', 'JM-Metabolic__GN0_2380',
                          'JM-Metabolic__GN0_2393', 'JM-Metabolic__GN0_2404',
                          '5B', '6A', 'BLANK_41_12G', 'BLANK_41_12H',
                          'Deoxyribose_PALE_ALE__MG1655_BOP27_4_14',
                          'Deoxyribose_PALE_ALE__MG1655_BOP27_4_23',
                          'Deoxyribose_PALE_ALE__MG1655_BOP27_4_48',
                          'Deoxyribose_PALE_ALE__MG1655_BOP27_6_21',
                          'Deoxyribose_PALE_ALE__MG1655_BOP27_6_35',
                          'Deoxyribose_PALE_ALE__MG1655_BOP27_10_13',
                          'Deoxyribose_PALE_ALE__MG1655_BOP27_10_28',
                          'Deoxyribose_PALE_ALE__MG1655_BOP27_10_51',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_18_19',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_18_59',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_18_35',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_20_16',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_20_43',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_20_71',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_22_16',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_22_28',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_22_52',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_24_9',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_24_24',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_24_52',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_26_6',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_26_27',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_26_69',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_28_13',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_28_28',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_28_53',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_30_7',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_30_22',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_30_60',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_32_6',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_32_20',
                          'Deoxyribose_PALE_ALE__MG1655_Lib4_32_56',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_1_24',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_1_57',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_1_69',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_3_23',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_3_50',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_3_61',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_5_22',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_5_36',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_5_46',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_7_23',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_7_41',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_7_51',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_17_25',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_17_58',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_17_64',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_19_25',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_19_55',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_19_63',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_21_23',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_21_46',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_21_51',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_29_25',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_29_49',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_29_57',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_31_24',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_31_42',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_31_62',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_33_21',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_33_41',
                          'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_33_50',
                          'JM-Metabolic__GN02514', 'JM-Metabolic__GN02529',
                          'JM-Metabolic__GN02531', 'JM-Metabolic__GN02567',
                          'JM-Metabolic__GN02590', 'JM-Metabolic__GN02657',
                          'JM-Metabolic__GN02748', 'JM-Metabolic__GN02766',
                          'JM-Metabolic__GN02769', 'JM-Metabolic__GN02787',
                          'JM-Metabolic__GN03132', 'JM-Metabolic__GN03218',
                          'JM-Metabolic__GN03252', 'JM-Metabolic__GN03409',
                          'JM-Metabolic__GN04014', 'JM-Metabolic__GN04094',
                          'JM-Metabolic__GN04255', 'JM-Metabolic__GN04306',
                          'JM-Metabolic__GN04428', 'JM-Metabolic__GN04488',
                          'JM-Metabolic__GN04540', 'JM-Metabolic__GN04563',
                          'JM-Metabolic__GN04612', 'JM-Metabolic__GN04665',
                          'JM-Metabolic__GN04682', 'JM-Metabolic__GN05002',
                          'JM-Metabolic__GN05109', 'JM-Metabolic__GN05128',
                          'JM-Metabolic__GN05367', 'JM-Metabolic__GN05377',
                          '7A', '8A', 'BLANK_42_12G', 'BLANK_42_12H',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0326',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0327',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0328',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0329',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0330',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0352',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0353',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0354',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0355',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0356',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0357',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0364',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0366',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0367',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0368',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0369',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0370',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0371',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0372',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0373',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0374',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0375',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0376',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0377',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0378',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0380',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0381',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0382',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0383',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0384',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0385',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0386',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0387',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0388',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0389',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0390',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0391',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0392',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0393',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0394',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0395',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0396',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0397',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0398',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0399',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0400',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0401',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0402',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0403',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0404',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0405',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0406',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0407',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0408',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0409',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0417',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0418',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0419',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0420',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0421',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0473',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0474',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0483',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0484',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0485',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0486',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0516',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0517',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0518',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0519',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0520',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0521',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0522',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0523',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0524',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0525',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R08624',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R08704',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R10727',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11044',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11078',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11101',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11102',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11103',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11135',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11153',
                          'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11154',
                          'JM-Metabolic__GN02424', 'JM-Metabolic__GN02446',
                          'JM-Metabolic__GN02449', 'JM-Metabolic__GN02487',
                          'JM-Metabolic__GN02501', 'ISB', 'GFR',
                          'BLANK_43_12G',
                          'BLANK_43_12H', 'RMA_KHP_rpoS_Mage_Q97D',
                          'RMA_KHP_rpoS_Mage_Q97L', 'RMA_KHP_rpoS_Mage_Q97N',
                          'RMA_KHP_rpoS_Mage_Q97E', 'JBI_KHP_HGL_021',
                          'JBI_KHP_HGL_022', 'JBI_KHP_HGL_023',
                          'JBI_KHP_HGL_024', 'JBI_KHP_HGL_025',
                          'JBI_KHP_HGL_026', 'JBI_KHP_HGL_027',
                          'JBI_KHP_HGL_028_Amitesh_soxR',
                          'JBI_KHP_HGL_029_Amitesh_oxyR',
                          'JBI_KHP_HGL_030_Amitesh_soxR_oxyR',
                          'JBI_KHP_HGL_031_Amitesh_rpoS', 'BLANK1_1A',
                          'BLANK1_1B', 'BLANK1_1C', 'BLANK1_1D', 'BLANK1_1E',
                          'BLANK1_1F', 'BLANK1_1G', 'BLANK1_1H', 'AP581451B02',
                          'EP256645B01', 'EP112567B02', 'EP337425B01',
                          'LP127890A01', 'EP159692B04', 'EP987683A01',
                          'AP959450A03', 'SP464350A04', 'C9', 'ep256643b01',
                          'EP121011B01', 'AP616837B04', 'SP506933A04',
                          'EP159695B01', 'EP256644B01', 'SP511289A02',
                          'EP305735B04', 'SP415030A01', 'AP549681B02',
                          'AP549678B01', 'EP260544B04', 'EP202452B01',
                          'EP282276B04', 'SP531696A04', 'SP515443A04',
                          'SP515763A04', 'EP184255B04', 'SP503615A02',
                          'EP260543B04', 'EP768748A04', 'AP309872B03',
                          'AP568785B04', 'EP721390A04', 'EP940013A01',
                          'EP291979B04', 'EP182065B04', 'EP128904B02',
                          'EP915769A04', 'SP464352A03', 'SP365864A04',
                          'SP511294A04', 'EP061002B01', 'SP410793A01',
                          'SP232077A04', 'EP128910B01', 'AP531397B04',
                          'EP043583B01', 'EP230245B01', 'EP606652B04',
                          'EP207041B01', 'EP727972A04', 'EP291980B04',
                          'EP087938B02', 'SP471496A04', 'SP573823A04',
                          'EP393718B01', 'SP612496A01', 'EP032410B02',
                          'EP073216B01', 'EP410046B01', 'SP561451A04',
                          'EP320438B01', 'SP612495A04', 'EP446604B03',
                          'EP446602B01', 'EP182243B02', 'EP333541B04',
                          'EP238034B01', 'AP298002B02', 'EP455759B04',
                          'EP207042B04', 'LP128479A01', 'LP128476A01',
                          'EP316863B03', 'C20', 'lp127896a01', 'SP491907A02',
                          'EP182060B03', 'EP422407B01', 'SP573859A04',
                          'SP584547A02', 'EP182346B04', 'AP668631B04',
                          'EP451428B04', 'LP128538A01', 'SP490298A02',
                          'SP573860A01', 'EP032412B02', 'EP163771B01',
                          'LP169879A01', 'EP729433A02', 'EP447940B04',
                          'SP584551A08', 'EP216516B04', 'EP023808B02',
                          'BLANK2_2A', 'BLANK2_2B', 'BLANK2_2C', 'BLANK2_2D',
                          'BLANK2_2E', 'BLANK2_2F', 'BLANK2_2G', 'BLANK2_2H',
                          'SP573843A04', 'EP683835A01', 'SP573824A04',
                          'SP335002A04', 'SP478193A02', 'SP232311A04',
                          'SP415021A02', 'SP231630A02', 'SP641029A02',
                          'SP232310A04', 'EP617442B01', 'EP587478B04',
                          'EP447928B04', 'EP587475B04', 'EP675042B01',
                          'EP554513B02', 'EP702221B04', 'AP568787B02',
                          'EP054632B01', 'EP121013B01', 'EP649418A02',
                          'EP573313B01', 'LP154981A01', 'AP470859B01',
                          'LP154986A01', 'AP732307B04', 'EP533426B03',
                          'EP587476B04', 'AP696363B02', 'EP587477B04',
                          'SP683466A02', 'EP554518B04', 'EP533429B04',
                          'EP431570B01', 'EP202095B04', 'EP504030B04',
                          'EP207036B01', 'EP393717B01', 'SP491898A02',
                          'EP484973B04', 'EP479794B02', 'EP554515B04',
                          'SP631994A04', 'EP921593A04', 'AP787247B04',
                          'EP090129B04', 'EP447975B02', 'EP212214B01',
                          'EP410042B01', 'SP404409A02', 'SP247340A04',
                          'AP029018B01', 'EP872341A01', 'AP062219B03',
                          'EP790020A02', 'EP808112A04', 'SP404403A02',
                          'EP073160B01', 'EP012991B03', 'SP317297A02',
                          'EP656055A04', 'EP649623A01', 'EP790019A01',
                          'SP257519A04', 'EP808104A01', 'EP808106A01',
                          'SP231629A02', 'EP675044A01', 'EP657260A01',
                          'EP808110A04', 'AP032413B04', 'EP843906A04',
                          'AP173305B04', 'SP231628A02', 'AP173301B04',
                          'SP404405A02', 'EP649653A04', 'EP718687A04',
                          'AP905750A02', 'EP738468A01', 'C6', 'EP890157A02',
                          'SP353893A02', 'EP944059A02', 'EP970005A01',
                          'EP927461A04', 'EP808111A03', 'EP927459A04',
                          'SP317293A02', 'SP235186A04', 'SP399724A04',
                          'EP738469A01', 'SP284095A03', 'C5', 'EP337325B04',
                          'EP759450A04', 'BLANK3_3A', 'BLANK3_3B', 'BLANK3_3C',
                          'BLANK3_3D', 'BLANK3_3E', 'BLANK3_3F', 'BLANK3_3G',
                          'BLANK3_3H', 'AP006367B02', 'EP929277A02',
                          'AP324642B04', 'EP786631A04', 'EP657385A04',
                          'SP235189A01', 'EP448041B04', 'SP231631A02',
                          'SP280481A02', 'AP032412B04', 'EP649737A03',
                          'AP967057A04', 'EP876243A04', 'SP229387A04',
                          'EP667743A04', 'SP246941A01', 'AP745799A04',
                          'SP205732A02', 'SP230382A04', 'SP230380A02',
                          'SP230381A01', 'SP205754A01', 'EP606662B04',
                          'AP780167B02', 'EP447927B04', 'C18', 'LP191039A01',
                          'EP606663B04', 'EP573296B01', 'EP447926B04',
                          'LP127767A01', 'EP479266B04', 'LP128543A01',
                          'EP479270B03', 'EP921594A04', 'EP554501B04',
                          'EP542577B04', 'EP487995B04', 'EP542578B04',
                          'EP573310B01', 'EP244366B01', 'EP533389B03',
                          'EP244360B01', 'AP911328B01', 'AP481403B02',
                          '22_001_801_552_503_00', 'EP372981B04',
                          'EP447929B04', 'SP573849A04', 'SP577399A02',
                          'EP606656B03', 'LP166715A01', 'AP668628B04',
                          'C14', 'EP446610B02', 'EP339061B02', 'SP681591A04',
                          'EP393712B02', 'EP410041B01', 'SP453872A01',
                          '22_001_710_503_791_00',
                          'LP128540A01', 'EP339053B02', 'EP617443B01',
                          'EP190307B01', 'AP795068B04', 'LP128541A01',
                          'EP584756B04', 'SP284096A02', 'EP431562B04',
                          'EP685640B01', 'EP339059B02', 'EP431575B01',
                          'EP379938B01', 'EP529635B02', 'EP554506B04',
                          'EP455757B04', 'SP491900A02', 'LP196272A01',
                          'SP704319A04', 'EP617441B01', 'AP687591B04',
                          'SP640978A02', 'EP981129A02', 'EP455763B04',
                          'EP339057B02', 'SP491897A02', 'EP980752B04',
                          'LP128539A01', 'EP996831B04', 'EP273332B04',
                          'EP483291B04', 'EP393715B01', 'EP617440B01',
                          'EP729434A01', 'SP645141A03', 'BLANK4_4A',
                          'BLANK4_4B', 'BLANK4_4C', 'BLANK4_4D', 'BLANK4_4E',
                          'BLANK4_4F', 'BLANK4_4G', 'BLANK4_4H', 'SP232114A04',
                          'EP393714B01', 'EP533388B01', 'EP724905B01',
                          'EP282108B01', 'EP282107B01', 'EP001625B01',
                          'EP073209B02', 'SP232079A01', 'EP772145A02',
                          'AP771472A04', 'AP223470B01', 'SP404412A02',
                          'EP772143A02', 'SP408629A01', 'EP749735A07',
                          'EP846485A01', 'EP808109A01', 'SP416130A04',
                          'EP882752A01', 'AP953594A02', 'AP046324B02',
                          'AP891020A04', 'EP790023A01', 'EP657386A01',
                          'EP805337A01', 'EP927458A04', 'AP173299B04',
                          'EP768164A02', 'EP886422A01', 'AP103463B01',
                          'AP744361A02', 'AP065292B01', 'SP257517A04',
                          'EP790021A04', 'EP675075A04', 'SP388683A02',
                          'SP232309A01', 'EP899038A04', 'EP636802A01',
                          'AP046327B02', 'EP905975A04', 'SP410796A02',
                          'EP784608A01', 'EP808105A01', 'SP331134A04',
                          'EP718688A01', 'SP232270A02', 'EP970001A01',
                          'EP001624B01', 'EP868682A01', 'EP927462A02', 'C3',
                          'EP890158A02', 'EP023801B04', 'EP400447B04',
                          'EP385379B01', 'EP385387B01', 'EP385384B01',
                          'SP754514A04', 'SP415025A01', 'SP415023A02',
                          'EP400448B04', 'EP479894B04']
        # test sample-information-file generation.
        pipeline = Pipeline(self.good_config_file, self.good_run_id,
                            self.good_sample_sheet_path,
                            self.output_file_path, self.qiita_id,
                            Pipeline.METAGENOMIC_PTYPE)

        obs = pipeline.get_sample_ids()
        self.assertEqual(sorted(obs), sorted(exp_sample_ids))

    def test_get_sample_names(self):
        exp_sample_ids = ['CDPH-SAL..Salmonella..Typhi..MDL-143',
                          'CDPH-SAL.Salmonella.Typhi.MDL-144',
                          'CDPH-SAL.Salmonella.Typhi.MDL-145',
                          'CDPH-SAL.Salmonella.Typhi.MDL-146',
                          'CDPH-SAL.Salmonella.Typhi.MDL-147',
                          'CDPH-SAL.Salmonella.Typhi.MDL-148',
                          'CDPH-SAL.Salmonella.Typhi.MDL-149',
                          'CDPH-SAL.Salmonella.Typhi.MDL-150',
                          'CDPH-SAL.Salmonella.Typhi.MDL-151',
                          'CDPH-SAL.Salmonella.Typhi.MDL-152',
                          'CDPH-SAL.Salmonella.Typhi.MDL-153',
                          'CDPH-SAL.Salmonella.Typhi.MDL-154',
                          'CDPH-SAL.Salmonella.Typhi.MDL-155',
                          'CDPH-SAL.Salmonella.Typhi.MDL-156',
                          'CDPH-SAL.Salmonella.Typhi.MDL-157',
                          'CDPH-SAL.Salmonella.Typhi.MDL-158',
                          'CDPH-SAL.Salmonella.Typhi.MDL-159',
                          'CDPH-SAL.Salmonella.Typhi.MDL-160',
                          'CDPH-SAL.Salmonella.Typhi.MDL-161',
                          'CDPH-SAL.Salmonella.Typhi.MDL-162',
                          'CDPH-SAL.Salmonella.Typhi.MDL-163',
                          'CDPH-SAL.Salmonella.Typhi.MDL-164',
                          'CDPH-SAL.Salmonella.Typhi.MDL-165',
                          'CDPH-SAL.Salmonella.Typhi.MDL-166',
                          'CDPH-SAL.Salmonella.Typhi.MDL-167',
                          'CDPH-SAL.Salmonella.Typhi.MDL-168',
                          'P21.E.coli.ELI344', 'P21.E.coli.ELI345',
                          'P21.E.coli.ELI347', 'P21.E.coli.ELI348',
                          'P21.E.coli.ELI349', 'P21.E.coli.ELI350',
                          'P21.E.coli.ELI351', 'P21.E.coli.ELI352',
                          'P21.E.coli.ELI353', 'P21.E.coli.ELI354',
                          'P21.E.coli.ELI355', 'P21.E.coli.ELI357',
                          'P21.E.coli.ELI358', 'P21.E.coli.ELI359',
                          'P21.E.coli.ELI361', 'P21.E.coli.ELI362',
                          'P21.E.coli.ELI363', 'P21.E.coli.ELI364',
                          'P21.E.coli.ELI365', 'P21.E.coli.ELI366',
                          'P21.E.coli.ELI367', 'P21.E.coli.ELI368',
                          'P21.E.coli.ELI369', 'stALE.E.coli.A1.F21.I1.R1',
                          'stALE.E.coli.A2.F21.I1.R1',
                          'stALE.E.coli.A3.F18.I1.R1',
                          'stALE.E.coli.A3.F40.I1.R1',
                          'stALE.E.coli.A4.F21.I1.R1',
                          'stALE.E.coli.A4.F21.I1.R2',
                          'stALE.E.coli.A4.F42.I1.R1',
                          'stALE.E.coli.A5.F21.I1.R1',
                          'stALE.E.coli.A5.F42.I1.R1',
                          'stALE.E.coli.A6.F21.I1.R1',
                          'stALE.E.coli.A6.F43.I1.R1',
                          'stALE.E.coli.A7.F21.I1.R1',
                          'stALE.E.coli.A7.F42.I1.R1',
                          'stALE.E.coli.A8.F20.I1.R1',
                          'stALE.E.coli.A8.F42.I1.R1',
                          'stALE.E.coli.A9.F21.I1.R1',
                          'stALE.E.coli.A9.F44.I1.R1',
                          'stALE.E.coli.A10.F21.I1.R1',
                          'stALE.E.coli.A10.F43.I1.R1',
                          'stALE.E.coli.A10.F131.I1.R1',
                          'stALE.E.coli.A11.F21.I1.R1',
                          'stALE.E.coli.A11.F43.I1.R1',
                          'stALE.E.coli.A11.F119.I1.R1',
                          'stALE.E.coli.A12.F21.I1.R1',
                          'stALE.E.coli.A12.F43.I1.R1',
                          'stALE.E.coli.A12.F136.I1.R1',
                          'stALE.E.coli.A13.F20.I1.R1',
                          'stALE.E.coli.A13.F42.I1.R1',
                          'stALE.E.coli.A13.F121.I1.R1',
                          'stALE.E.coli.A14.F20.I1.R1',
                          'stALE.E.coli.A14.F42.I1.R1',
                          'stALE.E.coli.A14.F133.I1.R1',
                          'stALE.E.coli.A15.F21.I1.R1',
                          'stALE.E.coli.A15.F42.I1.R1',
                          'stALE.E.coli.A15.F117.I1.R1',
                          'stALE.E.coli.A16.F20.I1.R1',
                          'stALE.E.coli.A16.F42.I1.R1',
                          'stALE.E.coli.A16.F134.I1.R1',
                          'stALE.E.coli.A17.F21.I1.R1',
                          'stALE.E.coli.A17.F118.I1.R1',
                          'stALE.E.coli.A18.F18.I1.R1',
                          'stALE.E.coli.A18.F39.I1.R1',
                          'stALE.E.coli.A18.F130.I1.R1', '3A', '4A',
                          'BLANK.40.12G', 'BLANK.40.12H',
                          'Pputida.JBEI.HGL.Pputida.107.BP6',
                          'Pputida.JBEI.HGL.Pputida.108.BP7',
                          'Pputida.JBEI.HGL.Pputida.109.BP8',
                          'Pputida.JBEI.HGL.Pputida.110.M2',
                          'Pputida.JBEI.HGL.Pputida.111.M5',
                          'Pputida.TALE.HGL.Pputida.112',
                          'Pputida.TALE.HGL.Pputida.113',
                          'Pputida.TALE.HGL.Pputida.114',
                          'Pputida.TALE.HGL.Pputida.115',
                          'Pputida.TALE.HGL.Pputida.116',
                          'Pputida.TALE.HGL.Pputida.117',
                          'Pputida.TALE.HGL.Pputida.118',
                          'Pputida.TALE.HGL.Pputida.119',
                          'Pputida.TALE.HGL.Pputida.120',
                          'Pputida.TALE.HGL.Pputida.121',
                          'Pputida.TALE.HGL.Pputida.122',
                          'Pputida.TALE.HGL.Pputida.123',
                          'Pputida.TALE.HGL.Pputida.124',
                          'Pputida.TALE.HGL.Pputida.125',
                          'Pputida.TALE.HGL.Pputida.126',
                          'Pputida.TALE.HGL.Pputida.127',
                          'Pputida.TALE.HGL.Pputida.128',
                          'Pputida.TALE.HGL.Pputida.129',
                          'Pputida.TALE.HGL.Pputida.130',
                          'Pputida.TALE.HGL.Pputida.131',
                          'Pputida.TALE.HGL.Pputida.132',
                          'Pputida.TALE.HGL.Pputida.133',
                          'Pputida.TALE.HGL.Pputida.134',
                          'Pputida.TALE.HGL.Pputida.135',
                          'Pputida.TALE.HGL.Pputida.136',
                          'Pputida.TALE.HGL.Pputida.137',
                          'Pputida.TALE.HGL.Pputida.138',
                          'Pputida.TALE.HGL.Pputida.139',
                          'Pputida.TALE.HGL.Pputida.140',
                          'Pputida.TALE.HGL.Pputida.141',
                          'Pputida.TALE.HGL.Pputida.142',
                          'Pputida.TALE.HGL.Pputida.143',
                          'Pputida.TALE.HGL.Pputida.144',
                          'Pputida.PALE.HGL.Pputida.145',
                          'Pputida.PALE.HGL.Pputida.146',
                          'Pputida.PALE.HGL.Pputida.147',
                          'Pputida.PALE.HGL.Pputida.148',
                          'Pputida.PALE.HGL.Pputida.149',
                          'Pputida.PALE.HGL.Pputida.150',
                          'Pputida.PALE.HGL.Pputida.151',
                          'Pputida.PALE.HGL.Pputida.152',
                          'Pputida.PALE.HGL.Pputida.153',
                          'Pputida.PALE.HGL.Pputida.154',
                          'Pputida.PALE.HGL.Pputida.155',
                          'Pputida.PALE.HGL.Pputida.156',
                          'Pputida.PALE.HGL.Pputida.157',
                          'Pputida.PALE.HGL.Pputida.158',
                          'Pputida.PALE.HGL.Pputida.159',
                          'Pputida.PALE.HGL.Pputida.160',
                          'Pputida.PALE.HGL.Pputida.161',
                          'Pputida.PALE.HGL.Pputida.162',
                          'Pputida.PALE.HGL.Pputida.163',
                          'Pputida.PALE.HGL.Pputida.164',
                          'Pputida.PALE.HGL.Pputida.165',
                          'Pputida.PALE.HGL.Pputida.166',
                          'Pputida.PALE.HGL.Pputida.167',
                          'Pputida.PALE.HGL.Pputida.168',
                          'Pputida.PALE.HGL.Pputida.169',
                          'Pputida.PALE.HGL.Pputida.170',
                          'Pputida.PALE.HGL.Pputida.171',
                          'Pputida.PALE.HGL.Pputida.172',
                          'Pputida.PALE.HGL.Pputida.173',
                          'Pputida.PALE.HGL.Pputida.174',
                          'Pputida.PALE.HGL.Pputida.175',
                          'Pputida.PALE.HGL.Pputida.176',
                          'JM-Metabolic.GN0.2005', 'JM-Metabolic.GN0.2007',
                          'JM-Metabolic.GN0.2009', 'JM-Metabolic.GN0.2094',
                          'JM-Metabolic.GN0.2099', 'JM-Metabolic.GN0.2148',
                          'JM-Metabolic.GN0.2165', 'JM-Metabolic.GN0.2169',
                          'JM-Metabolic.GN0.2172', 'JM-Metabolic.GN0.2175',
                          'JM-Metabolic.GN0.2183', 'JM-Metabolic.GN0.2215',
                          'JM-Metabolic.GN0.2254', 'JM-Metabolic.GN0.2277',
                          'JM-Metabolic.GN0.2290', 'JM-Metabolic.GN0.2337',
                          'JM-Metabolic.GN0.2317', 'JM-Metabolic.GN0.2354',
                          'JM-Metabolic.GN0.2375', 'JM-Metabolic.GN0.2380',
                          'JM-Metabolic.GN0.2393', 'JM-Metabolic.GN0.2404',
                          '5B', '6A', 'BLANK.41.12G', 'BLANK.41.12H',
                          'Deoxyribose.PALE.ALE.MG1655.BOP27.4.14',
                          'Deoxyribose.PALE.ALE.MG1655.BOP27.4.23',
                          'Deoxyribose.PALE.ALE.MG1655.BOP27.4.48',
                          'Deoxyribose.PALE.ALE.MG1655.BOP27.6.21',
                          'Deoxyribose.PALE.ALE.MG1655.BOP27.6.35',
                          'Deoxyribose.PALE.ALE.MG1655.BOP27.10.13',
                          'Deoxyribose.PALE.ALE.MG1655.BOP27.10.28',
                          'Deoxyribose.PALE.ALE.MG1655.BOP27.10.51',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.18.19',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.18.59',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.18.35',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.20.16',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.20.43',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.20.71',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.22.16',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.22.28',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.22.52',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.24.9',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.24.24',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.24.52',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.26.6',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.26.27',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.26.69',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.28.13',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.28.28',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.28.53',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.30.7',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.30.22',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.30.60',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.32.6',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.32.20',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.32.56',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.1.24',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.1.57',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.1.69',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.3.23',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.3.50',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.3.61',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.5.22',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.5.36',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.5.46',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.7.23',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.7.41',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.7.51',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.17.25',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.17.58',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.17.64',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.19.25',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.19.55',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.19.63',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.21.23',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.21.46',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.21.51',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.29.25',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.29.49',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.29.57',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.31.24',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.31.42',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.31.62',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.33.21',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.33.41',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.33.50',
                          'JM-Metabolic.GN02514', 'JM-Metabolic.GN02529',
                          'JM-Metabolic.GN02531', 'JM-Metabolic.GN02567',
                          'JM-Metabolic.GN02590', 'JM-Metabolic.GN02657',
                          'JM-Metabolic.GN02748', 'JM-Metabolic.GN02766',
                          'JM-Metabolic.GN02769', 'JM-Metabolic.GN02787',
                          'JM-Metabolic.GN03132', 'JM-Metabolic.GN03218',
                          'JM-Metabolic.GN03252', 'JM-Metabolic.GN03409',
                          'JM-Metabolic.GN04014', 'JM-Metabolic.GN04094',
                          'JM-Metabolic.GN04255', 'JM-Metabolic.GN04306',
                          'JM-Metabolic.GN04428', 'JM-Metabolic.GN04488',
                          'JM-Metabolic.GN04540', 'JM-Metabolic.GN04563',
                          'JM-Metabolic.GN04612', 'JM-Metabolic.GN04665',
                          'JM-Metabolic.GN04682', 'JM-Metabolic.GN05002',
                          'JM-Metabolic.GN05109', 'JM-Metabolic.GN05128',
                          'JM-Metabolic.GN05367', 'JM-Metabolic.GN05377',
                          '7A', '8A', 'BLANK.42.12G', 'BLANK.42.12H',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0326',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0327',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0328',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0329',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0330',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0352',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0353',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0354',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0355',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0356',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0357',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0364',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0366',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0367',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0368',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0369',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0370',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0371',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0372',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0373',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0374',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0375',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0376',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0377',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0378',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0380',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0381',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0382',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0383',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0384',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0385',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0386',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0387',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0388',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0389',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0390',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0391',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0392',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0393',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0394',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0395',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0396',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0397',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0398',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0399',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0400',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0401',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0402',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0403',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0404',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0405',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0406',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0407',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0408',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0409',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0417',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0418',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0419',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0420',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0421',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0473',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0474',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0483',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0484',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0485',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0486',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0516',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0517',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0518',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0519',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0520',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0521',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0522',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0523',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0524',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0525',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-R08624',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-R08704',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-R10727',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-R11044',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-R11078',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-R11101',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-R11102',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-R11103',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-R11135',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-R11153',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-R11154',
                          'JM-Metabolic.GN02424', 'JM-Metabolic.GN02446',
                          'JM-Metabolic.GN02449', 'JM-Metabolic.GN02487',
                          'JM-Metabolic.GN02501', 'ISB', 'GFR',
                          'BLANK.43.12G', 'BLANK.43.12H',
                          'RMA.KHP.rpoS.Mage.Q97D', 'RMA.KHP.rpoS.Mage.Q97L',
                          'RMA.KHP.rpoS.Mage.Q97N', 'RMA.KHP.rpoS.Mage.Q97E',
                          'JBI.KHP.HGL.021', 'JBI.KHP.HGL.022',
                          'JBI.KHP.HGL.023', 'JBI.KHP.HGL.024',
                          'JBI.KHP.HGL.025', 'JBI.KHP.HGL.026',
                          'JBI.KHP.HGL.027', 'JBI.KHP.HGL.028.Amitesh.soxR',
                          'JBI.KHP.HGL.029.Amitesh.oxyR',
                          'JBI.KHP.HGL.030.Amitesh.soxR.oxyR',
                          'JBI.KHP.HGL.031.Amitesh.rpoS', 'BLANK1.1A',
                          'BLANK1.1B', 'BLANK1.1C', 'BLANK1.1D', 'BLANK1.1E',
                          'BLANK1.1F', 'BLANK1.1G', 'BLANK1.1H', 'AP581451B02',
                          'EP256645B01', 'EP112567B02', 'EP337425B01',
                          'LP127890A01', 'EP159692B04', 'EP987683A01',
                          'AP959450A03', 'SP464350A04', 'C9', 'ep256643b01',
                          'EP121011B01', 'AP616837B04', 'SP506933A04',
                          'EP159695B01', 'EP256644B01', 'SP511289A02',
                          'EP305735B04', 'SP415030A01', 'AP549681B02',
                          'AP549678B01', 'EP260544B04', 'EP202452B01',
                          'EP282276B04', 'SP531696A04', 'SP515443A04',
                          'SP515763A04', 'EP184255B04', 'SP503615A02',
                          'EP260543B04', 'EP768748A04', 'AP309872B03',
                          'AP568785B04', 'EP721390A04', 'EP940013A01',
                          'EP291979B04', 'EP182065B04', 'EP128904B02',
                          'EP915769A04', 'SP464352A03', 'SP365864A04',
                          'SP511294A04', 'EP061002B01', 'SP410793A01',
                          'SP232077A04', 'EP128910B01', 'AP531397B04',
                          'EP043583B01', 'EP230245B01', 'EP606652B04',
                          'EP207041B01', 'EP727972A04', 'EP291980B04',
                          'EP087938B02', 'SP471496A04', 'SP573823A04',
                          'EP393718B01', 'SP612496A01', 'EP032410B02',
                          'EP073216B01', 'EP410046B01', 'SP561451A04',
                          'EP320438B01', 'SP612495A04', 'EP446604B03',
                          'EP446602B01', 'EP182243B02', 'EP333541B04',
                          'EP238034B01', 'AP298002B02', 'EP455759B04',
                          'EP207042B04', 'LP128479A01', 'LP128476A01',
                          'EP316863B03', 'C20', 'lp127896a01', 'SP491907A02',
                          'EP182060B03', 'EP422407B01', 'SP573859A04',
                          'SP584547A02', 'EP182346B04', 'AP668631B04',
                          'EP451428B04', 'LP128538A01', 'SP490298A02',
                          'SP573860A01', 'EP032412B02', 'EP163771B01',
                          'LP169879A01', 'EP729433A02', 'EP447940B04',
                          'SP584551A08', 'EP216516B04', 'EP023808B02',
                          'BLANK2.2A', 'BLANK2.2B', 'BLANK2.2C', 'BLANK2.2D',
                          'BLANK2.2E', 'BLANK2.2F', 'BLANK2.2G', 'BLANK2.2H',
                          'SP573843A04', 'EP683835A01', 'SP573824A04',
                          'SP335002A04', 'SP478193A02', 'SP232311A04',
                          'SP415021A02', 'SP231630A02', 'SP641029A02',
                          'SP232310A04', 'EP617442B01', 'EP587478B04',
                          'EP447928B04', 'EP587475B04', 'EP675042B01',
                          'EP554513B02', 'EP702221B04', 'AP568787B02',
                          'EP054632B01', 'EP121013B01', 'EP649418A02',
                          'EP573313B01', 'LP154981A01', 'AP470859B01',
                          'LP154986A01', 'AP732307B04', 'EP533426B03',
                          'EP587476B04', 'AP696363B02', 'EP587477B04',
                          'SP683466A02', 'EP554518B04', 'EP533429B04',
                          'EP431570B01', 'EP202095B04', 'EP504030B04',
                          'EP207036B01', 'EP393717B01', 'SP491898A02',
                          'EP484973B04', 'EP479794B02', 'EP554515B04',
                          'SP631994A04', 'EP921593A04', 'AP787247B04',
                          'EP090129B04', 'EP447975B02', 'EP212214B01',
                          'EP410042B01', 'SP404409A02', 'SP247340A04',
                          'AP029018B01', 'EP872341A01', 'AP062219B03',
                          'EP790020A02', 'EP808112A04', 'SP404403A02',
                          'EP073160B01', 'EP012991B03', 'SP317297A02',
                          'EP656055A04', 'EP649623A01', 'EP790019A01',
                          'SP257519A04', 'EP808104A01', 'EP808106A01',
                          'SP231629A02', 'EP675044A01', 'EP657260A01',
                          'EP808110A04', 'AP032413B04', 'EP843906A04',
                          'AP173305B04', 'SP231628A02', 'AP173301B04',
                          'SP404405A02', 'EP649653A04', 'EP718687A04',
                          'AP905750A02', 'EP738468A01', 'C6', 'EP890157A02',
                          'SP353893A02', 'EP944059A02', 'EP970005A01',
                          'EP927461A04', 'EP808111A03', 'EP927459A04',
                          'SP317293A02', 'SP235186A04', 'SP399724A04',
                          'EP738469A01', 'SP284095A03', 'C5', 'EP337325B04',
                          'EP759450A04', 'BLANK3.3A', 'BLANK3.3B', 'BLANK3.3C',
                          'BLANK3.3D', 'BLANK3.3E', 'BLANK3.3F', 'BLANK3.3G',
                          'BLANK3.3H', 'AP006367B02', 'EP929277A02',
                          'AP324642B04', 'EP786631A04', 'EP657385A04',
                          'SP235189A01', 'EP448041B04', 'SP231631A02',
                          'SP280481A02', 'AP032412B04', 'EP649737A03',
                          'AP967057A04', 'EP876243A04', 'SP229387A04',
                          'EP667743A04', 'SP246941A01', 'AP745799A04',
                          'SP205732A02', 'SP230382A04', 'SP230380A02',
                          'SP230381A01', 'SP205754A01', 'EP606662B04',
                          'AP780167B02', 'EP447927B04', 'C18', 'LP191039A01',
                          'EP606663B04', 'EP573296B01', 'EP447926B04',
                          'LP127767A01', 'EP479266B04', 'LP128543A01',
                          'EP479270B03', 'EP921594A04', 'EP554501B04',
                          'EP542577B04', 'EP487995B04', 'EP542578B04',
                          'EP573310B01', 'EP244366B01', 'EP533389B03',
                          'EP244360B01', 'AP911328B01', 'AP481403B02',
                          '22.001.801.552.503.00', 'EP372981B04',
                          'EP447929B04', 'SP573849A04', 'SP577399A02',
                          'EP606656B03', 'LP166715A01', 'AP668628B04', 'C14',
                          'EP446610B02', 'EP339061B02', 'SP681591A04',
                          'EP393712B02', 'EP410041B01', 'SP453872A01',
                          '22.001.710.503.791.00', 'LP128540A01',
                          'EP339053B02', 'EP617443B01', 'EP190307B01',
                          'AP795068B04', 'LP128541A01', 'EP584756B04',
                          'SP284096A02', 'EP431562B04', 'EP685640B01',
                          'EP339059B02', 'EP431575B01', 'EP379938B01',
                          'EP529635B02', 'EP554506B04', 'EP455757B04',
                          'SP491900A02', 'LP196272A01', 'SP704319A04',
                          'EP617441B01', 'AP687591B04', 'SP640978A02',
                          'EP981129A02', 'EP455763B04', 'EP339057B02',
                          'SP491897A02', 'EP980752B04', 'LP128539A01',
                          'EP996831B04', 'EP273332B04', 'EP483291B04',
                          'EP393715B01', 'EP617440B01', 'EP729434A01',
                          'SP645141A03', 'BLANK4.4A', 'BLANK4.4B', 'BLANK4.4C',
                          'BLANK4.4D', 'BLANK4.4E', 'BLANK4.4F', 'BLANK4.4G',
                          'BLANK4.4H', 'SP232114A04', 'EP393714B01',
                          'EP533388B01', 'EP724905B01', 'EP282108B01',
                          'EP282107B01', 'EP001625B01', 'EP073209B02',
                          'SP232079A01', 'EP772145A02', 'AP771472A04',
                          'AP223470B01', 'SP404412A02', 'EP772143A02',
                          'SP408629A01', 'EP749735A07', 'EP846485A01',
                          'EP808109A01', 'SP416130A04', 'EP882752A01',
                          'AP953594A02', 'AP046324B02', 'AP891020A04',
                          'EP790023A01', 'EP657386A01', 'EP805337A01',
                          'EP927458A04', 'AP173299B04', 'EP768164A02',
                          'EP886422A01', 'AP103463B01', 'AP744361A02',
                          'AP065292B01', 'SP257517A04', 'EP790021A04',
                          'EP675075A04', 'SP388683A02', 'SP232309A01',
                          'EP899038A04', 'EP636802A01', 'AP046327B02',
                          'EP905975A04', 'SP410796A02', 'EP784608A01',
                          'EP808105A01', 'SP331134A04', 'EP718688A01',
                          'SP232270A02', 'EP970001A01', 'EP001624B01',
                          'EP868682A01', 'EP927462A02', 'C3', 'EP890158A02',
                          'EP023801B04', 'EP400447B04', 'EP385379B01',
                          'EP385387B01', 'EP385384B01', 'SP754514A04',
                          'SP415025A01', 'SP415023A02', 'EP400448B04',
                          'EP479894B04']

        # test sample-information-file generation.
        pipeline = Pipeline(self.good_config_file, self.good_run_id,
                            self.good_sample_sheet_path,
                            self.output_file_path, self.qiita_id,
                            Pipeline.METAGENOMIC_PTYPE)

        obs = pipeline.get_sample_names()

        self.assertEqual(sorted(obs), sorted(exp_sample_ids))

        exp = {'3A', '4A', '5B', '6A', 'BLANK.41.12G', '7A', '8A', 'ISB',
               'GFR'}

        obs = set(pipeline.get_sample_names('Gerwick'))
        self.assertEqual(obs, exp)

    def test_get_project_info(self):
        exp_proj_info = [
            {'project_name': 'NYU_BMS_Melanoma_13059', 'qiita_id': '13059',
             'contains_replicates': False},
            {'project_name': 'Feist_11661', 'qiita_id': '11661',
             'contains_replicates': False},
            {'project_name': 'Gerwick_6123', 'qiita_id': '6123',
             'contains_replicates': False}]

        exp_project_names = ['NYU_BMS_Melanoma_13059', 'Feist_11661',
                             'Gerwick_6123']

        # test sample-information-file generation.
        pipeline = Pipeline(self.good_config_file, self.good_run_id,
                            self.good_sample_sheet_path,
                            self.output_file_path, self.qiita_id,
                            Pipeline.METAGENOMIC_PTYPE)

        obs_proj_info = pipeline.get_project_info()

        obs_project_names = []
        for d in obs_proj_info:
            obs_project_names.append(d['project_name'])

        self.assertEqual(sorted(obs_project_names), sorted(exp_project_names))

        for exp_d in exp_proj_info:
            for obs_d in obs_proj_info:
                if obs_d['project_name'] == exp_d['project_name']:
                    self.assertDictEqual(obs_d, exp_d)
                    break

        # repeat test, but set short_names to True and confirm that the Qiita
        # IDs are not part of the project_names.
        exp_project_names = ['NYU_BMS_Melanoma', 'Feist', 'Gerwick']

        obs_proj_info = pipeline.get_project_info(short_names=True)

        obs_project_names = []
        for d in obs_proj_info:
            obs_project_names.append(d['project_name'])

        self.assertEqual(sorted(obs_project_names), sorted(exp_project_names))

        pipeline = Pipeline(self.good_config_file, self.good_run_id,
                            self.good_sheet_w_replicates,
                            self.output_file_path, self.qiita_id,
                            Pipeline.METAGENOMIC_PTYPE)

        obs_proj_info = pipeline.get_project_info()

        # assert value is boolean type True, even though the method no longer
        # silently converts string values to bool.
        self.assertTrue(obs_proj_info[0]['contains_replicates'])

    def test_configuration_profiles(self):
        pipeline = Pipeline(self.good_config_file, self.good_run_id,
                            self.good_sample_sheet_path,
                            self.output_file_path, self.qiita_id,
                            Pipeline.METAGENOMIC_PTYPE)

        obs = pipeline.config_profile['profile']

        # assert a profile matching self.good_sample_sheet_path was found.
        self.assertEqual(obs['instrument_type'], "MiSeq")
        self.assertEqual(obs['assay_type'], "Metagenomic")

        obs = obs['configuration']

        # sample novaseq 6000/metagenomic profile doesn't contain settings for
        # bcl-convert and qc. Make sure that these settings exist in the final
        # output and make sure they are set to the default values found in
        # default.json.
        self.assertEqual(obs['bcl-convert']['nprocs'], 16)

        # assert increased values over default found in novaseq 6000/
        # metagenomic profile are found in the final configuration as well.
        self.assertEqual(obs['bcl2fastq']['nodes'], 2)
        self.assertEqual(obs['bcl2fastq']['nprocs'], 62)
        self.assertEqual(obs['nu-qc']['nodes'], 2)
        self.assertEqual(obs['nu-qc']['wallclock_time_in_minutes'], 2028)
        self.assertEqual(obs['nu-qc']['cpus_per_task'], 32)

    def test_parse_project_name(self):
        # test sample-information-file generation.
        pipeline = Pipeline(self.good_config_file, self.good_run_id,
                            self.good_sample_sheet_path,
                            self.output_file_path, self.qiita_id,
                            Pipeline.METAGENOMIC_PTYPE)

        tests = {
            'True': [('NYU_BMS_Melanoma_13059', ('NYU_BMS_Melanoma', '13059')),
                     ('Feist_11661', ('Feist', '11661')),
                     ('Gerwick_6123', ('Gerwick', '6123')),
                     ('bar.baz_123', ('bar.baz', '123')),
                     ('Foobar', None),
                     ('', None),
                     (None, None)],
            'False': [('NYU_BMS_Mel_13059', ('NYU_BMS_Mel_13059', '13059')),
                      ('Feist_11661', ('Feist_11661', '11661')),
                      ('Gerwick_6123', ('Gerwick_6123', '6123')),
                      ('bar.baz_123', ('bar.baz_123', '123')),
                      ('Foobar', None),
                      ('', None),
                      (None, None)]
            }

        for t_set in tests:
            for test, exp in tests[t_set]:
                if exp is None:
                    with self.assertRaises(ValueError):
                        pipeline._parse_project_name(test, t_set == 'True')
                else:
                    obs = pipeline._parse_project_name(test, t_set == 'True')
                    self.assertEqual(obs, exp)

    def test_identify_reserved_words(self):
        pipeline = Pipeline(self.good_config_file, self.good_run_id,
                            self.good_sample_sheet_path,
                            self.output_file_path, self.qiita_id,
                            Pipeline.METAGENOMIC_PTYPE)

        # assert that arbitrary strings are not reserved.
        obs = pipeline.identify_reserved_words(['NOT_A_RESERVED_WORD',
                                                'ANOTHER_WORD'])
        self.assertEqual(obs, [])

        # assert that 'well_id_384' is a reserved word.
        obs = pipeline.identify_reserved_words(['well_id_384',
                                                'NOT_A_RESERVED_WORD'])

        self.assertEqual(obs, ['well_id_384'])

        # create new pipeline using a/legacy (v90) metagenomic sample-sheet.
        pipeline = Pipeline(self.good_config_file, self.good_run_id,
                            self.good_legacy_sheet_path,
                            self.output_file_path, self.qiita_id,
                            Pipeline.METAGENOMIC_PTYPE)

        # assert that for legacy sample-sheets, well_id_384 is NOT a reserved
        # word and the appropriate reserved word is 'Sample_well'.
        obs = pipeline.identify_reserved_words(['well_id_384',
                                                'NOT_A_RESERVED_WORD',
                                                'Sample_well',
                                                'Sample_Well'])

        self.assertEqual(obs, ['sample_well'])


class TestAmpliconPipeline(unittest.TestCase):
    def setUp(self):
        package_root = abspath('./sequence_processing_pipeline')
        self.path = partial(join, package_root, 'tests', 'data')
        self.good_config_file = join(package_root, 'configuration.json')
        self.bad_config_file = self.path('bad_configuration.json')
        self.invalid_config_file = 'does/not/exist/configuration.json'
        self.good_run_id = '211021_A00000_0000_SAMPLE'
        # qiita_id is randomly-generated and does not match any known
        # existing qiita job_id.
        self.qiita_id = '077c4da8-74eb-4184-8860-0207f53623be'
        self.invalid_run_id = 'not-sample-sequence-directory'
        self.output_file_path = self.path('output_dir')
        makedirs(self.output_file_path, exist_ok=True)
        self.maxDiff = None
        self.good_mapping_file_path = self.path('good-mapping-file.txt')
        self.mf_missing_column = self.path('mf-missing-column.txt')
        self.mf_duplicate_sample = self.path('mf-duplicate-sample.txt')
        self.good_run_dir = self.path(self.good_run_id)
        self.runinfo_file = self.path(self.good_run_id, 'RunInfo.xml')
        self.rtacomplete_file = self.path(self.good_run_id, 'RTAComplete.txt')
        self.sample_sheet_path = self.path('good-sample-sheet.csv')

        # most of the tests here were written with the assumption that these
        # files already exist.
        self.create_runinfo_file()
        self.create_rtacomplete_file()

        # read good configuration file at initialization to avoid putting
        # most of the test code within 'with open()' expressions.
        with open(self.good_config_file, 'r') as f:
            self.good_config = json.load(f)

    def tearDown(self):
        # Pipeline is now the only class aware of these files, hence they
        # can be deleted at the end of testing.
        self.delete_runinfo_file()
        self.delete_rtacomplete_file()

    def make_runinfo_file_unreadable(self):
        os.chmod(self.runinfo_file, 0o000)

    def make_runinfo_file_readable(self):
        os.chmod(self.runinfo_file, 0o777)

    def create_runinfo_file(self, four_reads=False):
        # since good sample RunInfo.xml files already exist to support
        # other tests, reuse them here.

        f_name = 'RunInfo_Good2.xml' if four_reads else 'RunInfo_Good1.xml'
        copy(join('sequence_processing_pipeline/tests/data/', f_name),
             self.runinfo_file)

    def delete_runinfo_file(self):
        try:
            os.remove(self.runinfo_file)
        except FileNotFoundError:
            # make method idempotent
            pass

    def create_rtacomplete_file(self):
        with open(self.rtacomplete_file, 'w') as f:
            f.write("")

    def delete_rtacomplete_file(self):
        try:
            os.remove(self.rtacomplete_file)
        except FileNotFoundError:
            # make method idempotent
            pass

    def test_required_file_checks(self):
        # begin this test by deleting the RunInfo.txt file and verifying that
        # Pipeline object will raise an Error.
        self.delete_runinfo_file()

        with self.assertRaisesRegex(PipelineError, "required file 'RunInfo.xml"
                                                   "' is not present."):
            Pipeline(self.good_config_file, self.good_run_id,
                     self.good_mapping_file_path,
                     self.output_file_path,
                     self.qiita_id, Pipeline.AMPLICON_PTYPE)

        # delete RTAComplete.txt and recreate RunInfo.txt file to verify that
        # an Error is raised when only RTAComplete.txt is missing.
        self.delete_rtacomplete_file()
        self.create_runinfo_file()

        with self.assertRaisesRegex(PipelineError, "required file 'RTAComplete"
                                                   ".txt' is not present."):
            Pipeline(self.good_config_file, self.good_run_id,
                     self.good_mapping_file_path,
                     self.output_file_path,
                     self.qiita_id, Pipeline.AMPLICON_PTYPE)

        # make RunInfo.xml file unreadable and verify that Pipeline object
        # raises the expected Error.
        self.create_rtacomplete_file()
        self.make_runinfo_file_unreadable()

        with self.assertRaisesRegex(PipelineError, "RunInfo.xml is present, "
                                                   "but not readable"):
            Pipeline(self.good_config_file, self.good_run_id,
                     self.good_mapping_file_path, self.output_file_path,
                     self.qiita_id, Pipeline.AMPLICON_PTYPE)
            self.make_runinfo_file_readable()

    def test_creation(self):
        # Pipeline should assert due to config_file
        with self.assertRaises(PipelineError) as e:
            Pipeline(self.bad_config_file,
                     self.good_run_id,
                     self.good_mapping_file_path,
                     self.output_file_path,
                     self.qiita_id, Pipeline.AMPLICON_PTYPE)

        msg = re.sub(r'not a key in .*?/sequence_processing_pipeline',
                     r'not a key in sequence_processing_pipeline',
                     str(e.exception))
        self.assertEqual(msg, "'search_paths' is not a key in "
                              "sequence_processing_pipeline/tests"
                              "/data/bad_configuration.json")

        # Pipeline should assert due to an invalid config file path.
        with self.assertRaises(PipelineError) as e:
            Pipeline(self.invalid_config_file,
                     self.good_run_id,
                     self.good_mapping_file_path,
                     self.output_file_path,
                     self.qiita_id, Pipeline.AMPLICON_PTYPE)

        self.assertEqual(str(e.exception), 'does/not/exist/configuration.json '
                                           'does not exist.')

        # Pipeline should assert on config_file = None
        with self.assertRaises(PipelineError) as e:
            Pipeline(None,
                     self.good_run_id,
                     self.good_mapping_file_path,
                     self.output_file_path,
                     self.qiita_id, Pipeline.AMPLICON_PTYPE)

        self.assertEqual(str(e.exception), 'configuration_file_path cannot be '
                                           'None')

        # Pipeline should assert due to invalid_run_id
        with self.assertRaises(PipelineError) as e:
            Pipeline(self.good_config_file,
                     self.invalid_run_id,
                     self.good_mapping_file_path,
                     self.output_file_path,
                     self.qiita_id, Pipeline.AMPLICON_PTYPE)

        self.assertEqual(str(e.exception), "A run-dir for 'not-sample-sequence"
                                           "-directory' could not be found")

        # Pipeline should assert on run_id = None
        with self.assertRaises(PipelineError) as e:
            Pipeline(self.good_config_file,
                     None,
                     self.good_mapping_file_path,
                     self.output_file_path,
                     self.qiita_id, Pipeline.AMPLICON_PTYPE)

    def test_mapping_file_validation(self):
        # test successful validation of a good mapping-file.
        try:
            Pipeline(self.good_config_file, self.good_run_id,
                     self.good_mapping_file_path,
                     self.output_file_path,
                     self.qiita_id, Pipeline.AMPLICON_PTYPE)
        except PipelineError as e:
            self.fail(("test_filter_directories_for_time failed w/PipelineEr"
                       f"ror: {e.message}"))

        # test unsuccessful validation of a bad mapping-file.
        with self.assertRaises(PipelineError) as e:
            Pipeline(self.good_config_file, self.good_run_id,
                     self.mf_missing_column,
                     self.output_file_path,
                     self.qiita_id, Pipeline.AMPLICON_PTYPE)
        self.assertEqual(str(e.exception), ('Mapping-file is missing '
                                            'columns: tm10_8_tool, '
                                            'tm50_8_tool, well_id_384'))

        # test unsuccessful validation of a bad mapping-file.
        with self.assertRaises(PipelineError) as e:
            Pipeline(self.good_config_file, self.good_run_id,
                     self.mf_duplicate_sample,
                     self.output_file_path,
                     self.qiita_id, Pipeline.AMPLICON_PTYPE)
        self.assertEqual(str(e.exception), ("Mapping-file contains duplicate "
                                            "columns: ['barcode', 'Barcode']"))

    def test_is_mapping_file(self):
        # assert that a good sample-sheet is not a mapping-file
        self.assertFalse(Pipeline.is_mapping_file(self.sample_sheet_path))
        # assert that a good mapping-file returns correctly
        self.assertTrue(Pipeline.is_mapping_file(self.good_mapping_file_path))
        # assert that a mapping-file w/duplicate samples is still considered
        # a mapping file.
        self.assertTrue(Pipeline.is_mapping_file(self.mf_duplicate_sample))
        # assert that a mapping-file w/a missing columbn is still considered
        # a mapping file.
        self.assertTrue(Pipeline.is_mapping_file(self.mf_missing_column))

    def test_is_sample_sheet(self):
        # assert that a good sample-sheet returns correctly.
        self.assertTrue(Pipeline.is_sample_sheet(self.sample_sheet_path))
        # assert that a good mapping-file is not a sample-sheet.
        self.assertFalse(Pipeline.is_sample_sheet(self.good_mapping_file_path))

    def test_generate_sample_information_files(self):
        # test sample-information-file generation.
        pipeline = Pipeline(self.good_config_file, self.good_run_id,
                            self.good_mapping_file_path,
                            self.output_file_path,
                            self.qiita_id,
                            Pipeline.AMPLICON_PTYPE)
        paths = pipeline.generate_sample_info_files()

        # confirm file exists in the expected location and with the expected
        # filename.
        obs = [x.split('sequence_processing_pipeline/')[1] for x in paths]

        exp = [(f'tests/data/output_dir/{self.good_run_id}'
                '_ABTX_20230208_ABTX_11052_blanks.tsv')]

        # sort the lists to ensure both are in a fixed order.
        obs.sort()
        exp.sort()

        self.assertEqual(obs, exp)

        # generate_sample_information_file() remains 95% the same as before.
        # a few lines of code at the beginning of the method generate lists
        # of samples and projects from the mapping file versus a sample-
        # sheet. As long as the file is being generated, additional checks
        # aren't needed.

    def test_get_sample_ids(self):
        exp_sample_ids = ['11.1.21.RK.FH', '11.1.21.RK.LH',
                          '11.1.21.RK.RH', '11.10.21.RK.FH',
                          '11.10.21.RK.LH', '11.10.21.RK.RH',
                          '11.12.21.RK.FH', '11.12.21.RK.LH',
                          '11.12.21.RK.RH', '11.13.21.RK.FH',
                          '11.13.21.RK.LH', '11.13.21.RK.RH',
                          '11.17.21.RK.FH', '11.17.21.RK.LH',
                          '11.17.21.RK.RH', '11.2.21.RK.FH',
                          '11.2.21.RK.LH', '11.2.21.RK.RH',
                          '11.3.21.RK.FH', '11.3.21.RK.LH',
                          '11.3.21.RK.RH', '11.4.21.RK.FH',
                          '11.4.21.RK.LH', '11.4.21.RK.RH',
                          '11.5.21.RK.FH', '11.5.21.RK.LH',
                          '11.5.21.RK.RH', '11.6.21.RK.FH',
                          '11.6.21.RK.LH', '11.6.21.RK.RH',
                          '11.7.21.RK.FH', '11.7.21.RK.LH',
                          '11.7.21.RK.RH', '11.8.21.RK.FH',
                          '11.8.21.RK.LH', '5.1.22.RK.FH',
                          '5.1.22.RK.LH', '5.1.22.RK.RH',
                          '5.10.22.RK.RH', '5.11.22.RK.FH',
                          '5.11.22.RK.LH', '5.11.22.RK.RH',
                          '5.12.22.RK.FH', '5.12.22.RK.LH',
                          '5.12.22.RK.RH', '5.13.22.RK.FH',
                          '5.13.22.RK.LH', '5.13.22.RK.RH',
                          '5.14.22.RK.FH', '5.14.22.RK.LH',
                          '5.14.22.RK.RH', '5.15.22.RK.FH',
                          '5.15.22.RK.LH', '5.15.22.RK.RH',
                          '5.16.22.RK.FH', '5.16.22.RK.LH',
                          '5.16.22.RK.RH', '5.17.22.RK.FH',
                          '5.17.22.RK.LH', '5.17.22.RK.RH',
                          '5.18.22.RK.FH', '5.18.22.RK.LH',
                          '5.18.22.RK.RH', '5.19.22.RK.FH',
                          '5.19.22.RK.LH', '5.19.22.RK.RH',
                          '5.2.22.RK.FH', '5.2.22.RK.LH',
                          '5.2.22.RK.RH', '5.20.22.RK.FH',
                          '5.20.22.RK.LH', '5.20.22.RK.RH',
                          '5.21.22.RK.FH', '5.21.22.RK.LH',
                          '5.21.22.RK.RH', '5.22.22.RK.FH',
                          '5.22.22.RK.LH', '5.22.22.RK.RH',
                          '5.23.22.RK.FH', '5.23.22.RK.LH',
                          '5.23.22.RK.RH', '5.24.22.RK.FH',
                          '5.24.22.RK.LH', '5.24.22.RK.RH',
                          '5.27.22.RK.FH', '5.27.22.RK.LH',
                          '5.27.22.RK.RH', '5.29.22.RK.FH',
                          '5.29.22.RK.LH', '5.29.22.RK.RH',
                          '5.3.22.RK.FH', '5.3.22.RK.LH',
                          '5.3.22.RK.RH', '5.30.22.RK.FH',
                          '5.30.22.RK.LH', '5.30.22.RK.RH',
                          '5.31.22.RK.FH', '5.31.22.RK.LH',
                          '5.31.22.RK.RH', '5.4.22.RK.FH',
                          '5.4.22.RK.LH', '5.4.22.RK.RH',
                          '5.5.22.RK.FH', '5.5.22.RK.LH',
                          '5.5.22.RK.RH', '5.6.22.RK.FH',
                          '5.6.22.RK.LH', '5.6.22.RK.RH',
                          '5.7.22.RK.FH', '5.7.22.RK.LH',
                          '5.7.22.RK.RH', '5.8.22.RK.FH',
                          '5.8.22.RK.LH', '5.8.22.RK.RH',
                          '5.9.22.RK.FH', '5.9.22.RK.LH',
                          '5.9.22.RK.RH', '6.1.22.RK.FH',
                          '6.1.22.RK.LH', '6.1.22.RK.RH',
                          '6.10.22.RK.FH', '6.10.22.RK.LH',
                          '6.10.22.RK.RH', '6.11.22.RK.FH',
                          '6.11.22.RK.LH', '6.11.22.RK.RH',
                          '6.12.22.RK.FH', '6.12.22.RK.LH',
                          '6.12.22.RK.RH', '6.13.22.RK.FH',
                          '6.13.22.RK.LH', '6.13.22.RK.RH',
                          '6.14.22.RK.FH', '6.14.22.RK.LH',
                          '6.14.22.RK.RH', '6.15.22.RK.FH',
                          '6.15.22.RK.LH', '6.15.22.RK.RH',
                          '6.16.22.RK.FH', '6.16.22.RK.LH',
                          '6.16.22.RK.RH', '6.17.22.RK.FH',
                          '6.17.22.RK.LH', '6.17.22.RK.RH',
                          '6.18.22.RK.FH', '6.18.22.RK.LH',
                          '6.18.22.RK.RH', '6.19.22.RK.FH',
                          '6.19.22.RK.LH', '6.19.22.RK.RH',
                          '6.2.22.RK.FH', '6.2.22.RK.LH',
                          '6.2.22.RK.RH', '6.20.22.RK.FH',
                          '6.20.22.RK.LH', '6.20.22.RK.RH',
                          '6.21.22.RK.FH', '6.21.22.RK.LH',
                          '6.21.22.RK.RH', '6.22.22.RK.FH',
                          '6.22.22.RK.LH', '6.22.22.RK.RH',
                          '6.23.22.RK.FH', '6.23.22.RK.LH.A',
                          '6.23.22.RK.LH.B', '6.24.22.RK.FH',
                          '6.24.22.RK.RH', '6.25.22.RK.FH',
                          '6.25.22.RK.LH', '6.25.22.RK.RH',
                          '6.26.22.RK.FH', '6.26.22.RK.LH',
                          '6.26.22.RK.RH.A', '6.26.22.RK.RH.B',
                          '6.27.22.RK.FH', '6.27.22.RK.LH',
                          '6.27.22.RK.RH', '6.28.22.RK.FH',
                          '6.28.22.RK.LH', '6.28.22.RK.RH',
                          '6.29.22.RK.FH', '6.29.22.RK.LH',
                          '6.29.22.RK.RH', '6.3.22.RK.FH',
                          '6.3.22.RK.LH', '6.3.22.RK.RH',
                          '6.30.22.RK.FH', '6.30.22.RK.LH',
                          '6.30.22.RK.RH', '6.4.22.RK.FH',
                          '6.4.22.RK.LH', '6.4.22.RK.RH',
                          '6.5.22.RK.FH', '6.5.22.RK.LH',
                          '6.5.22.RK.RH', '6.6.22.RK.FH',
                          '6.6.22.RK.LH', '6.6.22.RK.RH',
                          '6.7.22.RK.FH', '6.7.22.RK.LH',
                          '6.7.22.RK.RH', '6.8.22.RK.FH',
                          '6.8.22.RK.LH', '6.8.22.RK.RH',
                          '6.9.22.RK.FH', '6.9.22.RK.LH',
                          '6.9.22.RK.RH', '9.1.22.RK.FH',
                          '9.1.22.RK.LH', '9.1.22.RK.RH',
                          '9.10.22.RK.FH', '9.10.22.RK.LH',
                          '9.10.22.RK.RH', '9.11.22.RK.FH',
                          '9.11.22.RK.LH', '9.11.22.RK.RH',
                          '9.12.22.RK.FH', '9.12.22.RK.LH',
                          '9.12.22.RK.RH', '9.13.22.RK.FH',
                          '9.13.22.RK.LH', '9.13.22.RK.RH',
                          '9.14.22.RK.FH', '9.14.22.RK.LH',
                          '9.14.22.RK.RH', '9.15.22.RK.FH',
                          '9.15.22.RK.LH', '9.15.22.RK.RH',
                          '9.16.22.RK.FH', '9.16.22.RK.LH',
                          '9.16.22.RK.RH', '9.17.22.RK.FH',
                          '9.17.22.RK.LH', '9.17.22.RK.RH',
                          '9.19.22.RK.FH', '9.19.22.RK.LH',
                          '9.19.22.RK.RH', '9.2.22.RK.FH',
                          '9.2.22.RK.LH', '9.2.22.RK.RH',
                          '9.20.22.RK.FH', '9.20.22.RK.LH',
                          '9.20.22.RK.RH', '9.21.22.RK.FH',
                          '9.21.22.RK.LH', '9.21.22.RK.RH',
                          '9.22.22.RK.FH', '9.22.22.RK.LH',
                          '9.22.22.RK.RH', '9.23.22.RK.FH',
                          '9.23.22.RK.LH', '9.23.22.RK.RH',
                          '9.24.22.RK.FH', '9.24.22.RK.LH',
                          '9.24.22.RK.RH', '9.25.22.RK.FH',
                          '9.25.22.RK.LH', '9.26.22.RK.FH',
                          '9.26.22.RK.LH', '9.26.22.RK.RH',
                          '9.27.22.RK.FH', '9.27.22.RK.LH',
                          '9.27.22.RK.RH', '9.29.22.RK.FH',
                          '9.29.22.RK.LH', '9.29.22.RK.RH',
                          '9.3.22.RK.FH', '9.3.22.RK.LH',
                          '9.3.22.RK.RH', '9.30.22.RK.FH',
                          '9.30.22.RK.LH', '9.30.22.RK.RH',
                          '9.4.22.RK.FH', '9.4.22.RK.LH',
                          '9.4.22.RK.RH', '9.5.22.RK.FH',
                          '9.5.22.RK.LH', '9.5.22.RK.RH',
                          '9.6.22.RK.FH', '9.6.22.RK.LH',
                          '9.6.22.RK.RH', '9.7.22.RK.FH',
                          '9.7.22.RK.LH', '9.7.22.RK.RH',
                          '9.8.22.RK.FH', '9.8.22.RK.LH',
                          '9.8.22.RK.RH', '9.9.22.RK.FH',
                          '9.9.22.RK.LH', '9.9.22.RK.RH',
                          'BLANK.242.4C', 'BLANK238.3A',
                          'BLANK238.3B', 'BLANK238.3C',
                          'BLANK238.3D', 'BLANK238.3E',
                          'BLANK238.3F', 'BLANK238.3G',
                          'BLANK238.3H', 'BLANK239.10A',
                          'BLANK239.10B', 'BLANK239.10C',
                          'BLANK239.10D', 'BLANK239.10E',
                          'BLANK239.10F', 'BLANK239.10G',
                          'BLANK239.10H', 'BLANK240.3A',
                          'BLANK240.3B', 'BLANK240.3C',
                          'BLANK240.3D', 'BLANK240.3E',
                          'BLANK240.3F', 'BLANK240.3G',
                          'BLANK240.3H', 'BLANK242.10A',
                          'BLANK242.10B', 'BLANK242.10C',
                          'BLANK242.10D', 'BLANK242.10E',
                          'BLANK242.10F', 'BLANK242.10G',
                          'BLANK242.10H', 'BLANK242.11A',
                          'BLANK242.11B', 'BLANK242.11C',
                          'BLANK242.11D', 'BLANK242.11E',
                          'BLANK242.11F', 'BLANK242.11G',
                          'BLANK242.11H', 'BLANK242.12A',
                          'BLANK242.12B', 'BLANK242.12C',
                          'BLANK242.12D', 'BLANK242.12E',
                          'BLANK242.12F', 'BLANK242.12G',
                          'BLANK242.12H', 'BLANK242.4D',
                          'BLANK242.4E', 'BLANK242.4F',
                          'BLANK242.4G', 'BLANK242.4H',
                          'BLANK242.5A', 'BLANK242.5B',
                          'BLANK242.5C', 'BLANK242.5D',
                          'BLANK242.5E', 'BLANK242.5F',
                          'BLANK242.5G', 'BLANK242.5H',
                          'BLANK242.6A', 'BLANK242.6B',
                          'BLANK242.6C', 'BLANK242.6D',
                          'BLANK242.6E', 'BLANK242.6F',
                          'BLANK242.6G', 'BLANK242.6H',
                          'BLANK242.7A', 'BLANK242.7B',
                          'BLANK242.7C', 'BLANK242.7D',
                          'BLANK242.7E', 'BLANK242.7F',
                          'BLANK242.7G', 'BLANK242.7H',
                          'BLANK242.8A', 'BLANK242.8B',
                          'BLANK242.8C', 'BLANK242.8D',
                          'BLANK242.8E', 'BLANK242.8F',
                          'BLANK242.8G', 'BLANK242.8H',
                          'BLANK242.9A', 'BLANK242.9B',
                          'BLANK242.9C', 'BLANK242.9D',
                          'BLANK242.9E', 'BLANK242.9F',
                          'BLANK242.9G', 'BLANK242.9H']
        # test sample-information-file generation.
        pipeline = Pipeline(self.good_config_file,
                            self.good_run_id,
                            self.good_mapping_file_path,
                            self.output_file_path,
                            self.qiita_id,
                            Pipeline.AMPLICON_PTYPE)

        obs = pipeline.get_sample_ids()
        self.assertEqual(sorted(obs), sorted(exp_sample_ids))

    def test_get_sample_names(self):
        exp = ['11.1.21.RK.FH', '11.1.21.RK.LH', '11.1.21.RK.RH',
               '11.10.21.RK.FH', '11.10.21.RK.LH', '11.10.21.RK.RH',
               '11.12.21.RK.FH', '11.12.21.RK.LH', '11.12.21.RK.RH',
               '11.13.21.RK.FH', '11.13.21.RK.LH', '11.13.21.RK.RH',
               '11.17.21.RK.FH', '11.17.21.RK.LH', '11.17.21.RK.RH',
               '11.2.21.RK.FH', '11.2.21.RK.LH', '11.2.21.RK.RH',
               '11.3.21.RK.FH', '11.3.21.RK.LH', '11.3.21.RK.RH',
               '11.4.21.RK.FH', '11.4.21.RK.LH', '11.4.21.RK.RH',
               '11.5.21.RK.FH', '11.5.21.RK.LH', '11.5.21.RK.RH',
               '11.6.21.RK.FH', '11.6.21.RK.LH', '11.6.21.RK.RH',
               '11.7.21.RK.FH', '11.7.21.RK.LH', '11.7.21.RK.RH',
               '11.8.21.RK.FH', '11.8.21.RK.LH', '5.1.22.RK.FH',
               '5.1.22.RK.LH', '5.1.22.RK.RH', '5.10.22.RK.RH',
               '5.11.22.RK.FH', '5.11.22.RK.LH', '5.11.22.RK.RH',
               '5.12.22.RK.FH', '5.12.22.RK.LH', '5.12.22.RK.RH',
               '5.13.22.RK.FH', '5.13.22.RK.LH', '5.13.22.RK.RH',
               '5.14.22.RK.FH', '5.14.22.RK.LH', '5.14.22.RK.RH',
               '5.15.22.RK.FH', '5.15.22.RK.LH', '5.15.22.RK.RH',
               '5.16.22.RK.FH', '5.16.22.RK.LH', '5.16.22.RK.RH',
               '5.17.22.RK.FH', '5.17.22.RK.LH', '5.17.22.RK.RH',
               '5.18.22.RK.FH', '5.18.22.RK.LH', '5.18.22.RK.RH',
               '5.19.22.RK.FH', '5.19.22.RK.LH', '5.19.22.RK.RH',
               '5.2.22.RK.FH', '5.2.22.RK.LH', '5.2.22.RK.RH', '5.20.22.RK.FH',
               '5.20.22.RK.LH', '5.20.22.RK.RH', '5.21.22.RK.FH',
               '5.21.22.RK.LH', '5.21.22.RK.RH', '5.22.22.RK.FH',
               '5.22.22.RK.LH', '5.22.22.RK.RH', '5.23.22.RK.FH',
               '5.23.22.RK.LH', '5.23.22.RK.RH', '5.24.22.RK.FH',
               '5.24.22.RK.LH', '5.24.22.RK.RH', '5.27.22.RK.FH',
               '5.27.22.RK.LH', '5.27.22.RK.RH', '5.29.22.RK.FH',
               '5.29.22.RK.LH', '5.29.22.RK.RH', '5.3.22.RK.FH',
               '5.3.22.RK.LH', '5.3.22.RK.RH', '5.30.22.RK.FH',
               '5.30.22.RK.LH', '5.30.22.RK.RH', '5.31.22.RK.FH',
               '5.31.22.RK.LH', '5.31.22.RK.RH', '5.4.22.RK.FH',
               '5.4.22.RK.LH', '5.4.22.RK.RH', '5.5.22.RK.FH', '5.5.22.RK.LH',
               '5.5.22.RK.RH', '5.6.22.RK.FH', '5.6.22.RK.LH', '5.6.22.RK.RH',
               '5.7.22.RK.FH', '5.7.22.RK.LH', '5.7.22.RK.RH', '5.8.22.RK.FH',
               '5.8.22.RK.LH', '5.8.22.RK.RH', '5.9.22.RK.FH', '5.9.22.RK.LH',
               '5.9.22.RK.RH', '6.1.22.RK.FH', '6.1.22.RK.LH', '6.1.22.RK.RH',
               '6.10.22.RK.FH', '6.10.22.RK.LH', '6.10.22.RK.RH',
               '6.11.22.RK.FH', '6.11.22.RK.LH', '6.11.22.RK.RH',
               '6.12.22.RK.FH', '6.12.22.RK.LH', '6.12.22.RK.RH',
               '6.13.22.RK.FH', '6.13.22.RK.LH', '6.13.22.RK.RH',
               '6.14.22.RK.FH', '6.14.22.RK.LH', '6.14.22.RK.RH',
               '6.15.22.RK.FH', '6.15.22.RK.LH', '6.15.22.RK.RH',
               '6.16.22.RK.FH', '6.16.22.RK.LH', '6.16.22.RK.RH',
               '6.17.22.RK.FH', '6.17.22.RK.LH', '6.17.22.RK.RH',
               '6.18.22.RK.FH', '6.18.22.RK.LH', '6.18.22.RK.RH',
               '6.19.22.RK.FH', '6.19.22.RK.LH', '6.19.22.RK.RH',
               '6.2.22.RK.FH', '6.2.22.RK.LH', '6.2.22.RK.RH', '6.20.22.RK.FH',
               '6.20.22.RK.LH', '6.20.22.RK.RH', '6.21.22.RK.FH',
               '6.21.22.RK.LH', '6.21.22.RK.RH', '6.22.22.RK.FH',
               '6.22.22.RK.LH', '6.22.22.RK.RH', '6.23.22.RK.FH',
               '6.23.22.RK.LH.A', '6.23.22.RK.LH.B', '6.24.22.RK.FH',
               '6.24.22.RK.RH', '6.25.22.RK.FH', '6.25.22.RK.LH',
               '6.25.22.RK.RH', '6.26.22.RK.FH', '6.26.22.RK.LH',
               '6.26.22.RK.RH.A', '6.26.22.RK.RH.B', '6.27.22.RK.FH',
               '6.27.22.RK.LH', '6.27.22.RK.RH', '6.28.22.RK.FH',
               '6.28.22.RK.LH', '6.28.22.RK.RH', '6.29.22.RK.FH',
               '6.29.22.RK.LH', '6.29.22.RK.RH', '6.3.22.RK.FH',
               '6.3.22.RK.LH', '6.3.22.RK.RH', '6.30.22.RK.FH',
               '6.30.22.RK.LH', '6.30.22.RK.RH', '6.4.22.RK.FH',
               '6.4.22.RK.LH', '6.4.22.RK.RH', '6.5.22.RK.FH', '6.5.22.RK.LH',
               '6.5.22.RK.RH', '6.6.22.RK.FH', '6.6.22.RK.LH', '6.6.22.RK.RH',
               '6.7.22.RK.FH', '6.7.22.RK.LH', '6.7.22.RK.RH', '6.8.22.RK.FH',
               '6.8.22.RK.LH', '6.8.22.RK.RH', '6.9.22.RK.FH', '6.9.22.RK.LH',
               '6.9.22.RK.RH', '9.1.22.RK.FH', '9.1.22.RK.LH', '9.1.22.RK.RH',
               '9.10.22.RK.FH', '9.10.22.RK.LH', '9.10.22.RK.RH',
               '9.11.22.RK.FH', '9.11.22.RK.LH', '9.11.22.RK.RH',
               '9.12.22.RK.FH', '9.12.22.RK.LH', '9.12.22.RK.RH',
               '9.13.22.RK.FH', '9.13.22.RK.LH', '9.13.22.RK.RH',
               '9.14.22.RK.FH', '9.14.22.RK.LH', '9.14.22.RK.RH',
               '9.15.22.RK.FH', '9.15.22.RK.LH', '9.15.22.RK.RH',
               '9.16.22.RK.FH', '9.16.22.RK.LH', '9.16.22.RK.RH',
               '9.17.22.RK.FH', '9.17.22.RK.LH', '9.17.22.RK.RH',
               '9.19.22.RK.FH', '9.19.22.RK.LH', '9.19.22.RK.RH',
               '9.2.22.RK.FH', '9.2.22.RK.LH', '9.2.22.RK.RH', '9.20.22.RK.FH',
               '9.20.22.RK.LH', '9.20.22.RK.RH', '9.21.22.RK.FH',
               '9.21.22.RK.LH', '9.21.22.RK.RH', '9.22.22.RK.FH',
               '9.22.22.RK.LH', '9.22.22.RK.RH', '9.23.22.RK.FH',
               '9.23.22.RK.LH', '9.23.22.RK.RH', '9.24.22.RK.FH',
               '9.24.22.RK.LH', '9.24.22.RK.RH', '9.25.22.RK.FH',
               '9.25.22.RK.LH', '9.26.22.RK.FH', '9.26.22.RK.LH',
               '9.26.22.RK.RH', '9.27.22.RK.FH', '9.27.22.RK.LH',
               '9.27.22.RK.RH', '9.29.22.RK.FH', '9.29.22.RK.LH',
               '9.29.22.RK.RH', '9.3.22.RK.FH', '9.3.22.RK.LH', '9.3.22.RK.RH',
               '9.30.22.RK.FH', '9.30.22.RK.LH', '9.30.22.RK.RH',
               '9.4.22.RK.FH', '9.4.22.RK.LH', '9.4.22.RK.RH', '9.5.22.RK.FH',
               '9.5.22.RK.LH', '9.5.22.RK.RH', '9.6.22.RK.FH', '9.6.22.RK.LH',
               '9.6.22.RK.RH', '9.7.22.RK.FH', '9.7.22.RK.LH', '9.7.22.RK.RH',
               '9.8.22.RK.FH', '9.8.22.RK.LH', '9.8.22.RK.RH', '9.9.22.RK.FH',
               '9.9.22.RK.LH', '9.9.22.RK.RH', 'BLANK.242.4C', 'BLANK238.3A',
               'BLANK238.3B', 'BLANK238.3C', 'BLANK238.3D', 'BLANK238.3E',
               'BLANK238.3F', 'BLANK238.3G', 'BLANK238.3H', 'BLANK239.10A',
               'BLANK239.10B', 'BLANK239.10C', 'BLANK239.10D', 'BLANK239.10E',
               'BLANK239.10F', 'BLANK239.10G', 'BLANK239.10H', 'BLANK240.3A',
               'BLANK240.3B', 'BLANK240.3C', 'BLANK240.3D', 'BLANK240.3E',
               'BLANK240.3F', 'BLANK240.3G', 'BLANK240.3H', 'BLANK242.10A',
               'BLANK242.10B', 'BLANK242.10C', 'BLANK242.10D', 'BLANK242.10E',
               'BLANK242.10F', 'BLANK242.10G', 'BLANK242.10H', 'BLANK242.11A',
               'BLANK242.11B', 'BLANK242.11C', 'BLANK242.11D', 'BLANK242.11E',
               'BLANK242.11F', 'BLANK242.11G', 'BLANK242.11H', 'BLANK242.12A',
               'BLANK242.12B', 'BLANK242.12C', 'BLANK242.12D', 'BLANK242.12E',
               'BLANK242.12F', 'BLANK242.12G', 'BLANK242.12H', 'BLANK242.4D',
               'BLANK242.4E', 'BLANK242.4F', 'BLANK242.4G', 'BLANK242.4H',
               'BLANK242.5A', 'BLANK242.5B', 'BLANK242.5C', 'BLANK242.5D',
               'BLANK242.5E', 'BLANK242.5F', 'BLANK242.5G', 'BLANK242.5H',
               'BLANK242.6A', 'BLANK242.6B', 'BLANK242.6C', 'BLANK242.6D',
               'BLANK242.6E', 'BLANK242.6F', 'BLANK242.6G', 'BLANK242.6H',
               'BLANK242.7A', 'BLANK242.7B', 'BLANK242.7C', 'BLANK242.7D',
               'BLANK242.7E', 'BLANK242.7F', 'BLANK242.7G', 'BLANK242.7H',
               'BLANK242.8A', 'BLANK242.8B', 'BLANK242.8C', 'BLANK242.8D',
               'BLANK242.8E', 'BLANK242.8F', 'BLANK242.8G', 'BLANK242.8H',
               'BLANK242.9A', 'BLANK242.9B', 'BLANK242.9C', 'BLANK242.9D',
               'BLANK242.9E', 'BLANK242.9F', 'BLANK242.9G', 'BLANK242.9H']

        # test sample-information-file generation.
        pipeline = Pipeline(self.good_config_file,
                            self.good_run_id,
                            self.good_mapping_file_path,
                            self.output_file_path,
                            self.qiita_id,
                            Pipeline.AMPLICON_PTYPE)

        obs = pipeline.get_sample_names()
        self.assertEqual(sorted(obs), sorted(exp))

        # mapping file only contains one project, but this test will
        # still exercise the correct code.
        obs = pipeline.get_sample_names('ABTX_20230208_ABTX_11052')
        self.assertEqual(sorted(obs), sorted(exp))

    def test_get_project_info(self):
        exp_proj_info = [
            {'project_name': 'ABTX_20230208_ABTX_11052',
             'qiita_id': '11052',
             'contains_replicates': False}]

        exp_project_names = ['ABTX_20230208_ABTX_11052']

        # test sample-information-file generation.
        pipeline = Pipeline(self.good_config_file, self.good_run_id,
                            self.good_mapping_file_path,
                            self.output_file_path, self.qiita_id,
                            Pipeline.AMPLICON_PTYPE)

        obs_proj_info = pipeline.get_project_info()
        obs_project_names = []
        for d in obs_proj_info:
            obs_project_names.append(d['project_name'])

        self.assertEqual(sorted(obs_project_names), sorted(exp_project_names))

        for exp_d in exp_proj_info:
            for obs_d in obs_proj_info:
                if obs_d['project_name'] == exp_d['project_name']:
                    self.assertDictEqual(obs_d, exp_d)
                    break

    def test_dummy_sheet_generation(self):
        # generate a RunInfo.xml file w/only one indexed read.
        self.create_runinfo_file(four_reads=False)

        _ = Pipeline(self.good_config_file,
                     self.good_run_id,
                     self.good_mapping_file_path,
                     self.output_file_path,
                     self.qiita_id,
                     Pipeline.AMPLICON_PTYPE)

        dummy_sheet_path = join(self.output_file_path,
                                'dummy_sample_sheet.csv')

        with open(dummy_sheet_path) as f:
            obs = f.readlines()
            obs = [x.strip() for x in obs]
            self.assertEqual(obs, good_dummy_sheet1)

        # generate a RunInfo.xml file w/two indexed reads.
        self.create_runinfo_file(four_reads=True)

        _ = Pipeline(self.good_config_file,
                     self.good_run_id,
                     self.good_mapping_file_path,
                     self.output_file_path,
                     self.qiita_id,
                     Pipeline.AMPLICON_PTYPE)

        dummy_sheet_path = join(self.output_file_path,
                                'dummy_sample_sheet.csv')

        with open(dummy_sheet_path) as f:
            obs = f.readlines()
            obs = [x.strip() for x in obs]
            self.assertEqual(obs, good_dummy_sheet2)

        # refrain from testing dummy_sheet_generation with bad inputs.
        # include those tests when testing process_run_info_file().

    def test_process_run_info_file(self):
        pipeline = Pipeline(self.good_config_file,
                            self.good_run_id,
                            self.good_mapping_file_path,
                            self.output_file_path,
                            self.qiita_id,
                            Pipeline.AMPLICON_PTYPE)

        obs = pipeline.process_run_info_file('sequence_processing_pipeline/'
                                             'tests/data/RunInfo_Good1.xml')

        exp = [{'NumCycles': 151, 'Number': 1, 'IsIndexedRead': False},
               {'NumCycles': 8, 'Number': 2, 'IsIndexedRead': True},
               {'NumCycles': 8, 'Number': 3, 'IsIndexedRead': True},
               {'NumCycles': 151, 'Number': 4, 'IsIndexedRead': False}]

        self.assertEqual(obs, exp)

        obs = pipeline.process_run_info_file('sequence_processing_pipeline/'
                                             'tests/data/RunInfo_Good2.xml')

        exp = [{'NumCycles': 151, 'Number': 1, 'IsIndexedRead': False},
               {'NumCycles': 8, 'Number': 2, 'IsIndexedRead': True},
               {'NumCycles': 151, 'Number': 3, 'IsIndexedRead': False}]

        self.assertEqual(obs, exp)

        # a ValueError should be raised when a file that is obviously not
        # a RunInfo.XML file is passed to the method.
        with self.assertRaisesRegex(ValueError, "Cannot extract read "
                                                "information"):
            pipeline.process_run_info_file('sequence_processing_pipeline/'
                                           'tests/data/good-sample-sheet.csv')

        # other errors like an improper list of results from
        # process_run_info_file() are handled by generate_dummy_sample_sheet().
        # These are indirectly tested as generate_dummy_sample_sheet() is
        # called by Pipeline's constructor.

    def test_identify_reserved_words(self):
        pipeline = Pipeline(self.good_config_file,
                            self.good_run_id,
                            self.good_mapping_file_path,
                            self.output_file_path,
                            self.qiita_id,
                            Pipeline.AMPLICON_PTYPE)

        # assert that arbitrary strings are not reserved.
        obs = pipeline.identify_reserved_words(['NOT_A_RESERVED_WORD',
                                                'ANOTHER_WORD'])
        self.assertEqual(obs, [])

        # assert that Sample_Well is okay for current pre-prep files but
        # well_id_384 is reserved. Show that all forms of tm300_8_tool are
        # also reserved.
        obs = pipeline.identify_reserved_words(['Sample_Well',
                                                'TM300_8_Tool',
                                                'tm300_8_tool',
                                                'well_id_384'])
        self.assertEqual(set(obs), {'tm300_8_tool', 'well_id_384'})


class TestInstrumentUtils(unittest.TestCase):
    def setUp(self):
        package_root = abspath('./sequence_processing_pipeline')
        self.path = partial(join, package_root, 'tests', 'data')

    def test_instrument_utils(self):
        iutils = InstrumentUtils()

        exp = {'231108_M04586_0992_000000000-L7342': {'id': 'M04586',
                                                      'type': 'MiSeq',
                                                      'date': '2023-11-08'},
               '200320_K00180_0957_AHCYKKBBXY_PE150_Knight': {'id': 'K00180',
                                                              'type': ('HiSeq '
                                                                       '4000'),
                                                              'date': ('2020-0'
                                                                       '3-20')
                                                              },
               '20220912_FS10001773_27_BSE39218-1017': {'id': 'FS10001773',
                                                        'type': 'iSeq',
                                                        'date': '2022-09-12'},
               '231215_LH00444_0031_B222WHFLT4': {'id': 'LH00444',
                                                  'type': 'NovaSeq X Plus',
                                                  'date': '2023-12-16'},
               '190809_D00611_0709_AH3CKJBCX3_RKL0040_Feist_36-39_2': {
                   'id': 'D00611',
                   'type': 'HiSeq 2500',
                   'date': '2019-08-09'},
               '231215_A01535_0435_BH23F5DSXC': {'id': 'A01535',
                                                 'type': 'NovaSeq 6000',
                                                 'date': '2023-12-15'},
               '150629_SN1001_0511_AH5L7GBCXX': {'id': 'SN1001',
                                                 'type': 'RapidRun',
                                                 'date': '2015-06-29'}}

        run_directories = []
        for root, dirs, files in walk(self.path('sample_run_directories')):
            for run_id in dirs:
                run_directories.append((run_id, join(root, run_id)))

            # don't walk recursively. stop after first level.
            break

        for run_id, run_dir in run_directories:
            self.assertEqual(iutils.get_instrument_id(run_dir),
                             exp[run_id]['id'])
            self.assertEqual(iutils.get_instrument_type(run_dir),
                             exp[run_id]['type'])
            self.assertEqual(iutils.get_date(run_dir),
                             exp[run_id]['date'])


good_dummy_sheet1 = [
    "[Header],,,,,,", "IEMFileVersion,4,,,,,", "Date,10/27/22,,,,,",
    "Workflow,GenerateFASTQ,,,,,", "Application,FASTQ Only,,,,,",
    "Assay,TruSeq HT,,,,,", "Description,test_run,,,,,",
    "Chemistry,Amplicon,,,,,", ",,,,,,", "[Reads],,,,,,", "151,,,,,,",
    "151,,,,,,", ",,,,,,", "[Settings],,,,,,",
    "OverrideCycles,Y151;N8;N8;Y151,,,,,", "MaskShortReads,1,,,,,",
    "CreateFastqForIndexReads,1,,,,,", ",,,,,,", "[Data],,,,,,",
    "Sample_ID,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2",
    "211021_A00000_0000_SAMPLE_SMPL1,,,,,,", ",,,,,,",
    "[Bioinformatics],,,,,,",
    ("Project,ForwardAdapter,ReverseAdapter,PolyGTrimming,HumanFiltering,"
     "QiitaID,"),
    "211021_A00000_0000_SAMPLE_SMPL1,NA,NA,FALSE,FALSE,14782,", ",,,,,,",
    "[Contact],,,,,,", "Email,Sample_Project,,,,,",
    "c2cowart@ucsd.edu,SomeProject,,,,,",
    "antgonza@gmail.com,AnotherProject,,,,,", ",,,,,,"]


good_dummy_sheet2 = [
    "[Header],,,,,,", "IEMFileVersion,4,,,,,", "Date,10/27/22,,,,,",
    "Workflow,GenerateFASTQ,,,,,", "Application,FASTQ Only,,,,,",
    "Assay,TruSeq HT,,,,,", "Description,test_run,,,,,",
    "Chemistry,Amplicon,,,,,", ",,,,,,", "[Reads],,,,,,", "151,,,,,,",
    "151,,,,,,", ",,,,,,", "[Settings],,,,,,",
    "OverrideCycles,Y151;N8;Y151,,,,,", "MaskShortReads,1,,,,,",
    "CreateFastqForIndexReads,1,,,,,", ",,,,,,", "[Data],,,,,,",
    "Sample_ID,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2",
    "211021_A00000_0000_SAMPLE_SMPL1,,,,,,", ",,,,,,",
    "[Bioinformatics],,,,,,",
    ("Project,ForwardAdapter,ReverseAdapter,PolyGTrimming,HumanFiltering,"
     "QiitaID,"),
    "211021_A00000_0000_SAMPLE_SMPL1,NA,NA,FALSE,FALSE,14782,", ",,,,,,",
    "[Contact],,,,,,", "Email,Sample_Project,,,,,",
    "c2cowart@ucsd.edu,SomeProject,,,,,",
    "antgonza@gmail.com,AnotherProject,,,,,", ",,,,,,"]


if __name__ == '__main__':
    unittest.main()
