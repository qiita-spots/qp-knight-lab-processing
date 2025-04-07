from json import load as json_load
from json import loads as json_loads
from json.decoder import JSONDecodeError
from os import makedirs, listdir, walk
from os.path import join, exists, isdir, basename
from metapool import (load_sample_sheet, AmpliconSampleSheet, is_blank,
                      parse_project_name, SAMPLE_NAME_KEY, QIITA_ID_KEY,
                      PROJECT_SHORT_NAME_KEY, PROJECT_FULL_NAME_KEY,
                      CONTAINS_REPLICATES_KEY)
from metapool.plate import ErrorMessage, WarningMessage
from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import PipelineError
import logging
from re import findall, search, match
import sample_sheet
import pandas as pd
from collections import defaultdict
from datetime import datetime
from xml.etree import ElementTree as ET
from metapool.prep import PREP_MF_COLUMNS


logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

_PROJECT_NAME_KEY = 'project_name'


class InstrumentUtils():
    types = {'A': 'NovaSeq 6000', 'D': 'HiSeq 2500', 'FS': 'iSeq',
             'K': 'HiSeq 4000', 'LH': 'NovaSeq X Plus', 'M': 'MiSeq',
             'MN': 'MiniSeq',
             # SN – RapidRun which is HiSeq 2500
             'SN': 'RapidRun'}

    @staticmethod
    def get_instrument_id(run_directory):
        run_info = join(run_directory, 'RunInfo.xml')

        if not exists(run_info):
            raise ValueError(f"'{run_info}' doesn't exist")

        with open(run_info) as f:
            # Instrument element should appear in all valid RunInfo.xml
            # files.
            return ET.fromstring(f.read()).find('Run/Instrument').text

    @staticmethod
    def get_instrument_type(run_directory):
        # extract all letters at the beginning of the string, stopping
        # at the first digit.
        code = match(r"^(.*?)\d.*",
                     InstrumentUtils.get_instrument_id(run_directory))

        if code is None:
            raise ValueError("Could not determine instrument code")
        else:
            code = code.group(1)

        # map instrument code to a name string and return it, if possible.
        try:
            return InstrumentUtils.types[code]
        except KeyError:
            raise ValueError(f"Instrument code '{code}' is of unknown type")

    @staticmethod
    def get_date(run_directory):
        run_info = join(run_directory, 'RunInfo.xml')

        if not exists(run_info):
            raise ValueError(f"'{run_info}' doesn't exist")

        with open(run_info) as f:
            # Date element should appear in all valid RunInfo.xml
            # files.
            date_string = ET.fromstring(f.read()).find('Run/Date').text

        # date is recorded in RunInfo.xml in various formats. Iterate
        # through all known formats until the date is properly parsed.

        # For now, assume timestamps w/Z are not actually in 'Zulu' or
        # UTC time w/out confirming the machines were actually set for
        # and/or reporting UTC time.
        formats = ["%y%m%d", "%Y-%m-%dT%H:%M:%s", "%Y-%m-%dT%H:%M:%SZ",
                   "%m/%d/%Y %I:%M:%S %p"]

        for format in formats:
            try:
                date = datetime.strptime(date_string, format)
                return str(date.date())
            except ValueError:
                # assume ValueErrors are due to incorrect format, rather than
                # incorrect value from the XML file.
                pass

        raise ValueError(f"'{date_string}' could not be parsed")


class Pipeline:
    _CONTROLS_SIF_SUFFIX = '_blanks.tsv'

    # TODO: replace these with spp_metadata package call based on qiimp2
    sif_header = [SAMPLE_NAME_KEY, 'collection_timestamp', 'elevation',
                  'empo_1', 'empo_2', 'empo_3',
                  'empo_4', 'env_biome', 'env_feature',
                  'env_material', 'env_package', 'geo_loc_name',
                  'host_subject_id', 'latitude', 'longitude',
                  'sample_type', 'scientific_name', 'taxon_id',
                  'description', 'title', 'dna_extracted',
                  'physical_specimen_location', 'physical_specimen_remaining']

    sif_defaults = [None, None, 193,
                    'Control', 'Negative', 'Sterile water blank',
                    'Sterile water blank', 'urban biome', 'research facility',
                    'sterile water', 'misc environment', 'USA:CA:San Diego',
                    None, 32.5, -117.25,
                    'control blank', 'metagenome', 256318,
                    None, None, 'TRUE',
                    'UCSD', 'FALSE']

    mapping_file_columns = {SAMPLE_NAME_KEY, 'barcode', 'center_name',
                            'center_project_name',
                            'experiment_design_description',
                            'instrument_model',
                            'library_construction_protocol',
                            'platform', 'run_center', 'run_date', 'run_prefix',
                            'runid', 'sample_plate', 'sequencing_meth',
                            'linker', 'primer', 'primer_plate', 'well_id_384',
                            'plating', 'extractionkit_lot', 'extraction_robot',
                            'tm1000_8_tool', 'primer_date', 'mastermix_lot',
                            'water_lot', 'processing_robot', 'tm300_8_tool',
                            'tm50_8_tool', _PROJECT_NAME_KEY, 'orig_name',
                            'well_description', 'pcr_primers', 'target_gene',
                            'tm10_8_tool', 'target_subfragment', 'well_id_96'}

    METAGENOMIC_PTYPE = 'Metagenomic'
    METATRANSCRIPTOMIC_PTYPE = 'Metatranscriptomic'
    AMPLICON_PTYPE = 'Amplicon'

    pipeline_types = {METAGENOMIC_PTYPE, AMPLICON_PTYPE,
                      METATRANSCRIPTOMIC_PTYPE}

    METAGENOMIC_ATYPE = 'Metagenomic'
    METATRANSCRIPTOMIC_ATYPE = 'Metatranscriptomic'
    AMPLICON_ATYPE = 'TruSeq HT'

    assay_types = [AMPLICON_ATYPE, METAGENOMIC_ATYPE, METATRANSCRIPTOMIC_ATYPE]

    @staticmethod
    def make_sif_fname(run_id, full_project_name):
        # TODO: the problem with this structure is that there's no clear way
        #  to figure out, for an arbitrary sif fname, where the run id ends and
        #  where the project name begins because single underscores are used as
        #  internal elements in both run ids and project names :(
        return f'{run_id}_{full_project_name}{Pipeline._CONTROLS_SIF_SUFFIX}'

    @staticmethod
    def is_sif_fp(fp):
        return fp.endswith(Pipeline._CONTROLS_SIF_SUFFIX)

    # get study_id from sif_file_name ...something_14385_blanks.tsv
    @staticmethod
    def get_qiita_id_from_sif_fp(fp):
        fname = basename(fp)
        temp_name = fname.replace(Pipeline._CONTROLS_SIF_SUFFIX, '')

        # This is a kind of hacky use of parse_project_name since
        # it is passing in something that *ends in* a project name but is not
        # all a project name (see above re no clear way to get just the
        # project name out of the sif name without knowing the internal
        # details of a project name format). So this is a bit of a hack, but at
        # least it is limiting the number of places in the code that know the
        # details of project name format.
        hacky_name_pieces_dict = parse_project_name(temp_name)
        return hacky_name_pieces_dict[QIITA_ID_KEY]

    def __init__(self, configuration_file_path, run_id, input_file_path,
                 output_path, qiita_job_id, pipeline_type, lane_number=None):
        """
        Initialize Pipeline object w/configuration information.
        :param configuration_file_path: Path to configuration.json file.
        :param run_id: Used w/search_paths to locate input run_directory.
        :param input_file_path: Path to sample-sheet or pre-prep file.
        :param output_path: Path where all pipeline-generated files live.
        :param qiita_job_id: Qiita Job ID creating this Pipeline.
        :param pipeline_type: Pipeline type ('Amplicon', 'Metagenomic', etc.)
        :param lane_number: (Optional) overwrite lane_number in input_file.
        """
        if input_file_path is None:
            raise PipelineError("user_input_file_path cannot be None")

        if pipeline_type not in Pipeline.pipeline_types:
            raise PipelineError(f"'{type}' is not a valid pipeline type.")

        self.pipeline_type = pipeline_type

        self.configuration_file_path = configuration_file_path

        # along with configuration profiles, a 'general' configuration file
        # is needed to provide the search paths to the run-directories. These
        # are used to locate the run-directory specified and parse the required
        # files to select the right configuration profile.
        try:
            f = open(configuration_file_path)
            self.configuration = json_load(f)
            f.close()
        except TypeError:
            raise PipelineError('configuration_file_path cannot be None')
        except FileNotFoundError:
            raise PipelineError(f'{configuration_file_path} does not '
                                'exist.')
        except JSONDecodeError:
            raise PipelineError(f'{configuration_file_path} is not a '
                                'valid json file')

        if run_id is None:
            raise PipelineError('run_id cannot be None')

        for key in ['search_paths', 'archive_path', 'amplicon_search_paths',
                    'profiles_path']:
            if key not in self.configuration:
                raise PipelineError(f"'{key}' is not a key in "
                                    f"{self.configuration_file_path}")

        # If our extended validate() method discovers any warnings or
        # Errors, it will raise a PipelineError and return them w/in the
        # error message as a single string separated by '\n'.
        self.warnings = []
        self._directory_check(output_path, create=False)
        self.output_path = output_path
        self.run_id = run_id
        self.qiita_job_id = qiita_job_id
        self.pipeline = []
        self.assay_type = None

        # this method will catch a run directory as well as its products
        # directory, which also has the same name. Hence, return the
        # shortest matching path as that will at least return the right
        # path between the two.
        results = []

        if pipeline_type == Pipeline.AMPLICON_PTYPE:
            self.search_paths = self.configuration['amplicon_search_paths']
            self.assay_type = Pipeline.AMPLICON_ATYPE
        else:
            self.search_paths = self.configuration['search_paths']

        for search_path in self.search_paths:
            logging.debug(f'Searching {search_path} for {self.run_id}')
            for entry in listdir(search_path):
                some_path = join(search_path, entry)
                # ensure some_path never ends in '/'
                some_path = some_path.rstrip('/')
                if isdir(some_path) and some_path.endswith(self.run_id):
                    logging.debug(f'Found {some_path}')
                    results.append(some_path)

        if results:
            results.sort(key=lambda s: len(s))
            self.run_dir = results[0]
        else:
            raise PipelineError(f"A run-dir for '{self.run_id}' could not be "
                                "found")

        # required files for successful operation
        # both RTAComplete.txt and RunInfo.xml should reside in the root of
        # the run directory.
        required_files = ['RTAComplete.txt', 'RunInfo.xml']
        for some_file in required_files:
            if not exists(join(self.run_dir, some_file)):
                raise PipelineError("required file '%s' is not present." %
                                    some_file)

        # verify that RunInfo.xml file is readable.
        try:
            fp = open(join(self.run_dir, 'RunInfo.xml'))
            fp.close()
        except PermissionError:
            raise PipelineError('RunInfo.xml is present, but not readable')

        self.input_file_path = input_file_path

        if pipeline_type == Pipeline.AMPLICON_PTYPE:
            # assume input_file_path references a pre-prep (mapping) file.

            self.mapping_file = self._validate_mapping_file(input_file_path)
            # unlike _validate_sample_sheet() which returns a SampleSheet
            # object that stores the path to the file it was created from,
            # _validate_mapping_file() just returns a DataFrame. Store the
            # path to the original mapping file itself as well.

            # create dummy sample-sheet
            output_fp = join(output_path, 'dummy_sample_sheet.csv')
            self.generate_dummy_sample_sheet(self.run_dir, output_fp)
            self.dummy_sheet_path = output_fp

            # Optional lane_number parameter is ignored for Amplicon
            # runs, as the only valid value is 1.
        else:
            if lane_number is not None:
                # confirm that the lane_number is a reasonable value.
                lane_number = int(lane_number)
                if lane_number < 1 or lane_number > 8:
                    raise ValueError(f"'{lane_number}' is not a valid name"
                                     " number")

                # overwrite sample-sheet w/DFSheets processed version
                # with overwritten Lane number.
                sheet = load_sample_sheet(input_file_path)
                with open(input_file_path, 'w') as f:
                    sheet.write(f, lane=lane_number)

            # assume user_input_file_path references a sample-sheet.
            self.sample_sheet = self._validate_sample_sheet(input_file_path)
            self.mapping_file = None

        if self.assay_type is None:
            # set self.assay_type for non-amplicon types.
            assay_type = self.sample_sheet.Header['Assay']
            if assay_type not in Pipeline.assay_types:
                raise ValueError(f"'{assay_type} is not a valid Assay type")
            self.assay_type = assay_type

        self._configure_profile()

    def get_sample_sheet_path(self):
        """
        Returns path to a sample-sheet or dummy sample-sheet for amplicon runs.
        """
        if self.assay_type == Pipeline.AMPLICON_ATYPE:
            # assume self.dummy_sheet_path has been created for amplicon runs.
            return self.dummy_sheet_path
        else:
            # assume input_file_path is a sample-sheet for non-amplicon runs.
            return self.input_file_path

    def get_software_configuration(self, software):
        if software is None or software == "":
            raise ValueError(f"'{software}' is not a valid value")

        key_order = ['profile', 'configuration', software]

        config = self.config_profile

        for key in key_order:
            if key in config:
                config = config[key]
            else:
                raise PipelineError(f"'{key}' is not defined in configuration")

        return config

    def identify_reserved_words(self, words):
        '''
        Returns a list of words that should not appear as column names in any
        project referenced in the Pipeline's sample-sheet/pre-prep file.
        :param words: A list of words that may include reserved words.
        :return: A list of words that are already reserved in upper, lower,
                 and mixed cases.
        '''

        # Only strings used as column names in pre-prep files are currently
        # considered 'reserved' as loading a pre-prep file containing these
        # column names will fail if one or more of the strings already appears
        # as a column name in a study's sample metadata table.

        # This implementation assumes some understanding of metapool's impl,
        # specifically how the proper set of prep-info file columns are
        # generated. For now the functionality will be defined here as this
        # area of metapool is currently in flux.
        if self.pipeline_type == Pipeline.AMPLICON_PTYPE:
            reserved = PREP_MF_COLUMNS
        else:
            # results will be dependent on SheetType and SheetVersion of
            # the sample-sheet. Since all columns in a prep-info file are
            # lower()ed before writing out to file, the word must be
            # reserved in all case forms. e.g.: 'Sample_Well' and 'Sample_well'
            # are both forms of 'sample_well'.
            reserved = [x.lower() for x in
                        self.sample_sheet.CARRIED_PREP_COLUMNS] + \
                        self.sample_sheet.GENERATED_PREP_COLUMNS

        return list(set([x.lower() for x in words]) & set(reserved))

    def _configure_profile(self):
        # extract the instrument type from self.run_dir and the assay type
        # from self.sample_sheet (or self.mapping_file).
        instr_type = InstrumentUtils.get_instrument_type(self.run_dir)

        # open the configuration profiles directory as specified by
        # profiles_path in the configuration.json file. parse each json into
        # a nested dictionary keyed by (instrument-type, assay-type) as
        # specified by the values inside each json.
        profile_dir = self.configuration['profiles_path']

        if not exists(profile_dir):
            raise ValueError(f"'{profile_dir}' doesn't exist")

        # profiles directory can be arbitrarily nested to help organize
        # profiles; profiles can also be named arbitrarily. Non-JSON files
        # such as notes can be in the directory as well. The only assertion
        # is that all JSON files found will be of the profile format, discussed
        # below.
        profile_paths = []
        for root, dirs, files in walk(profile_dir):
            for some_file in files:
                some_path = join(root, some_file)
                if some_path.endswith('.json'):
                    profile_paths.append(some_path)

        # There must be at least one valid profile for the Pipeline to
        # continue operation.
        if not profile_paths:
            raise ValueError(f"'{profile_dir}' doesn't contain profile files")

        profiles = []

        for profile_path in profile_paths:
            with open(profile_path, 'r') as f:
                # open each profile and perform minimum validation on its
                # contents.
                contents = json_load(f)

                # all files must contain a root element 'profile'. This helps
                # to identify it as a profile, rather than another type of
                # JSON file.
                if 'profile' not in contents:
                    raise ValueError("'profile' is not an attribute in "
                                     f"'{profile_path}'")

                # the 'profile' attribute must have a dictionary as its value.
                # all profiles must contain 'instrument_type' and 'assay_type',
                if 'instrument_type' not in contents['profile']:
                    raise ValueError("'instrument_type' is not an attribute "
                                     f"in '{profile_path}'.profile")

                if 'assay_type' not in contents['profile']:
                    raise ValueError("'assay_type' is not an attribute "
                                     f"in '{profile_path}'.profile")

                profiles.append(contents)

        selected_profile = None

        for profile in profiles:
            i_type = profile['profile']['instrument_type']
            a_type = profile['profile']['assay_type']

            if i_type == instr_type and a_type == self.assay_type:
                selected_profile = profile
                break

        if selected_profile is None:
            raise ValueError(f"a matching profile ({instr_type}, "
                             f"{self.assay_type}) was not found. Please notify"
                             " an administrator")

        self.config_profile = selected_profile

    def _directory_check(self, directory_path, create=False):
        if exists(directory_path):
            logging.debug("directory '%s' exists." % directory_path)
        else:
            if create:
                try:
                    makedirs(directory_path, exist_ok=True)
                except OSError as e:
                    # this is a known potential error. Re-raise it as a
                    # PipelineError, so it gets handled in the same location
                    # as the others.
                    raise PipelineError(str(e))
            else:
                raise PipelineError("directory_path '%s' does not exist." %
                                    directory_path)

    def run(self, callback=None):
        """
        Run all jobs added to Pipeline in the order they were added.
        :param callback: Optional function to call and upstate status with.
        :param callback(jid=): string identifying the current running process.
        :param callback(status=): a string message or description.
        :return:
        """
        for job in self.pipeline:
            job.run(callback=callback)

    def add(self, job):
        """
        Add a job to the Pipeline
        :param Job: A Job object
        :return: None
        """
        if isinstance(job, Job):
            self.pipeline.append(job)
        else:
            raise PipelineError("object is not a Job object.")

    def _validate_sample_sheet(self, sample_sheet_path):
        """
        Performs additional validation for sample-sheet on top of metapool.
        :return: If successful, a valid sample-sheet. Raises descriptive
                 PipelineError() on all failures. Warning messages are
                 appended to self.warnings.
        """
        # validate the sample-sheet using metapool package.
        sheet = load_sample_sheet(sample_sheet_path)

        msgs = sheet.quiet_validate_and_scrub_sample_sheet()

        if any([isinstance(m, ErrorMessage) for m in msgs]):
            # msgs will contain both ErrorMessages and WarningMessages.
            # we want to identify if there are any messages and if so, create
            # a separate list for them. An Error should only be raised on
            # Error messages and in this case, all error messages should be
            # concatenated.
            errors = [x for x in msgs if isinstance(x, ErrorMessage)]

            if errors:
                msgs = [str(x).replace('ErrorMessage: ', '') for x in msgs]
                msgs = 'Sample-sheet contains errors:\n' + '\n'.join(msgs)
                raise PipelineError(msgs)
            else:
                raise PipelineError('Cannot parse sample-sheet.')
        else:
            # perform extended validation based on required fields for
            # seqpro, and other issues encountered.
            bioinformatics = sheet.Bioinformatics
            if 'library_construction_protocol' not in bioinformatics:
                msgs.append(ErrorMessage("column 'library_construction_protoco"
                                         "l' not found in Bioinformatics secti"
                                         "on"))
            if 'experiment_design_description' not in bioinformatics:
                msgs.append(ErrorMessage("column 'experiment_design_descriptio"
                                         "n' not found in Bioinformatics secti"
                                         "on"))

            if sheet.Header['Assay'] not in Pipeline.assay_types:
                msgs.append(ErrorMessage("Valid Assay values are "
                                         f"{Pipeline.assay_types}"))

            # look for duplicate samples. metapool will allow two rows w/the
            # same lane and sample_id if one or more other columns are
            # different. However seqpro expects the tuple (lane, sample_id) to
            # be unique for indexing.
            unique_indexes = []
            for item in sheet.samples:
                unique_index = f'{item.lane}_{item.sample_id}'
                if unique_index in unique_indexes:
                    msgs.append(ErrorMessage("A sample already exists with la"
                                             f"ne {item.lane} and sample-id "
                                             f"{item.sample_id}"))
                else:
                    unique_indexes.append(unique_index)

            errors = [x for x in msgs if isinstance(x, ErrorMessage)]

            if errors:
                msgs = [str(x).replace('ErrorMessage: ', '') for x in msgs]
                msgs = 'Sample-sheet contains errors:\n' + '\n'.join(msgs)
                raise PipelineError(msgs)

            # return a valid sample-sheet, and preserve any warning
            # messages
            self.warnings += [str(x) for x in msgs if
                              isinstance(x, WarningMessage)]
            return sheet

    def _validate_mapping_file(self, mapping_file_path):
        """
        Performs validation for mapping-files.
        :return: If successful, a valid mapping-file. Raises descriptive
                 PipelineError() on all failures. Warning messages are
                 appended to self.warnings.
        """
        try:
            df = pd.read_csv(mapping_file_path, delimiter='\t', dtype=str)
        except pd.errors.ParserError:
            raise PipelineError('Cannot parse mapping-file.')

        # first, detect any duplicate column names, regardless of any mixed-
        # capitalization, and notify the user.
        d = defaultdict(list)
        for column in df.columns:
            d[column.lower()].append(column)

        # generate a list of all unique column names that appear more than
        # once, regardless of capitalization. Then generate a list containing
        # lists of duplicate column names in their original case to report to
        # the user.
        dupes = [d[column] for column in
                 [col for col in d.keys() if len(d[col]) > 1]]

        if dupes:
            # column-names are case-insensitive, and must be unique.
            # return groups of duplicate column names (differentiated only by
            # a different mixed-case) to the user.
            raise PipelineError("Mapping-file contains duplicate columns: "
                                "%s" % ', '.join([str(tpl) for tpl in dupes]))

        # if columns are unique, determine if any columns are missing and/or
        # unexpected and notify the user.
        obs = set(df.columns.str.lower())

        # Note that Pipeline.mapping_file_columns is expected to be all lower-
        # case.

        # if an expected column is missing in observed, that is an error.
        # Note that since a mapping-file is just a DataFrame, there isn't a
        # distinction between a mapping-file that is missing n columns and has
        # n additional columns and a dataframe that is not a mapping-file at
        # all. This method assumes an external test has determined that the
        # file is a mapping-file already.
        missing_columns = Pipeline.mapping_file_columns - obs
        if missing_columns:
            raise PipelineError("Mapping-file is missing columns: "
                                "%s" % ', '.join(sorted(missing_columns)))

        # if an observed column is unexpected, that is a warning.
        unexpected_columns = obs - Pipeline.mapping_file_columns
        if unexpected_columns:
            self.warnings += [("Mapping-file contains additional columns: "
                               "%s" % ', '.join(unexpected_columns))]

        # rename all columns to their lower-case versions.
        # we will want to return this version to the user.
        df.columns = df.columns.str.lower()

        return df

    def generate_sample_info_files(self, addl_info=None):
        """
        Generate sample-information files in self.output_path.
        :param addl_info: A df of (sample-name, project-name) pairs.
        :return: A list of paths to sample-information-files.
        """
        if self.pipeline_type == Pipeline.AMPLICON_PTYPE:
            # Generate a list of BLANKs for each project.
            temp_df = self.mapping_file[[SAMPLE_NAME_KEY, _PROJECT_NAME_KEY]]
            temp_df_as_dicts_list = temp_df.to_dict(orient='records')
            blanks_dicts_list = []
            for record in temp_df_as_dicts_list:
                if is_blank(record[SAMPLE_NAME_KEY]):
                    new_record = record.copy()
                    proj_info = parse_project_name(record[_PROJECT_NAME_KEY])
                    new_record.pop(_PROJECT_NAME_KEY)
                    new_record.update(proj_info)
                    blanks_dicts_list.append(new_record)
                # endif this is a blank
            # next record from mapping file df
            df = pd.DataFrame(blanks_dicts_list)
        else:
            controls = self.sample_sheet.get_denormalized_controls_list()
            df = pd.DataFrame(controls)

        projects = df[PROJECT_FULL_NAME_KEY].unique()

        paths = []
        for project in projects:
            project_info = parse_project_name(project)

            curr_fname = self.make_sif_fname(self.run_id, project)
            curr_fp = join(self.output_path, curr_fname)
            paths.append(curr_fp)

            controls_in_proj_df = \
                df.loc[df[PROJECT_FULL_NAME_KEY] == project].copy()

            # TODO: remove this loop and replace with spp_metadata call at end
            for column, default_value in zip(Pipeline.sif_header,
                                             Pipeline.sif_defaults):
                # ensure all defaults are converted to strings.
                if default_value is not None:
                    controls_in_proj_df[column] = str(default_value)
            # next metadata col/value

            # generate values for the four columns that must be
            # determined from sample-sheet information.
            TEMP_KEY = 'temp_name'
            controls_in_proj_df['title'] = project_info[PROJECT_SHORT_NAME_KEY]
            controls_in_proj_df[TEMP_KEY] = \
                controls_in_proj_df[SAMPLE_NAME_KEY].str.replace("_", ".")
            controls_in_proj_df['host_subject_id'] = \
                controls_in_proj_df[TEMP_KEY]
            controls_in_proj_df['description'] = controls_in_proj_df[TEMP_KEY]
            controls_in_proj_df.drop(columns=[TEMP_KEY], inplace=True)
            controls_in_proj_df['collection_timestamp'] = \
                self.get_date_from_run_id()

            controls_in_proj_df = controls_in_proj_df[Pipeline.sif_header]
            controls_in_proj_df.to_csv(curr_fp, sep='\t', index=False)

            # spp_metadata.write_extended_spp_metadata(
            #     controls_in_proj_df, self.output_path, curr_fname)

        return paths

    def get_date_from_run_id(self):
        # assume all run_ids begin with coded datestamp:
        # 210518_...
        # allow exception if substrings cannot convert to int
        # or if array indexes are out of bounds.
        year = int(self.run_id[0:2]) + 2000
        month = int(self.run_id[2:4])
        day = int(self.run_id[4:6])
        return f'{year}-{month}-{day}'

    def get_sample_ids(self):
        '''
        Returns list of sample-ids sourced from sample-sheet or pre-prep file
        :return: list of sample-ids
        '''

        # test for self.mapping_file, since self.sample_sheet will be
        # defined in both cases.
        if self.pipeline_type == Pipeline.AMPLICON_PTYPE:
            results = list(self.mapping_file.sample_name)
        else:
            results = [x.Sample_ID for x in self.sample_sheet.samples]

        return results

    def get_sample_names(self, project_name=None):
        '''
        Returns list of sample-names sourced from sample-sheet or pre-prep file
        :param project_name: If None, return all sample-names.
        :return: list of sample-names
        '''
        # test for self.mapping_file, since self.sample_sheet will be
        # defined in both cases.
        if self.pipeline_type == Pipeline.AMPLICON_PTYPE:
            return self._get_sample_names_from_mapping_file(project_name)
        else:
            return self._get_sample_names_from_sample_sheet(project_name)

    def _get_sample_names_from_sample_sheet(self, project_name):
        if project_name is None:
            return [x.Sample_Name for x in self.sample_sheet.samples]
        else:
            # Since the project-name is stored in an internal variable
            # in a third-party library, convert the data structure to
            # JSON using the exposed method and obtain from the result.
            jsn = json_loads(self.sample_sheet.to_json())

            results = []

            for sample in jsn['Data']:
                # handle case where project_name includes an appended qiita-id.
                if sample['Sample_Project'] == project_name:
                    results.append(sample['Sample_Name'])
                    continue

                # handle case where project_name does not include a qiita-id.
                # exact matching is required for cases where one project name
                # in a sheet is a superset of another project in the same
                # sheet.
                m = search(r'^(.+)_(\d+)$', sample['Sample_Project'])
                if m[1] == project_name:
                    results.append(sample['Sample_Name'])

            return results

    def get_orig_names_from_sheet(self, project_name):
        if project_name is None:
            results = [x.orig_name for x in self.sample_sheet.samples]
            # eliminate inavoidable duplicates and sort.
            return sorted(set(results))
        else:
            # Since the project-name is stored in an internal variable
            # in a third-party library, convert the data structure to
            # JSON using the exposed method and obtain from the result.
            jsn = json_loads(self.sample_sheet.to_json())

            results = []

            for sample in jsn['Data']:
                # handle case where project_name includes an appended qiita-id.
                if sample['Sample_Project'] == project_name:
                    results.append(sample['orig_name'])
                    continue

                # handle case where project_name does not include a qiita-id.
                # exact matching is required for cases where one project name
                # in a sheet is a superset of another project in the same
                # sheet.
                m = search(r'^(.+)_(\d+)$', sample['Sample_Project'])
                if m[1] == project_name:
                    results.append(sample['orig_name'])

            return sorted(set(results))

    def _get_sample_names_from_mapping_file(self, project_name):
        if project_name is None:
            return list(self.mapping_file.sample_name)
        else:
            df = self.mapping_file[self.mapping_file[_PROJECT_NAME_KEY] ==
                                   project_name]
            return list(df[SAMPLE_NAME_KEY])

    def _parse_project_name(self, project_name, short_names):
        """
        Split fully-qualified project_name into a project_name and a qiita-id
        if possible. Else return project_name and None.
        :param project_name: A fully-qualified project name e.g: Feist_1161.
        :param short_names: True returns orig. value. False returns name only.
        :return: Tuple (project-name, qiita-id)
        """
        # The main functionality of this method has been replaced by this call
        # to metapool's parse_project_name, but I can't guarantee this function
        # isn't used somewhere else in the codebase, so I'm not deleting it.
        proj_info = parse_project_name(project_name)

        if short_names is False:
            # return the fully-qualified project name w/Qiita ID.
            return project_name, proj_info[QIITA_ID_KEY]
        else:
            # return the project's name and qiita_id
            return proj_info[PROJECT_SHORT_NAME_KEY], proj_info[QIITA_ID_KEY]

    def get_project_info(self, short_names=False):
        results = []

        if self.pipeline_type == Pipeline.AMPLICON_PTYPE:
            if CONTAINS_REPLICATES_KEY in self.mapping_file:
                contains_replicates = True
            else:
                contains_replicates = False

            sample_project_map = {pn: _df.sample_name.values for pn, _df in
                                  self.mapping_file.groupby(_PROJECT_NAME_KEY)}
            projects_info = \
                {p: parse_project_name(p) for p in sample_project_map}
        else:
            projects_info = self.sample_sheet.get_projects_details()

        if short_names:
            proj_name_key = PROJECT_SHORT_NAME_KEY
        else:
            proj_name_key = PROJECT_FULL_NAME_KEY

        for curr_project_info in projects_info.values():
            curr_dict = {
                _PROJECT_NAME_KEY: curr_project_info[proj_name_key],
                QIITA_ID_KEY: curr_project_info[QIITA_ID_KEY]
            }

            if self.pipeline_type == Pipeline.AMPLICON_PTYPE:
                # this is a mapping file:
                curr_contains_reps = contains_replicates
            else:
                bi_df = self.sample_sheet.Bioinformatics
                if CONTAINS_REPLICATES_KEY in bi_df.columns.tolist():
                    # subselect rows in [Bioinformatics] based on whether they
                    # match the project name.

                    # whether short_names or full_names are requested in the
                    # results, the match will always need to be against the
                    # full project name, which is what's expected to be in
                    # the Sample_Project column.
                    sample_project = curr_project_info[PROJECT_FULL_NAME_KEY]
                    df = bi_df.loc[bi_df['Sample_Project'] == sample_project]
                    # since only one project can match by definition, convert
                    # to dict and extract the needed value.
                    curr_contains_reps = df.iloc[0].to_dict()[
                                         CONTAINS_REPLICATES_KEY]
                else:
                    curr_contains_reps = False

            curr_dict[CONTAINS_REPLICATES_KEY] = curr_contains_reps
            results.append(curr_dict)
        # next project

        return results

    @staticmethod
    def is_mapping_file(mapping_file_path):
        """
        Returns True if file follows basic mapping-file format.
        """
        try:
            df = pd.read_csv(mapping_file_path, delimiter='\t', dtype=str)
        except pd.errors.ParserError:
            return False

        # if the expected subset of columns required for a mapping-file
        # are present, then consider this a mapping file, even if it's
        # an invalid one.
        exp_columns = frozenset({'barcode', 'tm1000_8_tool',
                                 'extraction_robot', 'pcr_primers'})

        return set(df.columns.str.lower()).issuperset(exp_columns)

    @staticmethod
    def is_sample_sheet(sample_sheet_path):
        '''
        Returns True if file follows basic sample-sheet format.
        '''

        # Check to see if the file begins w/[Header].
        # Ignoring any legacy comments.
        with open(sample_sheet_path, 'r') as f:
            line = f.readline()
            while line:
                if line.startswith('#'):
                    line = f.readline()
                else:
                    break

            if line.startswith('[Header]'):
                return True

        return False

    def _generate_dummy_sample_sheet(self, first_read, last_read,
                                     indexed_reads, dummy_sample_id):
        # create object and initialize header
        sheet = AmpliconSampleSheet()
        sheet.Header['IEMFileVersion'] = '4'
        sheet.Header['Date'] = '10/27/22'
        sheet.Header['Workflow'] = 'GenerateFASTQ'
        sheet.Header['Application'] = 'FASTQ Only'
        sheet.Header['Assay'] = 'TruSeq HT'
        sheet.Header['Description'] = 'test_run'
        sheet.Header['Chemistry'] = 'Amplicon'

        # generate override_cycles string
        tmp = [f"N{x['NumCycles']}" for x in indexed_reads]
        tmp = ';'.join(tmp)
        override_cycles = f"Y{first_read};{tmp};Y{last_read}"

        # set Reads and Settings according to input values
        # we'll get this from the code on the server
        sheet.Reads = [first_read, last_read]
        sheet.Settings['OverrideCycles'] = override_cycles
        sheet.Settings['MaskShortReads'] = '1'
        sheet.Settings['CreateFastqForIndexReads'] = '1'

        dummy_samples = {'Sample_ID': dummy_sample_id,
                         'Sample_Plate': '',
                         'Sample_Well': '',
                         'I7_Index_ID': '',
                         'index': '',
                         'I5_Index_ID': '',
                         'index2': ''
                         }
        sheet.add_sample(sample_sheet.Sample(dummy_samples))

        # contacts won't matter for the dummy sample-sheet.
        contacts = [['c2cowart@ucsd.edu', 'SomeProject'],
                    ['antgonza@gmail.com', 'AnotherProject']]

        # we'll get these from input parameters as well.
        contacts = pd.DataFrame(columns=['Email', 'Sample_Project'],
                                data=contacts)
        sheet.Contact = contacts

        # add a dummy sample.
        samples = [[dummy_sample_id, 'NA', 'NA',
                    'FALSE', 'FALSE', '14782']]

        samples = pd.DataFrame(columns=['Project', 'ForwardAdapter',
                                        'ReverseAdapter', 'PolyGTrimming',
                                        'HumanFiltering', 'QiitaID'],
                               data=samples)

        sheet.Bioinformatics = samples

        return sheet

    def generate_dummy_sample_sheet(self, run_dir, output_fp):
        if exists(run_dir):
            reads = self.process_run_info_file(join(run_dir, 'RunInfo.xml'))
        else:
            raise ValueError("run_dir %s not found." % run_dir)

        # assumptions are first and last reads are non-indexed and there
        # are always two. Between them there is either 1 or 2 indexed
        # reads. If this is not true, raise an Error.

        if len(reads) < 3 or len(reads) > 4:
            # there must be a first and last read w/a minimum of one read
            # in the middle and maximum two in the middle.
            raise ValueError("RunInfo.xml contains abnormal reads.")

        first_read = reads.pop(0)
        last_read = reads.pop()

        if (first_read['IsIndexedRead'] is True or
                last_read['IsIndexedRead'] is True):
            raise ValueError("RunInfo.xml contains abnormal reads.")

        # confirm the interior read(s) are indexed ones.
        for read in reads:
            if read['IsIndexedRead'] is False:
                raise ValueError("RunInfo.xml contains abnormal reads.")

        dummy_sample_id = basename(run_dir) + '_SMPL1'

        sheet = self._generate_dummy_sample_sheet(first_read['NumCycles'],
                                                  last_read['NumCycles'],
                                                  reads, dummy_sample_id)

        with open(output_fp, 'w') as f:
            sheet.write(f, 1)

    def process_run_info_file(self, run_info_fp):
        def process_reads(reads):
            # extract all read elements as a list.
            # the contents of each Read element are highly regular.
            # for now, process w/out installing xml2dict or other
            # library into Qiita env.
            found = findall('<Read (.+?) />', reads)

            results = []
            for item in found:
                attributes = item.split(' ')
                d = {}
                for attribute in attributes:
                    k, v = attribute.split('=')
                    if k in ['NumCycles', 'Number']:
                        v = int(v.strip('"'))
                    elif k in ['IsIndexedRead']:
                        v = v.strip('"')
                        v = False if v == 'N' else True
                    else:
                        raise ValueError("Unknown key: %s" % k)
                    d[k] = v
                results.append(d)

            return results

        with open(run_info_fp, 'r') as f:
            s = f.read()
            reads = search('<Reads>(.+?)</Reads>', s.replace('\n', ''))
            if reads:
                result = reads.group(1)
            else:
                raise ValueError("Cannot extract read information")
            return process_reads(result)
