from .Instruments import Illumina, TellSeq
from os.path import join, abspath, exists, split
from os import walk, makedirs, listdir
import pandas as pd
from json import dumps
from subprocess import Popen, PIPE
from glob import glob
from shutil import copyfile
from sequence_processing_pipeline.Pipeline import Pipeline
from metapool import load_sample_sheet
import logging
from .Assays import Amplicon, Metagenomic
from .Assays import (METAOMIC_ASSAY_NAMES, ASSAY_NAME_AMPLICON,
                     ASSAY_NAME_METAGENOMIC)
from .Instruments import INSTRUMENT_NAME_ILLUMINA, INSTRUMENT_NAME_TELLSEQ


class WorkflowError(Exception):
    def __init__(self, message=None):
        self.message = message
        super().__init__(self.message)


class Workflow():
    def __init__(self, **kwargs):
        """
        base initializer allows WorkflowFactory to return the correct
        Workflow() type w/out having to include all possible starting
        parameters needed for each Workflow type.
        """
        self.kwargs = kwargs
        self.special_map = None
        self.prep_file_paths = None
        self.pipeline = None
        self.qclient = None
        self.sifs = None
        self.cmds_log_path = None
        self.cmds = None
        self.tube_id_map = None
        self.prep_copy_index = 0
        self.samples_in_qiita = None
        self.mandatory_attributes = []
        self.skip_steps = []

        if 'status_update_callback' in kwargs:
            self.status_update_callback = kwargs['status_update_callback']

    def confirm_mandatory_attributes(self):
        absent_list = []

        for attribute in self.mandatory_attributes:
            if attribute not in self.kwargs:
                absent_list.append(attribute)

        if absent_list:
            raise ValueError("The following values are not defined in kwargs:"
                             + " " + ', '.join(absent_list))

    def update_status(self, msg):
        if self.status_update_callback:
            self.status_update_callback(msg)

    def what_am_i(self):
        """
        Returns text description of Workflow's Instrument & Assay mixins.
        :return:
        """
        return (f"Instrument: {self.instrument_name()}" + " " +
                f"Assay: {self.assay_name()}")

    def generate_special_map(self):
        """
        Generates a list of tuples to support pipeline processing.
        :return: A list of triplets.
        """
        # this function should be able to be tested by passing in simulated =
        # results from qclient.

        # trimmed files are stored by qiita_id. Find the qiita_id
        # associated with each project and ensure a subdirectory exists
        # for when it comes time to move the trimmed files.

        special_map = []
        results = self.qclient.get("/qiita_db/artifacts/types/")
        projects = self.pipeline.get_project_info()
        for project in projects:
            upload_path = join(results['uploads'], project['qiita_id'])
            makedirs(upload_path, exist_ok=True)
            special_map.append((project['project_name'],
                                upload_path,
                                project['qiita_id']))

        self.special_map = special_map

    def generate_sifs(self):
        """
        TODO
        :return:
        """
        from_qiita = {}

        for study_id in self.prep_file_paths:
            url = f'/api/v1/study/{study_id}/samples'
            logging.debug(url)
            samples = list(self.qclient.get(url))
            from_qiita[study_id] = samples

        add_sif_info = []

        qid_pn_map = {proj['qiita_id']: proj['project_name'] for
                      proj in self.pipeline.get_project_info()}

        # in case we really do need to query for samples again:
        # assume set of valid study_ids can be determined from prep_file_paths.
        for study_id in from_qiita:
            samples = from_qiita[study_id]
            # generate a list of (sample-name, project-name) pairs.
            project_name = qid_pn_map[study_id]
            samples = [(x, project_name) for x in samples]
            add_sif_info.append(pd.DataFrame(data=samples,
                                             columns=['sample_name',
                                                      'project_name']))

        # convert the list of dataframes into a single dataframe.
        add_sif_info = pd.concat(add_sif_info,
                                 ignore_index=True).drop_duplicates()

        # generate SIF files with add_sif_info as additional metadata input.
        # duplicate sample-names and non-blanks will be handled properly.
        self.sifs = self.pipeline.generate_sample_info_files(add_sif_info)

        return self.sifs

    def update_blanks_in_qiita(self):
        """
        TODO
        :return:
        """
        for sif_path in self.sifs:
            # get study_id from sif_file_name ...something_14385_blanks.tsv
            study_id = sif_path.split('_')[-2]

            df = pd.read_csv(sif_path, delimiter='\t', dtype=str)

            # Prepend study_id to make them compatible w/list from Qiita.
            df['sample_name'] = f'{study_id}.' + df['sample_name'].astype(str)

            # SIFs only contain BLANKs. Get the list of potentially new BLANKs.
            blank_ids = [i for i in df['sample_name'] if 'blank' in i.lower()]
            blanks = df[df['sample_name'].isin(blank_ids)]['sample_name']
            if len(blanks) == 0:
                # we have nothing to do so let's return early
                return

            # Get list of BLANKs already registered in Qiita.
            from_qiita = self.qclient.get(f'/api/v1/study/{study_id}/samples')
            from_qiita = [x for x in from_qiita if
                          x.startswith(f'{study_id}.BLANK')]

            # Generate list of BLANKs that need to be ADDED to Qiita.
            new_blanks = (set(blanks) | set(from_qiita)) - set(from_qiita)

            if len(new_blanks):
                # Generate dummy entries for each new BLANK, if any.
                url = f'/api/v1/study/{study_id}/samples/info'
                logging.debug(url)
                categories = self.qclient.get(url)['categories']

                # initialize payload w/required dummy categories
                data = {i: {c: 'control sample' for c in categories} for i in
                        new_blanks}

                # populate payload w/additional columns and/or overwrite
                # existing columns w/metadata from SIF file.
                sif_data = df.set_index('sample_name').T.to_dict()
                for new_blank in new_blanks:
                    for column in sif_data[new_blank]:
                        data[new_blank][column] = sif_data[new_blank][column]

                # http_patch will raise Error if insert failed.
                self.qclient.http_patch(f'/api/v1/study/{study_id}/samples',
                                        data=dumps(data))

    def _helper_process_operations(self):
        """
        Helper method for generate_commands()
        :return:
        """
        RESULTS_DIR = 'final_results'
        TAR_CMD = 'tar zcvf'
        LOG_PREFIX = 'logs'
        REPORT_PREFIX = 'reports'
        PREP_PREFIX = 'prep-files'
        CONVERT_JOB = 'ConvertJob'
        QC_JOB = 'NuQCJob'
        FASTQC_JOB = 'FastQCJob'
        PREPFILE_JOB = 'GenPrepFileJob'
        TAR_EXT = 'tgz'

        op_meta = [(['ConvertJob/logs'], TAR_CMD,
                    f'{LOG_PREFIX}-{CONVERT_JOB}.{TAR_EXT}', 'OUTPUT_FIRST'),

                   (['ConvertJob/Reports', 'ConvertJob/logs'], TAR_CMD,
                    f'{REPORT_PREFIX}-{CONVERT_JOB}.{TAR_EXT}',
                    'OUTPUT_FIRST'),

                   (['NuQCJob/logs'], TAR_CMD,
                    f'{LOG_PREFIX}-{QC_JOB}.{TAR_EXT}', 'OUTPUT_FIRST'),

                   (['FastQCJob/logs'], TAR_CMD,
                    f'{LOG_PREFIX}-{FASTQC_JOB}.{TAR_EXT}', 'OUTPUT_FIRST'),

                   (['FastQCJob/fastqc'], TAR_CMD,
                    f'{REPORT_PREFIX}-{FASTQC_JOB}.{TAR_EXT}', 'OUTPUT_FIRST'),

                   (['GenPrepFileJob/logs'], TAR_CMD,
                    f'{LOG_PREFIX}-{PREPFILE_JOB}.{TAR_EXT}', 'OUTPUT_FIRST'),

                   (['GenPrepFileJob/PrepFiles'], TAR_CMD,
                    f'{PREP_PREFIX}.{TAR_EXT}', 'OUTPUT_FIRST'),

                   (['failed_samples.html', 'touched_studies.html'],
                    'mv', RESULTS_DIR, 'INPUTS_FIRST'),

                   (['FastQCJob/multiqc'], 'mv', RESULTS_DIR, 'INPUTS_FIRST')]

        cmds = []

        for inputs, action, output, order in op_meta:
            confirmed_inputs = []
            for input in inputs:
                if exists(join(self.pipeline.output_path, input)):
                    # it's expected that some inputs may not exist due to
                    # different pipeline types. If one or more inputs do not
                    # exist, do not include them in the command-line as they
                    # may cause an error.
                    confirmed_inputs.append(input)

            # do not add the command to the list unless at least one of
            # the inputs exists. It's okay for a command to go unprocessed.
            if confirmed_inputs:
                # convert to string form before using.
                confirmed_inputs = ' '.join(confirmed_inputs)
                if order == 'OUTPUT_FIRST':
                    cmds.append(f'{action} {output} {confirmed_inputs}')
                elif order == 'INPUTS_FIRST':
                    cmds.append(f'{action} {confirmed_inputs} {output}')
                else:
                    raise ValueError(f"'{order}' is not a defined order of "
                                     "operations")

        return cmds

    def _process_blanks(self):
        """
        Helper method for generate_commands().

        :return:
        """
        results = [x for x in listdir(self.pipeline.output_path) if
                   x.endswith('_blanks.tsv')]

        results.sort()

        if len(results) > 0:
            return 'tar zcvf sample-files.tgz' + ' ' + ' '.join(results)

    def _process_fastp_report_dirs(self):
        """
        Helper method for generate_commands().
        :return:
        """
        report_dirs = []

        for root, dirs, files in walk(self.pipeline.output_path):
            for dir_name in dirs:
                if dir_name == 'fastp_reports_dir':
                    # generate the full path for this directory before
                    # truncating everything up to the NuQCJob directory.
                    full_path = join(root, dir_name).split('NuQCJob/')
                    report_dirs.append(join('NuQCJob', full_path[1]))

        if report_dirs:
            report_dirs.sort()
            return 'tar zcvf reports-NuQCJob.tgz ' + ' '.join(report_dirs)
        else:
            # It is okay to return an empty list of commands if reports_dirs
            # is empty. Some pipelines do not generate fastp reports.
            return []

    def _write_commands_to_output_path(self):
        """
        Helper method for generate_commands().
        :return:
        """
        self.cmds_log_path = join(self.pipeline.output_path, 'cmds.log')
        with open(self.cmds_log_path, 'w') as f:
            for cmd in self.cmds:
                f.write(f'{cmd}\n')

    def generate_commands(self):
        """
        TODO
        :return:
        """
        cmds = self._helper_process_operations()

        result = self._process_fastp_report_dirs()

        if result:
            cmds.append(result)

        result = self._process_blanks()

        if result:
            cmds.append(result)

        # if one or more tar-gzip files are found (which we expect there to
        # be), move them into the 'final_results' directory. However, if none
        # are present, don't raise an error.
        cmds.append('(find *.tgz -maxdepth 1 -type f | xargs mv -t '
                    'final_results) || true')

        # prepend each command with a change-directory to the correct
        # location.
        cmds = [f'cd {self.pipeline.output_path}; {cmd}' for cmd in cmds]

        self.cmds = cmds

        self._write_commands_to_output_path()

    def execute_commands(self):
        """
        TODO
        :return:
        """
        # execute the list of commands in order
        for cmd in self.cmds:
            p = Popen(cmd,
                      universal_newlines=True,
                      shell=True,
                      stdout=PIPE,
                      stderr=PIPE)
            _, _ = p.communicate()
            return_code = p.returncode

            if return_code != 0:
                # during testing, ignore processes that fail and continue
                # to test other commands.
                raise WorkflowError(f"'{cmd}' returned {return_code}")

    def _project_metadata_check(self):
        """
        Helper method for pre_check()
        :return:
        """
        # Let Pipeline() retrieve the needed qiita study ids from the user
        # input while this plugin queries for the existing set of column
        # names in each project's sample metadata. We'll let Pipeline()
        # decide (using its metapool dependency) which column names are
        # reserved.
        qiita_ids = [x['qiita_id'] for x in self.pipeline.get_project_info()]

        results = []

        for qiita_id in qiita_ids:
            url = f"/api/v1/study/{qiita_id}/samples/info"
            # CHARLIE
            logging.debug(f"URL: {url}")
            categories = self.qclient.get(url)["categories"]

            res = self.pipeline.identify_reserved_words(categories)

            # if any reserved words were identified, generate an appropriate
            # error message for it and add it to the list of error messages
            # to return to the user.
            res = [f"'{x}' exists in Qiita study {qiita_id}'s sample metadata"
                   for x in res]

            results += res

        if results:
            # return any error messages generated across all the projects.
            raise WorkflowError("\n".join(results))

    def _process_tube_ids(self, qiita_id, samples):
        """
        Helper method for _compare_samples_against_qiita().
        :param qiita_id:
        :param samples:
        :return:
        """
        if qiita_id in self.tube_id_map:
            tids = [self.tube_id_map[qiita_id][sample] for sample in
                    self.tube_id_map[qiita_id]]

            not_in_qiita = samples - set(tids)

            if not_in_qiita:
                # strip any leading zeroes from the sample-ids. Note that
                # if a sample-id has more than one leading zero, all of
                # them will be removed.
                not_in_qiita = set([x.lstrip('0') for x in samples]) - \
                               set(tids)

            # convert examples to strings before returning
            examples = [str(example) for example in tids[:5]]

            number_in_project = len(set(tids))

            return not_in_qiita, examples, number_in_project

        # return None otherwise

    def _compare_samples_against_qiita(self):
        """
        Helper method for pre_check().
        :return:
        """
        projects = self.pipeline.get_project_info(short_names=True)

        results = []
        for project in projects:
            msgs = []
            self._get_tube_ids_from_qiita(self.qclient)
            p_name = project['project_name']
            qiita_id = str(project['qiita_id'])
            contains_replicates = project['contains_replicates']

            # get list of samples as presented by the sample-sheet or mapping
            # file and confirm that they are all registered in Qiita.
            if contains_replicates:
                # don't match against sample-names with a trailing well-id
                # if project contains replicates.
                msgs.append("This sample-sheet contains replicates. sample-"
                            "names will be sourced from orig_name column.")
                samples = set(self.pipeline.get_orig_names_from_sheet(p_name))
            else:
                samples = set(self.pipeline.get_sample_names(p_name))

            # do not include BLANKs. If they are unregistered, we will add
            # them downstream.
            samples = {smpl for smpl in samples
                       if not smpl.startswith('BLANK')}

            msgs.append(f"The total number of samples found in {p_name} that "
                        f"aren't BLANK is: {len(samples)}")

            results_sn = self._process_sample_names(p_name, qiita_id,
                                                    samples)

            msgs.append("Number of values in sheet that aren't sample-names in"
                        " Qiita: %s" % len(results_sn[0]))

            use_tids = False

            if len(results_sn[0]) == 0:
                msgs.append(f"All values in sheet matched sample-names "
                            f"registered with {p_name}")
            else:
                # not all values were matched to sample-names.
                # check for possible match w/tube-ids, if defined in project.
                results_tid = self._process_tube_ids(p_name, qiita_id,
                                                     samples)
                if results_tid:
                    msgs.append("Number of values in sheet that aren't "
                                "tube-ids in Qiita: %s" % len(results_tid[0]))

                    if len(results_tid[0]) == 0:
                        # all values were matched to tube-ids.
                        use_tids = True
                        msgs.append(f"All values in sheet matched tube-ids "
                                    f"registered with {p_name}")
                    else:
                        # we have sample-names and tube-ids and neither is
                        # a perfect match.
                        if len(results_tid[0]) < len(results_sn[0]):
                            # more tube-ids matched than sample-names.
                            use_tids = True
                            msgs.append(f"More values in sheet matched tube-"
                                        f"ids than sample-names with {p_name}")
                        elif len(results_tid[0]) == len(results_sn[0]):
                            msgs.append("Sample-names and tube-ids were "
                                        "equally non-represented in the "
                                        "sample-sheet")
                        else:
                            msgs.append(f"More values in sheet matched sample-"
                                        f"names than tube-ids with {p_name}")
                else:
                    msgs.append("there are no tube-ids registered with "
                                f"{p_name}")

            if use_tids:
                not_in_qiita = results_tid[0]
                examples = results_tid[1]
                total_in_qiita = results_tid[2]
            else:
                not_in_qiita = results_sn[0]
                examples = results_sn[1]
                total_in_qiita = results_sn[2]

            # return an entry for all projects, even when samples_not_in_qiita
            # is an empty list, as the information is still valuable.
            results.append({'samples_not_in_qiita': not_in_qiita,
                            'examples_in_qiita': examples,
                            'project_name': p_name,
                            'total_in_qiita': total_in_qiita,
                            'used_tids': use_tids,
                            'messages': msgs})

        return results

    @classmethod
    def update_sample_sheet(cls, sample_sheet_path, lane_number):
        # use KLSampleSheet functionality to add/overwrite lane number.
        sheet = load_sample_sheet(sample_sheet_path)
        for sample in sheet:
            sample['Lane'] = f'{lane_number}'

        with open(sample_sheet_path, 'w') as f:
            sheet.write(f)

    @classmethod
    def get_samples_in_qiita(cls, qclient, qiita_id):
        '''
        Obtain lists for sample-names and tube-ids registered in Qiita.
        :param qclient: QiitaClient object
        :param qiita_id: Qiita ID for the project in question.
        :return: a tuple of lists, one for sample-names, another for tube-ids.
        '''
        samples = qclient.get(f'/api/v1/study/{qiita_id}/samples')

        # remove Qiita ID as a prefix from the sample-names.
        samples = {x.replace(f'{qiita_id}.', '') for x in samples}

        # find out if tube-ids are registered in the study.
        categories = qclient.get(f'/api/v1/study/{qiita_id}'
                                 '/samples/info')['categories']

        if 'tube_id' in categories:
            tids = qclient.get(f'/api/v1/study/{qiita_id}/samples/'
                               'categories=tube_id')['samples']
        else:
            tids = None

        return (samples, tids)

    def _get_postqc_fastq_files(self, out_dir, project):
        af = None
        sub_folders = ['amplicon', 'filtered_sequences', 'trimmed_sequences']
        for sub_folder in sub_folders:
            sf = join(out_dir, 'NuQCJob', project, sub_folder)
            if exists(sf):
                af = [f for f in glob(join(sf, '*.fastq.gz'))]
                break
        if af is None or not af:
            raise WorkflowError("NuQCJob output not in expected location")

        files = {'raw_barcodes': [], 'raw_forward_seqs': [],
                 'raw_reverse_seqs': []}

        for fastq_file in af:
            if '_I1_' in fastq_file or '_I2_' in fastq_file:
                files['raw_barcodes'].append(fastq_file)
            elif '_R1_' in fastq_file:
                files['raw_forward_seqs'].append(fastq_file)
            elif '_R2_' in fastq_file:
                files['raw_reverse_seqs'].append(fastq_file)
            else:
                raise ValueError(f"Unrecognized file: {fastq_file}")

        files['raw_barcodes'].sort()
        files['raw_forward_seqs'].sort()
        files['raw_reverse_seqs'].sort()

        # Amplicon runs should contain raw_barcodes/I1 files.
        # Meta*omics files doesn't use them.
        if self.pipeline.pipeline_type != ASSAY_NAME_AMPLICON:
            del (files['raw_barcodes'])

        # confirm expected lists of reads are not empty.
        for f_type in files:
            if not files[f_type]:
                # if one or more of the expected list of reads is empty,
                # raise an Error.
                raise ValueError(f"'{f_type}' is empty")

        return files

    def _load_prep_into_qiita(self, qclient, prep_id, artifact_name,
                              qiita_id, project, fastq_files, atype):
        surl = f'{qclient._server_url}/study/description/{qiita_id}'
        prep_url = (f'{qclient._server_url}/study/description/'
                    f'{qiita_id}?prep_id={prep_id}')

        # ideally we would use the email of the user that started the SPP
        # run but at this point there is no easy way to retrieve it
        pdata = {'user_email': 'qiita.help@gmail.com',
                 'prep_id': prep_id,
                 'artifact_type': atype,
                 'command_artifact_name': artifact_name,
                 'add_default_workflow': True,
                 'files': dumps(fastq_files)}

        job_id = qclient.post('/qiita_db/artifact/', data=pdata)['job_id']

        return {'Project': project, 'Qiita Study ID': qiita_id,
                'Qiita Prep ID': prep_id, 'Qiita URL': surl,
                'Artifact Name': artifact_name,
                'Prep URL': prep_url, 'Linking JobID': job_id}

    def _copy_files(self, files):
        # increment the prep_copy_index before generating a new set of copies.
        self.prep_copy_index += 1
        new_files = {}
        for key in files:
            new_files[key] = []
            for some_path in files[key]:
                path_name, file_name = split(some_path)
                path_name = join(path_name, f'copy{self.prep_copy_index}')
                makedirs(path_name, exist_ok=True)
                new_files[key].append(join(path_name, file_name))

        for key in files:
            for src, dst in zip(files[key], new_files[key]):
                copyfile(src, dst)
        return new_files

    def get_prep_file_paths(self):
        return self.prep_file_paths

    def _get_tube_ids_from_qiita(self, qclient):
        # Update get_project_info() so that it can return a list of
        # samples in projects['samples']. Include blanks in projects['blanks']
        # just in case there are duplicate qiita_ids
        qiita_ids = [proj['qiita_id'] for proj in
                     self.pipeline.get_project_info(short_names=True)]

        tids_by_qiita_id = {}
        sample_names_by_qiita_id = {}

        for qiita_id in qiita_ids:
            # Qiita returns a set of sample-ids in qsam and a dictionary where
            # sample-names are used as keys and tube-ids are their values.
            qsam, tids = self.get_samples_in_qiita(qclient, qiita_id)

            sample_names_by_qiita_id[str(qiita_id)] = qsam

            if tids is not None:
                # fix values in tids to be a string instead of a list of one.
                # also, remove the qiita_id prepending each sample-name.
                tids = {k.replace(f'{qiita_id}.', ''): tids[k][0] for k in
                        tids}

                # the values Qiita returns for tids seems like it can include
                # empty strings if there is no tube-id associated with a
                # sample-name. For now assume it doesn't happen in production
                # and if prep-files have empty sample-names we'll know.
                tids_by_qiita_id[str(qiita_id)] = tids

        # use empty dict {} as an indication that get_tube_ids_from_qiita was
        # called but no tube-ids were found for any project.
        # to clarify, self.tube_id_map maps sample-names to tube-ids.
        self.tube_id_map = tids_by_qiita_id
        # should samples_in_qiita be none if tube_id_map is not?
        self.samples_in_qiita = sample_names_by_qiita_id

    def _process_sample_names(self, project_name, qiita_id, samples):
        not_in_qiita = samples - set(self.samples_in_qiita[qiita_id])
        examples = list(samples)[:5]

        # convert to strings before returning
        examples = [str(example) for example in examples]

        number_in_project = len(set(self.samples_in_qiita[qiita_id]))

        return not_in_qiita, examples, number_in_project

    def pre_check(self):
        # TODO: Note that amplicon/illumina workflows shouldn't call this
        #  method since demuxing samples is performed downstream, hence there
        #  isn't a list of real sample-names in the sample-sheet to check.
        #  This check could be performed using the pre-prep sample data
        #  however.

        # since one of the objectives of SPP is to generate prep-info files
        # and automatically load them into Qiita, confirm that all studies
        # mentioned in the sample-sheet/pre-prep do not contain sample
        # metadata that would cause an error in the pipeline after processing
        # has already completed but the results have not yet been loaded.
        self._project_metadata_check()

        # compare sample-ids/tube-ids in sample-sheet/mapping file
        # against what's in Qiita. Results are a list of dictionaries, one
        # per project.
        results = self._compare_samples_against_qiita()

        # obtain a list of non-zero counts of samples missing in Qiita, one
        # for each project. The names of the projects are unimportant. We
        # want to abort early if any project in the sample-sheet/pre-prep file
        # contains samples that aren't registered in Qiita.
        tmp = [len(project['samples_not_in_qiita']) for project in results]
        missing_counts = [count for count in tmp if count != 0]

        if missing_counts:
            msgs = []
            for result in results:
                msgs += result['messages']

            if msgs:
                raise WorkflowError('\n'.join(msgs))

    def scan_run_directory(self):
        # TODO: Impl
        pass


class StandardMetagenomicWorkflow(Workflow, Metagenomic, Illumina):
    def __init__(self, **kwargs):
        super().__init__(kwargs)

        self.mandatory_attributes = ['qclient', 'input_file_path',
                                     'lane_number', 'config_fp',
                                     'run_identifier', 'output_dir', 'job_id',
                                     'lane_number', 'is_restart']

        self.confirm_mandatory_attributes(kwargs)

        # second stage initializer that could conceivably be pushed down into
        # specific children requiring specific parameters.
        self.qclient = self.kwargs['qclient']

        self.pipeline = Pipeline(self.kwargs['config_fp'],
                                 self.kwargs['run_identifier'],
                                 self.kwargs['input_file_path'],
                                 self.kwargs['output_dir'],
                                 self.kwargs['job_id'],
                                 ASSAY_NAME_METAGENOMIC,
                                 lane_number=self.kwargs['lane_number'])

        self.master_qiita_job_id = None

        self.lane_number = self.kwargs['lane_number']
        self.is_restart = bool(self.kwargs['is_restart'])

        if self.is_restart is True:
            self.determine_steps_to_skip()

        self.update = True

        if 'update_qiita' in kwargs:
            if bool(kwargs['update_qiita']) is False:
                self.update = False

        """
        member variables from deprecated Step class. it may be good to
        add one or more of these here.

        self.lane_number = lane_number
        # self.generated_artifact_name =
        f'{self.pipeline.run_id}_{self.lane_number}'
        self.master_qiita_job_id = master_qiita_job_id

        self.is_restart = is_restart

        if status_update_callback is not None:
        self.update_callback = status_update_callback.update_job_status
        else:
        self.update_callback = None

        # for now, hardcode this at the legacy value, since we've never
        # changed it.
        self.job_pool_size = 30

        # initialize other member variables so that they're always present,
        # even when the step that populates them hasn't been run yet.
        self.project_names = None
        self.cmds = None
        self.cmds_log_path = None
        # set by child classes for use in parent class
        self.prep_file_paths = None
        # set by child classes for use in parent class
        self.has_replicates = None
        self.sifs = None
        self.tube_id_map = None
        self.samples_in_qiita = None
        self.output_path = None
        self.sample_state = None
        self.touched_studies_prep_info = None
        self.run_prefixes = {}
        self.prep_copy_index = 0
        self.use_tellread = False
        """

    def determine_steps_to_skip(self):
        out_dir = self.pipeline.output_path

        # figure out what actually needs to be skipped if restarting:
        if exists(join(out_dir, 'NuQCJob')):
            self.skip_steps.append('ConvertJob')

        if exists(join(out_dir, 'FastQCJob')):
            self.skip_steps.append('NuQCJob')

        if exists(join(out_dir, 'GenPrepFileJob')):
            self.skip_steps.append('FastQCJob')

        # it doesn't matter if cmds.log is a valid cmds.log or just
        # an empty file. The cmds.log will get overwritten downstream.
        if exists(join(out_dir, 'cmds.log')):
            self.skip_steps.append('GenPrepFileJob')

    def execute_pipeline(self):
        '''
        Executes steps of pipeline in proper sequence.
        :return: None
        '''
        if not self.is_restart:
            self.pre_check()

        # this is performed even in the event of a restart.
        self.generate_special_map()

        # even if a job is being skipped, it's being skipped because it was
        # determined that it already completed successfully. Hence,
        # increment the status because we are still iterating through them.

        self.update_status("Converting data", 1, 9)
        if "ConvertJob" not in self.skip_steps:
            # converting raw data to fastq depends heavily on the instrument
            # used to generate the run_directory. Hence this method is
            # supplied by the instrument mixin.
            results = self.convert_raw_to_fastq()
            self.fsr_write(results, 'ConvertJob')

        self.update_status("Performing quality control", 2, 9)
        if "NuQCJob" not in self.skip_steps:
            # amplicon runs do not currently perform qc as the demuxing of
            # samples is performed downstream of SPP. It also does not depend
            # on the instrument type since fastq files are by convention the
            # output in either case.
            #
            # Hence, quality control is associated w/the assay mixin (for now).
            self.quality_control(self.pipeline)

        self.update_status("Generating reports", 3, 9)
        if "FastQCJob" not in self.skip_steps:
            # reports are currently implemented by the assay mixin. This is
            # only because metagenomic runs currently require a failed-samples
            # report to be generated. This is not done for amplicon runs since
            # demultiplexing occurs downstream of SPP.
            self.generate_reports()

        self.update_status("Generating preps", 4, 9)
        if "GenPrepFileJob" not in self.skip_steps:
            # preps are currently associated with array mixin, but only
            # because there are currently some slight differences in how
            # FastQCJob gets instantiated(). This could get moved into a
            # shared method, but probably still in Assay.
            self.generate_prep_file()

        # moved final component of genprepfilejob outside of object.
        # obtain the paths to the prep-files generated by GenPrepFileJob
        # w/out having to recover full state.
        tmp = join(self.pipeline.output_path, 'GenPrepFileJob', 'PrepFiles')

        self.has_replicates = False

        prep_paths = []
        self.prep_file_paths = {}

        for root, dirs, files in walk(tmp):
            for _file in files:
                # breakup the prep-info-file into segments
                # (run-id, project_qid, other) and cleave
                # the qiita-id from the project_name.
                qid = _file.split('.')[1].split('_')[-1]

                if qid not in self.prep_file_paths:
                    self.prep_file_paths[qid] = []

                _path = abspath(join(root, _file))
                if _path.endswith('.tsv'):
                    prep_paths.append(_path)
                    self.prep_file_paths[qid].append(_path)

            for _dir in dirs:
                if _dir == '1':
                    # if PrepFiles contains the '1' directory, then it's a
                    # given that this sample-sheet contains replicates.
                    self.has_replicates = True

        # currently imported from Assay although it is a base method. it
        # could be imported into Workflows potentially, since it is a post-
        # processing step. All pairings of assay and instrument type need to
        # generate prep-info files in the same format.
        self.overwrite_prep_files(prep_paths)

        # for now, simply re-run any line below as if it was a new job, even
        # for a restart. functionality is idempotent, except for the
        # registration of new preps in Qiita. These will simply be removed
        # manually.

        # post-processing steps are by default associated with the Workflow
        # class, since they deal with fastq files and Qiita, and don't depend
        # on assay or instrument type.
        self.update_status("Generating sample information", 5, 9)
        self.sifs = self.generate_sifs()

        # post-processing step.
        self.update_status("Registering blanks in Qiita", 6, 9)
        if self.update:
            self.update_blanks_in_qiita()

        self.update_status("Loading preps into Qiita", 7, 9)
        if self.update:
            self.update_prep_templates()

        # before we load preps into Qiita we need to copy the fastq
        # files n times for n preps and correct the file-paths each
        # prep is pointing to.
        self.load_preps_into_qiita()

        self.update_status("Generating packaging commands", 8, 9)
        self.generate_commands()

        self.update_status("Packaging results", 9, 9)
        if self.update:
            self.execute_commands()


class StandardAmpliconWorkflow(Workflow, Amplicon, Illumina):
    @classmethod
    def generate_pipeline(cls, pipeline_type, input_file_path, lane_number,
                          config_fp,
                          run_identifier, out_dir, job_id):
        # TODO METAOMIC_ASSAY_NAMES could still fall out of sync with what;s
        # in mg-scripts.
        if pipeline_type in METAOMIC_ASSAY_NAMES:
            cls.update_sample_sheet(input_file_path, lane_number)
            return Pipeline(config_fp, run_identifier, input_file_path, None,
                            out_dir, job_id, pipeline_type)
        elif pipeline_type == ASSAY_NAME_AMPLICON:
            return Pipeline(config_fp, run_identifier, None, input_file_path,
                            out_dir, job_id, pipeline_type)
        else:
            raise WorkflowError(
                f"'{pipeline_type}' is not a valid Pipeline type.")


class TellSeqMetagenomicWorkflow(Workflow, Metagenomic, TellSeq):
    def __init__(self):
        pass

    def execute_pipeline(self, update_status, update=True, skip_steps=[]):
        '''
        Executes steps of pipeline in proper sequence.
        :param qclient: Qiita client library or equivalent.
        :param update_status: callback function to update status.
        :param update: Set False to prevent updates to Qiita.
        :return: None
        '''

        if not self.is_restart:
            self.pre_check()

        # this is performed even in the event of a restart.
        self.generate_special_map()

        # even if a job is being skipped, it's being skipped because it was
        # determined that it already completed successfully. Hence,
        # increment the status because we are still iterating through them.

        update_status("Converting data", 1, 9)
        if "ConvertJob" not in skip_steps:
            # converting raw data to fastq depends heavily on the instrument
            # used to generate the run_directory. Hence this method is
            # supplied by the instrument mixin.
            results = self.convert_raw_to_fastq()
            self.fsr_write(results, 'ConvertJob')

        update_status("Performing quality control", 2, 9)
        if "NuQCJob" not in skip_steps:
            # amplicon runs do not currently perform qc as the demuxing of
            # samples is performed downstream of SPP. It also does not depend
            # on the instrument type since fastq files are by convention the
            # output in either case.
            #
            # Hence, quality control is associated w/the assay mixin (for now).
            self.quality_control(self.pipeline)

        update_status("Generating reports", 3, 9)
        if "FastQCJob" not in skip_steps:
            # reports are currently implemented by the assay mixin. This is
            # only because metagenomic runs currently require a failed-samples
            # report to be generated. This is not done for amplicon runs since
            # demultiplexing occurs downstream of SPP.
            self.generate_reports()

        update_status("Generating preps", 4, 9)
        if "GenPrepFileJob" not in skip_steps:
            # preps are currently associated with array mixin, but only
            # because there are currently some slight differences in how
            # FastQCJob gets instantiated(). This could get moved into a
            # shared method, but probably still in Assay.
            self.generate_prep_file()

        # moved final component of genprepfilejob outside of object.
        # obtain the paths to the prep-files generated by GenPrepFileJob
        # w/out having to recover full state.
        tmp = join(self.pipeline.output_path, 'GenPrepFileJob', 'PrepFiles')

        self.has_replicates = False

        prep_paths = []
        self.prep_file_paths = {}

        for root, dirs, files in walk(tmp):
            for _file in files:
                # breakup the prep-info-file into segments
                # (run-id, project_qid, other) and cleave
                # the qiita-id from the project_name.
                qid = _file.split('.')[1].split('_')[-1]

                if qid not in self.prep_file_paths:
                    self.prep_file_paths[qid] = []

                _path = abspath(join(root, _file))
                if _path.endswith('.tsv'):
                    prep_paths.append(_path)
                    self.prep_file_paths[qid].append(_path)

            for _dir in dirs:
                if _dir == '1':
                    # if PrepFiles contains the '1' directory, then it's a
                    # given that this sample-sheet contains replicates.
                    self.has_replicates = True

        # currently imported from Assay although it is a base method. it
        # could be imported into Workflows potentially, since it is a post-
        # processing step. All pairings of assay and instrument type need to
        # generate prep-info files in the same format.
        self.overwrite_prep_files(prep_paths)

        # for now, simply re-run any line below as if it was a new job, even
        # for a restart. functionality is idempotent, except for the
        # registration of new preps in Qiita. These will simply be removed
        # manually.

        # post-processing steps are by default associated with the Workflow
        # class, since they deal with fastq files and Qiita, and don't depend
        # on assay or instrument type.
        update_status("Generating sample information", 5, 9)
        self.sifs = self.generate_sifs()

        # post-processing step.
        update_status("Registering blanks in Qiita", 6, 9)
        if update:
            self.update_blanks_in_qiita()

        update_status("Loading preps into Qiita", 7, 9)
        if update:
            self.update_prep_templates()

        # before we load preps into Qiita we need to copy the fastq
        # files n times for n preps and correct the file-paths each
        # prep is pointing to.
        self.load_preps_into_qiita()

        update_status("Generating packaging commands", 8, 9)
        self.generate_commands()

        update_status("Packaging results", 9, 9)
        if update:
            self.execute_commands()

        # DONE! :)


class WorkflowFactory():
    WORKFLOWS = [StandardMetagenomicWorkflow, TellSeqMetagenomicWorkflow]

    ST_TO_IN_MAP = {INSTRUMENT_NAME_ILLUMINA: ['standard_metag',
                                               'standard_metat',
                                               'absquant_metag',
                                               'absquant_metat'],
                    INSTRUMENT_NAME_TELLSEQ: ['tellseq_metag',
                                              'tellseq_absquant']}

    @classmethod
    def _get_instrument_type(cls, sheet):
        for instrument_type in cls.ST_TO_IN_MAP:
            if sheet.Header['SheetType'] in cls.ST_TO_IN_MAP[instrument_type]:
                return instrument_type

    @classmethod
    def generate_workflow(cls, **kwargs):
        if 'uif_path' not in kwargs:
            raise ValueError("The following values are not defined in "
                             "kwargs: 'uif_path'")

        # determine assay-type & instrument-type

        if Pipeline.is_sample_sheet(kwargs['uif_path']):
            sheet = load_sample_sheet(kwargs['uif_path'])
            assay_type = sheet.Header['Assay']
            if assay_type not in METAOMIC_ASSAY_NAMES:
                raise WorkflowError("Can't determine workflow from assay "
                                    "type: %s" % assay_type)
            instrument_type = cls._get_instrument_type(sheet)
        elif Pipeline.is_mapping_file(kwargs['uif_path']):
            # if file is readable as a basic TSV and contains all the required
            # headers, then treat this as a mapping file, even if it's an
            # invalid one.
            assay_type = ASSAY_NAME_AMPLICON
            # for Amplicon runs, the lane_number is always one, even if the
            # user supplies another value in the UI.
            kwargs['lane_number'] = 1
            # NB: For now, let's assume all Amplicon runs are Illumina, since
            # the entire Amplicon pipeline assumes as much.
            instrument_type = 'Illumina'
        else:
            raise WorkflowError("Your uploaded file doesn't appear to be a "
                                "sample-sheet or a mapping-file.")

        for workflow in WorkflowFactory.WORKFLOWS:
            if workflow.assay_name == assay_type:
                if workflow.instrument_name == instrument_type:
                    # return instantiated workflow object
                    kwargs
                    return workflow(**kwargs)

        raise ValueError(f"Assay type '{assay_type}' and Instrument type "
                         "'{instrument_type}' did not match any known workflow"
                         " configuration")
