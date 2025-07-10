from os.path import join, exists, split
from os import walk, makedirs, listdir, environ
import pandas as pd
from json import dumps
from subprocess import Popen, PIPE
from glob import glob
from shutil import copyfile
import logging
from shutil import rmtree
from .Assays import ASSAY_NAME_AMPLICON


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

        # initializing member variables known to be used in class and/or in
        # mixins.
        self.cmds_log_path = None
        self.cmds = None
        self.has_replicates = None
        self.job_pool_size = None
        self.mandatory_attributes = []
        self.master_qiita_job_id = None
        self.output_path = None
        self.pipeline = None
        self.prep_copy_index = 0
        self.prep_file_paths = None
        self.qclient = None
        self.run_prefixes = {}
        self.samples_in_qiita = None
        self.sample_state = None
        self.sifs = None
        self.skip_steps = []
        self.special_map = None
        self.status_msg = ''
        self.touched_studies_prep_info = None
        self.tube_id_map = None
        self.directories_to_check = [
            'ConvertJob', 'NuQCJob', 'FastQCJob', 'GenPrepFileJob']

        if 'status_update_callback' in kwargs:
            self.status_update_callback = kwargs['status_update_callback']
        else:
            self.status_update_callback = None

    def confirm_mandatory_attributes(self):
        """
        Confirms that all mandatory attributes are present in kwargs.
        """
        absent_list = []

        for attribute in self.mandatory_attributes:
            if attribute not in self.kwargs:
                absent_list.append(attribute)

        if absent_list:
            raise ValueError(f"The following values must also be defined in "
                             f"kwargs for {self.__class__.__name__} workflows"
                             + ": " + ', '.join(absent_list))

    def job_callback(self, jid, status):
        """
        Update main status message w/current child job status.
        """
        if self.status_update_callback:
            self.status_update_callback(self.status_msg +
                                        f" ({jid}: {status})")

    def update_status(self, msg, step_number, total_steps):
        """
        Prettify status message before updating.
        """

        # When this method is called, a new status message is created.
        # This is saved so that child jobs can use job_callback() to update
        # this message.

        # set self.status_msg even if self.status_update_callback() is None.
        msg = "Step %d of %d: %s" % (step_number, total_steps, msg)
        self.status_msg = msg

        if self.status_update_callback:
            self.status_update_callback(self.status_msg)

    def what_am_i(self):
        """
        Returns text description of Workflow's Instrument & Assay mixins.
        """
        return (f"Instrument: {self.protocol_type}" + "\t" +
                f"Assay: {self.assay_type}")

    def pre_check(self):
        if self.is_restart:
            self._get_tube_ids_from_qiita()
        else:
            # since one of the objectives of SPP is to generate prep-info files
            # and automatically load them into Qiita, confirm that all studies
            # mentioned in the sample-sheet/pre-prep do not contain sample
            # metadata that would cause an error in the pipeline after
            # processing has already completed but the results have not yet
            # been loaded.
            self._project_metadata_check()

            # compare sample-ids/tube-ids in sample-sheet/mapping file
            # against what's in Qiita. Results are a list of dictionaries, one
            # per project.
            results = self._compare_samples_against_qiita()

            # obtain a list of non-zero counts of samples missing in Qiita, one
            # for each project. The names of the projects are unimportant. We
            # want to abort early if any project in the sample-sheet/pre-prep
            # file contains samples that aren't registered in Qiita.
            tmp = [len(project['samples_not_in_qiita']) for project in results]
            missing_counts = [count for count in tmp if count != 0]

            if missing_counts:
                msgs = []
                for result in results:
                    msgs += result['messages']

                if msgs:
                    raise WorkflowError('\n'.join(msgs))

    def generate_special_map(self):
        """
        Generates a list of tuples to support pipeline processing.
        :return: A list of triplets.
        """
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
        Generates sample-info files for each project, containing
        metadata on blanks.
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
        Updates the blanks registered in a given project in Qiita.
        :return:
        """

        for sif_path in self.sifs:
            # get study_id from sif_file_name ...something_14385_blanks.tsv
            study_id = self.pipeline.get_qiita_id_from_sif_fp(sif_path)

            df = pd.read_csv(sif_path, delimiter='\t', dtype=str)

            # Prepend study_id to make them compatible w/list from Qiita.
            df['sample_name'] = f'{study_id}.' + df['sample_name'].astype(str)

            # SIFs only contain blanks. Get the list of potentially new blanks.
            blanks = df['sample_name'].tolist()
            if len(blanks) == 0:
                # we have nothing to do so let's return early
                return

            # Get list of samples already registered in Qiita
            # (will include any already-registered blanks)
            from_qiita = self.qclient.get(f'/api/v1/study/{study_id}/samples')

            # Generate list of blanks that need to be ADDED to Qiita.
            new_blanks = (set(blanks) | set(from_qiita)) - set(from_qiita)

            if len(new_blanks):
                # Generate dummy entries for each new blank, if any.
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

                   (['failed_samples.html',
                     'touched_studies.html',
                     'MultiQCJob/multiqc'],
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
                   self.pipeline.is_sif_fp(x)]

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
            self._get_tube_ids_from_qiita()
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

            # do not include blanks. If they are unregistered, we will add
            # them downstream.
            samples = {smpl for smpl in samples if
                       not self.pipeline.sample_sheet.sample_is_a_blank(smpl)}

            msgs.append(f"The total number of samples found in {p_name} that "
                        f"aren't BLANK is: {len(samples)}")

            results_sn = self._process_sample_names(p_name, qiita_id,
                                                    samples)
            rsn = results_sn[0]
            msgs.append('Number of sample-names not in Qiita: '
                        f'{len(rsn)}; {list(rsn)[:3]}')

            use_tids = False

            if len(results_sn[0]) == 0:
                msgs.append(f"All values in sheet matched sample-names "
                            f"registered with {p_name}")
            else:
                # not all values were matched to sample-names.
                # check for possible match w/tube-ids, if defined in project.
                results_tid = self._process_tube_ids(qiita_id, samples)
                if results_tid:
                    rtid = results_tid[0]
                    msgs.append('Number of tube-ids not in Qiita: '
                                f'{len(rtid)}; {list(rtid)[:3]}')

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

    @classmethod
    def _determine_orientation(cls, file_name):
        # aka forward, reverse, and indexed reads
        orientations = ['R1', 'R2', 'I1', 'I2']

        results = []

        # assume orientation is always present in the file's name.
        # assume that it is of one of the four forms above.
        # assume that it is always the right-most occurance of the four
        # orientations above.
        # assume that orientation is encapsulated with either '_' or '.'
        # e.g.: '_R1_', '.I2.'.
        # assume users can and will include any or all of the four
        # orientation as part of their filenames as well. e.g.:
        # ABC_7_04_1776_R1_SRE_S3_L007_R2_001.trimmed.fastq.gz
        for o in orientations:
            variations = [f"_{o}_", f".{o}."]
            for v in variations:
                # rfind searches from the end of the string, rather than
                # its beginning. It returns the position in the string
                # where the substring begins.
                results.append((file_name.rfind(v), o))

        # the orientation will be the substring found with the maximum
        # found value for pos. That is, it will be the substring that
        # begins at the rightest most position in the file name.
        results.sort(reverse=True)

        pos, orientation = results[0]

        # if no orientations were found, then return None.
        return None if pos == -1 else orientation

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
            _, file_name = split(fastq_file)
            orientation = self._determine_orientation(file_name)
            if orientation in ['I1', 'I2']:
                files['raw_barcodes'].append(fastq_file)
            elif orientation == 'R1':
                files['raw_forward_seqs'].append(fastq_file)
            elif orientation == 'R2':
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

    def _get_tube_ids_from_qiita(self):
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
            qsam, tids = self.get_samples_in_qiita(self.qclient, qiita_id)

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

    def determine_steps_to_skip(self):
        out_dir = self.pipeline.output_path

        # check if the test flag is on!
        test = False
        if 'PrepNuQCJob_TEST' in environ and \
                environ['PrepNuQCJob_TEST'] == 'true':
            test = True

        # Although amplicon runs don't perform host-filtering,
        # the output from ConvertJob is still copied and organized into
        # a form suitable for FastQCJob to process. Hence the presence or
        # absence of a 'NuQCJob' directory is still a thing (for now)
        for directory in self.directories_to_check:
            if exists(join(out_dir, directory)):
                if exists(join(out_dir, directory, 'job_completed')):
                    # this step completed successfully but
                    # TRIJ_Post_Processing is a special case and we
                    # need to look for post_processing_completed
                    if directory == 'TRIJ_Post_Processing':
                        if not exists(join(out_dir, directory,
                                      'post_processing_completed')):
                            rmtree(join(out_dir, directory))
                            break
                    self.skip_steps.append(directory)
                else:
                    if not test:
                        # work stopped before this job could be completed.
                        rmtree(join(out_dir, directory))
                    break
