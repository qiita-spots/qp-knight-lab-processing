from collections import defaultdict
from itertools import chain
from json import dumps
from metapool import KLSampleSheet
from os import makedirs, walk
from os.path import join, exists
from sequence_processing_pipeline.ConvertJob import ConvertJob
from sequence_processing_pipeline.FastQCJob import FastQCJob
from sequence_processing_pipeline.GenPrepFileJob import GenPrepFileJob
from sequence_processing_pipeline.PipelineError import PipelineError
from sequence_processing_pipeline.Pipeline import Pipeline
from sequence_processing_pipeline.QCJob import QCJob
from subprocess import Popen, PIPE
import pandas as pd


class FailedSamplesRecord:
    def __init__(self, output_dir, samples):
        # because we want to write out the list of samples that failed after
        # each Job is run, and we want to organize that output by project, we
        # need to keep a running state of failed samples, and reuse the method
        # to reorganize the running-results and write them out to disk.

        self.output_path = join(output_dir, 'failed_samples.html')

        # create an initial dictionary with sample-ids as keys and their
        # associated project-name and status as values. Afterwards, we'll
        # filter out the sample-ids w/no status (meaning they were
        # successfully processed) before writing the failed entries out to
        # file.
        self.sample_state = {x.Sample_ID: [x.Sample_Project, None] for x in
                             samples}

    def write(self, failed_ids, job_name):
        for failed_id in failed_ids:
            # as a rule, if a failed_id were to appear in more than one
            # audit(), preserve the earliest failure, rather than the
            # latest one.
            if self.sample_state[failed_id][1] is None:
                self.sample_state[failed_id][1] = job_name

        # filter out the sample-ids w/out a failure status
        filtered_fails = {x: self.sample_state[x] for x in self.sample_state if
                          self.sample_state[x][1] is not None}

        data = []
        for sample_id in filtered_fails:
            project_name = filtered_fails[sample_id][0]
            failed_at = filtered_fails[sample_id][1]
            data.append({'Project': project_name, 'Sample ID': sample_id,
                         'Failed at': failed_at})
        df = pd.DataFrame(data)

        with open(self.output_path, 'w') as f:
            f.write(df.to_html(border=2, index=False, justify="left",
                               render_links=True, escape=False))


class Step:
    '''
    The base Step class wraps the creation and running of the Job classes
    that are common to both Amplicon and Metagenomic Pipeline. Functionality
    specific to one pipeline or the other is handled in the appropriate
    subclass and makes calls to this base class as needed. In this way the
    codebase is kept DRY.
    '''

    AMPLICON_TYPE = 'Amplicon'
    METAGENOMIC_TYPE = 'Metagenomic'
    METATRANSCRIPTOMIC_TYPE = 'Metatranscriptomic'
    META_TYPES = {METAGENOMIC_TYPE, METATRANSCRIPTOMIC_TYPE}
    ALL_TYPES = META_TYPES.union(AMPLICON_TYPE)
    AMPLICON_SUB_TYPES = {'16S', '18S', 'ITS'}

    def __init__(self, pipeline, master_qiita_job_id,
                 status_update_callback=None,
                 lane_number=None):
        if pipeline is None:
            raise ValueError("A pipeline object is needed to initialize Step")

        if master_qiita_job_id is None:
            raise ValueError("A Qiita job-id is needed to initialize Step")

        self.pipeline = pipeline
        self.lane_number = lane_number
        self.generated_artifact_name = \
            f'{self.pipeline.run_id}_{self.lane_number}'
        self.master_qiita_job_id = master_qiita_job_id

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
        self.prep_file_paths = None
        self.sifs = None
        self.tube_id_map = None
        self.samples_in_qiita = None
        self.output_path = None
        self.prep_file_paths = None
        self.sample_state = None
        self.special_map = None
        self.touched_studies_prep_info = None

    @classmethod
    def generate_pipeline(cls, pipeline_type, input_file_path, lane_number,
                          config_fp,
                          run_identifier, out_dir, job_id):
        if pipeline_type in Step.META_TYPES:
            cls.update_sample_sheet(input_file_path, lane_number)
            return Pipeline(config_fp, run_identifier, input_file_path, None,
                            out_dir, job_id, pipeline_type)
        elif pipeline_type == Step.AMPLICON_TYPE:
            return Pipeline(config_fp, run_identifier, None, input_file_path,
                            out_dir, job_id, pipeline_type)
        else:
            raise PipelineError(
                f"'{pipeline_type}' is not a valid Pipeline type.")

    @classmethod
    def update_sample_sheet(cls, sample_sheet_path, lane_number):
        # use KLSampleSheet functionality to add/overwrite lane number.
        sheet = KLSampleSheet(sample_sheet_path)
        for sample in sheet:
            sample['Lane'] = f'{lane_number}'

        with open(sample_sheet_path, 'w') as f:
            sheet.write(f)

    @classmethod
    def parse_prep_file(cls, prep_file_path):
        metadata = pd.read_csv(prep_file_path,
                               dtype=str,
                               delimiter='\t',
                               # forces Pandas to not make the first column
                               # the index even when the values appear numeric.
                               index_col=False)

        if metadata is None:
            raise ValueError(f"{prep_file_path} does not exist.")

        metadata.set_index('sample_name', inplace=True)

        # convert to standard dictionary.
        return metadata.to_dict('index')

    def generate_special_map(self, qclient):
        # this function should be able to be tested by passing in simulated =
        # results from qclient.

        # trimmed files are stored by qiita_id. Find the qiita_id
        # associated with each project and ensure a subdirectory exists
        # for when it comes time to move the trimmed files.

        special_map = []
        results = qclient.get("/qiita_db/artifacts/types/")
        projects = self.pipeline.get_project_info()
        for project in projects:
            upload_path = join(results['uploads'], project['qiita_id'])
            makedirs(upload_path, exist_ok=True)
            special_map.append((project['project_name'], upload_path,
                                project['qiita_id']))

        self.special_map = special_map

    def update_prep_templates(self, qclient, prep_file_paths, pipeline_type):
        '''
        Update prep-template info in Qiita. Get dict of prep-ids by study-id.
        :param qclient:
        :param prep_file_paths:
        :return: A dict of lists of prep-ids, keyed by study-id.
        '''
        results = defaultdict(list)

        for study_id in prep_file_paths:
            for prep_file_path in prep_file_paths[study_id]:
                metadata = Step.parse_prep_file(prep_file_path)
                data = {'prep_info': dumps(metadata),
                        'study': study_id,
                        'data_type': None,
                        'job-id': self.master_qiita_job_id,
                        'name': self.generated_artifact_name}
                if pipeline_type in Step.META_TYPES:
                    data['data_type'] = pipeline_type
                elif pipeline_type == Step.AMPLICON_TYPE:
                    if 'target_gene' in metadata[list(metadata.keys())[0]]:
                        tg = metadata[list(metadata.keys())[0]]['target_gene']
                        for key in Step.AMPLICON_SUB_TYPES:
                            if key in tg:
                                data['data_type'] = key
                    else:
                        raise ValueError("target_gene must be specified for "
                                         "amplicon type")
                else:
                    raise ValueError(f"'{pipeline_type}' is not a valid "
                                     "pipeline type")

                reply = qclient.post('/qiita_db/prep_template/', data=data)
                prep_id = reply['prep']
                results[study_id].append(prep_id)

        self.touched_studies_prep_info = results
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

    def _convert_bcl_to_fastq(self, config, input_file_path):
        convert_job = ConvertJob(self.pipeline.run_dir,
                                 self.pipeline.output_path,
                                 input_file_path,
                                 config['queue'],
                                 config['nodes'],
                                 config['nprocs'],
                                 config['wallclock_time_in_hours'],
                                 config['per_process_memory_limit'],
                                 config['executable_path'],
                                 config['modules_to_load'],
                                 self.master_qiita_job_id)

        convert_job.run(callback=self.update_callback)

        return convert_job

    def _quality_control(self, config, input_file_path):
        qc_job = QCJob(join(self.pipeline.output_path, 'ConvertJob'),
                       self.pipeline.output_path,
                       input_file_path,
                       config['minimap_databases'],
                       config['kraken2_database'],
                       config['queue'],
                       config['nodes'],
                       config['nprocs'],
                       config['wallclock_time_in_hours'],
                       config['job_total_memory_limit'],
                       config['fastp_executable_path'],
                       config['minimap2_executable_path'],
                       config['samtools_executable_path'],
                       config['modules_to_load'],
                       self.master_qiita_job_id,
                       self.job_pool_size,
                       config['job_max_array_length'])

        qc_job.run(callback=self.update_callback)

        return qc_job

    def _generate_reports(self):
        config = self.pipeline.configuration['fastqc']
        is_amplicon = self.pipeline.pipeline_type == Step.AMPLICON_TYPE
        fastqc_job = FastQCJob(self.pipeline.run_dir,
                               self.pipeline.output_path,
                               join(self.pipeline.output_path, 'ConvertJob'),
                               join(self.pipeline.output_path, 'QCJob'),
                               config['nprocs'],
                               config['nthreads'],
                               config['fastqc_executable_path'],
                               config['modules_to_load'],
                               self.master_qiita_job_id,
                               config['queue'],
                               config['nodes'],
                               config['wallclock_time_in_hours'],
                               config['job_total_memory_limit'],
                               self.job_pool_size,
                               config['multiqc_config_file_path'],
                               config['job_max_array_length'],
                               is_amplicon)

        fastqc_job.run(callback=self.update_callback)

        return fastqc_job

    def _generate_prep_file(self, config, input_file_path, seqpro_path,
                            project_names):
        is_amplicon = self.pipeline.pipeline_type == Step.AMPLICON_TYPE

        gpf_job = GenPrepFileJob(
            self.pipeline.run_dir,
            join(self.pipeline.output_path, 'ConvertJob'),
            join(self.pipeline.output_path, 'QCJob'),
            self.pipeline.output_path,
            input_file_path,
            seqpro_path,
            project_names,
            config['modules_to_load'],
            self.master_qiita_job_id,
            is_amplicon=is_amplicon)

        gpf_job.run(callback=self.update_callback)

        # concatenate the lists of paths across all study_ids into a single
        # list. Replace sample-names w/tube-ids in all relevant prep-files.
        preps = list(chain.from_iterable(gpf_job.prep_file_paths.values()))
        self._overwrite_prep_files(preps)

        return gpf_job

    def _generate_commands(self):
        out_dir = self.pipeline.output_path
        qclient = self.qclient

        cmds = ['tar zcvf logs-ConvertJob.tgz ConvertJob/logs',
                'tar zcvf logs-FastQCJob.tgz FastQCJob/logs',
                'tar zcvf reports-FastQCJob.tgz FastQCJob/fastqc',
                'tar zcvf logs-GenPrepFileJob.tgz GenPrepFileJob/logs',
                'tar zcvf prep-files.tgz GenPrepFileJob/PrepFiles']
        self.cmds = [f'cd {out_dir}; {x}' for x in cmds]

        data = []
        for project, _, qiita_id in self.special_map:
            if self.pipeline.pipeline_type in Step.META_TYPES:
                self.cmds.append(f'cd {out_dir}; tar zcvf reports-QCJob.tgz '
                                 f'QCJob/{project}/fastp_reports_dir')

            if len(self.touched_studies_prep_info[qiita_id]) != 1:
                raise ValueError(
                    f"Too many preps for {qiita_id}: "
                    f"{self.touched_studies_prep_info[qiita_id]}")

            prep_id = self.touched_studies_prep_info[qiita_id][0]
            surl = f'{qclient._server_url}/study/description/{qiita_id}'
            prep_url = (f'{qclient._server_url}/study/description/'
                        f'{qiita_id}?prep_id={prep_id}')

            bd = f'{out_dir}/QCJob/{project}'
            if exists(f'{bd}/filtered_sequences'):
                atype = 'per_sample_FASTQ'
                af = [f for f in walk(f'{bd}/filtered_sequences/*.fastq.gz')]
                files = {'raw_forward_seqs': [], 'raw_reverse_seqs': []}
                for f in af:
                    if '.R1.' in f:
                        files['raw_forward_seqs'].append(f)
                    elif '.R2.' in f:
                        files['raw_reverse_seqs'].append(f)
                    else:
                        raise ValueError(f'Not recognized file: {f}')
            elif exists(f'{bd}/trimmed_sequences'):
                atype = 'per_sample_FASTQ'
                af = [f for f in walk(f'{bd}/trimmed_sequences/*.fastq.gz')]
                files = {'raw_forward_seqs': [], 'raw_reverse_seqs': []}
                for f in af:
                    if '.R1.' in f:
                        files['raw_forward_seqs'].append(f)
                    elif '.R2.' in f:
                        files['raw_reverse_seqs'].append(f)
                    else:
                        raise ValueError(f'Not recognized file: {f}')
            elif exists(f'{bd}/amplicon'):
                atype = 'FASTQ'
                af = sorted([f for f in walk(f'{bd}/amplicon/*.fastq.gz')])
                files = {'raw_barcoes': af[0],
                         'raw_forward_seqs': af[1],
                         'raw_reverse_seqs': af[2]}
            else:
                raise PipelineError("QCJob output not in expected location")
            # ideally we would use the email of the user that started the SPP
            # run but at this point there is no easy way to retrieve it
            data = {'user_email': 'qiita.help@gmail.com',
                    'prep_id': prep_id,
                    'artifact_type': atype,
                    'command_artifact_name': self.generated_artifact_name,
                    'files': files}
            job_id = qclient.post('/qiita_db/artifact/', data=data)

            data.append({'Project': project, 'Qiita Study ID': qiita_id,
                         'Qiita Prep ID': prep_id, 'Qiita URL': surl,
                         'Prep URL': prep_url, 'Linking JobID': job_id})

        df = pd.DataFrame(data)
        with open(join(out_dir, 'touched_studies.html'), 'w') as f:
            f.write(df.to_html(border=2, index=False, justify="left",
                               render_links=True, escape=False))

        if exists(join(out_dir, 'failed_samples.html')):
            tmp = f'cd {out_dir}; mv failed_samples.html final_results'
            self.cmds.append(tmp)

    def write_commands_to_output_path(self):
        self.cmds_log_path = join(self.pipeline.output_path, 'cmds.log')
        with open(self.cmds_log_path, 'w') as f:
            for cmd in self.cmds:
                f.write(f'{cmd}\n')

    def execute_commands(self):
        # execute the list of commands in order
        for cmd in self.cmds:
            p = Popen(cmd, universal_newlines=True, shell=True,
                      stdout=PIPE, stderr=PIPE)
            std_out, std_err = p.communicate()
            return_code = p.returncode

            if return_code != 0:
                # during testing, ignore processes that fail and continue
                # to test other commands.
                raise PipelineError(f"'{cmd}' returned {return_code}")

    def generate_sifs(self, qclient):
        from_qiita = {}

        for study_id in self.prep_file_paths:
            samples = list(qclient.get(f'/api/v1/study/{study_id}/samples'))
            from_qiita[study_id] = samples

        add_sif_info = []

        qid_pn_map = {x['qiita_id']: x['project_name'] for
                      x in self.pipeline.get_project_info()}

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

    def get_prep_file_paths(self):
        return self.prep_file_paths

    def get_tube_ids_from_qiita(self, qclient):
        # Update get_project_info() so that it can return a list of
        # samples in projects['samples']. Include blanks in projects['blanks']
        # just in case there are duplicate qiita_ids
        qiita_ids = [x['qiita_id'] for x in
                     self.pipeline.get_project_info(short_names=True)]

        tids_by_qiita_id = {}
        sample_names_by_qiita_id = {}

        for qiita_id in qiita_ids:
            # Qiita returns a set of sample-ids in qsam and a dictionary where
            # sample-names are used as keys and tube-ids are their values.
            qsam, tids = self.get_samples_in_qiita(qclient, qiita_id)

            if tids is None:
                sample_names_by_qiita_id[str(qiita_id)] = qsam
            else:
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
        self.tube_id_map = tids_by_qiita_id
        self.samples_in_qiita = sample_names_by_qiita_id

    def compare_samples_against_qiita(self):
        projects = self.pipeline.get_project_info(short_names=True)

        results = []
        for project in projects:
            project_name = project['project_name']
            qiita_id = str(project['qiita_id'])

            # get list of samples as presented by the sample-sheet or mapping
            # file and confirm that they are all registered in Qiita.
            samples = set(self.pipeline.get_sample_names())

            # strip any leading zeroes from the sample-ids. Note that
            # if a sample-id has more than one leading zero, all of
            # them will be removed.
            samples = {x.lstrip('0') for x in samples}

            # just get a list of the tube-ids themselves, not what they map
            # to.
            if qiita_id in self.tube_id_map:
                # if map is not empty
                tids = [self.tube_id_map[qiita_id][x] for x in
                        self.tube_id_map[qiita_id]]
                not_in_qiita = samples - set(tids)
                examples = tids[:5]
                used_tids = True
            else:
                # assume project is in samples_in_qiita
                not_in_qiita = samples - set(self.samples_in_qiita)
                examples = list(samples)[:5]
                used_tids = False

            # convert to strings before returning
            examples = [str(x) for x in examples]

            if not_in_qiita:
                # if there are actual differences
                results.append({'samples_not_in_qiita': not_in_qiita,
                                'examples_in_qiita': examples,
                                'project_name': project_name,
                                'tids': used_tids})
        return results

    @classmethod
    def _replace_with_tube_ids(cls, prep_file_path, qiita_id, tube_id_map):
        # passing tube_id_map as a parameter allows for easier testing.
        df = pd.read_csv(prep_file_path, sep='\t', dtype=str, index_col=False)
        # save copy of sample_name column as 'old_sample_name'
        df['old_sample_name'] = df['sample_name']
        for i in df.index:
            sample_name = df.at[i, "sample_name"]
            # blanks do not get their names swapped.
            if sample_name.startswith('BLANK'):
                continue

            # remove leading zeroes if they exist to match Qiita results.
            sample_name = sample_name.lstrip('0')
            if sample_name in tube_id_map[qiita_id]:
                df.at[i, "sample_name"] = tube_id_map[qiita_id][sample_name]

        df.to_csv(prep_file_path, index=False, sep="\t")

    def _overwrite_prep_files(self, prep_file_paths):
        # replaces sample-names in prep-files with tube-ids according to
        # a dict with project-names as keys and another dict as a value.
        # this dict uses sample-names as keys and tube-ids as values.
        if self.tube_id_map is None:
            raise ValueError("get_tube_ids_from_qiita() was not called")

        projects = self.pipeline.get_project_info(short_names=True)

        for project in projects:
            project_name = project['project_name']
            qiita_id = str(project['qiita_id'])

            if qiita_id not in self.tube_id_map:
                continue

            # prep files are named in the form:
            # 20220423_FS10001773_12_BRB11603-0615.Matrix_Tube_LBM_14332.1.tsv
            # search on project_name vs qiita_id since it's slightly more
            # unique.
            matching_files = [x for x in prep_file_paths if project_name in x]

            if len(matching_files) == 0:
                continue

            if len(matching_files) > 1:
                raise ValueError("More than one match found for project "
                                 f"'{project_name}': {str(matching_files)}")

            Step._replace_with_tube_ids(matching_files[0],
                                        qiita_id,
                                        self.tube_id_map)

    def update_blanks_in_qiita(self, qclient):
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
            from_qiita = qclient.get(f'/api/v1/study/{study_id}/samples')
            from_qiita = [x for x in from_qiita if
                          x.startswith(f'{study_id}.BLANK')]

            # Generate list of BLANKs that need to be ADDED to Qiita.
            new_blanks = (set(blanks) | set(from_qiita)) - set(from_qiita)

            if len(new_blanks):
                # Generate dummy entries for each new BLANK, if any.
                categories = qclient.get(f'/api/v1/study/{study_id}/samples/'
                                         'info')['categories']

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
                qclient.http_patch(f'/api/v1/study/{study_id}/samples',
                                   data=dumps(data))

                return data

    def execute_pipeline(self, qclient, increment_status, update=True):
        '''
        Executes steps of pipeline in proper sequence.
        :param qclient: Qiita client library or equivalent.
        :param increment_status: callback function to increment status.
        :param update: Set False to prevent updates to Qiita.
        :return: None
        '''
        self.generate_special_map(qclient)

        increment_status()
        self.convert_bcl_to_fastq()

        increment_status()
        self.quality_control()

        increment_status()
        self.generate_reports()

        increment_status()
        self.generate_prep_file()

        increment_status()
        self.sifs = self.generate_sifs(qclient)

        increment_status()

        if update:
            self.update_blanks_in_qiita(qclient)

        prep_file_paths = self.get_prep_file_paths()

        increment_status()
        ptype = self.pipeline.pipeline_type

        if update:
            self.update_prep_templates(qclient, prep_file_paths, ptype)

        self.generate_touched_studies(qclient)

        increment_status()
        self.generate_commands(qclient)

        increment_status()
        if update:
            self.execute_commands()
