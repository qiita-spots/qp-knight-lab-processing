from itertools import chain
from collections import defaultdict
from json import dumps
from metapool import load_sample_sheet
from os import makedirs, walk, listdir
from os.path import join, exists, split, basename, dirname
from sequence_processing_pipeline.ConvertJob import ConvertJob
from sequence_processing_pipeline.FastQCJob import FastQCJob
from sequence_processing_pipeline.GenPrepFileJob import GenPrepFileJob
from sequence_processing_pipeline.PipelineError import PipelineError
from sequence_processing_pipeline.Pipeline import Pipeline
from sequence_processing_pipeline.NuQCJob import NuQCJob
from subprocess import Popen, PIPE
import pandas as pd
from glob import glob
from shutil import copyfile


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
                 lane_number=None, is_restart=False):
        if pipeline is None:
            raise ValueError("A pipeline object is needed to initialize Step")

        if master_qiita_job_id is None:
            raise ValueError("A Qiita job-id is needed to initialize Step")

        self.pipeline = pipeline
        self.lane_number = lane_number
        self.generated_artifact_name = \
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
        self.special_map = None
        self.touched_studies_prep_info = None
        self.run_prefixes = {}
        self.prep_copy_index = 0

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
        sheet = load_sample_sheet(sample_sheet_path)
        for sample in sheet:
            sample['Lane'] = f'{lane_number}'

        with open(sample_sheet_path, 'w') as f:
            sheet.write(f)

    @classmethod
    def parse_prep_file(cls, prep_file_path, convert_to_dict=True):
        metadata = pd.read_csv(prep_file_path,
                               dtype=str,
                               delimiter='\t',
                               # forces Pandas to not make the first column
                               # the index even when the values appear numeric.
                               index_col=False)

        if metadata is None:
            raise ValueError(f"{prep_file_path} does not exist.")

        metadata.set_index('sample_name', inplace=True)

        if convert_to_dict:
            return metadata.to_dict('index')
        else:
            return metadata

    def generate_artifact_name(self, prep_file_path):
        a_name = f'{self.pipeline.run_id}_{self.lane_number}'
        repl_num = basename(dirname(prep_file_path))

        if self.has_replicates is True:
            # this is a replicate sheet file.
            # append a replication number to each name to
            # make it unique from other replicates.
            # return ('%s_r%s' % (a_name, result[1]), True)
            return ('%s_r%s' % (a_name, repl_num), True)
        else:
            # this is a normal pre-prep or sample-sheet.
            return (a_name, False)

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

    def get_data_type(self, prep_file_path):
        raise ValueError("get_data_type() not implemented for base-class.")

    def update_prep_templates(self, qclient):
        '''
        Update prep-template info in Qiita. Get dict of prep-ids by study-id.
        :param qclient:
        :return: A dict of lists of prep-ids, keyed by study-id.
        '''
        results = defaultdict(list)

        for study_id in self.prep_file_paths:
            for prep_fp in self.prep_file_paths[study_id]:
                metadata = Step.parse_prep_file(prep_fp)
                afact_name, is_repl = self.generate_artifact_name(prep_fp)
                data = {'prep_info': dumps(metadata),
                        'study': study_id,
                        'data_type': None,
                        'job-id': self.master_qiita_job_id,
                        'name': afact_name}
                if self.pipeline.pipeline_type in Step.META_TYPES:
                    data['data_type'] = self.pipeline.pipeline_type
                elif self.pipeline.pipeline_type == Step.AMPLICON_TYPE:
                    if 'target_gene' in metadata[list(metadata.keys())[0]]:
                        tg = metadata[list(metadata.keys())[0]]['target_gene']
                        for key in Step.AMPLICON_SUB_TYPES:
                            if key in tg:
                                data['data_type'] = key

                        if data['data_type'] is None:
                            raise ValueError("data_type could not be "
                                             "determined from target_gene "
                                             "column")
                    else:
                        raise ValueError("target_gene must be specified for "
                                         "amplicon type")
                else:
                    raise ValueError(f"'{self.pipeline.pipeline_type}' is not "
                                     " a valid pipeline type")

                reply = qclient.post('/qiita_db/prep_template/', data=data)
                prep_id = reply['prep']
                results[study_id].append((prep_id, afact_name, is_repl))
                self.run_prefixes[prep_id] = [metadata[sample]['run_prefix']
                                              for sample in metadata]

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
                                 config['wallclock_time_in_minutes'],
                                 config['per_process_memory_limit'],
                                 config['executable_path'],
                                 config['modules_to_load'],
                                 self.master_qiita_job_id)

        convert_job.run(callback=self.update_callback)

        return convert_job

    def _quality_control(self, config, input_file_path):
        nuqc_job = NuQCJob(join(self.pipeline.output_path, 'ConvertJob'),
                           self.pipeline.output_path,
                           input_file_path,
                           config['minimap2_databases'],
                           config['queue'],
                           config['nodes'],
                           config['wallclock_time_in_minutes'],
                           config['job_total_memory_limit'],
                           config['fastp_executable_path'],
                           config['minimap2_executable_path'],
                           config['samtools_executable_path'],
                           config['modules_to_load'],
                           self.master_qiita_job_id,
                           config['job_max_array_length'],
                           config['known_adapters_path'],
                           bucket_size=config['bucket_size'],
                           length_limit=config['length_limit'],
                           cores_per_task=config['cores_per_task'])

        nuqc_job.run(callback=self.update_callback)

        return nuqc_job

    def _generate_reports(self):
        config = self.pipeline.config_profile['profile']['configuration']
        is_amplicon = self.pipeline.pipeline_type == Step.AMPLICON_TYPE
        fastqc_job = FastQCJob(self.pipeline.run_dir,
                               self.pipeline.output_path,
                               join(self.pipeline.output_path, 'ConvertJob'),
                               join(self.pipeline.output_path, 'NuQCJob'),
                               config['fastqc']['nprocs'],
                               config['fastqc']['nthreads'],
                               config['fastqc']['fastqc_executable_path'],
                               config['fastqc']['modules_to_load'],
                               self.master_qiita_job_id,
                               config['fastqc']['queue'],
                               config['fastqc']['nodes'],
                               config['fastqc']['wallclock_time_in_minutes'],
                               config['fastqc']['job_total_memory_limit'],
                               self.job_pool_size,
                               config['fastqc']['multiqc_config_file_path'],
                               config['fastqc']['job_max_array_length'],
                               is_amplicon)

        fastqc_job.run(callback=self.update_callback)

        return fastqc_job

    def _generate_prep_file(self, config, input_file_path, seqpro_path):
        is_amplicon = self.pipeline.pipeline_type == Step.AMPLICON_TYPE

        gpf_job = GenPrepFileJob(
            self.pipeline.run_dir,
            join(self.pipeline.output_path, 'ConvertJob'),
            join(self.pipeline.output_path, 'NuQCJob'),
            self.pipeline.output_path,
            input_file_path,
            seqpro_path,
            config['modules_to_load'],
            self.master_qiita_job_id,
            is_amplicon=is_amplicon)

        gpf_job.run(callback=self.update_callback)

        # concatenate the lists of paths across all study_ids into a single
        # list. Replace sample-names w/tube-ids in all relevant prep-files.
        preps = list(chain.from_iterable(gpf_job.prep_file_paths.values()))
        self._overwrite_prep_files(preps)

        return gpf_job

    def _helper_process_fastp_report_dirs(self):
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

    def _helper_process_blanks(self):
        results = [x for x in listdir(self.pipeline.output_path) if
                   x.endswith('_blanks.tsv')]

        results.sort()

        if len(results) > 0:
            return 'tar zcvf sample-files.tgz' + ' ' + ' '.join(results)

    def _helper_process_operations(self):
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

    def generate_commands(self):
        cmds = self._helper_process_operations()

        result = self._helper_process_fastp_report_dirs()

        if result:
            cmds.append(result)

        result = self._helper_process_blanks()

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

        self.write_commands_to_output_path()

    def _get_fastq_files(self, out_dir, project):
        af = None
        sub_folders = ['amplicon', 'filtered_sequences', 'trimmed_sequences']
        for sub_folder in sub_folders:
            sf = f'{out_dir}/NuQCJob/{project}/{sub_folder}'
            if exists(sf):
                af = [f for f in glob(f'{sf}/*.fastq.gz')]
                break
        if af is None or not af:
            raise PipelineError("NuQCJob output not in expected location")

        files = {'raw_barcodes': [], 'raw_forward_seqs': [],
                 'raw_reverse_seqs': []}

        for fastq_file in af:
            if '_I1_' in fastq_file:
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
        if self.pipeline.pipeline_type != Step.AMPLICON_TYPE:
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

    def load_preps_into_qiita(self, qclient):
        atype = 'per_sample_FASTQ'
        if self.pipeline.pipeline_type == Step.AMPLICON_TYPE:
            atype = 'FASTQ'

        data = []
        for project, _, qiita_id in self.special_map:
            fastq_files = self._get_fastq_files(
                self.pipeline.output_path, project)

            for vals in self.touched_studies_prep_info[qiita_id]:
                prep_id, artifact_name, is_repl = vals
                if self.pipeline.pipeline_type == Step.AMPLICON_TYPE:
                    if is_repl:
                        # for Amplicon runs, each prep needs a copy of the
                        # entire set of fastq files, because demuxing samples
                        # happens downstream. If we don't make copies of the
                        # files, Qiita will move the files when loading the
                        # first prep and they won't be available for the
                        # second prep and after.
                        # Note that this will leave the original files present
                        # in the working directory after processing instead of
                        # being moved.
                        working_set = self._copy_files(fastq_files)
                    else:
                        working_set = fastq_files
                else:
                    # for meta*omics, generate the subset of files used by
                    # this prep only.
                    working_set = {}
                    for key in fastq_files:
                        working_set[key] = []
                        for run_prefix in self.run_prefixes[prep_id]:
                            working_set[key] += [fastq for fastq in
                                                 fastq_files[key] if
                                                 run_prefix in fastq]

                    if is_repl:
                        working_set = self._copy_files(working_set)

                data.append(self._load_prep_into_qiita(
                    qclient, prep_id, artifact_name, qiita_id, project,
                    working_set, atype))

        df = pd.DataFrame(data)
        opath = join(self.pipeline.output_path, 'touched_studies.html')
        with open(opath, 'w') as f:
            f.write(df.to_html(border=2, index=False, justify="left",
                               render_links=True, escape=False))

        return df

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
        # should samples_in_qiita be none if tube_id_map is not?
        self.samples_in_qiita = sample_names_by_qiita_id

    def _compare_samples_against_qiita(self, qclient):
        projects = self.pipeline.get_project_info(short_names=True)
        self._get_tube_ids_from_qiita(qclient)

        results = []
        for project in projects:
            project_name = project['project_name']
            qiita_id = str(project['qiita_id'])

            # get list of samples as presented by the sample-sheet or mapping
            # file and confirm that they are all registered in Qiita.
            samples = set(self.pipeline.get_sample_names(project_name))

            # do not include BLANKs. If they are unregistered, we will add
            # them downstream.
            samples = {smpl for smpl in samples
                       if not smpl.startswith('BLANK')}

            # just get a list of the tube-ids themselves, not what they map
            # to.
            if qiita_id in self.tube_id_map:
                # if map is not empty
                tids = [self.tube_id_map[qiita_id][sample] for sample in
                        self.tube_id_map[qiita_id]]

                not_in_qiita = samples - set(tids)

                if not_in_qiita:
                    # strip any leading zeroes from the sample-ids. Note that
                    # if a sample-id has more than one leading zero, all of
                    # them will be removed.
                    not_in_qiita = set([x.lstrip('0') for x in samples]) - \
                                   set(tids)

                examples = tids[:5]
                used_tids = True
            else:
                # assume project is in samples_in_qiita
                not_in_qiita = samples - set(self.samples_in_qiita[qiita_id])
                examples = list(samples)[:5]
                used_tids = False

            # convert to strings before returning
            examples = [str(example) for example in examples]

            # return an entry for all projects, even when samples_not_in_qiita
            # is an empty list, as the information is still valuable.

            results.append({'samples_not_in_qiita': not_in_qiita,
                            'examples_in_qiita': examples,
                            'project_name': project_name,
                            'tids': used_tids})

        return results

    @classmethod
    def _replace_with_tube_ids(cls, prep_file_path, tube_id_map):
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

            reversed_map = {tube_id_map[k]: k for k in tube_id_map}
            if sample_name in reversed_map:
                df.at[i, "sample_name"] = reversed_map[sample_name]

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
            matching_files = [prep_file for prep_file in prep_file_paths if
                              project_name in prep_file]

            if len(matching_files) == 0:
                continue

            if len(matching_files) > 1:
                raise ValueError("More than one match found for project "
                                 f"'{project_name}': {str(matching_files)}")

            Step._replace_with_tube_ids(matching_files[0],
                                        self.tube_id_map[qiita_id])

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

    def precheck(self, qclient):
        # compare sample-ids/tube-ids in sample-sheet/mapping file
        # against what's in Qiita. Results are a list of dictionaries, one
        # per project.
        results = self._compare_samples_against_qiita(qclient)

        # obtain a list of non-zero counts of samples missing in Qiita, one
        # for each project. The names of the projects are unimportant. We
        # want to abort early if any project in the sample-sheet/pre-prep file
        # contains samples that aren't registered in Qiita.
        tmp = [len(project['samples_not_in_qiita']) for project in results]
        missing_counts = [count for count in tmp if count != 0]

        if missing_counts:
            msgs = []
            for comparison in results:
                not_in_qiita = list(comparison['samples_not_in_qiita'])
                not_in_qiita_count = len(not_in_qiita)
                examples_in_qiita = ', '.join(comparison['examples_in_qiita'])
                p_name = comparison['project_name']
                uses_tids = comparison['tids']

                msgs.append(
                    f"<br/><b>Project '{p_name}'</b> has {not_in_qiita_count} "
                    f"samples not registered in Qiita: {not_in_qiita[:5]}")

                msgs.append(f"Some registered samples in Project '{p_name}'"
                            f" include: {examples_in_qiita}")

                if uses_tids:
                    msgs.append(f"Project '{p_name}' is using tube-ids. You "
                                "may be using sample names in your file.")

            if msgs:
                raise PipelineError('\n'.join(msgs))

    def execute_pipeline(self, qclient, increment_status, update=True,
                         skip_steps=[]):
        '''
        Executes steps of pipeline in proper sequence.
        :param qclient: Qiita client library or equivalent.
        :param increment_status: callback function to increment status.
        :param update: Set False to prevent updates to Qiita.
        :return: None
        '''
        # this is performed even in the event of a restart.
        self.generate_special_map(qclient)

        # even if a job is being skipped, it's being skipped because it was
        # determined that it already completed successfully. Hence,
        # increment the status because we are still iterating through them.

        increment_status()
        if "ConvertJob" not in skip_steps:
            self.convert_bcl_to_fastq()

        increment_status()
        if "NuQCJob" not in skip_steps:
            self.quality_control()

        increment_status()
        if "FastQCJob" not in skip_steps:
            self.generate_reports()

        increment_status()
        if "GenPrepFileJob" not in skip_steps:
            self.generate_prep_file()

        # for now, simply re-run any line below as if it was a new job, even
        # for a restart. functionality is idempotent, except for the
        # registration of new preps in Qiita. These will simply be removed
        # manually.
        increment_status()
        self.sifs = self.generate_sifs(qclient)

        increment_status()  # increment status regardless of update
        if update:
            self.update_blanks_in_qiita(qclient)

        increment_status()  # status encompasses multiple operations
        if update:
            self.update_prep_templates(qclient)

        # before we load preps into Qiita we need to copy the fastq
        # files n times for n preps and correct the file-paths each
        # prep is pointing to.
        self.load_preps_into_qiita(qclient)

        increment_status()
        self.generate_commands()

        increment_status()
        if update:
            self.execute_commands()
