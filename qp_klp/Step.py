from collections import defaultdict
from itertools import chain
from json import dumps
from metapool import KLSampleSheet
from os import makedirs
from os.path import join
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
        self.failed = {x.Sample_ID: [x.Sample_Project, None] for x in samples}

    def write(self, failed_ids, job_name):
        for failed_id in failed_ids:
            # as a rule, if a failed_id were to appear in more than one
            # audit(), preserve the earliest failure, rather than the
            # latest one.
            if self.failed[failed_id][1] is None:
                self.failed[failed_id][1] = job_name

        # filter out the sample-ids w/out a failure status
        filtered_fails = {x: self.failed[x] for x in self.failed if
                          self.failed[x][1] is not None}

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
    def __init__(self, pipeline, master_qiita_job_id, sn_tid_map_by_project,
                 status_update_callback=None):

        if pipeline is None:
            raise ValueError("A pipeline object is needed to initialize Step")

        if master_qiita_job_id is None:
            raise ValueError("A Qiita job-id is needed to initialize Step")

        if sn_tid_map_by_project is None:
            raise ValueError("sn_tid_map_by_project is needed to initialize"
                             " Step")

        self.pipeline = pipeline
        self.master_qiita_job_id = master_qiita_job_id

        if status_update_callback is not None:
            self.update_callback = status_update_callback.update_job_status
        else:
            self.update_callback = None

        # for now, hardcode this at the legacy value, since we've never
        # changed it.
        self.job_pool_size = 30
        self.project_names = None
        self.cmds = None
        self.cmds_log_path = None
        self.sn_tid_map_by_project = sn_tid_map_by_project
        self.prep_file_paths = None
        self.sifs = None

    @classmethod
    def generate_pipeline(cls, pipeline_type, input_file_path, lane_number,
                          config_fp,
                          run_identifier, out_dir, job_id):
        if pipeline_type in ['metagenomic', 'metatranscriptomic']:
            cls.update_sample_sheet(input_file_path, lane_number)
            return Pipeline(config_fp, run_identifier, input_file_path, None,
                            out_dir, job_id, pipeline_type)
        elif pipeline_type == 'amplicon':
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

    @classmethod
    def generate_special_map(cls, results, projects):
        # this function should be able to be tested by passing in simulated =
        # results from qclient.

        # trimmed files are stored by qiita_id. Find the qiita_id
        # associated with each project and ensure a subdirectory exists
        # for when it comes time to move the trimmed files.
        special_map = []
        for project in projects:
            upload_path = join(results['uploads'], project['qiita_id'])
            makedirs(upload_path, exist_ok=True)
            special_map.append((project['project_name'], upload_path,
                                project['qiita_id']))

        return special_map

    @classmethod
    def update_prep_templates(cls, qclient, prep_file_paths, pipeline_type):
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
                        'data_type': None}
                if pipeline_type in ['metagenomic', 'metatranscriptomic']:
                    data['data_type'] = pipeline_type
                elif pipeline_type == 'amplicon':
                    if 'target_gene' in metadata[list(metadata.keys())[0]]:
                        tg = metadata[list(metadata.keys())[0]]['target_gene']
                        for key in {'16S', '18S', 'ITS'}:
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

    def _generate_reports(self, input_file_path):
        config = self.pipeline.configuration['fastqc']
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
                               self.pipeline.pipeline_type == 'amplicon')

        fastqc_job.run(callback=self.update_callback)

        return fastqc_job

    def _generate_prep_file(self, config, input_file_path, seqpro_path,
                            project_names):
        is_amplicon = self.pipeline.pipeline_type == 'amplicon'

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
        Step.map_sample_names_to_tube_ids(preps, self.sn_tid_map_by_project)

        return gpf_job

    def _generate_commands(self):
        cmds = ['tar zcvf logs-ConvertJob.tgz ConvertJob/logs',
                'tar zcvf logs-FastQCJob.tgz FastQCJob/logs',
                'tar zcvf reports-FastQCJob.tgz FastQCJob/fastqc',
                'tar zcvf logs-GenPrepFileJob.tgz GenPrepFileJob/logs',
                'tar zcvf prep-files.tgz GenPrepFileJob/PrepFiles']
        self.cmds = [f'cd {self.pipeline.output_path}; {x}' for x in cmds]

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

    def generate_sifs(self, from_qiita):
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

    def map_sample_names_to_tube_ids(self, prep_info_file_paths,
                                     sn_tid_map_by_proj):
        for proj in sn_tid_map_by_proj:
            if sn_tid_map_by_proj[proj] is not None:
                # this project has tube-ids registered in Qiita.
                # find the prep-file associated with this project.
                for prep_file in prep_info_file_paths:
                    # not the best check but good enough for now.
                    if proj in prep_file:
                        df = pd.read_csv(prep_file, sep='\t',
                                         dtype=str, index_col=False)
                        # save copy of sample_name column as 'old_sample_name'
                        df['old_sample_name'] = df['sample_name']
                        for i in df.index:
                            smpl_name = df.at[i, "sample_name"]
                            if not smpl_name.startswith('BLANK'):
                                # remove any leading zeroes if they exist
                                smpl_name = smpl_name.lstrip('0')
                                if smpl_name in sn_tid_map_by_proj[proj]:
                                    tube_id = sn_tid_map_by_proj[proj][
                                        smpl_name]
                                    df.at[i, "sample_name"] = tube_id
                        df.to_csv(prep_file, index=False, sep="\t")

    def update_blanks_in_qiita(self, sifs, qclient):
        for sif_path in sifs:
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
                data = {i: {c: 1 for c in categories} for i in new_blanks}

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
