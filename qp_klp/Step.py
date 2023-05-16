from itertools import chain
from os.path import join
from qp_klp.klp_util import map_sample_names_to_tube_ids
from sequence_processing_pipeline.ConvertJob import ConvertJob
from sequence_processing_pipeline.FastQCJob import FastQCJob
from sequence_processing_pipeline.GenPrepFileJob import GenPrepFileJob
from sequence_processing_pipeline.PipelineError import PipelineError
from sequence_processing_pipeline.QCJob import QCJob
from subprocess import Popen, PIPE
import pandas as pd


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
        self.status_update_callback = status_update_callback
        # for now, hardcode this at the legacy value, since we've never
        # changed it.
        self.job_pool_size = 30
        self.project_names = None
        self.cmds = None
        self.cmds_log_path = None
        self.sn_tid_map_by_project = sn_tid_map_by_project
        self.prep_file_paths = None
        self.sifs = None

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

        convert_job.run(callback=self.status_update_callback.update_job_status)

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

        qc_job.run(callback=self.status_update_callback.update_job_status)

        return qc_job

    def _generate_reports(self, config, input_file_path):
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
                               self.pipeline.type == 'amplicon')

        fastqc_job.run(callback=self.status_update_callback.update_job_status)

        return fastqc_job

    def _generate_prep_file(self, config, input_file_path, seqpro_path,
                            project_names):
        is_amplicon = self.pipeline.type == 'amplicon'

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

        gpf_job.run(callback=self.status_update_callback.update_job_status)

        # concatenate the lists of paths across all study_ids into a single
        # list. Replace sample-names w/tube-ids in all relevant prep-files.
        preps = list(chain.from_iterable(gpf_job.prep_file_paths.values()))
        map_sample_names_to_tube_ids(preps, self.sn_tid_map_by_project)

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
