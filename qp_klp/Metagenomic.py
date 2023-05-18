from os import walk
from os.path import exists, join, basename
from qp_klp.klp_util import FailedSamplesRecord
from sequence_processing_pipeline.PipelineError import PipelineError
import pandas as pd
from qp_klp.Step import Step


class Metagenomic(Step):
    def __init__(self, pipeline, master_qiita_job_id, sn_tid_map_by_project,
                 status_update_callback=None):
        super().__init__(pipeline,
                         master_qiita_job_id,
                         sn_tid_map_by_project,
                         status_update_callback)

        if pipeline.pipeline_type != 'metagenomic':
            raise ValueError("Cannot instantiate Metagenomic object from "
                             f"pipeline of type '{pipeline.pipeline_type}'")

        # Note: FailedSamplesRecord is not used when processing amplicon as the
        # samples are processed as a single fastq file and hence that info
        # is not available.
        self.fsr = FailedSamplesRecord(self.pipeline.output_path,
                                       pipeline.sample_sheet.samples)
        self.project_names = None

    def convert_bcl_to_fastq(self):
        # The 'bcl-convert' key is a convention hard-coded into mg-scripts and
        # qp-klp projects. Currently metagenomic jobs use bcl-convert for its
        # improved performance over bcl2fastq. The name and path of the
        # executable, the resource requirements to instantiate a SLURM job
        # with, etc. are stored in configuration['bcl-convert''].
        config = self.pipeline.configuration['bcl-convert']
        job = super()._convert_bcl_to_fastq(config,
                                            self.pipeline.sample_sheet.path)
        self.fsr.write(job.audit(self.pipeline.get_sample_ids()), 'ConvertJob')

    def quality_control(self):
        config = self.pipeline.configuration['qc']
        job = super()._quality_control(config, self.pipeline.sample_sheet.path)
        self.fsr.write(job.audit(self.pipeline.get_sample_ids()), 'QCJob')

    def generate_reports(self):
        config = self.pipeline.configuration['fastqc']
        job = super()._generate_reports(config,
                                        self.pipeline.sample_sheet.path)
        self.fsr.write(job.audit(self.pipeline.get_sample_ids()), 'FastQCJob')

        self.project_names = job.project_names

    def generate_prep_file(self):
        config = self.pipeline.configuration['seqpro']

        if self.project_names is None:
            raise ValueError("reports not yet generated")

        job = super()._generate_prep_file(config,
                                          self.pipeline.sample_sheet,
                                          config['seqpro_path'],
                                          self.project_names)

        self.prep_file_paths = job.prep_file_paths

    def generate_commands(self, special_map, server_url,
                          touched_studies_prep_info):
        super()._generate_commands()

        out_dir = self.pipeline.output_path
        output_path = self.pipeline.output_path

        self.cmds.append(f'cd {self.pipeline.output_path}; '
                         'tar zcvf logs-QCJob.tgz QCJob/logs')

        self.cmds.append(f'cd {self.pipeline.output_path}; '
                         'tar zcvf reports-ConvertJob.tgz ConvertJob/Reports '
                         'ConvertJob/Logs')

        self.write_commands_to_output_path()

        if self.sifs:
            # just use the filenames for tarballing the sifs.
            # the sifs should all be stored in the {out_dir} by default.
            tmp = [basename(x) for x in self.sifs]
            # convert sifs into a list of filenames.
            tmp = ' '.join(tmp)
            self.cmds.append(f'cd {out_dir}; tar zcvf sample-files.tgz {tmp}')

        csv_fps = []
        for root, dirs, files in walk(join(output_path, 'PrepFiles')):
            for csv_file in files:
                csv_fps.append(join(root, csv_file))

        touched_studies = []

        for project, upload_dir, qiita_id in special_map:
            # sif filenames are of the form:
            blanks_file = f'{self.pipeline.run_id}_{project}_blanks.tsv'
            if self.sifs and [x for x in self.sifs if blanks_file in x]:
                # move uncompressed sifs to upload_dir.
                tmp = f'cd {out_dir}; mv {blanks_file} {upload_dir}'
                self.cmds.append(tmp)

            # record that something is being moved into a Qiita Study.
            # this will allow us to notify the user which Studies to
            # review upon completion.
            touched_studies.append((qiita_id, project))

            if self.pipeline.pipeline_type == 'metagenomic':
                self.cmds.append(f'cd {out_dir}; tar zcvf reports-QCJob.tgz '
                                 f'QCJob/{project}/fastp_reports_dir')

            if exists(f'{out_dir}/QCJob/{project}/filtered_sequences'):
                self.cmds.append(f'cd {out_dir}; mv '
                                 f'QCJob/{project}/filtered_sequences/* '
                                 f'{upload_dir}')
            elif exists(f'{out_dir}/QCJob/{project}/trimmed_sequences'):
                self.cmds.append(f'cd {out_dir}; mv '
                                 f'QCJob/{project}/trimmed_sequences/* '
                                 f'{upload_dir}')
            elif exists(f'{out_dir}/QCJob/{project}/amplicon'):
                self.cmds.append(f'cd {out_dir}; mv '
                                 f'QCJob/{project}/amplicon/* '
                                 f'{upload_dir}')
            else:
                raise PipelineError("QCJob output not in expected location")

            for csv_file in csv_fps:
                if project in csv_file:
                    tmp = f'cd {out_dir}; mv {csv_file} {upload_dir}'
                    self.cmds.append(tmp)
                    break

        # create a set of unique study-ids that were touched by the Pipeline
        # and return this information to the user.
        touched_studies = sorted(list(set(touched_studies)))

        data = []
        for qiita_id, project in touched_studies:
            for prep_id in touched_studies_prep_info[qiita_id]:
                study_url = f'{server_url}/study/description/{qiita_id}'
                prep_url = (f'{server_url}/study/description/'
                            f'{qiita_id}?prep_id={prep_id}')
                data.append({'Project': project, 'Qiita Study ID': qiita_id,
                             'Qiita Prep ID': prep_id, 'Qiita URL': study_url,
                             'Prep URL': prep_url})

        df = pd.DataFrame(data)

        with open(join(out_dir, 'touched_studies.html'), 'w') as f:
            f.write(df.to_html(border=2, index=False, justify="left",
                               render_links=True, escape=False))

        # copy all tgz files, including sample-files.tgz, to final_results.
        self.cmds.append(f'cd {out_dir}; mv *.tgz final_results')
        self.cmds.append(f'cd {out_dir}; mv FastQCJob/multiqc final_results')

        if exists(join(out_dir, 'touched_studies.html')):
            tmp = f'cd {out_dir}; mv touched_studies.html final_results'
            self.cmds.append(tmp)

        if exists(join(out_dir, 'failed_samples.html')):
            tmp = f'cd {out_dir}; mv failed_samples.html final_results'
            self.cmds.append(tmp)
