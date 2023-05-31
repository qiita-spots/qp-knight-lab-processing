from os import walk
from os.path import exists, join
from sequence_processing_pipeline.PipelineError import PipelineError
import pandas as pd
from qp_klp.Step import FailedSamplesRecord, Step


class Metagenomic(Step):
    def __init__(self, pipeline, master_qiita_job_id,
                 status_update_callback=None, lane_number=None):
        super().__init__(pipeline,
                         master_qiita_job_id,
                         status_update_callback,
                         lane_number)

        if pipeline.pipeline_type not in Step.META_TYPES:
            raise ValueError("Cannot instantiate Metagenomic object from "
                             f"pipeline of type '{pipeline.pipeline_type}'")

        # Note: FailedSamplesRecord is not used when processing amplicon as the
        # samples are processed as a single fastq file and hence that info
        # is not available.
        self.fsr = FailedSamplesRecord(self.pipeline.output_path,
                                       pipeline.sample_sheet.samples)

    def convert_bcl_to_fastq(self):
        # The 'bcl-convert' key is a convention hard-coded into mg-scripts and
        # qp-klp projects. Currently meta*omic jobs use bcl-convert for its
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
        job = super()._generate_reports()
        self.fsr.write(job.audit(self.pipeline.get_sample_ids()), 'FastQCJob')

        self.project_names = job.project_names

    def generate_prep_file(self):
        config = self.pipeline.configuration['seqpro']

        if self.project_names is None:
            raise ValueError("reports not yet generated")

        job = super()._generate_prep_file(config,
                                          self.pipeline.sample_sheet.path,
                                          config['seqpro_path'],
                                          self.project_names)

        self.prep_file_paths = job.prep_file_paths

    def generate_touched_studies(self, qclient):
        results = {}

        for study_id, pf_paths in self.prep_file_paths.items():
            for pf_path in pf_paths:
                # record the data-type as either metagenomic or
                # metatranscriptomic, according to what's stored in the
                # pipeline.
                results[pf_path] = self.pipeline.pipeline_type

    def generate_commands(self, qclient):
        super()._generate_commands()

        out_dir = self.pipeline.output_path

        self.cmds.append(f'cd {self.pipeline.output_path}; '
                         'tar zcvf logs-QCJob.tgz QCJob/logs')

        self.cmds.append(f'cd {self.pipeline.output_path}; '
                         'tar zcvf reports-ConvertJob.tgz ConvertJob/Reports '
                         'ConvertJob/Logs')

        self.write_commands_to_output_path()

        data = []
        for project, _, qiita_id in self.special_map:
            if self.pipeline.pipeline_type in Step.META_TYPES:
                self.cmds.append(f'cd {out_dir}; tar zcvf reports-QCJob.tgz '
                                 f'QCJob/{project}/fastp_reports_dir')

            # AFAIK there can only be 1 prep per project
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
            reply = qclient.post('/qiita_db/artifact/', data=data)

            data.append({'Project': project, 'Qiita Study ID': qiita_id,
                         'Qiita Prep ID': prep_id, 'Qiita URL': surl,
                         'Prep URL': prep_url, 'Linking Job': reply})

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
