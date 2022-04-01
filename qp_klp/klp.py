# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from functools import partial
from inspect import stack
from os import environ, walk
from os import makedirs
from os.path import basename, join, exists
from qiita_client import ArtifactInfo
from sequence_processing_pipeline.ConvertJob import ConvertJob
from sequence_processing_pipeline.FastQCJob import FastQCJob
from sequence_processing_pipeline.GenPrepFileJob import GenPrepFileJob
from sequence_processing_pipeline.Pipeline import Pipeline
from sequence_processing_pipeline.PipelineError import PipelineError
from sequence_processing_pipeline.QCJob import QCJob
from subprocess import Popen, PIPE
from metapool import KLSampleSheet
from metapool.sample_sheet import sample_sheet_to_dataframe
from metapool.prep import remove_qiita_id
from random import choices
import pandas as pd


CONFIG_FP = environ["QP_KLP_CONFIG_FP"]


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


def sequence_processing_pipeline(qclient, job_id, parameters, out_dir):
    """Sequence Processing Pipeline command

    Parameters
    ----------
    qclient : tgp.qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    parameters : dict
        The parameter values for this job
    out_dir : str
        The path to the job's output directory

    Returns
    -------
    bool, list, str
        The results of the job
    """
    run_identifier = parameters.pop('run_identifier')
    sample_sheet = parameters.pop('sample_sheet')
    lane_number = parameters.pop('lane_number')
    job_pool_size = 30

    # checking if this is running as part of the unittest
    # https://stackoverflow.com/a/25025987
    skip_exec = True if [x for x in stack() if
                         'unittest' in x.filename] else False

    success = True
    ainfo = None
    msg = None

    # maintains state of the current message, minus additional updates from
    # a callback function. E.g. "Step 1 of 6: Setting up pipeline"

    def _update_job_step(id, status):
        # internal function implements a callback function for Pipeline.run().
        # :param id: PBS/Torque/or some other informative and current job id.
        # :param status: status message
        qclient.update_job_step(job_id,
                                _update_job_step.msg + f" ({id}: {status})")

    # initialize static variable to maintain current message
    _update_job_step.msg = ""

    def _update_current_message(msg):
        # internal function that sets current_message to the new value before
        # updating the job step in the UI.
        _update_job_step.msg = msg
        qclient.update_job_step(job_id, msg)

    _update_current_message("Step 1 of 6: Setting up pipeline")

    if {'body', 'content_type', 'filename'} == set(sample_sheet):
        # Pipeline now takes the path to a sample-sheet as a parameter.
        # We must now save the sample-sheet to disk before creating a
        # Pipeline object.
        outpath = partial(join, out_dir)
        final_results_path = outpath('final_results')
        makedirs(final_results_path, exist_ok=True)
        # replace any whitespace in the filename with underscores
        sample_sheet_path = outpath(sample_sheet['filename'].replace(' ',
                                                                     '_'))
        # save raw data to file
        with open(sample_sheet_path, 'w') as f:
            f.write(sample_sheet['body'])

        # open new file as a KLSampleSheet
        # use KLSampleSheet functionality to add/overwrite lane number.
        sheet = KLSampleSheet(sample_sheet_path)
        for sample in sheet:
            sample['Lane'] = f'{lane_number}'

        sheet_df = sample_sheet_to_dataframe(sheet)
        errors = []
        for project, _df in sheet_df.groupby('sample_project'):
            project_name = remove_qiita_id(project)
            qiita_id = project.replace(f'{project_name}_', '')
            qurl = f'/api/v1/study/{qiita_id}/samples'
            sheet_samples = {
                s for s in _df['sample_name'] if not s.startswith('BLANK')}
            qsamples = {
                s.replace(f'{qiita_id}.', '') for s in qclient.get(qurl)}
            sample_name_diff = sheet_samples - qsamples
            if sample_name_diff:
                # before we report as an error, check tube_id
                error_tube_id = 'No tube_id column in Qiita.'
                if 'tube_id' in qclient.get(f'{qurl}/info')['categories']:
                    tids = qclient.get(f'{qurl}/categories=tube_id')['samples']
                    tids = {tid[0] for _, tid in tids.items()}
                    tube_id_diff = sheet_samples - tids
                    if not tube_id_diff:
                        continue
                    len_tube_id_overlap = len(tube_id_diff)
                    tids_example = ', '.join(choices(list(tids), k=5))
                    error_tube_id = (
                        f'tube_id in Qiita but {len_tube_id_overlap} missing '
                        f'samples. Some samples from tube_id: {tids_example}.')

                len_overlap = len(sample_name_diff)
                # selecting at random k=5 samples to minimize space in display
                samples_example = ', '.join(choices(list(qsamples), k=5))
                # selecting the up to 4 first samples to minimize space in
                # display
                missing = ', '.join(sorted(sample_name_diff)[:4])
                errors.append(
                    f'{project} has {len_overlap} missing samples (i.e. '
                    f'{missing}). Some samples from Qiita: {samples_example}. '
                    f'{error_tube_id}')

        if errors:
            return False, None, '\n'.join(errors)

        with open(sample_sheet_path, 'w') as f:
            sheet.write(f)

        # Create a Pipeline object
        try:
            pipeline = Pipeline(CONFIG_FP, run_identifier, sample_sheet_path,
                                out_dir, job_id)
        except PipelineError as e:
            # Pipeline is the object that finds the input fp, based on
            # a search directory set in configuration.json and a run_id.
            if str(e).endswith("could not be found"):
                msg = f"A path for {run_identifier} could not be found."
                return False, None, msg
            elif str(e).startswith("Sample-sheet has the following errors:"):
                msg = str(e)
                _update_current_message(msg)
                raise ValueError(msg)
            else:
                raise e

        # if an Error has not been raised, assume there are no errors
        # in the sample-sheet. A successfully-created Pipeline object can
        # still contain a list of warnings about sample-sheet validation.
        # If any warnings are present, they should be reported to the user.
        if pipeline.warnings:
            msg = '\n'.join(pipeline.warnings)
            _update_current_message('Sample-sheet has been flagged with the '
                                    'following warnings: {msg}')

        sifs = pipeline.generate_sample_information_files()

        # get a list of unique sample ids from the sample-sheet.
        # comparing them against the list of samples found in the results
        # of each stage of the pipeline will let us know when and where
        # samples failed to process.
        samples = pipeline.get_sample_ids()

        fsr = FailedSamplesRecord(out_dir, pipeline.sample_sheet.samples)

        # find the uploads directory all trimmed files will need to be
        # moved to.
        results = qclient.get("/qiita_db/artifacts/types/")

        # trimmed files are stored by qiita_id. Find the qiita_id
        # associated with each project and ensure a subdirectory exists
        # for when it comes time to move the trimmed files.
        special_map = []
        for project in pipeline.get_project_info():
            upload_path = join(results['uploads'], project['qiita_id'])
            makedirs(upload_path, exist_ok=True)
            special_map.append((project['project_name'], upload_path,
                                project['qiita_id']))

        _update_current_message("Step 2 of 6: Converting BCL to fastq")

        config = pipeline.configuration['bcl-convert']
        convert_job = ConvertJob(pipeline.run_dir,
                                 pipeline.output_path,
                                 sample_sheet_path,
                                 config['queue'],
                                 config['nodes'],
                                 config['nprocs'],
                                 config['wallclock_time_in_hours'],
                                 config['per_process_memory_limit'],
                                 config['executable_path'],
                                 config['modules_to_load'],
                                 job_id)

        # if skip_execution is True, then each Pipeline object will be
        # initialized, their assertions tested, and an ainfo will be
        # returned to the caller. However the Jobs will not actually
        # be executed. This is useful for testing.
        if not skip_exec:
            convert_job.run(callback=_update_job_step)
            fsr.write(convert_job.audit(samples), 'ConvertJob')

        _update_current_message("Step 3 of 6: Adaptor & Host [optional] "
                                "trimming")

        raw_fastq_files_path = join(pipeline.output_path, 'ConvertJob')

        config = pipeline.configuration['qc']
        qc_job = QCJob(raw_fastq_files_path,
                       pipeline.output_path,
                       sample_sheet_path,
                       config['mmi_db'],
                       config['queue'],
                       config['nodes'],
                       config['nprocs'],
                       config['wallclock_time_in_hours'],
                       config['job_total_memory_limit'],
                       config['fastp_executable_path'],
                       config['minimap2_executable_path'],
                       config['samtools_executable_path'],
                       config['modules_to_load'],
                       job_id,
                       job_pool_size,
                       config['job_max_array_length'])

        if not skip_exec:
            qc_job.run(callback=_update_job_step)
            fsr.write(qc_job.audit(samples), 'QCJob')

        _update_current_message("Step 4 of 6: Generating FastQC & "
                                "MultiQC reports")

        config = pipeline.configuration['fastqc']

        raw_fastq_files_path = join(pipeline.output_path, 'ConvertJob')
        processed_fastq_files_path = join(pipeline.output_path, 'QCJob')

        fastqc_job = FastQCJob(pipeline.run_dir,
                               pipeline.output_path,
                               raw_fastq_files_path,
                               processed_fastq_files_path,
                               config['nprocs'],
                               config['nthreads'],
                               config['fastqc_executable_path'],
                               config['modules_to_load'],
                               job_id,
                               config['queue'],
                               config['nodes'],
                               config['wallclock_time_in_hours'],
                               config['job_total_memory_limit'],
                               job_pool_size,
                               config['multiqc_config_file_path'],
                               config['job_max_array_length'])

        if not skip_exec:
            fastqc_job.run(callback=_update_job_step)
            fsr.write(fastqc_job.audit(samples), 'FastQCJob')

        project_list = fastqc_job.project_names

        _update_current_message("Step 5 of 6: Generating Prep "
                                "Information Files")

        config = pipeline.configuration['seqpro']
        gpf_job = GenPrepFileJob(
            pipeline.run_dir,
            raw_fastq_files_path,
            processed_fastq_files_path,
            pipeline.output_path,
            sample_sheet_path,
            config['seqpro_path'],
            project_list,
            config['modules_to_load'],
            job_id)

        if not skip_exec:
            gpf_job.run(callback=_update_job_step)

        _update_current_message("Step 6 of 6: Copying results to archive")

        cmds = [f'cd {out_dir}; tar zcvf logs-ConvertJob.tgz ConvertJob/logs',
                f'cd {out_dir}; tar zcvf reports-ConvertJob.tgz '
                'ConvertJob/Reports ConvertJob/Logs',
                f'cd {out_dir}; tar zcvf logs-QCJob.tgz QCJob/logs',
                f'cd {out_dir}; tar zcvf logs-FastQCJob.tgz '
                'FastQCJob/logs',
                f'cd {out_dir}; tar zcvf reports-FastQCJob.tgz '
                'FastQCJob/fastqc',
                f'cd {out_dir}; tar zcvf logs-GenPrepFileJob.tgz '
                'GenPrepFileJob/logs',
                f'cd {out_dir}; tar zcvf prep-files.tgz '
                'GenPrepFileJob/PrepFiles']

        # just use the filenames for tarballing the sifs.
        # the sifs should all be stored in the {out_dir} by default.
        if sifs:
            tmp = [basename(x) for x in sifs]
            # convert sifs into a list of filenames.
            tmp = ' '.join(tmp)
            cmds.append(f'cd {out_dir}; tar zcvf sample-files.tgz {tmp}')

        csv_fps = []
        for root, dirs, files in walk(join(gpf_job.output_path, 'PrepFiles')):
            for csv_file in files:
                csv_fps.append(join(root, csv_file))

        touched_studies = []

        for project, upload_dir, qiita_id in special_map:
            if sifs and [x for x in sifs if f'{project}_blanks.tsv' in x]:
                # move uncompressed sifs to upload_dir.
                # sif filenames are of the form: '{project}_blanks.tsv'
                cmds.append(f'cd {out_dir}; mv {project}_blanks.tsv'
                            f' {upload_dir}')

            # record that something is being moved into a Qiita Study.
            # this will allow us to notify the user which Studies to
            # review upon completion.
            touched_studies.append((qiita_id, project))

            cmds.append(f'cd {out_dir}; tar zcvf reports-QCJob.tgz '
                        f'QCJob/{project}/fastp_reports_dir')

            if exists(f'{out_dir}/QCJob/{project}/filtered_sequences'):
                cmds.append(f'cd {out_dir}; mv '
                            f'QCJob/{project}/filtered_sequences/* '
                            f'{upload_dir}')
            else:
                cmds.append(f'cd {out_dir}; mv '
                            f'QCJob/{project}/trimmed_sequences/* '
                            f'{upload_dir}')

            for csv_file in csv_fps:
                if project in csv_file:
                    cmds.append(f'cd {out_dir}; mv {csv_file} {upload_dir}')
                    break

        # create a set of unique study-ids that were touched by the Pipeline
        # and return this information to the user.
        touched_studies = sorted(list(set(touched_studies)))

        data = []
        for qiita_id, project in touched_studies:
            url = f'https://{qclient._server_url}/study/description/{qiita_id}'
            data.append({'Project': project, 'Qiita Study ID': qiita_id,
                         'Qiita URL': url})
        df = pd.DataFrame(data)

        with open(join(out_dir, 'touched_studies.html'), 'w') as f:
            f.write(df.to_html(border=2, index=False, justify="left",
                               render_links=True, escape=False))

        # copy all tgz files, including sample-files.tgz, to final_results.
        cmds.append(f'cd {out_dir}; mv *.tgz final_results')
        cmds.append(f'cd {out_dir}; mv FastQCJob/multiqc final_results')

        if exists(join(out_dir, 'touched_studies.html')):
            cmds.append(f'cd {out_dir}; mv touched_studies.html final_results')

        if exists(join(out_dir, 'failed_samples.html')):
            cmds.append(f'cd {out_dir}; mv failed_samples.html final_results')

        # allow the writing of commands out to cmds.log, even if skip_exec
        # is True. This allows for unit-testing of cmds generation.
        cmd_log_fp = join(out_dir, 'cmds.log')
        with open(cmd_log_fp, 'w') as cmd_log_f:
            for cmd in cmds:
                cmd_log_f.write(f'{cmd}\n')

        # if execution was skipped, reinitialze the cmds list to empty after
        # writing to log and before actually executing commands.
        if skip_exec:
            cmds = []

        for cmd in cmds:
            p = Popen(cmd, universal_newlines=True, shell=True,
                      stdout=PIPE, stderr=PIPE)
            std_out, std_err = p.communicate()
            return_code = p.returncode
            if return_code != 0:
                raise PipelineError(f"'{cmd}' returned {return_code}")

        ainfo = [
            ArtifactInfo('output', 'job-output-folder',
                         [(f'{final_results_path}/', 'directory')])
        ]
    else:
        success = False
        msg = "This doesn't appear to be a valid sample sheet; please review."

    _update_current_message("Main Pipeline Finished, processing results")

    return success, ainfo, msg
