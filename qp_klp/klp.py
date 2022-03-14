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


CONFIG_FP = environ["QP_KLP_CONFIG_FP"]


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

    qclient.update_job_step(job_id, "Step 1 of 6: Setting up pipeline")

    if {'body', 'content_type', 'filename'} == set(sample_sheet):
        # Pipeline now takes the path to a sample-sheet as a parameter.
        # We must now save the sample-sheet to disk before creating a
        # Pipeline object.
        outpath = partial(join, out_dir)
        final_results_path = outpath('final_results')
        makedirs(final_results_path, exist_ok=True)
        # replace any whitespace in the filename with underscores
        sample_sheet_path = outpath(sample_sheet['filename']).replace(' ',
                                                                      '_')

        # save raw data to file
        with open(sample_sheet_path, 'w') as f:
            f.write(sample_sheet)

        # open new file as a KLSampleSheet
        # use KLSampleSheet functionality to add/overwrite lane number.
        sheet = KLSampleSheet(sample_sheet_path)
        for sample in sheet['body']:
            sample['Lane'] = '%d' % lane_number
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
                qclient.update_job_step(job_id, msg)
                raise ValueError(msg)
            else:
                raise e

        # if an Error has not been raised, assume there are no errors
        # in the sample-sheet. A successfully-created Pipeline object can
        # still contain a list of warnings about sample-sheet validation.
        # If any warnings are present, they should be reported to the user.
        if pipeline.warnings:
            msg = '\n'.join(pipeline.warnings)
            qclient.update_job_step(job_id, 'Sample-sheet has been flagged '
                                            'with the following warnings: '
                                            f'{msg}')

        sifs = pipeline.generate_sample_information_files()

        # get a list of unique sample ids from the sample-sheet.
        # comparing them against the list of samples found in the results
        # of each stage of the pipeline will let us know when and where
        # samples failed to process.
        samples = pipeline.get_sample_ids()
        # remove the BLANKs. We won't need them for our purposes.
        samples = [x for x in samples if 'BLANK' not in x]

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

        qclient.update_job_step(job_id,
                                "Step 2 of 6: Converting BCL to fastq")

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
            convert_job.run()
            failed_samples = convert_job.audit(samples)
            # write list of failed samples out to file and update after
            # each Job completes.
            with open(join(out_dir, 'failed_samples.txt'), 'a') as f:
                f.write(f"ConvertJob: {failed_samples}\n")

        qclient.update_job_step(job_id,
                                "Step 3 of 6: Adaptor & Host [optional] "
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
            qc_job.run()
            failed_samples = qc_job.audit(samples)
            with open(join(out_dir, 'failed_samples.txt'), 'a') as f:
                f.write(f"QCJob: {failed_samples}\n")

        qclient.update_job_step(job_id, "Step 4 of 6: Generating FastQC & "
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
            fastqc_job.run()
            failed_samples = fastqc_job.audit(samples)
            with open(join(out_dir, 'failed_samples.txt'), 'a') as f:
                f.write(f"FastQCJob: {failed_samples}\n")

        project_list = fastqc_job.project_names

        qclient.update_job_step(job_id, "Step 5 of 6: Generating Prep "
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
            gpf_job.run()

        qclient.update_job_step(job_id, "Step 6 of 6: Copying results to "
                                        "archive")

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
                touched_studies.append((qiita_id, project))
            else:
                cmds.append(f'cd {out_dir}; mv '
                            f'QCJob/{project}/trimmed_sequences/* '
                            f'{upload_dir}')
                touched_studies.append((qiita_id, project))

            for csv_file in csv_fps:
                if project in csv_file:
                    cmds.append(f'cd {out_dir}; mv {csv_file} {upload_dir}')
                    touched_studies.append((qiita_id, project))
                    break

        # create a set of unique study-ids that were touched by the Pipeline
        # and return this information to the user.
        touched_studies = sorted(list(set(touched_studies)))
        with open(join(out_dir, 'touched_studies.tsv'), 'a') as f:
            f.write('Project\tQiita Study ID\tQiita URL\n')
            for qiita_id, project in touched_studies:
                f.write((f'{project}\t{qiita_id}\thttps://'
                         f'{qclient.server_url}/study/description/'
                         f'{qiita_id}\n'))

        # copy all tgz files, including sample-files.tgz, to final_results.
        cmds.append(f'cd {out_dir}; mv *.tgz final_results')
        cmds.append(f'cd {out_dir}; mv FastQCJob/multiqc final_results')

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

    qclient.update_job_step(job_id, "Main Pipeline Finished, processing "
                                    "results")

    return success, ainfo, msg
