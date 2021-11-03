# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from functools import partial
from qiita_client import ArtifactInfo
from sequence_processing_pipeline.Pipeline import Pipeline
from os import makedirs
from os.path import join
from sequence_processing_pipeline.ConvertJob import ConvertJob
from sequence_processing_pipeline.QCJob import QCJob
from sequence_processing_pipeline.FastQCJob import FastQCJob
from sequence_processing_pipeline.GenPrepFileJob import GenPrepFileJob
from sequence_processing_pipeline.SequenceDirectory import SequenceDirectory
from sequence_processing_pipeline.PipelineError import PipelineError
from metapool import KLSampleSheet, validate_and_scrub_sample_sheet
from subprocess import Popen, PIPE
from os import environ, walk
from inspect import stack


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
        outpath = partial(join, out_dir)
        final_results_path = outpath('final_results')
        makedirs(final_results_path, exist_ok=True)

        # the user is uploading a sample-sheet to us, but we need the
        # sample-sheet as a file to pass to the Pipeline().
        sample_sheet_path = outpath(sample_sheet['filename'])
        with open(sample_sheet_path, 'w') as f:
            f.write(sample_sheet['body'])

        # validate the sample-sheet using metapool package.
        sheet = KLSampleSheet(sample_sheet_path)
        val_sheet = validate_and_scrub_sample_sheet(sheet)
        if not val_sheet:
            qclient.update_job_step(job_id,
                                    "Sample sheet failed validation.")
            raise ValueError("Sample sheet failed validiation")
        else:
            # get project names and their associated qiita ids
            bioinformatics = val_sheet.Bioinformatics
            lst = bioinformatics.to_dict('records')

        # find the uploads directory all trimmed files will need to be
        # moved to.
        results = qclient.get("/qiita_db/artifacts/types/")

        # trimmed files are stored by qiita_id. Find the qiita_id
        # associated with each project and ensure a subdirectory exists
        # for when it comes time to move the trimmed files.
        special_map = []
        for result in lst:
            project_name = result['Sample_Project']
            qiita_id = result['QiitaID']
            upload_path = join(results['uploads'], qiita_id)
            makedirs(upload_path, exist_ok=True)
            special_map.append((project_name, upload_path))

        try:
            pipeline = Pipeline(CONFIG_FP, run_identifier, out_dir, job_id)
        except PipelineError as e:
            # Pipeline is the object that finds the input fp, based on
            # a search directory set in configuration.json and a run_id.
            if str(e).endswith("could not be found"):
                msg = f"A path for {run_identifier} could not be found."
                return False, None, msg
            else:
                raise e

        makedirs(join(pipeline.run_dir, run_identifier), exist_ok=True)

        sdo = SequenceDirectory(pipeline.run_dir, sample_sheet_path)

        qclient.update_job_step(job_id,
                                "Step 2 of 6: Converting BCL to fastq")

        config = pipeline.configuration['bcl-convert']
        convert_job = ConvertJob(pipeline.run_dir,
                                 pipeline.output_path,
                                 sdo.sample_sheet_path,
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

        qclient.update_job_step(job_id,
                                "Step 3 of 6: Adaptor & Host [optional] "
                                "trimming")

        raw_fastq_files_path = join(pipeline.output_path, 'ConvertJob')

        config = pipeline.configuration['qc']
        qc_job = QCJob(raw_fastq_files_path,
                       pipeline.output_path,
                       sdo.sample_sheet_path,
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
                       job_pool_size)

        if not skip_exec:
            qc_job.run()

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
                               config['multiqc_config_file_path'])

        if not skip_exec:
            fastqc_job.run()

        project_list = fastqc_job.project_names

        qclient.update_job_step(job_id, "Step 5 of 6: Generating Prep "
                                        "Information Files")

        config = pipeline.configuration['seqpro']
        gpf_job = GenPrepFileJob(
            pipeline.run_dir,
            raw_fastq_files_path,
            processed_fastq_files_path,
            pipeline.output_path,
            sdo.sample_sheet_path,
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

        csv_fps = []
        for root, dirs, files in walk(join(gpf_job.output_path, 'PrepFiles')):
            for csv_file in files:
                csv_fps.append(join(root, csv_file))

        for project, upload_dir in special_map:
            cmds.append(f'cd {out_dir}; tar zcvf reports-QCJob.tgz '
                        'QCJob/{project}/fastp_reports_dir')
            cmds.append(f'cd {out_dir}; mv '
                        'QCJob/{project}/filtered_sequences/* '
                        '{upload_dir}')
            for csv_file in csv_fps:
                if project in csv_file:
                    cmds.append(f'cd {out_dir}; mv {csv_file} {upload_dir}')
                    break

        cmds.append(f'cd {out_dir}; mv *.tgz final_results')
        cmds.append(f'cd {out_dir}; mv FastQCJob/multiqc final_results')

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
