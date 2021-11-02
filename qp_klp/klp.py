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
    config_fp = parameters.pop('config_filepath')
    skip_execution = parameters.pop('skip_execution')
    job_pool_size = 30

    try:
        qclient.update_job_step(job_id, "Step 1 of 6: Setting up pipeline")

        # Although KLSampleSheet thoroughly validates the sample-sheeet,
        # verify that we at least have a body to write to file, etc.
        if {'body', 'content_type', 'filename'} != set(sample_sheet):
            raise PipelineError("This doesn't appear to be a valid sample "
                                "sheet; please review.")

        # the final_results directory is the subdirectory of out_dir where
        # all final products are moved to. This is distinct from out_dir,
        # which is used as a working directory.
        outpath = partial(join, out_dir)
        final_results_path = outpath('final_results')
        makedirs(final_results_path, exist_ok=True)

        # the user is uploading a sample-sheet to us, but we need the sample-
        # sheet as a file to pass to the Pipeline().
        sample_sheet_path = outpath(sample_sheet['filename'])
        with open(sample_sheet_path, 'w') as f:
            f.write(sample_sheet['body'])

        # validate the sample-sheet using metapool package.
        sheet = KLSampleSheet(sample_sheet_path)
        val_sheet = validate_and_scrub_sample_sheet(sheet)
        if not val_sheet:
            qclient.update_job_step(job_id,
                                    "Sample sheet failed validation.")
            raise PipelineError("Sample sheet failed validiation")
        else:
            # get project names and their associated qiita ids
            bioinformatics = val_sheet.Bioinformatics
            lst = bioinformatics.to_dict('records')

        # find the uploads directory all trimmed files will need to be moved
        # to.
        results = qclient.get("/qiita_db/artifacts/types/")

        # trimmed files are stored by qiita_id. Find the qiita_id associated
        # with each project and ensure a subdirectory exists for when it comes
        # time to move the trimmed files.
        special_map = []
        for result in lst:
            project_name = result['Sample_Project']
            qiita_id = result['QiitaID']
            upload_path = join(results['uploads'], qiita_id)
            makedirs(upload_path, exist_ok=True)
            special_map.append((project_name, upload_path))

        # Create a Pipeline() object. Pipeline() will raise PipelineErrors()
        # if files and directories don't exist as they should, etc. These are
        # all caught below. SequenceDirectory() and Jobs() also raise
        # PipelineErrors().
        # Note that the location of the BCL files (input data) is discovered
        # by Pipeline(), based on the run_identifier given, and the search
        # paths stored in the file pointed to by config_fp. The input path
        # to the BCL files is not passed directly to Pipeline by the user.

        pipeline = Pipeline(config_fp, run_identifier, out_dir, job_id)

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
        if not skip_execution:
            convert_job.run()

        qclient.update_job_step(job_id,
                                "Step 3 of 6: Adaptor & Host [optional]"
                                " trimming")

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

        if not skip_execution:
            qc_job.run()

        qclient.update_job_step(job_id, "Step 4 of 6: Generating FastQC &"
                                        " MultiQC reports")

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

        if not skip_execution:
            fastqc_job.run()

        project_list = fastqc_job.project_names

        qclient.update_job_step(job_id, "Step 5 of 6: Generating Prep Info"
                                        "rmation Files")

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

        if not skip_execution:
            gpf_job.run()

        qclient.update_job_step(job_id, "Step 6 of 6: Copying results to "
                                        "archive")

        cmds = [f'cd {out_dir}; tar zcvf logs-ConvertJob.tgz ConvertJob/'
                'logs',
                f'cd {out_dir}; tar zcvf reports-ConvertJob.tgz '
                'ConvertJob/Reports ConvertJob/Logs',
                f'cd {out_dir}; tar zcvf logs-QCJob.tgz QCJob/logs',
                f'cd {out_dir}; tar zcvf logs-FastQCJob.tgz '
                'FastQCJob/logs',
                f'cd {out_dir}; tar zcvf reports-FastQCJob.tgz '
                'FastQCJob/fastqc',
                f'cd {out_dir}; tar zcvf logs-GenPrepFileJob.tgz '
                'GenPrepFileJob/logs']

        for project, upload_dir in special_map:
            cmds.append(f'cd {out_dir}; tar zcvf reports-QCJob.tgz '
                        'QCJob/{project}/fastp_reports_dir')
            cmds.append(f'cd {out_dir}; mv '
                        'QCJob/{project}/filtered_sequences/* '
                        '{upload_dir}')

        cmds.append(f'cd {out_dir}; mv *.tgz final_results')
        cmds.append(f'cd {out_dir}; mv FastQCJob/multiqc final_results')
        cmds.append(f'cd {out_dir}; mv GenPrepFileJob/PrepFiles/* {upload_dir}')

        if skip_execution:
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

        qclient.update_job_step(job_id, "Main Pipeline Finished, processing "
                                        "results")
        return True, ainfo, None

    except PipelineError as e:
        # Pipeline() will identify directories and files that don't exist.
        # Catch all PipelineErrors() here and handle them by returning
        # success=False, None (ainfo), and a string representation of the
        # error.
        return False, None, str(e)

    # sequence_processing_pipeline() should not arrive at this point.
    raise PipelineError("sequence_processing_pipeline() did not complete.")
