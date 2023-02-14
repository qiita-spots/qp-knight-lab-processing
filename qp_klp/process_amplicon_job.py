from os import listdir, makedirs
from os.path import exists, join, isfile, basename
from qiita_client import ArtifactInfo
from sequence_processing_pipeline.ConvertJob import ConvertJob
from sequence_processing_pipeline.FastQCJob import FastQCJob
# from sequence_processing_pipeline.GenPrepFileJob import GenPrepFileJob
from sequence_processing_pipeline.Pipeline import Pipeline
from sequence_processing_pipeline.PipelineError import PipelineError
import shutil
import pandas as pd
from subprocess import Popen, PIPE


def process_amplicon(mapping_file_path, qclient,
                     run_identifier, out_dir, job_id,
                     _update_current_message, skip_exec,
                     _update_job_step, job_pool_size,
                     final_results_path, success, msg, config_fp):

    # TODO: get sample_ids from mapping file for each project and confirm
    #  the values are present in Qiita.

    # Create a Pipeline object
    try:
        pipeline = Pipeline(config_fp, run_identifier, None,
                            mapping_file_path, out_dir, job_id)
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

    # TODO: Confirm that we need SIFs for amplicon
    # sifs = pipeline.generate_sample_information_files()

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
                             pipeline.sample_sheet,
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

    _update_current_message("Step 3 of 6: Generating FastQC & "
                            "MultiQC reports")

    config = pipeline.configuration['fastqc']

    raw_fastq_files_path = join(pipeline.output_path, 'ConvertJob')

    # Fake QCJob output directory
    my_projects = pipeline.get_project_info()
    my_projects = [x['project_name'] for x in my_projects]

    for some_project in my_projects:
        # copy the files from ConvertJob output to QCJob 'processed' output
        some_path = join(pipeline.output_path, 'QCJob', some_project,
                         'amplicon')
        # making some_path: $WKDIR/$RUN_ID/QCJob/$PROJ_NAME/amplicon
        makedirs(some_path)

        # build file-paths
        job_output = [join(raw_fastq_files_path, x) for x in
                      listdir(raw_fastq_files_path)]
        # filter for only the files
        job_output = [x for x in job_output if isfile(x)]
        # filter for only fastq files
        job_output = [x for x in job_output if x.endswith('fastq.gz')]

        for a_file in job_output:
            # assume Undetermined files are very small and should be
            # filtered out.
            file_name = basename(a_file)
            if not file_name.startswith('Undetermined'):
                new_path = join(some_path, file_name)
                shutil.copyfile(a_file, new_path)

        # Now do the same for ConvertJob - make copies of the fastq files in
        # project directories for FastQC.
        some_path = join(raw_fastq_files_path, some_project)
        makedirs(some_path)
        job_output = [join(raw_fastq_files_path, x) for x in
                      listdir(raw_fastq_files_path)]
        job_output = [x for x in job_output if isfile(x)]
        job_output = [x for x in job_output if x.endswith('fastq.gz')]

        for a_file in job_output:
            file_name = basename(a_file)
            if not file_name.startswith('Undetermined'):
                new_path = join(some_path, file_name)
                shutil.copyfile(a_file, new_path)

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
                           config['job_max_array_length'],
                           True)

    if not skip_exec:
        fastqc_job.run(callback=_update_job_step)

    _update_current_message("Step 5 of 6: Generating Prep "
                            "Information Files")

    """
    project_list = fastqc_job.project_names
    config = pipeline.configuration['seqpro']
    gpf_job = GenPrepFileJob(
        pipeline.run_dir,
        raw_fastq_files_path,
        processed_fastq_files_path,
        pipeline.output_path,
        mapping_file_path,
        config['seqpro_path'],
        project_list,
        config['modules_to_load'],
        job_id)

    if not skip_exec:
        gpf_job.run(callback=_update_job_step)
    """

    _update_current_message("Step 6 of 6: Copying results to archive")

    cmds = [f'cd {out_dir}; tar zcvf logs-ConvertJob.tgz ConvertJob/logs',
            f'cd {out_dir}; tar zcvf reports-ConvertJob.tgz '
            'ConvertJob/Reports ConvertJob/Logs',
            f'cd {out_dir}; tar zcvf logs-FastQCJob.tgz '
            'FastQCJob/logs',
            f'cd {out_dir}; tar zcvf reports-FastQCJob.tgz '
            'FastQCJob/fastqc',
            # f 'cd {out_dir}; tar zcvf logs-GenPrepFileJob.tgz '
            # 'GenPrepFileJob/logs',
            # f'cd {out_dir}; tar zcvf prep-files.tgz '
            # 'GenPrepFileJob/PrepFiles'
            ]

    '''
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
    '''

    touched_studies = []

    for project, upload_dir, qiita_id in special_map:
        # sif filenames are of the form:
        # blanks_file = f'{run_identifier}_{project}_blanks.tsv'
        # if sifs and [x for x in sifs if blanks_file in x]:
        #     # move uncompressed sifs to upload_dir.
        #     cmds.append(f'cd {out_dir}; mv {blanks_file} {upload_dir}')

        # record that something is being moved into a Qiita Study.
        # this will allow us to notify the user which Studies to
        # review upon completion.
        touched_studies.append((qiita_id, project))

        if exists(f'{out_dir}/QCJob/{project}/filtered_sequences'):
            cmds.append(f'cd {out_dir}; mv '
                        f'QCJob/{project}/filtered_sequences/* '
                        f'{upload_dir}')
        elif exists(f'{out_dir}/QCJob/{project}/trimmed_sequences'):
            cmds.append(f'cd {out_dir}; mv '
                        f'QCJob/{project}/trimmed_sequences/* '
                        f'{upload_dir}')
        elif exists(f'{out_dir}/QCJob/{project}/amplicon'):
            cmds.append(f'cd {out_dir}; mv '
                        f'QCJob/{project}/amplicon/* '
                        f'{upload_dir}')
        else:
            raise PipelineError("QCJob output not in expected location")

        '''
        # requires GenPrepFileJob
        for csv_file in csv_fps:
            if project in csv_file:
                cmds.append(f'cd {out_dir}; mv {csv_file} {upload_dir}')
                break
        '''

    # create a set of unique study-ids that were touched by the Pipeline
    # and return this information to the user.
    touched_studies = sorted(list(set(touched_studies)))

    data = []
    for qiita_id, project in touched_studies:
        url = f'{qclient._server_url}/study/description/{qiita_id}'
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

    return success, ainfo, msg
