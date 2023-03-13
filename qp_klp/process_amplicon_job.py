from metapool.prep import remove_qiita_id
from os import listdir, makedirs
from os import walk
from os.path import exists, join, isfile, basename
from qiita_client import ArtifactInfo
from qp_klp.klp_util import map_sample_names_to_tube_ids
from random import choices
from sequence_processing_pipeline.ConvertJob import ConvertJob
from sequence_processing_pipeline.FastQCJob import FastQCJob
from sequence_processing_pipeline.GenPrepFileJob import GenPrepFileJob
from sequence_processing_pipeline.Pipeline import Pipeline
from sequence_processing_pipeline.PipelineError import PipelineError
from subprocess import Popen, PIPE
import pandas as pd
import shutil


def process_amplicon(mapping_file_path, qclient, run_identifier, out_dir,
                     job_id, skip_exec, job_pool_size, final_results_path,
                     config_fp, status_line):

    # Note: FailedSamplesRecord is not used by process_amplicon as the
    # samples are processed as a single fastq file and hence that info
    # is not available.

    # Create a Pipeline object
    try:
        pipeline = Pipeline(config_fp, run_identifier, None,
                            mapping_file_path, out_dir, job_id)
    except PipelineError as e:
        # Pipeline is the object that finds the input fp, based on
        # a search directory set in configuration.json and a run_id.
        if str(e).endswith("could not be found"):
            msg = f"A path for {run_identifier} could not be found."
            raise PipelineError(msg)
        elif str(e).startswith("Sample-sheet has the following errors:"):
            status_line.update_current_message(str(e))
            raise ValueError(str(e))
        else:
            raise e

    # perform sample-id validation against Qiita
    projects = pipeline.get_project_info()

    errors = []
    sn_tid_map_by_project = {}

    for project in projects:
        project_name_with_qid = project['project_name']

        # remove the qiita-id prepending the project_name
        project_name = remove_qiita_id(project['project_name'])
        qiita_id = project['qiita_id']

        if qiita_id == project_name:
            raise PipelineError("Values in the project_name column must "
                                "be appended with a Qiita ID.")

        # assume the BLANKS in the mapping-file are not prepended w/qiita-id
        # or some other value. Confirmed w/wet-lab.
        df = pipeline.mapping_file
        df = df[df['project_name'] == project_name_with_qid]

        mf_samples = {s for s in df['sample_name']
                      if not s.startswith('BLANK')}

        # collect needed info from Qiita here.
        url = f'/api/v1/study/{qiita_id}/samples'
        qsamples = qclient.get(url)
        qsamples = {x.replace(f'{qiita_id}.', '') for x in qsamples}
        tube_id_present = 'tube_id' in qclient.get(f'{url}/info')['categories']
        if tube_id_present:
            tids = qclient.get(f'{url}/categories=tube_id')['samples']

        # compare the list of samples from the mapping file against the
        # list of samples from Qiita.
        sample_name_diff = mf_samples - qsamples
        sn_tid_map_by_project[project_name] = None

        if sample_name_diff:
            # if tube_id is defined in the Qiita study, then any sample_names
            # missing from the mapping-file may simply have a leading zero
            # present.
            if tube_id_present:
                # strip any leading zeroes from the sample-ids. Note that
                # if a sample-id has more than one leading zero, all of
                # them will be removed.
                mf_samples = {x.lstrip('0') for x in mf_samples}
                # once any leading zeros have been removed, recalculate
                # sample_name_diff before continuing processing.
                sample_name_diff = mf_samples - qsamples

                # before we report as an error, check tube_id.

                # generate a map of sample_names to tube_ids for
                # GenPrepFileJob.
                sn_tid_map_by_project[project_name] = {
                    y[0]: x.replace(f'{qiita_id}.', '') for x, y in
                    tids.items()}
                tids = set(sn_tid_map_by_project[project_name].keys())
                tube_id_diff = mf_samples - tids
                if not tube_id_diff:
                    continue
                len_tube_id_overlap = len(tube_id_diff)
                tids_example = ', '.join(choices(list(tids), k=5))
                error_tube_id = (
                    f'tube_id in Qiita but {len_tube_id_overlap} missing '
                    f'samples. Some samples from tube_id: {tids_example}.')
            else:
                error_tube_id = 'No tube_id column in Qiita.'

            len_overlap = len(sample_name_diff)
            # selecting at random k=5 samples to minimize space in display
            samples_example = ', '.join(choices(list(qsamples), k=5))
            # selecting the up to 4 first samples to minimize space in
            # display
            missing = ', '.join(sorted(sample_name_diff)[:4])
            errors.append(
                f'{project_name} has {len_overlap} missing samples (i.e. '
                f'{missing}). Some samples from Qiita: {samples_example}. '
                f'{error_tube_id}')

    if errors:
        raise PipelineError('\n'.join(errors))

    sifs = pipeline.generate_sample_information_files()

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

    status_line.update_current_message("Step 2 of 5: Converting BCL to fastq")

    config = pipeline.configuration['bcl2fastq']
    convert_job = ConvertJob(pipeline.run_dir,
                             pipeline.output_path,
                             # note that pipeline.sample_sheet in this case
                             # is the dummy sample-sheet created by Pipeline.
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
        convert_job.run(callback=status_line.update_job_step)

    status_line.update_current_message("Step 3 of 5: Generating FastQC & "
                                       "MultiQC reports")

    config = pipeline.configuration['fastqc']

    raw_fastq_files_path = join(pipeline.output_path, 'ConvertJob')

    # Since amplicon jobs do not require human-filtering or adapter-trimming,
    # QCJob is not used. Simulate QCJob's output directory for use as input
    # into FastQCJob.
    projects = pipeline.get_project_info()
    projects = [x['project_name'] for x in projects]

    for project_name in projects:
        # copy the files from ConvertJob output to faked QCJob output folder.
        # $WKDIR/$RUN_ID/QCJob/$PROJ_NAME/amplicon
        faked_output_folder = join(pipeline.output_path, 'QCJob',
                                   project_name, 'amplicon')
        makedirs(faked_output_folder)

        # get list of all raw output files to be copied.
        job_output = [join(raw_fastq_files_path, x) for x in
                      listdir(raw_fastq_files_path)]
        job_output = [x for x in job_output if isfile(x)]
        job_output = [x for x in job_output if x.endswith('fastq.gz')]
        # Undetermined files are very small and should be filtered from
        # results.
        job_output = [x for x in job_output if
                      not basename(x).startswith('Undetermined')]

        # copy the file
        for fastq_file in job_output:
            new_path = join(faked_output_folder, basename(fastq_file))
            shutil.copyfile(fastq_file, new_path)

        # FastQC expects the ConvertJob output to also be organized by
        # project. Since this would entail running the same ConvertJob
        # multiple times on the same input with just a name-change in
        # the dummy sample-sheet, we instead perform ConvertJob once
        # and copy the output from ConvertJob's output folder into
        # ConvertJob's output folder/project1, project2 ... projectN.
        faked_output_folder = join(raw_fastq_files_path, project_name)
        makedirs(faked_output_folder)
        job_output = [join(raw_fastq_files_path, x) for x in
                      listdir(raw_fastq_files_path)]
        job_output = [x for x in job_output if isfile(x)]
        job_output = [x for x in job_output if x.endswith('fastq.gz')]
        job_output = [x for x in job_output if
                      not basename(x).startswith('Undetermined')]

        for raw_fastq_file in job_output:
            new_path = join(faked_output_folder, basename(raw_fastq_file))
            shutil.copyfile(raw_fastq_file, new_path)

    # Run FastQCJob on faked directories
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
        fastqc_job.run(callback=status_line.update_job_step)

    status_line.update_current_message("Step 4 of 5: Generating Prep "
                                       "Information Files")

    config = pipeline.configuration['seqpro']
    # until seqpro_mf is merged into seqpro, modify seqpro_path to
    # reference the former rather than the latter.
    seqpro_path = config['seqpro_path'].replace('seqpro', 'seqpro_mf')
    project_list = fastqc_job.project_names
    gpf_job = GenPrepFileJob(
        pipeline.run_dir,
        raw_fastq_files_path,
        processed_fastq_files_path,
        pipeline.output_path,
        mapping_file_path,
        seqpro_path,
        project_list,
        config['modules_to_load'],
        job_id,
        is_amplicon=True)

    if not skip_exec:
        gpf_job.run(callback=status_line.update_job_step)

    prep_file_paths = gpf_job._get_prep_file_paths
    map_sample_names_to_tube_ids(prep_file_paths, sn_tid_map_by_project)

    status_line.update_current_message("Step 5 of 5: Copying results to "
                                       "archive")

    cmds = [f'cd {out_dir}; tar zcvf logs-ConvertJob.tgz ConvertJob/logs',
            f'cd {out_dir}; tar zcvf reports-ConvertJob.tgz '
            'ConvertJob/Reports',
            f'cd {out_dir}; tar zcvf logs-FastQCJob.tgz '
            'FastQCJob/logs',
            f'cd {out_dir}; tar zcvf reports-FastQCJob.tgz '
            'FastQCJob/fastqc',
            f'cd {out_dir}; tar zcvf logs-GenPrepFileJob.tgz '
            'GenPrepFileJob/logs',
            f'cd {out_dir}; tar zcvf prep-files.tgz '
            'GenPrepFileJob/PrepFiles'
            ]

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
        # sif filenames are of the form:
        blanks_file = f'{run_identifier}_{project}_blanks.tsv'
        if sifs and [x for x in sifs if blanks_file in x]:
            # move uncompressed sifs to upload_dir.
            cmds.append(f'cd {out_dir}; mv {blanks_file} {upload_dir}')

        # record that something is being moved into a Qiita Study.
        # this will allow us to notify the user which Studies to
        # review upon completion.
        touched_studies.append((qiita_id, project))

        # Note that even though QCJob was invoked, the output directories
        # were faked and hence existing tested code is reused here.
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

        for csv_file in csv_fps:
            if project in csv_file:
                cmds.append(f'cd {out_dir}; mv {csv_file} {upload_dir}')
                break

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

    for cmd in cmds:
        p = Popen(cmd, universal_newlines=True, shell=True,
                  stdout=PIPE, stderr=PIPE)
        std_out, std_err = p.communicate()
        return_code = p.returncode

        if return_code != 0 and not skip_exec:
            # during testing, ignore processes that fail and continue
            # to test other commands.
            raise PipelineError(f"'{cmd}' returned {return_code}")

    return [ArtifactInfo('output', 'job-output-folder',
                         [(f'{final_results_path}/', 'directory')])]
