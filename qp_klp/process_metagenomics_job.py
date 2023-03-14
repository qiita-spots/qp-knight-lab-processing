from os import walk
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
from qp_klp.klp_util import map_sample_names_to_tube_ids, FailedSamplesRecord
from json import dumps
from itertools import chain


def process_metagenomics(sample_sheet_path, lane_number, qclient,
                         run_identifier, out_dir, job_id, skip_exec,
                         job_pool_size, final_results_path, config_fp,
                         status_line):

    # open new file as a KLSampleSheet
    # use KLSampleSheet functionality to add/overwrite lane number.
    sheet = KLSampleSheet(sample_sheet_path)
    for sample in sheet:
        sample['Lane'] = f'{lane_number}'

    assay_type = sheet.Header.Assay

    sheet_df = sample_sheet_to_dataframe(sheet)
    errors = []
    sn_tid_map_by_project = {}
    for project, _df in sheet_df.groupby('sample_project'):
        project_name = remove_qiita_id(project)
        qiita_id = project.replace(f'{project_name}_', '')
        qurl = f'/api/v1/study/{qiita_id}/samples'

        sheet_samples = {
            s for s in _df['sample_name'] if not s.startswith('BLANK')}
        qsamples = {
            s.replace(f'{qiita_id}.', '') for s in qclient.get(qurl)}
        sample_name_diff = sheet_samples - qsamples
        sn_tid_map_by_project[project_name] = None

        # check that tube_id is defined in the Qiita study. If so,
        # then any sample_names missing from the sample-sheet may simply
        # have a leading zero present.
        tube_id_present = 'tube_id' in qclient.get(f'{qurl}/info')[
            'categories']

        if sample_name_diff:
            # before we report as an error, check tube_id
            error_tube_id = 'No tube_id column in Qiita.'

            if tube_id_present:
                # strip any leading zeroes from the sample-ids. Note that
                # if a sample-id has more than one leading zero, all of
                # them will be removed.
                sheet_samples = {x.lstrip('0') for x in sheet_samples}
                # once any leading zeros have been removed, recalculate
                # sample_name_diff before continuing processing.
                sample_name_diff = sheet_samples - qsamples

                tids = qclient.get(f'{qurl}/categories=tube_id')['samples']
                # generate a map of sample_names to tube_ids for
                # GenPrepFileJob.
                sn_tid_map_by_project[project_name] = {
                    y[0]: x.replace(f'{qiita_id}.', '') for x, y in
                    tids.items()}
                tids = set(sn_tid_map_by_project[project_name].keys())
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
        raise PipelineError('\n'.join(errors))

    with open(sample_sheet_path, 'w') as f:
        sheet.write(f)

    # Create a Pipeline object
    try:
        pipeline = Pipeline(config_fp, run_identifier, sample_sheet_path,
                            None, out_dir, job_id)

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

    # if an Error has not been raised, assume there are no errors
    # in the sample-sheet. A successfully-created Pipeline object can
    # still contain a list of warnings about sample-sheet validation.
    # If any warnings are present, they should be reported to the user.
    if pipeline.warnings:
        msg = '\n'.join(pipeline.warnings)
        status_line.update_current_message('Sample-sheet has been flagged with'
                                           f' the following warnings: {msg}')

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

    status_line.update_current_message("Step 2 of 6: Converting BCL to fastq")

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
        convert_job.run(callback=status_line.update_job_step)
        fsr.write(convert_job.audit(samples), 'ConvertJob')

    status_line.update_current_message("Step 3 of 6: Adaptor & Host "
                                       "[optional] trimming")

    raw_fastq_files_path = join(pipeline.output_path, 'ConvertJob')

    config = pipeline.configuration['qc']
    qc_job = QCJob(raw_fastq_files_path,
                   pipeline.output_path,
                   sample_sheet_path,
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
                   job_id,
                   job_pool_size,
                   config['job_max_array_length'])

    if not skip_exec:
        qc_job.run(callback=status_line.update_job_step)
        fsr.write(qc_job.audit(samples), 'QCJob')

    status_line.update_current_message("Step 4 of 6: Generating FastQC & "
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
                           config['job_max_array_length'],
                           False)

    if not skip_exec:
        fastqc_job.run(callback=status_line.update_job_step)
        fsr.write(fastqc_job.audit(samples), 'FastQCJob')

    project_list = fastqc_job.project_names

    status_line.update_current_message("Step 5 of 6: Generating Prep "
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

    touched_studies_prep_info = {}

    if not skip_exec:
        gpf_job.run(callback=status_line.update_job_step)

        # concatenate the lists of paths across all study_ids into a single
        # list. Replace sample-names w/tube-ids in all relevant prep-files.
        preps = list(chain.from_iterable(gpf_job.prep_file_paths.values()))
        map_sample_names_to_tube_ids(preps, sn_tid_map_by_project)

        for study_id in gpf_job.prep_file_paths:
            for prep_file_path in gpf_job.prep_file_paths[study_id]:
                metadata_dict = pd.read_csv(prep_file_path,
                                            delimiter='\t').to_dict('index')

                # determine data_type based on sample-sheet
                # value will be from the Assay field
                data = {'prep_info': dumps(metadata_dict),
                        'study': study_id,
                        'data_type': assay_type}

                reply = qclient.post('/apitest/prep_template/', data=data)
                prep_id = reply['prep']

                if study_id not in touched_studies_prep_info:
                    touched_studies_prep_info[study_id] = []
                touched_studies_prep_info[study_id].append(prep_id)
    else:
        # replace sample-names w/tube-ids in all relevant prep-files.
        map_sample_names_to_tube_ids(join(pipeline.output_path,
                                          'GenPrepFileJob', 'PrepFiles',
                                          ('good-prep-file.txt')),
                                     sn_tid_map_by_project)

        # assume testing conditions and assign preps to study 1.
        metadata_dict = {
            'SKB8.640193': {'primer': 'GTGCCAGCMGCCGCGGTAA',
                            'barcode': 'GTCCGCAAGTTA',
                            'platform': 'Illumina',
                            'instrument_model': 'Illumina MiSeq'},
            'SKD8.640184': {'primer': 'GTGCCAGCMGCCGCGGTAA',
                            'barcode': 'GTCCGCAAGTTA',
                            'platform': 'Illumina',
                            'instrument_model': 'Illumina MiSeq'}}

        # valid entries include:
        # ['16S', '18S', 'ITS', 'Proteomic', 'Metabolomic', 'Metagenomic',
        #  'Multiomic', 'Metatranscriptomics', 'Viromics', 'Genomics',
        #  'Transcriptomics']

        data = {'prep_info': dumps(metadata_dict),
                'study': '1',
                'data_type': 'Metagenomic'}

        reply = qclient.post('/apitest/prep_template/', data=data)
        prep_id = reply['prep']
        touched_studies_prep_info['1'] = [prep_id]

    status_line.update_current_message("Step 6 of 6: Copying results to "
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
        # sif filenames are of the form:
        blanks_file = f'{run_identifier}_{project}_blanks.tsv'
        if sifs and [x for x in sifs if blanks_file in x]:
            # move uncompressed sifs to upload_dir.
            cmds.append(f'cd {out_dir}; mv {blanks_file} {upload_dir}')

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
