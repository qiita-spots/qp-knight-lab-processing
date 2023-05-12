# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from functools import partial
from os import environ
from qp_klp.klp_util import StatusUpdate
from os import makedirs
from os.path import join
from qiita_client import ArtifactInfo
from qp_klp.klp_util import (update_blanks_in_qiita, generate_special_map,
                             generate_pipeline, update_prep_templates,
                             get_registered_samples_in_qiita)
from random import sample as rsampl
from sequence_processing_pipeline.Pipeline import Pipeline
from sequence_processing_pipeline.PipelineError import PipelineError
from qp_klp.AmpliconStep import AmpliconStep
from qp_klp.MetagenomicStep import MetagenomicStep


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
    user_input_file = parameters.pop('sample_sheet')
    lane_number = parameters.pop('lane_number')

    if {'body', 'content_type', 'filename'} != set(user_input_file):
        return False, None, ("This doesn't appear to be a valid sample sheet "
                             "or mapping file; please review.")

    out_path = partial(join, out_dir)
    final_results_path = out_path('final_results')
    makedirs(final_results_path, exist_ok=True)
    # replace any whitespace in the filename with underscores
    uif_path = out_path(user_input_file['filename'].replace(' ', '_'))

    if Pipeline.is_sample_sheet(uif_path):
        # if file follows basic sample-sheet format, then it is most likely
        # a sample-sheet, even if it's an invalid one.
        pipeline_type = 'metagenomic'
        step_count = 6
    elif Pipeline.is_mapping_file(uif_path):
        # if file is readable as a basic TSV and contains all of the required
        # headers, then treat this as a mapping file, even if it's an invalid
        # one.
        pipeline_type = 'amplicon'
        lane_number = -1
        step_count = 5
    else:
        # file doesn't look like a sample-sheet, or a valid mapping file.
        return False, None, ("Your uploaded file doesn't appear to be a sample"
                             "-sheet or a mapping-file.")

    status_line = StatusUpdate(qclient, job_id, step_count)
    status_line.update_current_message("Setting up pipeline")

    # save raw data to file
    with open(uif_path, 'w') as f:
        f.write(user_input_file['body'])

    try:
        pipeline = generate_pipeline(pipeline_type,
                                     uif_path,
                                     lane_number,
                                     CONFIG_FP,
                                     run_identifier,
                                     out_dir,
                                     job_id)

        errors = []
        sn_tid_map_by_project = {}

        # TODO: Update get_project_info() so that it can return a list of
        #  samples in projects['samples']. Include blanks in projects['blanks']
        projects = pipeline.get_project_info(short_names=True)

        for project in projects:
            qsam, tids = get_registered_samples_in_qiita(qclient,
                                                         project['qiita_id'])

            # compare the list of samples from the user against the list of
            # registered samples from Qiita. If the project has tube-ids, then
            # compare against tids rather than qsam.
            samples = project['samples']
            sample_name_diff = samples - tids if tids else samples - qsam

            if tids:
                # strip any leading zeroes from the sample-ids. Note that
                # if a sample-id has more than one leading zero, all of
                # them will be removed.
                my_samples = {x.lstrip('0') for x in my_samples}
                # once any leading zeros have been removed, recalculate
                # sample_name_diff before continuing processing.
                sample_name_diff = my_samples - qsam

                # before we report as an error, check tube_id.
                if pipeline_type == 'metagenomic':
                    tids = qclient.get(f'{qurl}/categories=tube_id')[
                        'samples']

                # generate a map of sample_names to tube_ids for
                # GenPrepFileJob.
                sn_tid_map_by_project[project_name] = {
                    y[0]: x.replace(f'{qiita_id}.', '') for x, y in
                    tids.items()}
                tids = set(sn_tid_map_by_project[project_name].keys())
                tube_id_diff = my_samples - tids
                if not tube_id_diff:
                    continue
                len_tube_id_overlap = len(tube_id_diff)
                tidsx = ', '.join(
                    tids if len(tids) < 6 else rsampl(tids, k=5))

                error_tube_id = (
                    f'tube_id in Qiita but {len_tube_id_overlap} missing '
                    f'samples. Some samples from tube_id: {tidsx}.')

            len_overlap = len(sample_name_diff)
            # selecting at random k=5 samples to minimize space in display
            samx = ', '.join(qsam if len(qsam) < 6 else rsampl(qsam, k=5))

            # selecting the up to 4 first samples to minimize space in
            # display
            missing = ', '.join(sorted(sample_name_diff)[:4])
            errors.append(
                f'{project_name} has {len_overlap} missing samples (i.e. '
                f'{missing}). Some samples from Qiita: {samx}. '
                f'{error_tube_id}')

        if errors:
            raise PipelineError('\n'.join(errors))

        # find the uploads directory all trimmed files will need to be
        # moved to and generate a map.
        special_map = generate_special_map(
            qclient.get("/qiita_db/artifacts/types/"),
            pipeline.get_project_info())

        if pipeline.type == 'metagenomic':
            step = MetagenomicStep(pipeline, job_id, status_line,
                                   sn_tid_map_by_project)
        else:
            # pipeline.type == 'amplicon':
            step = AmpliconStep(pipeline, job_id, status_line,
                                sn_tid_map_by_project)

        step.convert_bcl_to_fastq()

        step.quality_control()

        step.generate_reports()

        step.generate_prep_file()

        sifs = step.work_in_progress()

        update_blanks_in_qiita(sifs, qclient)

        prep_file_paths = step.get_prep_file_paths()

        update_prep_templates(qclient, prep_file_paths)

        # generate commands to execute
        step.generate_commands()

        step.execute_commands()

    except PipelineError as e:
        return False, None, str(e)

    # return success, ainfo, and the last status message.
    paths = [(f'{final_results_path}/', 'directory')]
    return (True, [ArtifactInfo('output', 'job-output-folder', paths)],
            status_line.msg)
