# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from functools import partial
from os import environ
from qiita_client import ArtifactInfo
from random import sample as rsampl
from qp_klp.Amplicon import Amplicon
from qp_klp.Metagenomic import Metagenomic
from qp_klp.Step import Step
from json import dumps
from os import makedirs
from os.path import join
from sequence_processing_pipeline.Pipeline import Pipeline
from sequence_processing_pipeline.PipelineError import PipelineError
from collections import defaultdict


CONFIG_FP = environ["QP_KLP_CONFIG_FP"]


class StatusUpdate():
    def __init__(self, qclient, job_id, msgs):
        self.qclient = qclient
        self.job_id = job_id
        self.msg = ''
        self.current_step = 1
        self.step_count = len(msgs)
        self.msgs = msgs

    def update_job_status(self, status, id):
        # internal function implements a callback function for Pipeline.run().
        # :param id: PBS/Torque/or some other informative and current job id.
        # :param status: status message
        self.qclient.update_job_step(self.job_id,
                                     self.msg + f" ({id}: {status})")

    def update_current_message(self, include_step=True):
        # internal function that sets current_message to the new value before
        # updating the job step in the UI.
        if include_step:
            self.msg = (f"Step {self.current_step} of {self.step_count}: "
                        f"{self.msgs[self.step_count]}")
        else:
            self.msg = (f"{self.msgs[self.step_count]}")

        self.current_step += 1

        self.qclient.update_job_step(self.job_id, self.msg)


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

    # available fields for parameters are:
    #   run_identifier, sample_sheet, content_type, filename, lane_number
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

    # save raw data to file
    with open(uif_path, 'w') as f:
        f.write(user_input_file['body'])

    if Pipeline.is_sample_sheet(uif_path):
        # if file follows basic sample-sheet format, then it is most likely
        # a sample-sheet, even if it's an invalid one.

        # a valid sample-sheet is going to have one and only one occurance of
        # 'Assay,Metagenomic' or 'Assay,Metatranscriptomic'. Anything else is
        # an error.

        tmp = [x for x in user_input_file['body'] if 'Assay' in x]

        if len(tmp) != 1:
            return False, None, ("Your uploaded file doesn't appear to be a"
                                 "valid sample-sheet.")

        pipeline_type = None
        for p_type in Step.META_TYPES:
            if p_type in tmp[0].lower():
                pipeline_type = p_type

        if pipeline_type is None:
            return False, None, ("Your uploaded file doesn't appear to be a"
                                 "valid sample-sheet.")
    elif Pipeline.is_mapping_file(uif_path):
        # if file is readable as a basic TSV and contains all the required
        # headers, then treat this as a mapping file, even if it's an invalid
        # one.
        pipeline_type = Step.AMPLICON_TYPE
        lane_number = -1
    else:
        # file doesn't look like a sample-sheet, or a valid mapping file.
        return False, None, ("Your uploaded file doesn't appear to be a sample"
                             "-sheet or a mapping-file.")

    msgs = ["Setting up pipeline", "Getting project information",
            "Converting BCL to Fastq", "Performing Quality Control",
            "Generating FastQC Reports", "Generating Prep Files",
            "Generating Sample Information ", "Updating Blanks in Qiita",
            "Adding Prep Templates to Qiita", "Packaging results",
            "Finishing up"]

    status_line = StatusUpdate(qclient, job_id, msgs)
    status_line.update_current_message()

    try:
        pipeline = Step.generate_pipeline(pipeline_type,
                                          uif_path,
                                          lane_number,
                                          CONFIG_FP,
                                          run_identifier,
                                          out_dir,
                                          job_id)
    except PipelineError as e:
        # assume AttributeErrors are issues w/bad sample-sheets or
        # mapping-files.
        return False, None, str(e)

    try:
        errors = []
        sn_tid_map_by_project = {}

        status_line.update_current_message()

        # Update get_project_info() so that it can return a list of
        # samples in projects['samples']. Include blanks in projects['blanks']
        projects = pipeline.get_project_info(short_names=True)

        for project in projects:
            qsam, tids = Step.get_samples_in_qiita(qclient,
                                                   project['qiita_id'])
            qsam = set(qsam)
            tids = {tids[k][0] for k in tids} if tids else None

            # compare the list of samples from the user against the list of
            # registered samples from Qiita. If the project has tube-ids, then
            # compare against tids rather than qsam.
            samples = set(pipeline.get_sample_ids())

            # strip any leading zeroes from the sample-ids. Note that
            # if a sample-id has more than one leading zero, all of
            # them will be removed.
            # my_samples = {x.lstrip('0') for x in samples}

            sample_name_diff = samples - tids if tids else samples - qsam

            project_name = project['project_name']

            if sample_name_diff:
                len_overlap = len(sample_name_diff)
                # selecting at random k=5 samples to minimize space in display
                samx = ', '.join(qsam if len(qsam) < 6 else rsampl(qsam, k=5))

                # selecting the up to 4 first samples to minimize space in
                # display
                missing = ', '.join(sorted(sample_name_diff)[:4])
                errors.append(
                    f'{project_name} has {len_overlap} missing samples (i.e. '
                    f'{missing}). Some samples from Qiita: {samx}. ')

        if errors:
            raise PipelineError('\n'.join(errors))

        # find the uploads directory all trimmed files will need to be
        # moved to and generate a map.
        special_map = Step.generate_special_map(
            qclient.get("/qiita_db/artifacts/types/"),
            pipeline.get_project_info())

        if pipeline.pipeline_type in Step.META_TYPES:
            step = Metagenomic(pipeline, job_id, sn_tid_map_by_project,
                               status_line)
        else:
            # pipeline.pipeline_type == Step.AMPLICON_TYPE:
            step = Amplicon(pipeline, job_id, sn_tid_map_by_project,
                            status_line)

        status_line.update_current_message()
        step.convert_bcl_to_fastq()

        status_line.update_current_message()
        step.quality_control()

        status_line.update_current_message()
        step.generate_reports()

        status_line.update_current_message()
        step.generate_prep_file()

        from_qiita = {}

        for study_id in step.prep_file_paths:
            samples = list(qclient.get(f'/api/v1/study/{study_id}/samples'))
            from_qiita[study_id] = samples

        status_line.update_current_message()
        sifs = step.generate_sifs(from_qiita)

        status_line.update_current_message()
        Step.update_blanks_in_qiita(sifs, qclient)

        prep_file_paths = step.get_prep_file_paths()

        status_line.update_current_message()
        ptype = step.pipeline.pipeline_type
        Step.update_prep_templates(qclient, prep_file_paths, ptype)

        touched_studies_prep_info = defaultdict(list)

        for study_id in step.prep_file_paths:
            for prep_file_path in step.prep_file_paths[study_id]:
                metadata = Step.parse_prep_file(prep_file_path)
                data = {'prep_info': dumps(metadata),
                        'study': study_id,
                        'data_type': step.pipeline.pipeline_type}

                reply = qclient.post('/qiita_db/prep_template/', data=data)
                prep_id = reply['prep']

                touched_studies_prep_info[study_id].append(prep_id)

        # generate commands to execute
        status_line.update_current_message()
        step.generate_commands(special_map, touched_studies_prep_info)

        status_line.update_current_message()
        step.execute_commands()

    except PipelineError as e:
        return False, None, str(e)

    # return success, ainfo, and the last status message.
    paths = [(f'{final_results_path}/', 'directory')]
    return (True, [ArtifactInfo('output', 'job-output-folder', paths)],
            status_line.msg)
