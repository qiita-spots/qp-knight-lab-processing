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
from qp_klp.Amplicon import Amplicon
from qp_klp.Metagenomic import Metagenomic
from qp_klp.Step import Step
from os import makedirs
from os.path import join
from sequence_processing_pipeline.Pipeline import Pipeline
from sequence_processing_pipeline.PipelineError import PipelineError


CONFIG_FP = environ["QP_KLP_CONFIG_FP"]


class StatusUpdate():
    def __init__(self, qclient, job_id, msgs):
        self.qclient = qclient
        self.job_id = job_id
        self.msg = ''
        self.current_step = 0
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
            msg = "Step %d of %d: " % (self.current_step + 1, self.step_count)
        else:
            msg = ""

        self.msg = msg + self.msgs[self.current_step]

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

        # a valid sample-sheet is going to have one and only one occurrence of
        # 'Assay,Metagenomic' or 'Assay,Metatranscriptomic'. Anything else is
        # an error.

        # works best from file
        with open(uif_path, 'r') as f:
            assay = [x for x in f.readlines() if 'Assay' in x]

            if len(assay) == 0:
                return False, None, ("Assay type is not defined in the "
                                     "sample-sheet")

            if len(assay) > 1:
                return False, None, ("Assay type is defined multiple times "
                                     "within the sample-sheet")

        pipeline_type = None
        for p_type in Step.META_TYPES:
            if p_type in assay[0]:
                pipeline_type = p_type

        if pipeline_type is None:
            msg = [f"'{x}'" for x in Step.META_TYPES]
            return False, None, ("The following Assay types are valid within"
                                 " a sample-sheet: %s" % ', '.join(msg))

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
        if pipeline.pipeline_type in Step.META_TYPES:
            step = Metagenomic(
                pipeline, job_id, status_line, lane_number)
        else:
            step = Amplicon(
                pipeline, job_id, status_line, lane_number)

        step.get_tube_ids_from_qiita(qclient)

        status_line.update_current_message()

        # compare sample-ids/tube-ids in sample-sheet/mapping file
        # against what's in Qiita.
        results = step.compare_samples_against_qiita()

        if results is not None:
            msgs = []
            for comparison in results:
                not_in_qiita_count = len(comparison['samples_not_in_qiita'])
                examples_in_qiita = ', '.join(comparison['examples_in_qiita'])
                p_name = comparison['project_name']
                uses_tids = comparison['tids']

                msgs.append(f"Project '{p_name}' has {not_in_qiita_count} "
                            "samples not registered in Qiita.")

                msgs.append(f"Some registered samples in Project '{p_name}'"
                            f" include: {examples_in_qiita}")

                if uses_tids:
                    msgs.append(f"Project '{p_name}' is using tube-ids. You "
                                "may be using sample names in your file.")

            if msgs:
                raise PipelineError('\n'.join(msgs))

        # set update=False to prevent updating Qiita database and copying
        # files into uploads directory. Useful for testing.
        step.execute_pipeline(qclient,
                              status_line.update_current_message,
                              update=True)

    except PipelineError as e:
        return False, None, str(e)

    # return success, ainfo, and the last status message.
    paths = [(f'{final_results_path}/', 'directory')]
    return (True, [ArtifactInfo('output', 'job-output-folder', paths)],
            status_line.msg)
