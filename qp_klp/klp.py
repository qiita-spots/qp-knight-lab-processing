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
from os.path import join, split, exists
from sequence_processing_pipeline.Pipeline import Pipeline
from sequence_processing_pipeline.PipelineError import PipelineError
from sequence_processing_pipeline.ConvertJob import ConvertJob
from metapool import load_sample_sheet


CONFIG_FP = environ["QP_KLP_CONFIG_FP"]


class StatusUpdate():
    def __init__(self, qclient, job_id, msgs):
        self.qclient = qclient
        self.job_id = job_id
        self.msg = ''
        self.current_step = 0
        self.step_count = len(msgs)
        self.msgs = msgs

    def update_job_status(self, status, jid):
        # internal function implements a callback function for Pipeline.run().
        # :param id: PBS/Torque/or some other informative and current job id.
        # :param status: status message
        self.qclient.update_job_step(self.job_id,
                                     self.msg + f" ({jid}: {status})")

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
    # Assume that for a job to be considered a restart, there must be work
    # performed worth re-starting for. Since the working directory for each
    # step is created only if the previous steps were successful, testing
    # for the presence of them ensures that n-1 steps exist and were
    # successful.

    # at minimum, ConvertJob needs to have been successful.
    is_restart = True if exists(join(out_dir, 'NuQCJob')) else False

    if is_restart:
        # Assume ConvertJob directory exists and parse the job-script found
        # there. If this is a restart, we won't be given the run-identifier,
        # the lane number, and the sample-sheet as input parameters.
        some_path = join(out_dir, 'ConvertJob', 'ConvertJob.sh')
        result = ConvertJob.parse_job_script(some_path)
        run_identifier = split(result['out_directory'])[-1]
        user_input_file = result['sample_sheet_path']
        sheet = load_sample_sheet(user_input_file)
        # on Amplicon runs, lane_number is always 1, and this will be
        # properly reflected in the dummy sample-sheet as well.
        lane_number = sheet.get_lane_number()

        # check if sample-sheet is a dummy-sample-sheet. If this is an
        # Amplicon run, then Assay type will be 'TruSeq HT' and Chemistry
        # will be 'Amplicon'. For now, raise Error on restarting an
        # Amplicon run so we don't have to search for the pre-prep file.
        if sheet.Header['Assay'] == 'TruSeq HT' and \
           sheet.Header['Chemistry'] == 'Amplicon':
            raise ValueError("Restarting Amplicon jobs currently unsupported")
    else:
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

    if is_restart:
        pass
    else:
        # save raw data to file
        with open(uif_path, 'w') as f:
            f.write(user_input_file['body'])

    if Pipeline.is_sample_sheet(uif_path):
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
        lane_number = 1
    else:
        # file doesn't look like a sample-sheet, or a valid mapping file.
        return False, None, ("Your uploaded file doesn't appear to be a sample"
                             "-sheet or a mapping-file.")

    msgs = ["Setting up pipeline", "Getting project information",
            "Converting data", "Performing quality control",
            "Generating reports", "Generating preps",
            "Generating sample information ", "Registering blanks in Qiita",
            "Loading preps into Qiita", "Generating packaging commands",
            "Packaging results", "SPP finished"]

    status_line = StatusUpdate(qclient, job_id, msgs)
    status_line.update_current_message()

    skip_steps = []
    if is_restart:
        # figure out what actually needs to be skipped if restarting:
        if exists(join(out_dir, 'NuQCJob')):
            skip_steps.append('ConvertJob')

        if exists(join(out_dir, 'FastQCJob')):
            skip_steps.append('NuQCJob')

        if exists(join(out_dir, 'GenPrepFileJob')):
            skip_steps.append('FastQCJob')

        if exists(join(out_dir, 'cmds.log')):
            skip_steps.append('GenPrepFileJob')

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
                pipeline, job_id, status_line, lane_number,
                is_restart=is_restart)
        else:
            step = Amplicon(
                pipeline, job_id, status_line, lane_number,
                is_restart=is_restart)

        status_line.update_current_message()

        step.precheck(qclient)

        # set update=False to prevent updating Qiita database and copying
        # files into uploads directory. Useful for testing.
        step.execute_pipeline(qclient,
                              status_line.update_current_message,
                              update=True,
                              skip_steps=skip_steps)

        status_line.update_current_message()

    except PipelineError as e:
        return False, None, str(e)

    # return success, ainfo, and the last status message.
    paths = [(f'{final_results_path}/', 'directory')]
    return (True, [ArtifactInfo('output', 'job-output-folder', paths)],
            status_line.msg)
