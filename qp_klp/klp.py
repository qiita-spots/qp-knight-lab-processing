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
from os import makedirs
from os.path import join, split, exists
from sequence_processing_pipeline.PipelineError import PipelineError
from sequence_processing_pipeline.ConvertJob import ConvertJob
from metapool import load_sample_sheet
from .Workflows import WorkflowFactory, WorkflowError


CONFIG_FP = environ["QP_KLP_CONFIG_FP"]


class StatusUpdate():
    def __init__(self, qclient, job_id):
        self.qclient = qclient
        self.job_id = job_id
        self.msg = ''

    def update_job_status(self, status, jid):
        # internal function implements a callback function for Pipeline.run().
        # :param id: PBS/Torque/or some other informative and current job id.
        # :param status: status message
        self.qclient.update_job_step(self.job_id,
                                     self.msg + f" ({jid}: {status})")

    def update_current_message(self, step=None):
        # internal function that sets current_message to the new value before
        # updating the job step in the UI.
        if step:
            msg = "Step %d of %d: " % (step[0], step[1])
        else:
            msg = ""

        self.msg = msg + self.msgs[self.current_step]

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

    # if a file exists in the working (out) directory named, 'restart_me'
    # then consider this a job that executed previously but failed.
    out_path = partial(join, out_dir)
    is_restart = True if exists(out_path('restart_me')) else False

    status_line = StatusUpdate(qclient, job_id)

    if is_restart:
        status_line.update_current_message("Restarting pipeline")
        # Assume ConvertJob directory exists and parse the job-script found
        # there. If this is a restart, we won't be given the run-identifier,
        # the lane number, and the sample-sheet as input parameters.
        some_path = join(out_dir, 'ConvertJob', 'ConvertJob.sh')
        result = ConvertJob.parse_job_script(some_path)
        run_identifier = split(result['run_directory'])[-1]
        uif_path = result['sample_sheet_path']
        sheet = load_sample_sheet(uif_path)
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

        # add a note for the wetlab that this job was restarted.
        with open(join(out_dir, 'notes.txt'), 'w') as f:
            f.write("This job was restarted.\n"
                    "failed_samples.html may contain incorrect data.\n")
    else:
        status_line.update_current_message("Setting up pipeline")
        # available fields for parameters are:
        #   run_identifier, sample_sheet, content_type, filename, lane_number
        run_identifier = parameters.pop('run_identifier')
        user_input_file = parameters.pop('sample_sheet')
        lane_number = parameters.pop('lane_number')

        if {'body', 'content_type', 'filename'} != set(user_input_file):
            return False, None, ("This doesn't appear to be a valid sample "
                                 "sheet or mapping file; please review.")
        uif_path = out_path(user_input_file['filename'].replace(' ', '_'))
        # save raw data to file
        with open(uif_path, 'w') as f:
            f.write(user_input_file['body'])

    final_results_path = out_path('final_results')
    makedirs(final_results_path, exist_ok=True)

    try:
        kwargs = {'qclient': qclient,
                  'uif_path': uif_path,
                  'lane_number': lane_number,
                  'config_path': CONFIG_FP,
                  'run_identifier': run_identifier,
                  'output_dir': out_dir,
                  'job_id': job_id,
                  'status_update_callback': status_line.update_current_message,
                  # set 'update_qiita' to False to avoid updating Qiita DB
                  # and copying files into uploads dir. Useful for testing.
                  'update_qiita': True,
                  'is_restart': is_restart}

        workflow = WorkflowFactory().generate_workflow(**kwargs)

        status_line.update_current_message("Getting project information")

        workflow.execute_pipeline()

        status_line.update_current_message("SPP finished")

    except (PipelineError, WorkflowError) as e:
        # assume AttributeErrors are issues w/bad sample-sheets or
        # mapping-files.
        return False, None, str(e)

    # return success, ainfo, and the last status message.
    paths = [(f'{final_results_path}/', 'directory')]
    return (True, [ArtifactInfo('output', 'job-output-folder', paths)],
            status_line.msg)
