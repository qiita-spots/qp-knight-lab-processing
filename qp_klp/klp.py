# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from functools import partial
from inspect import stack
from os import environ
from os import makedirs
from os.path import join
from qp_klp.process_amplicon_job import process_amplicon
from qp_klp.process_metagenomics_job import process_metagenomics
from sequence_processing_pipeline.Pipeline import Pipeline


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
    job_pool_size = 30

    # checking if this is running as part of the unittest
    # https://stackoverflow.com/a/25025987
    skip_exec = True if [x for x in stack() if
                         'unittest' in x.filename] else False

    success = True
    ainfo = None
    msg = None

    # maintains state of the current message, minus additional updates from
    # a callback function. E.g. "Step 1 of 6: Setting up pipeline"

    def _update_job_step(id, status):
        # internal function implements a callback function for Pipeline.run().
        # :param id: PBS/Torque/or some other informative and current job id.
        # :param status: status message
        qclient.update_job_step(job_id,
                                _update_job_step.msg + f" ({id}: {status})")

    # initialize static variable to maintain current message
    _update_job_step.msg = ""

    def _update_current_message(msg):
        # internal function that sets current_message to the new value before
        # updating the job step in the UI.
        _update_job_step.msg = msg
        qclient.update_job_step(job_id, msg)

    _update_current_message("Step 1 of 6: Setting up pipeline")

    if {'body', 'content_type', 'filename'} == set(user_input_file):
        outpath = partial(join, out_dir)
        final_results_path = outpath('final_results')
        makedirs(final_results_path, exist_ok=True)
        # replace any whitespace in the filename with underscores
        uif_path = outpath(user_input_file['filename'].replace(' ', '_'))
        # save raw data to file
        with open(uif_path, 'w') as f:
            f.write(user_input_file['body'])

        if Pipeline.is_mapping_file(uif_path):
            success, ainfo, msg = process_amplicon(uif_path,
                                                   qclient,
                                                   run_identifier,
                                                   out_dir,
                                                   job_id,
                                                   _update_current_message,
                                                   skip_exec,
                                                   _update_job_step,
                                                   job_pool_size,
                                                   final_results_path,
                                                   success, msg,
                                                   CONFIG_FP)
        else:
            success, ainfo, msg = process_metagenomics(uif_path,
                                                       lane_number,
                                                       qclient,
                                                       run_identifier,
                                                       out_dir,
                                                       job_id,
                                                       _update_current_message,
                                                       skip_exec,
                                                       _update_job_step,
                                                       job_pool_size,
                                                       final_results_path,
                                                       success, msg,
                                                       CONFIG_FP)
    else:
        success = False
        msg = "This doesn't appear to be a valid sample sheet; please review."

    _update_current_message("Main Pipeline Finished, processing results")

    return success, ainfo, msg
