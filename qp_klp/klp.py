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
from qp_klp.process_metagenomics_job import process_metagenomics
from sequence_processing_pipeline.Pipeline import Pipeline
from sequence_processing_pipeline.PipelineError import PipelineError
from qp_klp.klp_util import StatusUpdate


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

    # StatusUpdate uses qclient to easily update the status-line and maintain
    # a certain amount of state needed to provide the updates users are
    # accustomed to. E.g. "Step 1 of 6: Setting up pipeline"
    status_line = StatusUpdate(qclient, job_id)

    status_line.update_current_message("Step 1 of 6: Setting up pipeline")

    if {'body', 'content_type', 'filename'} != set(user_input_file):
        return False, None, ("This doesn't appear to be a valid sample sheet "
                             "or mapping file; please review.")

    outpath = partial(join, out_dir)
    final_results_path = outpath('final_results')
    makedirs(final_results_path, exist_ok=True)
    # replace any whitespace in the filename with underscores
    uif_path = outpath(user_input_file['filename'].replace(' ', '_'))
    # save raw data to file
    with open(uif_path, 'w') as f:
        f.write(user_input_file['body'])

    if Pipeline.is_mapping_file(uif_path):
        return False, None, "Not Implemented"
    else:
        try:
            ainfo = process_metagenomics(uif_path, lane_number, qclient,
                                         run_identifier, out_dir, job_id,
                                         skip_exec, job_pool_size,
                                         final_results_path, CONFIG_FP,
                                         status_line)
        except PipelineError as e:
            return False, None, str(e)
                
    status_line.update_current_message("Main Pipeline Finished, processing "
                                       "results")
    
    # return success, ainfo, and the last status message.
    return True, ainfo, status_line.msg
