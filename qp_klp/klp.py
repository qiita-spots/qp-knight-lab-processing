# -----------------------------------------------------------------------------
# Copyright (c) 2020--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from os.path import exists, isdir
# from glob import glob


def list_folder(qclient, job_id, parameters, out_dir):
    """Run Woltka with the given parameters

    Parameters
    ----------
    qclient : tgp.qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    parameters : dict
        The parameter values to run split libraries
    out_dir : str
        The path to the job's output directory

    Returns
    -------
    bool, list, str
        The results of the job
    """
    qclient.update_job_step(job_id, "Step 1 of 3: Collecting information")
    input_folder = parameters.pop('input_folder')

    success = True
    ainfo = None
    msg = None
    if exists(input_folder) and isdir(input_folder):
        ainfo = ''
        # [
        #     ArtifactInfo('output', 'raw_job_folder',
        #                  [(f'{out_dir}/', 'string')])]
    else:
        success = False
        msg = "The path doesn't exist or is not a folder"

    return success, ainfo, msg
