# -----------------------------------------------------------------------------
# Copyright (c) 2020--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from os.path import exists, isdir, join
from glob import glob

from qiita_client import ArtifactInfo


def list_folder(qclient, job_id, parameters, out_dir):
    """Create file listing in output directory

    Parameters
    ----------
    qclient : tgp.qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    parameters : dict
        The parameter values for input folder to use
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
        with open(join(out_dir, 'listing.txt'), 'w') as f:
            f.write('\n'.join(glob(join(input_folder, '*'))))

        ainfo = [
            ArtifactInfo('output', 'job-output-folder',
                         [(f'{out_dir}/', 'directory')])
        ]
    else:
        success = False
        msg = "The path doesn't exist or is not a folder"

    return success, ainfo, msg
