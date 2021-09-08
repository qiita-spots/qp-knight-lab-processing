# -----------------------------------------------------------------------------
# Copyright (c) 2020--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from os.path import exists, isdir, join
from glob import glob
from functools import partial

from qiita_client import ArtifactInfo


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
    qclient.update_job_step(job_id, "Step 1 of 3: Collecting information")
    run_identifier = parameters.pop('run_identifier')

    success = True
    ainfo = None
    msg = None
    if exists(run_identifier) and isdir(run_identifier):
        outpath = partial(join, out_dir)

        # copytree('metat-pilot', outpath('multiqc-report'))

        with open(outpath('listing.txt'), 'w') as f:
            f.write('\n'.join(glob(join(run_identifier, '*'))))

        file = parameters['sample_sheet']
        with open(outpath(file['filename']), 'w') as f:
            f.write(file['body'])

        ainfo = [
            ArtifactInfo('output', 'job-output-folder',
                         [(f'{out_dir}/', 'directory')])
        ]
    else:
        success = False
        msg = "The path doesn't exist or is not a folder"

    return success, ainfo, msg
