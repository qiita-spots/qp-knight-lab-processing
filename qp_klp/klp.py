# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from functools import partial
from os import environ, symlink
from qiita_client import ArtifactInfo
from os import makedirs
from os.path import join, exists, basename
from datetime import datetime
from sequence_processing_pipeline.PipelineError import PipelineError
from metapool import load_sample_sheet, MetagenomicSampleSheetv90
import sample_sheet
import pandas as pd
from .Workflows import WorkflowError
from .WorkflowFactory import WorkflowFactory


CONFIG_FP = environ["QP_KLP_CONFIG_FP"]


class StatusUpdate():
    def __init__(self, qclient, job_id):
        self.qclient = qclient
        self.job_id = job_id
        self.msg = ''

    def update_job_status(self, msg):
        self.qclient.update_job_step(self.job_id, msg)


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
    status_line = StatusUpdate(qclient, job_id)

    # NB: RESTART FUNCTIONALITY:
    # To restart a job that's failed, simply create a ProcessingJob() object
    # in Qiita's interpreter using the Qiita Job ID aka the name of the
    # working directory. _set_status("in construction") and then call
    # submit(). New SPP jobs will record all the information needed to perform
    # a restart should it fail.
    out_path = partial(join, out_dir)
    is_restart = True if exists(out_path('restart_me')) else False

    if is_restart:
        status_line.update_job_status("Restarting pipeline")
        with open(out_path('restart_me'), 'r') as f:
            lines = f.readlines()
            lines = [x.strip() for x in lines]
            lines = [x for x in lines if x != '']
            if len(lines) != 2:
                raise ValueError(f"{out_path('restart_me')} can only contain "
                                 "the run-identifier on line 1 and the name "
                                 "of the user input file on line 2")
        run_identifier = lines[0]
        uif_path = out_path(lines[1])

        if not exists(uif_path):
            raise ValueError(f"{uif_path} does not exist")

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
        # This is a new, fresh SPP job. Get the parameters from the user and
        # create a 'restart_me' file in the working directory.

        status_line.update_job_status("Setting up pipeline")

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

        # the run_identifier must be saved because it is not always preserved
        # in a dependable location downstream. The user input file must be
        # saved because it is always a unique name and it cannot be guaranteed
        # to be the only .csv or .txt file at the root level of the working
        # directory. The Lane number is not needed because it is written into
        # the user_input file on the first run.
        restart_file_path = out_path('restart_me')
        with open(restart_file_path, 'w') as f:
            f.write(f"{run_identifier}\n{uif_path}")

    final_results_path = out_path('final_results')
    makedirs(final_results_path, exist_ok=True)

    try:
        kwargs = {'qclient': qclient,
                  'uif_path': uif_path,
                  'lane_number': lane_number,
                  'config_fp': CONFIG_FP,
                  'run_identifier': run_identifier,
                  'output_dir': out_dir,
                  'job_id': job_id,
                  'status_update_callback': status_line.update_job_status,
                  # set 'update_qiita' to False to avoid updating Qiita DB
                  # and copying files into uploads dir. Useful for testing.
                  'update_qiita': True,
                  'is_restart': is_restart}

        workflow = WorkflowFactory().generate_workflow(**kwargs)

        status_line.update_job_status("Getting project information")

        workflow.execute_pipeline()

        status_line.update_job_status("SPP finished")

    except (PipelineError, WorkflowError) as e:
        # assume AttributeErrors are issues w/bad sample-sheets or
        # mapping-files.
        return False, None, str(e)

    # return success, ainfo, and the last status message.
    paths = [(f'{final_results_path}/', 'directory')]
    return (True, [ArtifactInfo('output', 'job-output-folder', paths)],
            status_line.msg)


def prep_NuQCJob(qclient, job_id, parameters, out_dir):
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
    status_line = StatusUpdate(qclient, job_id)

    status_line.update_job_status("Setting up pipeline")

    pid = parameters.pop('prep_id')

    prep_info = qclient.get(f'/qiita_db/prep_template/{pid}/')
    dt = prep_info['data_type']
    sid = prep_info['study']
    if dt not in {'Metagenomic', 'Metatranscriptomic'}:
        return (False, [], f'Prep {pid} has a not valid data type: {dt}')
    aid = prep_info['artifact']
    if not str(aid).isnumeric():
        return (False, [], f'Prep {pid} has a not valid artifact: {aid}')

    files, pt = qclient.artifact_and_preparation_files(aid)
    pt.set_index('sample_name', inplace=True)

    out_path = partial(join, out_dir)
    final_results_path = out_path('final_results')
    makedirs(final_results_path, exist_ok=True)

    project_name = f'qiita-{pid}-{aid}_{sid}'

    sheet = MetagenomicSampleSheetv90()
    sheet.Header['IEMFileVersion'] = '4'
    sheet.Header['Date'] = datetime.today().strftime('%m/%d/%y')
    sheet.Header['Workflow'] = 'GenerateFASTQ'
    sheet.Header['Application'] = 'FASTQ Only'
    sheet.Header['Assay'] = prep_info['data_type']
    sheet.Header['Description'] = f'prep_NuQCJob - {pid}'
    sheet.Header['Chemistry'] = 'Default'
    sheet.Header['SheetType'] = 'standard_metag'
    sheet.Header['SheetVersion'] = '90'
    sheet.Header['Investigator Name'] = 'Qiita'
    sheet.Header['Experiment Name'] = project_name

    sheet.Bioinformatics = pd.DataFrame(
        columns=['Sample_Project', 'ForwardAdapter', 'ReverseAdapter',
                 'library_construction_protocol',
                 'experiment_design_description',
                 'PolyGTrimming', 'HumanFiltering', 'QiitaID'],
        data=[[project_name, 'NA', 'NA', 'NA', 'NA', 'FALSE', 'TRUE', sid]])

    for k, vals in pt.iterrows():
        k = k.split('.', 1)[-1]
        sample = {
            'Sample_Name': k,
            'Sample_ID': k.replace('.', '_'),
            'Sample_Plate': '',
            'well_id_384': '',
            'I7_Index_ID': '',
            'index': vals['index'],
            'I5_Index_ID': '',
            'index2': vals['index2'],
            'Sample_Project': project_name,
            'Well_description': '',
            'Sample_Well': '',
            'Lane': '1'}
        sheet.add_sample(sample_sheet.Sample(sample))

    sheet.Contact = pd.DataFrame(
        columns=['Email', 'Sample_Project'],
        data=[['qiita.help@gmail.com', project_name]])

    new_sample_sheet = out_path('sample-sheet.csv')
    with open(new_sample_sheet, 'w') as f:
        sheet.write(f, 1)

    # now that we have a sample_sheet we can fake the
    # ConvertJob folder so we are ready for the restart
    convert_path = out_path('ConvertJob')
    project_folder = out_path('ConvertJob', project_name)
    makedirs(project_folder, exist_ok=True)

    for _, fs in files.items():
        for f in fs:
            bn = basename(f['filepath'])
            symlink(f['filepath'], f'{project_folder}/{bn}')

    # create job_completed file to skip this step
    with open(f'{convert_path}/job_completed', 'w') as f:
        f.write('')

    try:
        kwargs = {'qclient': qclient,
                  'uif_path': new_sample_sheet,
                  'lane_number': "1",
                  'config_fp': CONFIG_FP,
                  'run_identifier': '211021_A00000_0000_SAMPLE',
                  'output_dir': out_dir,
                  'job_id': job_id,
                  'status_update_callback': status_line.update_job_status,
                  # set 'update_qiita' to False to avoid updating Qiita DB
                  # and copying files into uploads dir. Useful for testing.
                  'update_qiita': True,
                  'is_restart': True}

        workflow = WorkflowFactory().generate_workflow(**kwargs)

        status_line.update_job_status("Getting project information")

        workflow.execute_pipeline()

        status_line.update_job_status("SPP finished")

    except (PipelineError, WorkflowError) as e:
        # assume AttributeErrors are issues w/bad sample-sheets or
        # mapping-files.
        return False, None, str(e)

    # # return success, ainfo, and the last status message.
    paths = [(f'{final_results_path}/', 'directory')]
    return (True, [ArtifactInfo('output', 'job-output-folder', paths)], '')
