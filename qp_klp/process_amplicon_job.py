from os import makedirs
from os.path import join
from qiita_client import ArtifactInfo
from sequence_processing_pipeline.ConvertJob import ConvertJob
from sequence_processing_pipeline.FastQCJob import FastQCJob
from sequence_processing_pipeline.GenPrepFileJob import GenPrepFileJob
from sequence_processing_pipeline.AmpliconPipeline import AmpliconPipeline
from sequence_processing_pipeline.PipelineError import PipelineError
from sequence_processing_pipeline.QCJob import QCJob
import shutil


def process_amplicon_job(mapping_file_path, lane_number, qclient,
                         run_identifier, out_dir, job_id,
                         _update_current_message, skip_exec,
                         _update_job_step, job_pool_size,
                         final_results_path, success, msg, config_fp):

    # TODO: get sample_ids from mapping file for each project and confirm
    #  the values are present in Qiita.

    # Create a Pipeline object
    try:
        pipeline = AmpliconPipeline(config_fp, run_identifier,
                                    mapping_file_path, out_dir, job_id)
    except PipelineError as e:
        # Pipeline is the object that finds the input fp, based on
        # a search directory set in configuration.json and a run_id.
        if str(e).endswith("could not be found"):
            msg = f"A path for {run_identifier} could not be found."
            return False, None, msg
        elif str(e).startswith("Sample-sheet has the following errors:"):
            msg = str(e)
            _update_current_message(msg)
            raise ValueError(msg)
        else:
            raise e

    # find the uploads directory all trimmed files will need to be
    # moved to.
    results = qclient.get("/qiita_db/artifacts/types/")

    # trimmed files are stored by qiita_id. Find the qiita_id
    # associated with each project and ensure a subdirectory exists
    # for when it comes time to move the trimmed files.
    special_map = []
    for project in pipeline.get_project_info():
        upload_path = join(results['uploads'], project['qiita_id'])
        makedirs(upload_path, exist_ok=True)
        special_map.append((project['project_name'], upload_path,
                            project['qiita_id']))

    _update_current_message("Step 2 of 6: Converting BCL to fastq")

    config = pipeline.configuration['bcl-convert']
    convert_job = ConvertJob(pipeline.run_dir,
                             pipeline.output_path,
                             pipeline.sample_sheet,
                             config['queue'],
                             config['nodes'],
                             config['nprocs'],
                             config['wallclock_time_in_hours'],
                             config['per_process_memory_limit'],
                             config['executable_path'],
                             config['modules_to_load'],
                             job_id)

    # if skip_execution is True, then each Pipeline object will be
    # initialized, their assertions tested, and an ainfo will be
    # returned to the caller. However the Jobs will not actually
    # be executed. This is useful for testing.
    if not skip_exec:
        convert_job.run(callback=_update_job_step)

    _update_current_message("Step 3 of 6: Generating FastQC & "
                            "MultiQC reports")

    config = pipeline.configuration['fastqc']

    raw_fastq_files_path = join(pipeline.output_path, 'ConvertJob')

    # Fake QCJob output directory
    my_projects = AmpliconPipeline.get_project_info()
    my_projects = [x['project_name'] for x in my_projects]

    # ConvertJob should have files organized by project:
    # /$WKDIR/$QIITA_JOB_ID/ConvertJob/$PROJECTNAME_QIITAID/BLANK_XXXXX_2_1A_S97_L001_R1_001.fastq.gz

    # outputs to QCJob should be at:
    # /$WKDIR/$QIITA_JOB_ID/QCJob/$PROJECTNAME_QIITAID/fastp_reports_dir/json/BLANK_XXXXX_2_1A_S97_L001_R1_001.json
    # /$WKDIR/$QIITA_JOB_ID/QCJob/$PROJECTNAME_QIITAID/fastp_reports_dir/html/BLANK_XXXXX_2_1A_S97_L001_R1_001.html
    # /$WKDIR/$QIITA_JOB_ID/QCJob/$PROJECTNAME_QIITAID/filtered_sequences/BLANK_XXXXX_2_1A_S97_L001_R1_001.trimmed.fastq.gz
    # /$WKDIR/$QIITA_JOB_ID/QCJob/$PROJECTNAME_QIITAID/filtered_sequences/BLANK_XXXXX_2_1A_S97_L001_R2_001.trimmed.fastq.gz'

    # however, for our dummy sample-sheet, we will have just one project and
    # the results will be in a project directory but just the one. We need to
    # clone the one project directory n times, and name them after the n
    # projects in the mapping file. This is for the convertjob directory.

    # I think we may need to run QCJob on Amplicon after all, because why shouldn't
    # we be able to run fastp only on the sole fastq file?

    # we'll just want to copy into a faked QCJob folder that's 'amplicon' instead of
    # filtered_sequences or trimmed_sequences

    for some_project in my_projects:
        # copy the files from ConvertJob output to QCJob 'processed' output
        some_path = join(pipeline.output_path, 'QCJob', some_project,
                         'amplicon')
        makedirs(some_path)

        # COPY here
        src_folder = r"E:\demos\files\reports"
        dst_folder = r"E:\demos\files\account"

        # file names
        src_file = src_folder + "\profit.txt"
        dst_file = dst_folder + "\profit.txt"

        shutil.copyfile(src_file, dst_file)
        print('Copied')

    processed_fastq_files_path = join(pipeline.output_path, 'QCJob')

    fastqc_job = FastQCJob(pipeline.run_dir,
                           pipeline.output_path,
                           raw_fastq_files_path,
                           processed_fastq_files_path,
                           config['nprocs'],
                           config['nthreads'],
                           config['fastqc_executable_path'],
                           config['modules_to_load'],
                           job_id,
                           config['queue'],
                           config['nodes'],
                           config['wallclock_time_in_hours'],
                           config['job_total_memory_limit'],
                           job_pool_size,
                           config['multiqc_config_file_path'],
                           config['job_max_array_length'])

    if not skip_exec:
        fastqc_job.run(callback=_update_job_step)

    project_list = fastqc_job.project_names

    """
    _update_current_message("Step 5 of 6: Generating Prep "
                            "Information Files")

    config = pipeline.configuration['seqpro']
    gpf_job = GenPrepFileJob(
        pipeline.run_dir,
        raw_fastq_files_path,
        processed_fastq_files_path,
        pipeline.output_path,
        mapping_file_path,
        config['seqpro_path'],
        project_list,
        config['modules_to_load'],
        job_id)

    if not skip_exec:
        gpf_job.run(callback=_update_job_step)
    """

    ainfo = [
        ArtifactInfo('output', 'job-output-folder',
                     [(f'{final_results_path}/', 'directory')])
    ]

    return success, ainfo, msg