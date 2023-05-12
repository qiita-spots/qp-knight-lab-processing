import pandas as pd
from json import dumps
from metapool import KLSampleSheet
from os import makedirs
from os.path import join
from sequence_processing_pipeline.Pipeline import Pipeline
from sequence_processing_pipeline.PipelineError import PipelineError
from collections import defaultdict


def update_blanks_in_qiita(sifs, qclient):
    for sif_path in sifs:
        # get study_id from sif_file_name ...something_14385_blanks.tsv
        study_id = sif_path.split('_')[-2]

        df = pd.read_csv(sif_path, delimiter='\t', dtype=str)

        # Prepend study_id to make them compatible w/list from Qiita.
        df['sample_name'] = f'{study_id}.' + df['sample_name'].astype(str)

        # SIFs only contain BLANKs. Get the list of potentially new BLANKs.
        blank_ids = [i for i in df['sample_name'] if 'blank' in i.lower()]
        blanks = df[df['sample_name'].isin(blank_ids)]['sample_name']
        if len(blanks) == 0:
            # we have nothing to do so let's return early
            return

        # Get list of BLANKs already registered in Qiita.
        from_qiita = qclient.get(f'/api/v1/study/{study_id}/samples')
        from_qiita = [x for x in from_qiita if
                      x.startswith(f'{study_id}.BLANK')]

        # Generate list of BLANKs that need to be ADDED to Qiita.
        new_blanks = (set(blanks) | set(from_qiita)) - set(from_qiita)

        if len(new_blanks):
            # Generate dummy entries for each new BLANK, if any.
            categories = qclient.get(f'/api/v1/study/{study_id}/samples/'
                                     'info')['categories']

            # initialize payload w/required dummy categories
            data = {i: {c: 1 for c in categories} for i in new_blanks}

            # populate payload w/additional columns and/or overwrite existing
            # columns w/metadata from SIF file.
            sif_data = df.set_index('sample_name').T.to_dict()
            for new_blank in new_blanks:
                for column in sif_data[new_blank]:
                    data[new_blank][column] = sif_data[new_blank][column]

            # http_patch will raise Error if insert failed.
            qclient.http_patch(f'/api/v1/study/{study_id}/samples',
                               data=dumps(data))

            return data


def map_sample_names_to_tube_ids(prep_info_file_paths, sn_tid_map_by_proj):
    for proj in sn_tid_map_by_proj:
        if sn_tid_map_by_proj[proj] is not None:
            # this project has tube-ids registered in Qiita.
            # find the prep-file associated with this project.
            for prep_file in prep_info_file_paths:
                # not the best check but good enough for now.
                if proj in prep_file:
                    df = pd.read_csv(prep_file, sep='\t',
                                     dtype=str, index_col=False)
                    # save a copy of sample_name column as 'old_sample_name'
                    df['old_sample_name'] = df['sample_name']
                    for i in df.index:
                        smpl_name = df.at[i, "sample_name"]
                        if not smpl_name.startswith('BLANK'):
                            # remove any leading zeroes if they exist
                            smpl_name = smpl_name.lstrip('0')
                            if smpl_name in sn_tid_map_by_proj[proj]:
                                tube_id = sn_tid_map_by_proj[proj][smpl_name]
                                df.at[i, "sample_name"] = tube_id
                    df.to_csv(prep_file, index=False, sep="\t")


class StatusUpdate():
    def __init__(self, qclient, job_id, step_count):
        self.qclient = qclient
        self.job_id = job_id
        self.msg = ''
        self.current_step = 0
        self.step_count = step_count

    def update_job_status(self, status, id):
        # internal function implements a callback function for Pipeline.run().
        # :param id: PBS/Torque/or some other informative and current job id.
        # :param status: status message
        self.qclient.update_job_step(self.job_id,
                                       self.msg + f" ({id}: {status})")

    def update_current_message(self, msg, include_step=True):
        # internal function that sets current_message to the new value before
        # updating the job step in the UI.
        if include_step:
            self.current_step += 1
            self.msg = f"Step {self.current_step} of {self.step_count}: {msg}"
        else:
            self.msg = msg

        self.qclient.update_job_step(self.job_id, self.msg)


class FailedSamplesRecord:
    def __init__(self, output_dir, samples):
        # because we want to write out the list of samples that failed after
        # each Job is run, and we want to organize that output by project, we
        # need to keep a running state of failed samples, and reuse the method
        # to reorganize the running-results and write them out to disk.

        self.output_path = join(output_dir, 'failed_samples.html')

        # create an initial dictionary with sample-ids as keys and their
        # associated project-name and status as values. Afterwards, we'll
        # filter out the sample-ids w/no status (meaning they were
        # successfully processed) before writing the failed entries out to
        # file.
        self.failed = {x.Sample_ID: [x.Sample_Project, None] for x in samples}

    def write(self, failed_ids, job_name):
        for failed_id in failed_ids:
            # as a rule, if a failed_id were to appear in more than one
            # audit(), preserve the earliest failure, rather than the
            # latest one.
            if self.failed[failed_id][1] is None:
                self.failed[failed_id][1] = job_name

        # filter out the sample-ids w/out a failure status
        filtered_fails = {x: self.failed[x] for x in self.failed if
                          self.failed[x][1] is not None}

        data = []
        for sample_id in filtered_fails:
            project_name = filtered_fails[sample_id][0]
            failed_at = filtered_fails[sample_id][1]
            data.append({'Project': project_name, 'Sample ID': sample_id,
                         'Failed at': failed_at})
        df = pd.DataFrame(data)

        with open(self.output_path, 'w') as f:
            f.write(df.to_html(border=2, index=False, justify="left",
                               render_links=True, escape=False))


def parse_prep_file(prep_file_path):
    metadata = pd.read_csv(prep_file_path,
                           dtype=str,
                           delimiter='\t',
                           # forces Pandas to not make the first column the
                           # index even when the values appear numeric.
                           index_col=False)

    if metadata is None:
        raise ValueError(f"{prep_file_path} does not exist.")

    metadata.set_index('sample_name', inplace=True)

    # convert to standard dictionary.
    return metadata.to_dict('index')


def update_sample_sheet(sample_sheet_path, lane_number):
    # use KLSampleSheet functionality to add/overwrite lane number.
    sheet = KLSampleSheet(sample_sheet_path)
    for sample in sheet:
        sample['Lane'] = f'{lane_number}'

    with open(sample_sheet_path, 'w') as f:
        sheet.write(f)

def generate_special_map(results, projects):
    # this function should be able to be tested by passing in simulated =
    # results from qclient.

    # trimmed files are stored by qiita_id. Find the qiita_id
    # associated with each project and ensure a subdirectory exists
    # for when it comes time to move the trimmed files.
    special_map = []
    for project in projects:
        upload_path = join(results['uploads'], project['qiita_id'])
        makedirs(upload_path, exist_ok=True)
        special_map.append((project['project_name'], upload_path,
                            project['qiita_id']))

    return special_map

def generate_pipeline(pipeline_type, input_file_path, lane_number,config_fp,
                      run_identifier, out_dir, job_id):
    if pipeline_type in ['metagenomic', 'metatranscriptomic']:
        update_sample_sheet(input_file_path, lane_number)
        return Pipeline(config_fp, run_identifier, input_file_path, None,
                        out_dir, job_id, pipeline_type)
    elif pipeline_type == 'amplicon':
        return Pipeline(config_fp, run_identifier, None, input_file_path,
                        out_dir, job_id, pipeline_type)
    else:
        raise PipelineError(f"'{pipeline_type}' is not a valid Pipeline type.")

def get_data_type(pipeline_type, target_gene=None):
    if pipeline_type in ['metagenomic', 'metatranscriptomic']:
        return pipeline_type
    elif pipeline_type == 'amplicon':
        if target_gene:
            for key in {'16S', '18S', 'ITS'}:
                if key in target_gene:
                    return key
        else:
            raise ValueError("target_gene must be specified for amplicon type")
    else:
        raise ValueError(f"'{pipeline_type}' is not a valid pipeline type")


def update_prep_templates(qclient, prep_file_paths):
    '''
    Update prep-template info in Qiita. Get breakdown of prep-ids by study-id.
    :param qclient:
    :param prep_file_paths:
    :return: A dict of lists of prep-ids, keyed by study-id.
    '''
    results = defaultdict(list)

    for study_id in prep_file_paths:
        for prep_file_path in prep_file_paths[study_id]:
            metadata = parse_prep_file(prep_file_path)
            target_gene = metadata[list(metadata.keys())[0]]['target_gene']

            data = {'prep_info': dumps(metadata),
                    'study': study_id,
                    'data_type': get_data_type(target_gene)}

            reply = qclient.post('/qiita_db/prep_template/', data=data)
            prep_id = reply['prep']
            results[study_id].append(prep_id)

    return results

def get_registered_samples_in_qiita(qclient, qiita_id):
    '''
    Obtain lists for sample-names and tube-ids registered in Qiita.
    :param qclient: QiitaClient object
    :param qiita_id: Qiita ID for the project in question.
    :return: a tuple of lists, one for sample-names and another for tube-ids.
    '''
    samples = qclient.get(f'/api/v1/study/{qiita_id}/samples')

    # remove Qiita ID as a prefix from the sample-names.
    samples = {x.replace(f'{qiita_id}.', '') for x in samples}

    # find out if tube-ids are registered in the study.
    categories = qclient.get(f'/api/v1/study/{qiita_id}'
                             '/samples/info')['categories']

    if 'tube_id' in categories:
        tids = qclient.get(f'/api/v1/study/{qiita_id}/samples/'
                           'categories=tube_id')['samples']
    else:
        tids = None

    return (samples, tids)
