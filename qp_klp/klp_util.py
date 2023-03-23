from os.path import join
import pandas as pd
from json import dumps


def update_blanks_in_qiita(sifs, qclient):
    results = {}
    for sif_path in sifs:
        results[sif_path] = {'sif_path': sif_path}

        # get study_id from sif_file_name ...something_14385_blanks.tsv
        study_id = sif_path.split('_')[-2]
        results[sif_path]['study_id'] = study_id

        # SIFs only contain BLANKs. Get the list of potentially new BLANKs.
        blanks = pd.read_csv(sif_path, delimiter='\t')['sample_name']

        # Prepend study_id to make them compatible w/list from Qiita.
        blanks = [f'{study_id}.{x}' for x in blanks]

        results[sif_path]['blanks'] = blanks

        # Get list of BLANKs already registered in Qiita.
        from_qiita = qclient.get(f'/api/v1/study/{study_id}/samples')
        from_qiita = [x for x in from_qiita if
                      x.startswith(f'{study_id}.BLANK')]

        results[sif_path]['from_qiita'] = from_qiita

        # Generate list of BLANKs that need to be ADDED to Qiita.
        new_blanks = (set(blanks) | set(from_qiita)) - set(from_qiita)

        results[sif_path]['new_blanks'] = new_blanks

        if len(new_blanks):
            # Generate dummy entries for each new BLANK, if any.
            categories = qclient.get(f'/api/v1/study/{study_id}/samples/'
                                     'info')['categories']
            results[sif_path]['categories'] = categories
            data = {i: {c: 1 for c in categories} for i in new_blanks}
            results[sif_path]['data'] = data
            # http_patch will raise Error if insert failed.
            qclient.http_patch(f'/api/v1/study/{study_id}/samples',
                               data=dumps(data))
    return results


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
    def __init__(self, qclient, job_id):
        self.qclient = qclient
        self.job_id = job_id
        self.msg = ''

    def update_job_step(self, status, id):
        # internal function implements a callback function for Pipeline.run().
        # :param id: PBS/Torque/or some other informative and current job id.
        # :param status: status message
        self.qclient.update_job_step(self.job_id,
                                     self.msg + f" ({id}: {status})")

    def update_current_message(self, msg):
        # internal function that sets current_message to the new value before
        # updating the job step in the UI.
        self.msg = msg
        self.qclient.update_job_step(self.job_id, msg)


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
