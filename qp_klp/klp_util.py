from os import walk
from os.path import join
import pandas as pd


def map_sample_names_to_tube_ids(sn_tid_map_by_proj, output_dir):
    prep_files = []
    for root, dirs, files in walk(output_dir):
        for prep_file in files:
            # sanity check the os.walk results
            if prep_file.endswith('.tsv'):
                # store the full path to the prep-file.
                prep_files.append(join(root, prep_file))

    results = {}

    for project in sn_tid_map_by_proj:
        if sn_tid_map_by_proj[project] is not None:
            # this project has tube-ids registered in Qiita.
            # find the prep-file associated with this project.
            results[project] = {}
            for prep_file in prep_files:
                # not the best check but good enough for now.
                if project in prep_file:
                    df = pd.read_csv(prep_file, sep='\t',
                                     dtype=str, index_col=False)
                    # save a copy of sample_name column as 'old_sample_name'
                    df['old_sample_name'] = df['sample_name']
                    for i in df.index:
                        sample_name = df.at[i, "sample_name"]
                        if not sample_name.startswith('BLANK'):
                            # remove any leading zeroes if they exist
                            sample_name = sample_name.lstrip('0')
                            sample_name = sn_tid_map_by_proj[project][
                                sample_name]
                            df.at[i, "sample_name"] = sample_name
                    results[project][prep_file] = df
                    break

    return results


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
