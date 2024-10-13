from json import dumps, load
from os.path import join, exists
import pandas as pd


class FailedSamplesRecord:
    def __init__(self, output_dir, samples):
        # because we want to write out the list of samples that failed after
        # each Job is run, and we want to organize that output by project, we
        # need to keep a running state of failed samples, and reuse the method
        # to reorganize the running-results and write them out to disk.
        self.output_path = join(output_dir, 'failed_samples.json')
        self.report_path = join(output_dir, 'failed_samples.html')

        # create an initial dictionary with sample-ids as keys and their
        # associated project-name and status as values. Afterwards, we'll
        # filter out the sample-ids w/no status (meaning they were
        # successfully processed) before writing the failed entries out to
        # file.
        self.sample_state = {x.Sample_ID: None for x in samples}
        self.project_map = {x.Sample_ID: x.Sample_Project for x in samples}

    def dump(self):
        output = {'sample_state': self.sample_state,
                  'project_map': self.project_map}

        with open(self.output_path, 'w') as f:
            f.write(dumps(output, indent=2, sort_keys=True))

    def load(self):
        # if recorded state exists, overwrite initial state.
        if exists(self.output_path):
            with open(self.output_path, 'r') as f:
                state = load(f)

            self.sample_state = state['sample_state']
            self.project_map = state['project_map']

    def update(self, failed_ids, job_name):
        # as a rule, if a failed_id were to appear in more than one
        # audit(), preserve the earliest failure, rather than the
        # latest one.
        for failed_id in failed_ids:
            if self.sample_state[failed_id] is None:
                self.sample_state[failed_id] = job_name

    def write(self, failed_ids, job_name):
        # a convenience method to support legacy behavior.
        # specifically, reload recorded state, if it exists.
        # then update state before recording to file.
        self.load()
        self.update(failed_ids, job_name)
        self.dump()

    def generate_report(self):
        # filter out the sample-ids w/out a failure status
        filtered_fails = {x: self.sample_state[x] for x in self.sample_state if
                          self.sample_state[x] is not None}

        data = []
        for sample_id in filtered_fails:
            data.append({'Project': filtered_fails[sample_id],
                         'Sample ID': sample_id,
                         'Failed at': self.project_map[sample_id]
                         })
        df = pd.DataFrame(data)

        with open(self.report_path, 'w') as f:
            f.write(df.to_html(border=2, index=False, justify="left",
                               render_links=True, escape=False))
