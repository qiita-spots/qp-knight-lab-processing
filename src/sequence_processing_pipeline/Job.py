from jinja2 import BaseLoader, TemplateNotFound
from os.path import getmtime
import pathlib
from itertools import zip_longest
from os import makedirs, walk
from os.path import basename, exists, split, join
from sequence_processing_pipeline.PipelineError import (PipelineError,
                                                        JobFailedError,
                                                        ExecFailedError)
from subprocess import Popen, PIPE
from time import sleep
import logging
from inspect import stack
import re
from collections import Counter
from glob import glob


# taken from https://jinja.palletsprojects.com/en/3.0.x/api/#jinja2.BaseLoader
class KISSLoader(BaseLoader):
    def __init__(self, path):
        # pin the path for loader to the location sequence_processing_pipeline
        # (the location of this file), along w/the relative path to the
        # templates directory.
        self.path = join(pathlib.Path(__file__).parent.resolve(), path)

    def get_source(self, environment, template):
        path = join(self.path, template)
        if not exists(path):
            raise TemplateNotFound(template)
        mtime = getmtime(path)
        with open(path) as f:
            source = f.read()
        return source, path, lambda: mtime == getmtime(path)


class Job:
    slurm_status_terminated = ['BOOT_FAIL', 'CANCELLED', 'DEADLINE', 'FAILED',
                               'NODE_FAIL', 'OUT_OF_MEMORY', 'PREEMPTED',
                               'REVOKED', 'TIMEOUT']

    slurm_status_successful = ['COMPLETED']

    slurm_status_running = ['COMPLETING', 'CONFIGURING', 'PENDING', 'REQUEUED',
                            'REQUEUE_FED', 'REQUEUE_HOLD', 'RESIZING',
                            'RESV_DEL_HOLD', 'RUNNING', 'SIGNALING',
                            'SPECIAL_EXIT', 'STAGE_OUT', 'STOPPED',
                            'SUSPENDED']

    slurm_status_not_running = (slurm_status_terminated +
                                slurm_status_successful)

    slurm_status_all_states = (slurm_status_terminated +
                               slurm_status_successful +
                               slurm_status_running)

    polling_interval_in_seconds = 60
    squeue_retry_in_seconds = 10

    def __init__(self, root_dir, output_path, job_name, executable_paths,
                 max_array_length, modules_to_load=None):
        """
        Base-class to implement Jobs from.
        :param job_name: A name for the job. Used to create log files.
        :param root_dir: The path to a Job's root input directory.
        :param output_path: The root path to store all job products.
        :param executable_paths: A list of executables to validate.
        :param modules_to_load: A list of modules to load before validation.
        """
        self.job_name = job_name
        self.root_dir = root_dir
        self._directory_check(self.root_dir, create=False)
        self.force_job_fail = False

        self.output_path = join(output_path, self.job_name)
        self._directory_check(self.output_path, create=True)

        self.log_path = join(self.output_path, 'logs')
        self._directory_check(self.log_path, create=True)

        self.modules_to_load = modules_to_load
        self.max_array_length = max_array_length

        self.script_count = 0

        self.bypass_exec_check = ['bcl-convert']

        self.suffix = None

        # checking if this is running as part of the unittest
        # https://stackoverflow.com/a/25025987
        self.is_test = True if [
            x for x in stack() if 'unittest' in x.filename] else False

        self.audit_folders = None

        # For each executable in the list, get its filename and use _which()
        # to see if it can be found. Directly pass an optional list of modules
        # to load before-hand, so that the binary can be found.
        # If the executable can't be found or doesn't have the same path as
        # the version given, raise a PipelineError.
        for executable_path in executable_paths:
            file_path, file_name = split(executable_path)

            # bcl-convert module is not installed on the node this test gets
            # run on, so forego it entirely.

            # Some modules such as bcl-convert are not available to the
            # servers running this check. However they are still available
            # on Barnacle nodes and this is what's important. For now simply
            # bypass the check for known situations.
            for name in self.bypass_exec_check:
                if name in file_name:
                    continue

            # No need to test results. _which() will raise a PipelineError if
            # file_name is a path and the path found does not match. It will
            # also raise a PipelineError if the file could not be found.
            if not self.is_test:
                self._which(file_name, modules_to_load=self.modules_to_load)

    def run(self):
        """
        Since a Job object can encapsulate one or more submit_job() or system()
        calls, the base run() method remains unimplemented. It is the job of
        each sub-class to define what needs to be run in order to generate the
        expected output.
        """
        raise PipelineError("Base class run() method not implemented.")

    def mark_job_completed(self):
        with open(join(self.output_path, 'job_completed'), 'w') as f:
            f.write("job_completed")

    def mark_post_processing_completed(self):
        with open(join(self.output_path,
                       'post_processing_completed'), 'w') as f:
            f.write("post_processing_completed")

    def parse_logs(self):
        # by default, look for anything to parse in the logs directory.
        log_path = join(self.output_path, 'logs')
        files = sorted(glob(join(log_path, '*')))
        msgs = []

        for some_file in files:
            with open(some_file, 'r') as f:
                msgs += [line for line in f.readlines()
                         if 'error:' in line.lower()]

        return [msg.strip() for msg in msgs]

    def _which(self, file_path, modules_to_load=None):
        """
        Returns file_path if file_path exists and file_path is a full path.
        Otherwise returns the path to a file named 'file_path' found in PATH.
        :param file_path: The path of the executable to find.
        :param modules_to_load: A list of Linux module names to load.
        :return: A path to 'file_name'.
        """
        tmp = split(file_path)
        # remove any elements that are empty string.
        tmp = [x for x in tmp if x]

        isPath = True if len(tmp) > 1 else False

        cmd = 'which ' + file_path

        if modules_to_load:
            cmd = 'module load ' + ' '.join(modules_to_load) + ';' + cmd

        results = self._system_call(cmd)
        result = results['stdout'].strip()

        if not result:
            raise PipelineError("File '%s' does not exist." % file_path)

        if isPath is True and result != file_path:
            raise PipelineError(f"Found path '{result} does not match "
                                f"{file_path}")

        return result

    def _file_check(self, file_path):
        if exists(file_path):
            logging.debug("file '%s' exists." % file_path)
            return True
        else:
            raise PipelineError("file '%s' does not exist." % file_path)

    def _find_files(self, search_path):
        lst = []
        for root, dirs, files in walk(search_path):
            lst += [join(root, x) for x in files]
        return lst

    def _directory_check(self, directory_path, create=False):
        if exists(directory_path):
            logging.debug("directory '%s' exists." % directory_path)
        else:
            if create:
                try:
                    makedirs(directory_path, exist_ok=True)
                except OSError as e:
                    # this is a known potential error. Re-raise it as a
                    # PipelineError, so it gets handled in the same location
                    # as the others.
                    raise PipelineError(str(e))
            else:
                raise PipelineError(
                    "directory_path '%s' does not exist." % directory_path)

    def _system_call(self, cmd, allow_return_codes=[], callback=None):
        """
        Call command and return (stdout, stderr, return_value)
        :param cmd: The string containing the command to be run, or a sequence
                    of strings that are the tokens of the command.
        :param allow_return_codes: optional user-defined list of successful
                    return codes in addition to 0.
        :param callback: optional function taking two parameters (id, status)
                         that is called when a running process's status is
                         changed.
        :return: a dictionary containing stdout, stderr, and return_code as
                 key/value pairs.
        """
        proc = Popen(cmd, universal_newlines=True, shell=True,
                     stdout=PIPE, stderr=PIPE)

        if callback is not None:
            callback(jid=proc.pid, status='RUNNING')

        # Communicate pulls all stdout/stderr from the PIPEs
        # This call blocks until the command is done
        stdout, stderr = proc.communicate()
        return_code = proc.returncode

        logging.debug("stdout: %s" % stdout)
        logging.debug("stderr: %s" % stderr)
        logging.debug("return code: %s" % return_code)

        acceptable_return_codes = [0] + allow_return_codes

        if return_code not in acceptable_return_codes:
            if callback is not None:
                callback(jid=proc.pid, status='ERROR')
            msg = (
                'Execute command-line statement failure:\n'
                f'Command: {cmd}\n'
                f'return code: {return_code}\n'
                f'stdout: {stdout}\n'
                f'stderr: {stderr}\n')
            logging.error(msg)
            raise ExecFailedError(message=msg)

        if callback is not None:
            callback(jid=proc.pid, status='COMPLETED')

        return {'stdout': stdout, 'stderr': stderr, 'return_code': return_code}

    def _query_slurm(self, job_ids):
        # query_slurm encapsulates the handling of squeue.
        count = 0
        while True:
            result = self._system_call("squeue -t all -j "
                                       f"{','.join(job_ids)} "
                                       "-o '%i,%T'")

            if result['return_code'] == 0:
                # there was no issue w/squeue, break this loop and
                # continue.
                break
            else:
                # there was likely an intermittent issue w/squeue. Pause
                # and wait before trying a few more times. If the problem
                # persists then report the error and exit.
                count += 1

                if count > 3:
                    raise ExecFailedError(result['stderr'])

                sleep(Job.squeue_retry_in_seconds)

        lines = result['stdout'].split('\n')
        lines.pop(0)  # remove header
        lines = [x.split(',') for x in lines if x != '']

        jobs = {}
        for job_id, state in lines:
            # ensure unique_id is of type string for downstream use.
            job_id = str(job_id)
            jobs[job_id] = state

        return jobs

    def wait_on_job_ids(self, job_ids, callback=None):
        '''
        Wait for the given job-ids to finish running before returning.
        :param job_ids: A list of Slurm job-ids
        :param callback: Set callback function that receives status updates.
        :return: A dictionary of job-ids and their current statuses.
        '''

        # wait_on_job_ids was broken out of submit_job() and updated to monitor
        # multiple job ids. This will allow multiple jobs to be submitted to
        # Slurm in parallel and a single wait_on_job_ids() can wait on all of
        # them before returning, optionally submitting callbacks for each
        # job-id.

        # ensure all ids are strings to ensure proper working w/join().
        job_ids = [str(x) for x in job_ids]

        while True:
            # Because query_slurm only returns state on the job-ids we specify,
            # the wait process is a simple check to see whether any of the
            # states are 'running' states or not.
            jobs = self._query_slurm(job_ids)

            # jobs will be a dict of job-ids or array-ids for jobs that
            # are array-jobs. the value of jobs[id] will be a state e.g.:
            # 'RUNNING', 'FAILED', 'COMPLETED'.
            states = [jobs[x] in Job.slurm_status_not_running for x in jobs]

            if set(states) == {True}:
                # if all the states are either FAILED or COMPLETED
                # then the set of those states no matter how many
                # array-jobs there were will ultimately be the set of
                # {True}. If not then that means there are still jobs
                # that are running.
                break

            logging.debug(f"sleeping {Job.polling_interval_in_seconds} "
                          "seconds...")
            sleep(Job.polling_interval_in_seconds)

        return jobs

    def submit_job(self, script_path, job_parameters=None,
                   script_parameters=None, wait=True,
                   exec_from=None, callback=None):
        """
        Submit a Slurm job script and optionally wait for it to finish.
        :param script_path: The path to a Slurm job (bash) script.
        :param job_parameters: Optional parameters for scheduler submission.
        :param script_parameters: Optional parameters for your job script.
        :param wait: Set to False to submit job and not wait.
        :param exec_from: Set working directory to execute command from.
        :param callback: Set callback function that receives status updates.
        :return: If wait is True, a dictionary containing the job's id and
                 status. If wait is False, the Slurm job-id of the submitted
                 job. Raises PipelineError if job could not be submitted or if
                 job was unsuccessful.
        """
        if job_parameters:
            cmd = 'sbatch %s %s' % (job_parameters, script_path)
        else:
            cmd = 'sbatch %s' % (script_path)

        if script_parameters:
            cmd += ' %s' % script_parameters

        if exec_from:
            cmd = f'cd {exec_from};' + cmd

        logging.debug("job scheduler call: %s" % cmd)

        if self.force_job_fail:
            raise JobFailedError("This job died.")

        # if system_call does not raise a PipelineError(), then the scheduler
        # successfully submitted the job. In this case, it should return
        # the id of the job in stdout.
        results = self._system_call(cmd)
        stdout = results['stdout']

        job_id = stdout.strip().split()[-1]

        # Just to give some time for everything to be set up properly
        sleep(10)

        if wait is False:
            # return job_id since that is the only information for this new
            # job that we have available. User should expect that this is
            # not a dict if they explicitly set wait=False.
            return job_id

        # the user is expecting a dict with 'job_id' and 'job_state'
        # attributes. This method will return a dict w/job_ids as keys and
        # their job status as values. This must be munged before returning
        # to the user.
        results = self.wait_on_job_ids([job_id], callback=callback)

        if job_id in results:
            # job is a non-array job
            job_result = {'job_id': job_id, 'job_state': results[job_id]}
        else:
            # job is an array job
            # assume all array jobs in this case will be associated w/job_id.
            counts = Counter()
            for array_id in results:
                counts[results[array_id]] += 1

            # for array jobs we won't be returning a string representing the
            # state of a single job. Instead we're returning a dictionary of
            # the number of unique states the set of array-jobs ended up in and
            # the number for each one.
            job_result = {'job_id': job_id, 'job_state': dict(counts)}

        if callback is not None:
            if isinstance(job_result['job_state'], dict):
                # this is an array job
                states = []
                for key in counts:
                    states.append(f"{key}: {counts[key]}")

                callback(jid=job_id, status=", ".join(states))

            else:
                # this is a standard job
                callback(jid=job_id, status=job_result['job_state'])

        if isinstance(job_result['job_state'], dict):
            states = list(job_result['job_state'].keys())
            if states == ['COMPLETED']:
                return job_result
            else:
                raise JobFailedError(f"job {job_id} exited with jobs in the "
                                     f"following states: {', '.join(states)}")
        else:
            if job_result['job_state'] == 'COMPLETED':
                return job_result
            else:
                raise JobFailedError(f"job {job_id} exited with status "
                                     f"{job_result['job_state']}")

    def _group_commands(self, cmds):
        # break list of commands into chunks of max_array_length (Typically
        # 1000 for Slurm job arrays). To ensure job arrays are never more
        # than 1000 jobs long, we'll chain additional commands together, and
        # evenly distribute them amongst the first 1000.
        cmds.sort()
        chunks = [cmds[i:i + self.max_array_length] for i in
                  range(0, len(cmds), self.max_array_length)]

        results = []

        # create a chained command by taking one command from each list.
        # zip_longest() allows us to handle lists of different lengths, as the
        # last chunk will always be of different length than 1000.
        for tuple in zip_longest(*chunks):
            # zip_longest() pads shorter lists with None. In our case, we
            # don't want an additional command named 'None'.
            chained_cmd = [x for x in list(tuple) if x is not None]
            chained_cmd = ';'.join(chained_cmd)
            results.append(chained_cmd)

        return results

    # assume for now that a corresponding zip file exists for each html
    # file found. Assume for now that all html files will be found in a
    # 'filtered_sequences' or 'trimmed_sequences' subdirectory.
    #
    # verify that the entire list of sample-ids found match what's
    # expected. Since the list of expected ids is very small, we'll
    # perform an exact comparison.

    def audit(self, sample_ids):
        """
        Audit the results of a run.
        :param sample_ids: A list of sample-ids that require results.
        :return: A list of sample-ids that were not found.
        """
        files_found = []

        if self.suffix is None:
            raise PipelineError("Audit() method called on base Job object.")

        for root, dirs, files in walk(self.output_path):
            if 'zero_files' in root:
                continue
            if self.audit_folders is not None:
                # let's check that any of the audit_folders is in root
                if not [f for f in self.audit_folders if f in root]:
                    continue
            files_found += [join(root, x) for x in files if
                            x.endswith(self.suffix)]

        found = []
        for sample_id in sample_ids:
            for found_file in files_found:
                # the trailing underscore is important as it can be assumed
                # that all fastq.gz files will begin with sample_id followed
                # by an '_', and then one or more additional parameters
                # separated by underscores. This substring is unlikely to be
                if basename(found_file).startswith('%s_' % sample_id):
                    found.append(sample_id)
                    break

        return sorted(list(set(found) ^ set(sample_ids)))

    def _toggle_force_job_fail(self):
        if self.force_job_fail is True:
            self.force_job_fail = False
        else:
            self.force_job_fail = True

        return self.force_job_fail

    def extract_project_names_from_fastq_dir(self, path_to_fastq_dir):
        '''
        Returns list of <project_qiita-ids> based on a fastq path.
        :param path_to_fastq_dir: Path to a nested dir containing fastq files.
        :return: list of strings of the form PROJECT-NAME_QIITA-ID.
        '''
        tmp = []

        # given a nested directory containing fastq.gz files from any part
        # of the SPP, get a list of paths to just the fastq files.
        for root, dirs, files in walk(path_to_fastq_dir):
            for file_name in files:
                if file_name.endswith('.fastq.gz'):
                    if not file_name.startswith('Undetermined'):
                        # do not include the Undetermined files found in
                        # ConvertJob, because they do not contain a project
                        # name in their path, by definition.
                        # don't record the full_path w/filename, as we're
                        # not interested in the files themselves, just their
                        # basedirs().
                        tmp.append(root)

        # break up the path into a set of unique directory names. at least
        # some of these will be of the form PROJECT-NAME_QIITA-ID. Flatten
        # out the list of lists into a single list.
        tmp = [some_root.split('/') for some_root in tmp]
        tmp = [dir_name for sub_list in tmp for dir_name in sub_list]

        # remove any duplicate entries from the list.
        tmp = list(set(tmp))

        # file-paths for fastq.gz files are going to include the names of
        # standard directories e.g. 'ConvertJob', 'trimmed_sequences' as well
        # as directories named PROJECT-NAME_QIITA-ID. Base Job() class is not
        # aware of all the possibilities in the path supplied by the user.
        # Hence, this method filters for names that fit the above format.
        results = []
        for dir_name in tmp:
            m = re.match(r'^(\w[\w\-_]*_\d{5})$', dir_name)
            if m:
                results.append(m.groups(0)[0])

        return sorted(results)
