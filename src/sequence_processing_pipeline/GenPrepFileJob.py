from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import PipelineError
from os import makedirs, symlink
from os.path import isdir, join, exists, basename
from shutil import copy, copytree
from functools import partial
from collections import defaultdict
from metapool import (demux_sample_sheet, parse_prep,
                      demux_pre_prep, pre_prep_needs_demuxing,
                      sheet_needs_demuxing, load_sample_sheet)


class GenPrepFileJob(Job):
    def __init__(self, run_dir, convert_job_path, qc_job_path, output_path,
                 input_file_path, seqpro_path, modules_to_load,
                 qiita_job_id, reports_path, is_amplicon=False):

        super().__init__(run_dir,
                         output_path,
                         'GenPrepFileJob',
                         [seqpro_path],
                         1000,
                         modules_to_load=modules_to_load)

        self.run_id = basename(run_dir)
        self.input_file_path = input_file_path
        self.seqpro_path = seqpro_path
        self.qiita_job_id = qiita_job_id
        self.is_amplicon = is_amplicon
        self.prep_file_paths = None
        self.commands = []
        self.has_replicates = False
        self.replicate_count = 0

        self.reports_path = reports_path

        # make the 'root' of the 'run_directory'. on restarts it will exist
        # already.
        makedirs(join(self.output_path, self.run_id), exist_ok=True)

        # This directory will already exist on restarts, hence avoid
        # copying. To support legacy seqpro, We will copy the single file
        # seqpro needs into a clean sub-directory named 'Reports'. This can
        # be fixed when seqpro is refactored.
        reports_dir = join(self.output_path, self.run_id, 'Reports')

        # handle reports_path being either a directory or a file.
        if exists(reports_dir):
            self.is_restart = True
        else:
            self.is_restart = False

            if isdir(self.reports_path):
                copytree(self.reports_path, reports_dir)
            elif not self.is_amplicon:
                # assume self.reports_path is a file.
                makedirs(reports_dir)
                copy(self.reports_path, reports_dir)

        # extracting from either convert_job_path or qc_job_path should
        # produce equal results.
        projects = self.extract_project_names_from_fastq_dir(qc_job_path)

        for project in projects:
            src_path = partial(join, qc_job_path, project)
            filtered_seq_dir = src_path('filtered_sequences')
            trimmed_seq_dir = src_path('trimmed_sequences')
            fastp_rept_dir = src_path('fastp_reports_dir', 'json')
            amplicon_seq_dir = join(convert_job_path, project)

            dst = join(self.output_path, self.run_id, project)

            if not self.is_restart:
                # these will already be created if restarted.
                if self.is_amplicon:
                    if exists(amplicon_seq_dir):
                        makedirs(dst, exist_ok=True)
                        symlink(amplicon_seq_dir, join(dst, 'amplicon'))
                else:
                    if exists(filtered_seq_dir):
                        makedirs(dst, exist_ok=True)
                        symlink(filtered_seq_dir, join(dst,
                                                       'filtered_sequences'))

                    if exists(trimmed_seq_dir):
                        makedirs(dst, exist_ok=True)
                        symlink(trimmed_seq_dir, join(dst,
                                                      'trimmed_sequences'))

                    if exists(fastp_rept_dir):
                        makedirs(dst, exist_ok=True)
                        symlink(fastp_rept_dir, join(dst, 'json'))

        # seqpro usage:
        # seqpro path/to/run_dir path/to/sample/sheet /path/to/fresh/output_dir

        # note that seqpro takes in a sample-sheet as input. Unlike other
        # Jobs that process the sample-sheet, validate parameters, and ensure
        # that output is segregated by project name, seqpro will be doing that
        # for us. A single call to seqpro will generate n output files, one
        # for each project described in the sample-sheet's Bioinformatics
        # heading.

        # by default, set file_paths to the default:
        file_paths = [self.input_file_path]

        if self.is_amplicon:
            # parse_prep extended to support parsing pre-prep files as well.
            fp = parse_prep(self.input_file_path)
            if pre_prep_needs_demuxing(fp):
                self.has_replicates = True

                # overwrite default setting
                file_paths = self._write_to_file(demux_pre_prep(fp))
        else:
            fp = load_sample_sheet(self.input_file_path)
            if sheet_needs_demuxing(fp):
                self.has_replicates = True

                # overwrite default setting
                file_paths = self._write_to_file(demux_sample_sheet(fp))

        for fp in file_paths:
            # generate a seqpro command-line using the new sample-sheet.
            if self.has_replicates:
                self.replicate_count += 1
                out_path = join(self.output_path,
                                'PrepFiles',
                                str(self.replicate_count))
            else:
                out_path = join(self.output_path, 'PrepFiles')

            self.commands.append([self.seqpro_path, '--verbose',
                                  join(self.output_path, self.run_id),
                                  f'"{fp}"',
                                  out_path])

    def _write_to_file(self, demuxed):
        '''
        Saves the new plate-replicate-specific sample-sheet or pre-prep file
        w/a unique name. Assume demuxed is a list of DataFrames originating
        from a single sample-sheet or pre-prep file.
        :param demuxed:
        :return:
        '''
        results = []
        for count, replicate in enumerate(demuxed, 1):
            if self.is_amplicon:
                replicate['sample_name'] = replicate['orig_name']
                fp = join(self.output_path, f"replicate_sheet_{count}.txt")
                replicate.to_csv(fp, sep='\t', index=False, header=True)
                results.append(fp)
            else:
                fp = join(self.output_path, f"replicate_sheet_{count}.csv")
                with open(fp, 'w') as f:
                    replicate.write(f)
                results.append(fp)

        return results

    def _get_prep_file_paths(self, stdout):
        # Strip UserWarnings and empty lines that appear on stdout.
        tmp = [x for x in stdout.split('\n') if x != ''
               and 'UserWarning' not in x]

        results = defaultdict(list)

        for line in tmp:
            qiita_id, prep_file_fp = line.strip().split('\t')
            results[qiita_id].append(prep_file_fp)

        return results

    def run(self, callback=None):
        results = defaultdict(list)

        for command in self.commands:
            # note that if GenPrepFileJob will be run after QCJob in a
            # Pipeline, and QCJob currently moves its products to the final
            # location. It would be cleaner if it did not do this, but
            # currently that is how it's done. Hence, self.output_directory
            # and the path to run_dir might be different locations than the
            # others.
            res = self._system_call(' '.join(command), callback=callback)

            if res['return_code'] != 0:
                raise PipelineError("Seqpro encountered an error")

            # if successful, store results.
            cmd_results = self._get_prep_file_paths(res['stdout'])

            for qiita_id in cmd_results:
                results[qiita_id] += cmd_results[qiita_id]

        self.prep_file_paths = results

        self.mark_job_completed()
