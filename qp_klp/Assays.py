from os import listdir, makedirs
from os.path import isfile
from shutil import copyfile
from .FailedSamplesRecord import FailedSamplesRecord
from sequence_processing_pipeline import NuQCJob, FastQCJob, GenPrepFileJob
from os.path import join
import pandas as pd
from json import dumps
from collections import defaultdict
from os.path import basename, dirname


ASSAY_NAME_NONE = "Assay"
ASSAY_NAME_AMPLICON = "Amplicon"
ASSAY_NAME_METAGENOMIC = "Metagenomic"
ASSAY_NAME_METATRANSCRIPTOMIC = "Metatranscriptomic"
METAOMIC_ASSAY_NAMES = [ASSAY_NAME_METAGENOMIC, ASSAY_NAME_METATRANSCRIPTOMIC]


class Assay():
    """
        Assay encapsulate Job()s and other functionality that varies on the
         assay-type of the run. All Assays are mixins for Workflow() classes
         and shouldn't define their own initialization.

        Functionality specific to one assay-type is assigned to a specific
         type. Functionality used by more than one Assay is defined in the base
         Assay() class with a unique name and wrapped by the child classes that
         use it. Functionality used by all children is defined in Assay() class
         w/the name meant to be shared by all and used by the user. Helper
         methods used by other functions in Assay() or its children begin
         w/'_'.

        Assume the following members are present in Workflow()
        self._copy_files()
        self._generate_artifact_name()
        self._get_postqc_fastq_files()
        self._load_prep_into_qiita()
        self.has_replicates
        self.job_pool_size
        self.master_qiita_job_id
        self.pipeline
        self.prep_file_paths
        self.qclient
        self.run_prefixes
        self.special_map
        self.touched_studies_prep_info
        self.update_callback (can be None)
    """
    assay_type = ASSAY_NAME_NONE

    @classmethod
    def _replace_tube_ids_w_sample_names(cls, prep_file_path, tube_id_map):
        """
        Helper method for overwrite_prep_files().
        :param prep_file_path: The path to a generated prep-info file.
        :param tube_id_map: A mapping (dict) of sample-names to tube-ids .
        :return: None
        """
        # reversed_map maps tube-ids to sample-names
        reversed_map = {tube_id_map[k]: k for k in tube_id_map}

        # passing tube_id_map as a parameter allows for easier testing.
        df = pd.read_csv(prep_file_path, sep='\t', dtype=str, index_col=False)
        # save copy of sample_name column as 'old_sample_name'
        df['old_sample_name'] = df['sample_name']
        for i in df.index:
            sample_name = df.at[i, "sample_name"]
            # blanks do not get their names swapped.
            if sample_name.startswith('BLANK'):
                continue

            # remove leading zeroes if they exist to match Qiita results.
            sample_name = sample_name.lstrip('0')

            if sample_name in reversed_map:
                df.at[i, "sample_name"] = reversed_map[sample_name]

        df.to_csv(prep_file_path, index=False, sep="\t")

    def overwrite_prep_files(self, prep_file_paths):
        """
        Replace tube-ids in prep-info files w/sample-names.
        :param prep_file_paths: A list of generated prep-info files.
        :return: None
        """
        # replace tube-ids in prep-info files w/sample-names.
        if self.tube_id_map is None:
            raise ValueError("get_tube_ids_from_qiita() was not called")

        projects = self.pipeline.get_project_info(short_names=True)

        for project in projects:
            project_name = project['project_name']
            qiita_id = str(project['qiita_id'])

            if qiita_id not in self.tube_id_map:
                continue

            # prep files are named in the following form:
            # 20220423_FS10001773_12_BRB11603-0615.Matrix_Tube_LBM_14332.1.tsv
            fqp_name = "%s_%s" % (project_name, qiita_id)
            matching_files = [prep_file for prep_file in prep_file_paths if
                              fqp_name in prep_file]

            if len(matching_files) == 0:
                continue

            for matching_file in matching_files:
                Assay._replace_tube_ids_w_sample_names(matching_file,
                                                       self.tube_id_map[
                                                           qiita_id])

    @classmethod
    def _parse_prep_file(cls, prep_file_path, convert_to_dict=True):
        """
        Helper method for update_prep_templates().
        :param prep_file_path: The path to a generated prep-info file.
        :param convert_to_dict: If True, a dict() is returned.
        :return: A DataFrame() is returned, unless convert_to_dict is True.
        """
        metadata = pd.read_csv(prep_file_path,
                               dtype=str,
                               delimiter='\t',
                               # forces Pandas to not make the first column
                               # the index even when the values appear numeric.
                               index_col=False)

        if metadata is None:
            raise ValueError(f"{prep_file_path} does not exist.")

        metadata.set_index('sample_name', inplace=True)

        if convert_to_dict:
            return metadata.to_dict('index')
        else:
            return metadata

    def _generate_artifact_name(self, prep_file_path):
        """
        Helper method for update_prep_templates().
        :param prep_file_path: The path to a generated prep-info file.
        :return: If prep is a replicate, returns artifact-name, repl-number,
         and True. Otherwise, returns artifact-name and False.
        """
        a_name = f'{self.pipeline.run_id}_{self.lane_number}'
        repl_num = basename(dirname(prep_file_path))

        if self.has_replicates is True:
            # this is a replicate sheet file.
            # append a replication number to each name to
            # make it unique from other replicates.
            # return ('%s_r%s' % (a_name, result[1]), True)
            return ('%s_r%s' % (a_name, repl_num), True)
        else:
            # this is a normal pre-prep or sample-sheet.
            return (a_name, False)


class Amplicon(Assay):
    AMPLICON_TYPE = 'Amplicon'
    AMPLICON_SUB_TYPES = {'16S', '18S', 'ITS'}
    assay_type = ASSAY_NAME_AMPLICON

    def quality_control(self):
        """
        TODO
        :return:
        """
        # Since demuxing and thus quality control occurs downstream of SPP
        # for amplicon runs, there is no QC to be performed at this time.
        # However, we do want to fake 'QCJob' output so that downstream Jobs()
        # can act on the results of ConvertJob() w/out modification.

        # Simulate NuQCJob's output directory for use as input into FastQCJob.
        projects = self.pipeline.get_project_info()
        projects = [x['project_name'] for x in projects]

        for project_name in projects:
            # copy the files from ConvertJob output to faked NuQCJob output
            # folder: $WKDIR/$RUN_ID/NuQCJob/$PROJ_NAME/amplicon
            output_folder = join(self.pipeline.output_path,
                                 'NuQCJob',
                                 project_name,
                                 # for legacy purposes, output folders are
                                 # either 'trimmed_sequences', 'amplicon', or
                                 # 'filtered_sequences'. Hence, this folder
                                 # is not defined using AMPLICON_TYPE as that
                                 # value may or may not equal the needed value.
                                 'amplicon')

            makedirs(output_folder)

            raw_fastq_files_path = join(self.pipeline.output_path,
                                        'ConvertJob')

            # get list of all raw output files to be copied.
            job_output = [join(raw_fastq_files_path, x) for x in
                          listdir(raw_fastq_files_path)]
            job_output = [x for x in job_output if isfile(x)]
            job_output = [x for x in job_output if x.endswith('fastq.gz')]
            # Undetermined files are very small and should be filtered from
            # results.
            job_output = [x for x in job_output if
                          not basename(x).startswith('Undetermined')]

            # copy the file
            for fastq_file in job_output:
                new_path = join(output_folder, basename(fastq_file))
                copyfile(fastq_file, new_path)

            # FastQC expects the ConvertJob output to also be organized by
            # project. Since this would entail running the same ConvertJob
            # multiple times on the same input with just a name-change in
            # the dummy sample-sheet, we instead perform ConvertJob once
            # and copy the output from ConvertJob's output folder into
            # ConvertJob's output folder/project1, project2 ... projectN.
            output_folder = join(raw_fastq_files_path, project_name)
            makedirs(output_folder)

            job_output = [join(raw_fastq_files_path, x) for x in
                          listdir(raw_fastq_files_path)]
            job_output = [x for x in job_output if isfile(x) and x.endswith(
                'fastq.gz') and not basename(x).startswith('Undetermined')]

            for raw_fastq_file in job_output:
                new_path = join(output_folder, basename(raw_fastq_file))
                copyfile(raw_fastq_file, new_path)

    def generate_reports(self):
        """
        TODO
        :return:
        """
        config = self.pipeline.get_software_configuration('fastqc')
        job = FastQCJob(self.pipeline.run_dir,
                        self.pipeline.output_path,
                        join(self.pipeline.output_path, 'ConvertJob'),
                        join(self.pipeline.output_path, 'NuQCJob'),
                        config['nprocs'],
                        config['nthreads'],
                        config['fastqc_executable_path'],
                        config['modules_to_load'],
                        self.master_qiita_job_id,
                        config['queue'],
                        config['nodes'],
                        config['wallclock_time_in_minutes'],
                        config['job_total_memory_limit'],
                        self.job_pool_size,
                        config['multiqc_config_file_path'],
                        config['job_max_array_length'],
                        True)

        job.run(callback=self.update_callback)

    def generate_prep_file(self):
        """
        TODO
        :return:
        """
        config = self.pipeline.get_software_configuration('seqpro')

        # NB: For amplicon runs, the executable used is currently a variant
        # of seqpro called 'seqpro_mf'. It is stored in the same location as
        # 'seqpro'.
        seqpro_path = config['seqpro_path'].replace('seqpro', 'seqpro_mf')

        job = GenPrepFileJob(self.pipeline.run_dir,
                             join(self.pipeline.output_path, 'ConvertJob'),
                             join(self.pipeline.output_path, 'NuQCJob'),
                             self.pipeline.output_path,
                             self.pipeline.input_file_path,
                             seqpro_path,
                             config['modules_to_load'],
                             self.master_qiita_job_id,
                             is_amplicon=True)

        job.run(callback=self.update_callback)

        self.prep_file_paths = job.prep_file_paths
        self.has_replicates = job.has_replicates

    def update_prep_templates(self):
        """
        Update prep-template info in Qiita. Get dict of prep-ids by study-id.
        :return: A dict of lists of prep-ids, keyed by study-id.
        """
        results = defaultdict(list)

        for study_id in self.prep_file_paths:
            for prep_fp in self.prep_file_paths[study_id]:
                metadata = Assay._parse_prep_file(prep_fp)
                afact_name, is_repl = self._generate_artifact_name(prep_fp)
                data = {'prep_info': dumps(metadata),
                        'study': study_id,
                        'data_type': None,
                        'job-id': self.master_qiita_job_id,
                        'name': afact_name}

                if 'target_gene' in metadata[list(metadata.keys())[0]]:
                    tg = metadata[list(metadata.keys())[0]]['target_gene']
                    for key in Amplicon.AMPLICON_SUB_TYPES:
                        if key in tg:
                            data['data_type'] = key

                    if data['data_type'] is None:
                        raise ValueError("data_type could not be "
                                         "determined from target_gene "
                                         "column")
                else:
                    raise ValueError("target_gene must be specified for "
                                     "amplicon type")

                reply = self.qclient.post('/qiita_db/prep_template/',
                                          data=data)
                prep_id = reply['prep']
                results[study_id].append((prep_id, afact_name, is_repl))
                self.run_prefixes[prep_id] = [metadata[sample]['run_prefix']
                                              for sample in metadata]

        self.touched_studies_prep_info = results
        return results

    def load_preps_into_qiita(self):
        """
        TODO
        :return:
        """
        atype = 'FASTQ'

        data = []
        for project, _, qiita_id in self.special_map:
            fastq_files = self._get_postqc_fastq_files(
                self.pipeline.output_path, project)

            for vals in self.touched_studies_prep_info[qiita_id]:
                prep_id, artifact_name, is_repl = vals
                if is_repl:
                    # for Amplicon runs, each prep needs a copy of the
                    # entire set of fastq files, because demuxing samples
                    # happens downstream. If we don't make copies of the
                    # files, Qiita will move the files when loading the
                    # first prep and they won't be available for the
                    # second prep and after.
                    # Note that this will leave the original files present
                    # in the working directory after processing instead of
                    # being moved.
                    working_set = self._copy_files(fastq_files)
                else:
                    working_set = fastq_files

                data.append(self._load_prep_into_qiita(
                    self.qclient, prep_id, artifact_name, qiita_id, project,
                    working_set, atype))

        df = pd.DataFrame(data)
        opath = join(self.pipeline.output_path, 'touched_studies.html')
        with open(opath, 'w') as f:
            f.write(df.to_html(border=2, index=False, justify="left",
                               render_links=True, escape=False))

        return df


class MetaOmic(Assay, FailedSamplesRecord):
    """
    MetaOmic() is a base class for Metagenomic() and Metatranscriptomic(),
    which are currently identical in functionality.
    """
    assay_type = ASSAY_NAME_NONE

    def quality_control(self):
        # because this is a mixin, assume containing object will contain
        # a pipeline object.
        config = self.pipeline.get_software_configuration('nu-qc')

        # base quality control used by multiple Assay types.
        job = NuQCJob(join(self.pipeline.output_path, 'ConvertJob'),
                      self.pipeline.output_path,
                      self.pipeline.sample_sheet.path,
                      config['minimap2_databases'],
                      config['queue'],
                      config['nodes'],
                      config['wallclock_time_in_minutes'],
                      config['job_total_memory_limit'],
                      config['fastp_executable_path'],
                      config['minimap2_executable_path'],
                      config['samtools_executable_path'],
                      config['modules_to_load'],
                      self.master_qiita_job_id,
                      config['job_max_array_length'],
                      config['known_adapters_path'],
                      bucket_size=config['bucket_size'],
                      length_limit=config['length_limit'],
                      cores_per_task=config['cores_per_task'],
                      movi_path=config['movi_executable_path'],
                      gres_value=config['gres_value'],
                      pmls_path=config['pmls_path'])

        job.run(callback=self.update_callback)

        self.fsr_write(job.audit(self.pipeline.get_sample_ids()), 'NuQCJob')

    def generate_reports(self):
        config = self.pipeline.get_software_configuration('fastqc')
        job = FastQCJob(self.pipeline.run_dir,
                        self.pipeline.output_path,
                        join(self.pipeline.output_path, 'ConvertJob'),
                        join(self.pipeline.output_path, 'NuQCJob'),
                        config['nprocs'],
                        config['nthreads'],
                        config['fastqc_executable_path'],
                        config['modules_to_load'],
                        self.master_qiita_job_id,
                        config['queue'],
                        config['nodes'],
                        config['wallclock_time_in_minutes'],
                        config['job_total_memory_limit'],
                        self.job_pool_size,
                        config['multiqc_config_file_path'],
                        config['job_max_array_length'],
                        False)

        job.run(callback=self.update_callback)

        self.fsr_write(job.audit(self.pipeline.get_sample_ids()), 'FastQCJob')

        # as generate_reports is the final step that updates self.fsr,
        # we can generate the final human-readable report following this step.
        self.fsr_generate_report()

    def generate_prep_file(self):
        config = self.pipeline.get_software_configuration('seqpro')

        job = GenPrepFileJob(self.pipeline.run_dir,
                             join(self.pipeline.output_path, 'ConvertJob'),
                             join(self.pipeline.output_path, 'NuQCJob'),
                             self.pipeline.output_path,
                             self.pipeline.input_file_path,
                             config['seqpro_path'],
                             config['modules_to_load'],
                             self.master_qiita_job_id,
                             is_amplicon=False)

        job.run(callback=self.update_callback)

        self.prep_file_paths = job.prep_file_paths
        self.has_replicates = job.has_replicates

    def update_prep_templates(self):
        """
        Update prep-template info in Qiita. Get dict of prep-ids by study-id.
        :return: A dict of lists of prep-ids, keyed by study-id.
        """
        results = defaultdict(list)

        for study_id in self.prep_file_paths:
            for prep_fp in self.prep_file_paths[study_id]:
                metadata = Assay._parse_prep_file(prep_fp)
                afact_name, is_repl = self._generate_artifact_name(prep_fp)
                data = {'prep_info': dumps(metadata),
                        'study': study_id,
                        'job-id': self.master_qiita_job_id,
                        'name': afact_name,
                        'data_type': self.pipeline.pipeline_type}

                # since all Assays are mixins for Workflows, assume
                # self.qclient exists and available.
                reply = self.qclient.post('/qiita_db/prep_template/',
                                          data=data)
                prep_id = reply['prep']
                results[study_id].append((prep_id, afact_name, is_repl))
                self.run_prefixes[prep_id] = [metadata[sample]['run_prefix']
                                              for sample in metadata]

        self.touched_studies_prep_info = results
        return results

    def load_preps_into_qiita(self):
        atype = 'per_sample_FASTQ'

        data = []
        for project, _, qiita_id in self.special_map:
            fastq_files = self._get_postqc_fastq_files(
                self.pipeline.output_path, project)

            for vals in self.touched_studies_prep_info[qiita_id]:
                prep_id, artifact_name, is_repl = vals
                # for meta*omics, generate the subset of files used by
                # this prep only.
                working_set = {}
                for key in fastq_files:
                    working_set[key] = []
                    for run_prefix in self.run_prefixes[prep_id]:
                        working_set[key] += [fastq for fastq in
                                             fastq_files[key] if
                                             run_prefix in fastq]

                if is_repl:
                    working_set = self._copy_files(working_set)

                data.append(self._load_prep_into_qiita(
                    self.qclient, prep_id, artifact_name, qiita_id, project,
                    working_set, atype))

        df = pd.DataFrame(data)
        opath = join(self.pipeline.output_path, 'touched_studies.html')
        with open(opath, 'w') as f:
            f.write(df.to_html(border=2, index=False, justify="left",
                               render_links=True, escape=False))

        return df


class Metagenomic(MetaOmic):
    METAGENOMIC_TYPE = 'Metagenomic'
    assay_type = ASSAY_NAME_METAGENOMIC


class Metatranscriptomic(MetaOmic):
    METATRANSCRIPTOMIC_TYPE = 'Metatranscriptomic'

    assay_type = ASSAY_NAME_METATRANSCRIPTOMIC
