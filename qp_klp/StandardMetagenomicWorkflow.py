from functools import partial
import sample_sheet
import pandas as pd
from os.path import basename, join
from os import symlink, makedirs
from datetime import datetime
from metapool import MetagenomicSampleSheetv90
from pathlib import Path

from .Protocol import Illumina
from sequence_processing_pipeline.Pipeline import Pipeline
from .Assays import Metagenomic
from .Assays import ASSAY_NAME_METAGENOMIC
from .FailedSamplesRecord import FailedSamplesRecord
from .Workflows import Workflow, WorkflowError


class StandardMetagenomicWorkflow(Workflow, Metagenomic, Illumina):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.mandatory_attributes = ['qclient', 'uif_path',
                                     'lane_number', 'config_fp',
                                     'run_identifier', 'output_dir', 'job_id',
                                     'lane_number', 'is_restart']

        self.confirm_mandatory_attributes()

        # second stage initializer that could conceivably be pushed down into
        # specific children requiring specific parameters.
        self.qclient = self.kwargs['qclient']

        self.pipeline = Pipeline(self.kwargs['config_fp'],
                                 self.kwargs['run_identifier'],
                                 self.kwargs['uif_path'],
                                 self.kwargs['output_dir'],
                                 self.kwargs['job_id'],
                                 ASSAY_NAME_METAGENOMIC,
                                 lane_number=self.kwargs['lane_number'])

        self.fsr = FailedSamplesRecord(self.kwargs['output_dir'],
                                       self.pipeline.sample_sheet.samples)

        self.master_qiita_job_id = self.kwargs['job_id']

        self.lane_number = self.kwargs['lane_number']
        self.is_restart = bool(self.kwargs['is_restart'])

        if self.is_restart is True:
            self.determine_steps_to_skip()

        # this is a convenience member to allow testing w/out updating Qiita.
        self.update = True

        if 'update_qiita' in kwargs:
            if not isinstance(kwargs['update_qiita'], bool):
                raise ValueError("value for 'update_qiita' must be of "
                                 "type bool")

            self.update = kwargs['update_qiita']


class PrepNuQC(StandardMetagenomicWorkflow):
    def __init__(self, **kwargs):
        qclient = kwargs['qclient']
        job_id = kwargs['job_id']
        parameters = kwargs['parameters']
        out_dir = kwargs['out_dir']
        config_fp = kwargs['config_fp']
        status_line = kwargs['status_line']

        out_path = partial(join, out_dir)
        self.final_results_path = out_path('final_results')
        makedirs(self.final_results_path, exist_ok=True)

        pid = parameters.pop('prep_id')

        prep_info = qclient.get(f'/qiita_db/prep_template/{pid}/')
        dt = prep_info['data_type']
        sid = prep_info['study']
        if dt not in {'Metagenomic', 'Metatranscriptomic'}:
            raise WorkflowError(f'Prep {pid} has a not valid data type: {dt}')
        aid = prep_info['artifact']
        if not str(aid).isnumeric():
            raise WorkflowError(f'Prep {pid} has a not valid artifact: {aid}')

        files, pt = qclient.artifact_and_preparation_files(aid)
        html_summary = qclient.get_artifact_html_summary(aid)
        if html_summary is None:
            raise WorkflowError(f'Artifact {aid} doesnot have a summary, '
                                'please generate one.')
        df_summary = pd.read_html(html_summary)[0]
        pt.set_index('sample_name', inplace=True)

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
            data=[[project_name, 'NA', 'NA', 'NA', 'NA',
                   'FALSE', 'TRUE', sid]])

        df_summary = df_summary[df_summary.file_type == 'raw_forward_seqs']
        data = []
        for k, vals in pt.iterrows():
            k = k.split('.', 1)[-1]
            rp = vals['run_prefix']
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
            _d = df_summary[
                df_summary.filename.str.startswith(rp)]
            if _d.shape[0] != 1:
                ValueError(f'The run_prefix {rp} from {k} has {_d.shape[0]} '
                           'matches with files')
            data.append({
                'Lane': '1', 'SampleID': rp, 'Sample_Project': project_name,
                'Index': vals['index'], '# Reads': _d.reads.values[0]})

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
        # creating Demultiplex_Stats.csv
        reports_folder = out_path('ConvertJob', 'Reports')
        makedirs(reports_folder, exist_ok=True)
        pd.DataFrame(data).set_index('SampleID').to_csv(
            f'{reports_folder}/Demultiplex_Stats.csv')

        for fs in files.values():
            for f in fs:
                bn = basename(f['filepath']).replace(
                    '.trimmed.fastq.gz', '.fastq.gz')
                symlink(f['filepath'], f'{project_folder}/{bn}')

        # create job_completed file to skip this step
        Path(f'{convert_path}/job_completed').touch()

        kwargs = {'qclient': qclient,
                  'uif_path': new_sample_sheet,
                  'lane_number': "1",
                  'config_fp': config_fp,
                  'run_identifier': '250225_LH00444_0301_B22N7T2LT4',
                  'output_dir': out_dir,
                  'job_id': job_id,
                  'status_update_callback': status_line.update_job_status,
                  # set 'update_qiita' to False to avoid updating Qiita DB
                  # and copying files into uploads dir. Useful for testing.
                  'update_qiita': True,
                  'is_restart': True}

        super().__init__(**kwargs)
