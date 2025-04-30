# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from qp_klp.WorkflowFactory import WorkflowFactory
from unittest import main
from os import makedirs
from os.path import dirname, abspath, join
from shutil import copyfile
from pathlib import Path
import pandas as pd
from qp_klp.Protocol import (PROTOCOL_NAME_ILLUMINA,
                             PROTOCOL_NAME_TELLSEQ)
from qp_klp.Assays import (ASSAY_NAME_METAGENOMIC,
                           ASSAY_NAME_METATRANSCRIPTOMIC,
                           ASSAY_NAME_AMPLICON)
from shutil import rmtree
from qiita_client.testing import PluginTestCase
from qiita_client.exceptions import NotFoundError
from sequence_processing_pipeline.PipelineError import PipelineError


class WorkflowFactoryTests(PluginTestCase):
    def setUp(self):
        self.remove_these = []
        self.base_dir = dirname(abspath(__file__))
        self.output_dir = join(self.base_dir, 'tests', 'test_output')

    def tearDown(self):
        for fp in self.remove_these:
            rmtree(fp)

    def _create_directory(self, fp):
        makedirs(fp, exist_ok=True)
        self.remove_these.append(fp)

    def test_no_parameters(self):
        # confirm factory doesn't create workflows when given no parameters.
        msg = "kwargs must not be None and must define 'uif_path'"
        with self.assertRaisesRegex(ValueError, msg):
            WorkflowFactory.generate_workflow()

        with self.assertRaisesRegex(ValueError, msg):
            WorkflowFactory.generate_workflow(**{})

        with self.assertRaisesRegex(ValueError, msg):
            WorkflowFactory.generate_workflow(**{"foo": "bar"})

    def test_workflow_performs_parameter_checking(self):
        # confirm that once factory is given enough information to select a
        # Workflow() class, the class itself confirms that it has all the
        # parameters it needs.
        kwargs = {"uif_path": "tests/data/sample-sheets/metagenomic/"
                  "illumina/good_sheet1.csv"}

        msg = ("The following values must also be defined in kwargs for "
               "StandardMetagenomicWorkflow workflows: qclient, lane_number,"
               " config_fp, run_identifier, output_dir, job_id, lane_number,"
               " is_restart")

        with self.assertRaisesRegex(ValueError, msg):
            WorkflowFactory.generate_workflow(**kwargs)

        kwargs = {"uif_path": "tests/data/pre-preps/good_pre_prep1.txt"}

        msg = ("The following values must also be defined in kwargs for "
               "StandardAmpliconWorkflow workflows: qclient, config_fp, "
               "run_identifier, output_dir, job_id, is_restart")

        with self.assertRaisesRegex(ValueError, msg):
            WorkflowFactory.generate_workflow(**kwargs)

    def test_invalid_sample_sheets(self):
        # confirm an otherwise good-sample-sheet w/a bad SheetType is going
        # fail because it doesn't pass validation. SheetType directly
        # determines the Instrument Mixin to be used.
        kwargs = {"uif_path": "tests/data/sample-sheets/metagenomic/"
                  "illumina/bad_sheet1.csv"}

        msg = ("'tests/data/sample-sheets/metagenomic/illumina/"
               "bad_sheet1.csv' does not appear to be a valid sample-sheet.")

        with self.assertRaisesRegex(ValueError, msg):
            WorkflowFactory.generate_workflow(**kwargs)

        # confirm an otherwise good-sample-sheet w/a bad Assay value is going
        # to fail because it doesn't pass validation.
        kwargs = {"uif_path": "tests/data/sample-sheets/metagenomic/"
                  "illumina/bad_sheet2.csv"}

        msg = ("'tests/data/sample-sheets/metagenomic/illumina/"
               "bad_sheet2.csv' does not appear to be a valid sample-sheet.")

        with self.assertRaisesRegex(ValueError, msg):
            WorkflowFactory.generate_workflow(**kwargs)

        # confirm a file that is obviously not a mapping-file or sample-sheet
        # file will not follow through.
        kwargs = {"uif_path": "tests/data/Demultiplex_Stats.csv"}

        msg = ("Your uploaded file doesn't appear to be a sample-sheet or a "
               "mapping-file.")
        with self.assertRaisesRegex(ValueError, msg):
            WorkflowFactory.generate_workflow(**kwargs)

    def _inject_data(self, wf):
        '''This is a helper method for testing that all the steps are run
           when testing wf.execute_pipeline()
        '''
        # inject Convert/NuQC/FastQC/MultiQC/GenPrepFileJob files so we can
        # move down the pipeline; first let's create the base folders
        gz_source = f'{self.base_dir}/data/dummy.fastq.gz'
        convert_dir = f'{self.output_dir}/ConvertJob'
        reports_dir = f'{self.output_dir}/ConvertJob/Reports'
        nuqc_dir = f'{self.output_dir}/NuQCJob'
        fastqc_dir = f'{self.output_dir}/FastQCJob/logs/'
        multiqc_dir = f'{self.output_dir}/MultiQCJob/logs/'
        genprep_dir = (f'{self.output_dir}/GenPrepFileJob/'
                       '211021_A00000_0000_SAMPLE/')
        makedirs(convert_dir, exist_ok=True)
        Path(f'{convert_dir}/job_completed').touch()
        makedirs(reports_dir, exist_ok=True)
        makedirs(nuqc_dir, exist_ok=True)
        makedirs(fastqc_dir, exist_ok=True)
        makedirs(multiqc_dir, exist_ok=True)
        makedirs(genprep_dir, exist_ok=True)
        # now let's create the required project folders
        for project in wf.pipeline.get_project_info():
            sp = project['project_name']
            makedirs(f'{convert_dir}/{sp}', exist_ok=True)
            makedirs(f'{nuqc_dir}/filtered_sequences/{sp}', exist_ok=True)
            makedirs(f'{genprep_dir}/{sp}/filtered_sequences/', exist_ok=True)

        # then loop over samples and stage all fastq.gz files
        demux_stats = []
        samples = wf.pipeline.sample_sheet.samples
        for i, sample in enumerate(samples):
            rp = sample["Sample_ID"]
            sp = sample["Sample_Project"]

            # ConvertJob
            demux_stats.append(
                {'Lane': sample['Lane'], 'SampleID': rp,
                    'Sample_Project': sp, 'Index': sample['index'],
                 '# Reads': 2})
            dname = f'{convert_dir}/{sp}'
            copyfile(gz_source, f'{dname}/{rp}_L001_R1_001.fastq.gz')
            copyfile(gz_source, f'{dname}/{rp}_L001_R2_001.fastq.gz')

            # NuQCJob
            dname = f'{nuqc_dir}/filtered_sequences/{sp}'
            copyfile(gz_source, f'{dname}/{rp}_L001_R1_001.fastq.gz')
            copyfile(gz_source, f'{dname}/{rp}_L001_R2_001.fastq.gz')

            # GenPrepFileJob
            gprep_base = f'{genprep_dir}/{sp}/filtered_sequences/{rp}'
            copyfile(gz_source, f'{gprep_base}_L001_R1_001.fastq.gz')
            copyfile(gz_source, f'{gprep_base}_L001_R2_001.fastq.gz')

        # this is required by the Convert step
        pd.DataFrame(demux_stats).set_index('SampleID').to_csv(
            f'{reports_dir}/Demultiplex_Stats.csv')

        # generating the "*.completed" files
        for i in range(len(samples)*3):
            Path(f'{fastqc_dir}/FastQCJob_{i}.completed').touch()
            Path(f'{multiqc_dir}/MultiQCJob_{i}.completed').touch()


    def test_metagenomic_workflow_creation(self):
        kwargs = {"uif_path": "tests/data/sample-sheets/metagenomic/"
                  "illumina/good_sheet1.csv",
                  "qclient": self.qclient,
                  "lane_number": "1",
                  "config_fp": "tests/configuration.json",
                  "run_identifier": "211021_A00000_0000_SAMPLE",
                  "output_dir": self.output_dir,
                  "job_id": "78901",
                  "is_restart": False
                  }

        self._create_directory(kwargs['output_dir'])

        wf = WorkflowFactory.generate_workflow(**kwargs)

        # confirm that the proper type of workflow was generated.
        self.assertEqual(wf.protocol_type, PROTOCOL_NAME_ILLUMINA)
        self.assertEqual(wf.assay_type, ASSAY_NAME_METAGENOMIC)

        with self.assertRaisesRegex(
                PipelineError, 'There are no fastq files for FastQCJob'):
            wf.execute_pipeline()

        self._inject_data(wf)
        wf.execute_pipeline()

    def test_metatranscriptomic_workflow_creation(self):
        kwargs = {"uif_path": "tests/data/sample-sheets/"
                  "metatranscriptomic/illumina/good_sheet1.csv",
                  "qclient": self.qclient,
                  "lane_number": "1",
                  "config_fp": "tests/configuration.json",
                  "run_identifier": "211021_A00000_0000_SAMPLE",
                  "output_dir": self.output_dir,
                  "job_id": "78901",
                  "is_restart": False
                  }

        self._create_directory(kwargs['output_dir'])

        wf = WorkflowFactory.generate_workflow(**kwargs)

        # confirm that the proper type of workflow was generated.
        self.assertEqual(wf.protocol_type, PROTOCOL_NAME_ILLUMINA)
        self.assertEqual(wf.assay_type, ASSAY_NAME_METATRANSCRIPTOMIC)

        with self.assertRaisesRegex(
                PipelineError, 'There are no fastq files for FastQCJob'):
            wf.execute_pipeline()

        self._inject_data(wf)
        wf.execute_pipeline()

    def test_amplicon_workflow_creation(self):
        kwargs = {"uif_path": "tests/data/pre-preps/good_pre_prep1.txt",
                  "qclient": self.qclient,
                  "config_fp": "tests/configuration.json",
                  "run_identifier": "211021_A00000_0000_SAMPLE",
                  "output_dir": self.output_dir,
                  "job_id": "78901",
                  "is_restart": False
                  }

        self._create_directory(kwargs['output_dir'])

        wf = WorkflowFactory.generate_workflow(**kwargs)

        # confirm that the proper type of workflow was generated.
        self.assertEqual(wf.protocol_type, PROTOCOL_NAME_ILLUMINA)
        self.assertEqual(wf.assay_type, ASSAY_NAME_AMPLICON)

        with self.assertRaisesRegex(
                NotFoundError, '{"message": "Study not found"}'):
            wf.execute_pipeline()

        self._inject_data(wf)
        wf.execute_pipeline()

    def test_tellseq_workflow_creation(self):
        kwargs = {"uif_path": "tests/data/sample-sheets/metagenomic/"
                  "tellseq/good_sheet1.csv",
                  "qclient": self.qclient,
                  "config_fp": "tests/configuration.json",
                  "run_identifier": "211021_A00000_0000_SAMPLE",
                  "output_dir": self.output_dir,
                  "job_id": "78901",
                  "lane_number": "1",
                  "is_restart": False
                  }

        self._create_directory(kwargs['output_dir'])

        wf = WorkflowFactory.generate_workflow(**kwargs)

        # confirm that the proper type of workflow was generated.
        self.assertEqual(wf.protocol_type, PROTOCOL_NAME_TELLSEQ)
        self.assertEqual(wf.assay_type, ASSAY_NAME_METAGENOMIC)

        with self.assertRaisesRegex(
                FileNotFoundError,
                'TellReadJob/sample_index_list_TellReadJob.txt'):
            wf.execute_pipeline()

        self._inject_data(wf)
        wf.execute_pipeline()



if __name__ == '__main__':
    main()
