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
from os.path import dirname, abspath, join, exists
from shutil import copyfile
from shutil import rmtree
from pathlib import Path
import pandas as pd
from metapool.sample_sheet import (PROTOCOL_NAME_PACBIO_SMRT)
from sequence_processing_pipeline.PipelineError import PipelineError
from qp_klp.Assays import ASSAY_NAME_METAGENOMIC
from qiita_client.testing import PluginTestCase


class WorkflowFactoryTests(PluginTestCase):
    def setUp(self):
        self.base_dir = dirname(abspath(__file__))
        self.output_dir = join(self.base_dir, 'test_output')
        self.remove_these = [self.output_dir]

        self.gz_source = f'{self.base_dir}/data/dummy.fastq.gz'

    def tearDown(self):
        for fp in self.remove_these:
            if exists(fp):
                rmtree(fp)

    def _create_directory(self, fp):
        # just making sure that we always start with a clean folder:
        if exists(fp):
            rmtree(fp)
        makedirs(fp)
        self.remove_these.append(fp)

    def _inject_data(self, wf):
        '''This is a helper method for testing that all the steps are run
           when testing wf.execute_pipeline()
        '''
        samples = wf.pipeline.sample_sheet.samples

        convert_dir = f'{self.output_dir}/ConvertJob'
        reports_dir = f'{self.output_dir}/ConvertJob/Reports'
        makedirs(convert_dir, exist_ok=True)
        makedirs(reports_dir, exist_ok=True)
        Path(f'{convert_dir}/job_completed').touch()
        nuqc_dir = f'{self.output_dir}/NuQCJob'
        fastqc_dir = f'{self.output_dir}/FastQCJob/logs/'
        multiqc_dir = f'{self.output_dir}/MultiQCJob/logs/'
        genprep_dir = (f'{self.output_dir}/GenPrepFileJob/'
                       '211021_A00000_0000_SAMPLE/')
        makedirs(nuqc_dir, exist_ok=True)
        makedirs(fastqc_dir, exist_ok=True)
        makedirs(multiqc_dir, exist_ok=True)
        makedirs(genprep_dir, exist_ok=True)
        # now let's create the required project folders
        for project in wf.pipeline.get_project_info():
            sp = project['project_name']
            makedirs(f'{convert_dir}/{sp}', exist_ok=True)
            makedirs(f'{nuqc_dir}/filtered_sequences/{sp}', exist_ok=True)
            makedirs(f'{genprep_dir}/{sp}/filtered_sequences/',
                     exist_ok=True)

        # # then loop over samples and stage all fastq.gz files
        dstats = []
        for sample in samples:
            sn = sample['Sample_Name']
            sp = sample["Sample_Project"]
            dstats.append({'SampleID': sn, '# Reads': 2})
            dname = f'{convert_dir}/{sp}'
            Path(f'{dname}/{sn}_R1.fastq.gz').touch()
            with open(f'{dname}/{sn}_R1.counts.txt', 'w') as f:
                f.write("2")

            # NuQCJob
            dname = f'{nuqc_dir}/filtered_sequences/{sp}'
            copyfile(self.gz_source, f'{dname}/{sn}_R1.fastq.gz')

            # GenPrepFileJob
            gprep_base = f'{genprep_dir}/{sp}/filtered_sequences/{sn}'
            Path(f'{gprep_base}_R1.fastq.gz').touch()

        pd.DataFrame(dstats).set_index('SampleID').to_csv(
            f'{reports_dir}/Demultiplex_Stats.csv')

        # generating the "*.completed" files
        for i in range(len(samples)*3):
            Path(f'{fastqc_dir}/FastQCJob_{i}.completed').touch()
            Path(f'{multiqc_dir}/MultiQCJob_{i}.completed').touch()

    def test_pacbio_metagenomic_workflow_creation(self):
        kwargs = {"uif_path": "tests/data/sample-sheets/metagenomic/"
                  "pacbio/good_pacbio_metagv10.csv",
                  "qclient": self.qclient,
                  "lane_number": "1",
                  "config_fp": "tests/configuration.json",
                  "run_identifier": "r11111_20250101_111111",
                  "output_dir": self.output_dir,
                  "job_id": "78901",
                  "is_restart": False
                  }

        self._create_directory(kwargs['output_dir'])

        wf = WorkflowFactory.generate_workflow(**kwargs)

        # confirm that the proper type of workflow was generated.
        self.assertEqual(wf.protocol_type, PROTOCOL_NAME_PACBIO_SMRT)
        self.assertEqual(wf.assay_type, ASSAY_NAME_METAGENOMIC)

        self._inject_data(wf)

        # we can add some checks/tests of the initial scripts (mainly Convert)
        # but not doing now as is not required

        # Metagenomic is a valid data type in the default qiita test
        # database but job-id: 78901 doesn't exist; however, if we get
        # to here, it means that all steps have run to completion
        # and the system is trying to create the preps.
        # Note: Qiita job_id's are UUID in the database and this tests
        # uses 78901 as the job_id so the db complains about the format
        # with self.assertRaisesRegex(RuntimeError, 'invalid input '
        #                             'syntax for type uuid: "78901"'):
        with self.assertRaisesRegex(PipelineError, 'There are no fastq '
                                    'files for FastQCJob'):
            wf.execute_pipeline()


if __name__ == '__main__':
    main()
