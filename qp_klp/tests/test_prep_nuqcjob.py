# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from unittest import main
from os.path import join, exists, isdir
from os import remove
from shutil import rmtree, copyfile
from tempfile import mkdtemp
from qp_klp import plugin
from qiita_client.testing import PluginTestCase
from qp_klp.klp import prep_NuQCJob
from qp_klp import __version__ as version_qp_klp
from json import dumps


class Test_prep_NuQCJob(PluginTestCase):
    def setUp(self):
        plugin("https://localhost:21174", 'register', 'ignored')

        out_dir = mkdtemp()
        self.maxDiff = None
        self.out_dir = out_dir
        self._clean_up_files = []
        self._clean_up_files.append(out_dir)

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def _setup_test(self):
        prep_info_dict = {
            'SKB8.640193': {'run_prefix': 'S22205_S104', 'index': 'GTCACTGT',
                            'index2': 'GCAAGATC'},
            'SKD8.640184': {'run_prefix': 'S22282_S102', 'index': 'GCCATAAC',
                            'index2': 'AGTCTCAC'}}
        data = {'prep_info': dumps(prep_info_dict),
                # magic #1 = testing study
                'study': 1,
                'data_type': 'Metagenomic'}
        pid = self.qclient.post('/apitest/prep_template/', data=data)['prep']

        # inserting artifacts
        in_dir = mkdtemp()
        self._clean_up_files.append(in_dir)

        source_gz = 'qp_klp/tests/data/dummy.fastq.gz'
        fps = [
            join(in_dir, 'S22205_S104_L001_R1_001.fastq.gz'),
            join(in_dir, 'S22205_S104_L001_R2_001.fastq.gz'),
            join(in_dir, 'S22282_S102_L001_R1_001.fastq.gz'),
            join(in_dir, 'S22282_S102_L001_R2_001.fastq.gz')
        ]
        for fp in fps:
            copyfile(source_gz, fp)
        fp_summary = join(in_dir, 'summary.html')
        copyfile('qp_klp/tests/data/summary.html', fp_summary)

        data = {
            'filepaths': dumps([
                (fps[0], 'raw_forward_seqs'),
                (fps[1], 'raw_reverse_seqs'),
                (fps[2], 'raw_forward_seqs'),
                (fps[3], 'raw_reverse_seqs'),
                (fp_summary, 'html_summary')]),
            'type': "per_sample_FASTQ",
            'name': "Test Woltka artifact",
            'prep': pid}
        self.qclient.post('/apitest/artifact/', data=data)['artifact']

        data = {'user': 'demo@microbio.me',
                'command': dumps(['qp-klp', version_qp_klp, 'prep_NuQCJob']),
                'status': 'running',
                'parameters': dumps({'prep_id': pid})}
        job_id = self.qclient.post(
            '/apitest/processing_job/', data=data)['job']

        return pid, job_id

    def test_prep_nuqcjob(self):
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # test failure
        jid = 'bcc7ebcd-39c1-43e4-af2d-822e3589f14d'
        success, ainfo, msg = prep_NuQCJob(
            self.qclient, jid, {'prep_id': '1'}, out_dir)
        self.assertEqual(msg, "Prep 1 has a not valid data type: 18S")
        self.assertFalse(success)

        # test success; note that success at this point is that
        # the job fails when submitting via sbatch
        pid, job_id = self._setup_test()
        success, ainfo, msg = prep_NuQCJob(
            self.qclient, job_id, {'prep_id': pid}, out_dir)
        self.assertRegex(msg, r'(?s)sbatch:(?s)found(?s)')
        self.assertFalse(success)


if __name__ == '__main__':
    main()
