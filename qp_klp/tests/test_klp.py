# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import main
from os import remove
from shutil import rmtree
from json import dumps
from tempfile import mkdtemp
from os.path import exists, isdir, realpath, dirname

from qiita_client.testing import PluginTestCase
from qiita_client import ArtifactInfo

from qp_klp import __version__, plugin

from qp_klp.klp import sequence_processing_pipeline


class KLPTests(PluginTestCase):
    def setUp(self):
        # this will allow us to see the full errors
        self.maxDiff = None

        plugin("https://localhost:8383", 'register', 'ignored')
        self._clean_up_files = []

        self.basedir = dirname(realpath(__file__))

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def test_sequence_processing_pipeline(self):
        # not a valid run_identifier folder and sample_sheet
        params = {
            'run_identifier': '/this/path/doesnt/exist',
            'sample_sheet': 'NA'}

        data = {
            'user': 'demo@microbio.me',
            'command': dumps(['qp-klp', __version__,
                              'Sequence Processing Pipeline']),
            'status': 'running',
            'parameters': dumps(params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        success, ainfo, msg = sequence_processing_pipeline(
            self.qclient, jid, params, out_dir)
        self.assertFalse(success)
        self.assertEqual(msg, "The path doesn't exist or is not a folder")

        # valid run_identifier folder but not sample_sheet
        # NOTE: we are no creating a new job for this test, which is fine
        params = {
            'run_identifier': out_dir,
            'sample_sheet': 'NA'}

        success, ainfo, msg = sequence_processing_pipeline(
            self.qclient, jid, params, out_dir)
        self.assertFalse(success)
        self.assertEqual(
            msg, "Doesn't look like a valid uploaded file; please review.")

        # test success
        # both valid run_identifier and sample_sheet
        # NOTE: we are no creating a new job for this test, which is fine
        params = {
            'run_identifier': out_dir,
            'sample_sheet': {
                'body': 'sample_name\trun_prefix\n1.SKB1.640202\tSKB1.640202',
                'content_type': 'text/plain',
                'filename': 'prep_16S.txt'}}

        success, ainfo, msg = sequence_processing_pipeline(
            self.qclient, jid, params, out_dir)
        self.assertTrue(success)
        exp = [
            ArtifactInfo(
                'output', 'job-output-folder', [(f'{out_dir}/', 'directory')])
        ]

        self.assertEqual(ainfo, exp)


if __name__ == '__main__':
    main()
