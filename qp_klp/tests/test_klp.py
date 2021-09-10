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
        params = {'input_folder': '/this/path/doesnt/exist'}

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

        params = {'input_folder': out_dir}
        data['parameters'] = dumps(params)
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']
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
