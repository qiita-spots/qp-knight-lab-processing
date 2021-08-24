# -----------------------------------------------------------------------------
# Copyright (c) 2020--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaPlugin, QiitaCommand

from .klp import list_folder


__version__ = '2021.08'

plugin = QiitaPlugin('qp-klp', __version__, 'Knight Lab Processing')

req_params = {'input_folder': ('string', [''])}
opt_params = dict()
outputs = {'output': 'job-output-folder'}
dflt_param_set = dict()

list_folder_cmd = QiitaCommand(
    'List Folders', 'List the files (no subfolders) available', list_folder,
    req_params, opt_params, outputs, dflt_param_set)

plugin.register_command(list_folder_cmd)
