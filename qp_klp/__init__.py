# -----------------------------------------------------------------------------
# Copyright (c) 2020--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaPlugin, QiitaCommand

from .klp import list_folder


__version__ = '08.2021'


# Initialize the plugin
class QiitaPluginAdmin(QiitaPlugin):
    _plugin_type = "private"


plugin = QiitaPluginAdmin('qp-klp', __version__, 'Knight Lab Processing')

req_params = {'input_folder': ('string', [''])}
opt_params = dict()
outputs = None
# {'output': 'raw_job_folder'}
dflt_param_set = dict()

list_folder_cmd = QiitaCommand(
    'List Folders', 'List the files (no subfolders) available', list_folder,
    req_params, opt_params, outputs, dflt_param_set)

plugin.register_command(list_folder_cmd)
