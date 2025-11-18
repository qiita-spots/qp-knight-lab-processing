# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from qiita_client import QiitaCommand, QiitaPlugin

from .klp import PrepNuQCJob, sequence_processing_pipeline


class QiitaPluginAdmin(QiitaPlugin):
    _plugin_type = "private"


__version__ = "2024.11"

# Adding SPP job

plugin = QiitaPluginAdmin("qp-klp", __version__, "Knight Lab Processing")

req_params = {
    "run_identifier": ("string", [""]),
    "sample_sheet": ("prep_template", [""]),
    "lane_number": ("integer", [None]),
}
opt_params = dict()
outputs = {"output": "job-output-folder"}
dflt_param_set = dict()

sequence_processing_pipeline_cmd = QiitaCommand(
    "Sequence Processing Pipeline",
    "Associate a sample sheet with a "
    "sequencing run and generate the corresponding sequence files",
    sequence_processing_pipeline,
    req_params,
    opt_params,
    outputs,
    dflt_param_set,
)

plugin.register_command(sequence_processing_pipeline_cmd)

# adding PrepNuQC job

req_params = {"prep_id": ("integer", [None])}
outputs = {"output": "job-output-folder"}

prepNuQC_cmd = QiitaCommand(
    "Human Filter & QC existing Prep",
    "Apply NuQCJob to an existing Qiita prep raw data",
    PrepNuQCJob,
    req_params,
    dict(),
    outputs,
    dict(),
)

plugin.register_command(prepNuQC_cmd)
