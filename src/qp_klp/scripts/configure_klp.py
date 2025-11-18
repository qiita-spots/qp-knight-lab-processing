#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

import click

from qp_klp import plugin


@click.command()
@click.option(
    "--env-script",
    prompt="Environment script",
    required=True,
    default="source activate qp-klp",
)
@click.option(
    "--ca-cert",
    prompt="Server certificate",
    required=False,
    default="None",
    show_default=True,
)
def config(env_script, ca_cert):
    """Generates the Qiita configuration files"""
    if ca_cert == "None":
        ca_cert = None
    plugin.generate_config(env_script, "start_klp", ca_cert)


if __name__ == "__main__":
    config()
