#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
import click
import pandas as pd
from glob import glob
from os import makedirs


@click.command()
@click.argument('sample_list', required=True)
@click.argument('run_folder', required=True)
@click.argument('outdir', required=True)
@click.argument('threads', required=True, default=1)
def generate_bam2fastq_commands(sample_list, run_folder, outdir, threads):
    """Generates the bam2fastq commands"""
    df = pd.read_csv(sample_list, sep='\t')
    files = {f.split('.')[-2]: f
             for f in glob(f'{run_folder}/*/hifi_reads/*.bam')}

    makedirs(outdir, exist_ok=True)

    commands, missing_files = [], []
    for _, row in df.iterrows():
        bc = row['barcode']
        sn = row['sample_name']
        pn = row['project_name']
        lane = row['lane']
        if bc not in files:
            missing_files.append(bc)
            continue
        od = f'{outdir}/{pn}'
        makedirs(od, exist_ok=True)
        fn = f'{od}/{sn}_L00{lane}_R1_001'
        cmd = (f'bam2fastq -j {threads} -o {fn} -c 9 '
               f'{files[bc]}; '
               f'fqtools count {fn}.fastq.gz > '
               f'{fn}.counts.txt')
        commands.append(cmd)

    if missing_files:
        raise ValueError(
            f'{run_folder} is missing barcodes: {missing_files}')

    for cmd in commands:
        print(cmd)
