#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from glob import glob
from os import makedirs

import click
import numpy as np
import pandas as pd


@click.command()
@click.argument("sample_list", required=True)
@click.argument("run_folder", required=True)
@click.argument("outdir", required=True)
@click.argument("threads", required=True, default=1)
def generate_bam2fastq_commands(sample_list, run_folder, outdir, threads):
    """Generates the bam2fastq commands"""
    df = pd.read_csv(sample_list, sep="\t", dtype=str)

    # pacbio raw files are in a hifi_reads folder, wihtin multiple folders, which represent
    # smartcells (1_A01, 2_A02, ect), within the run-id folder; and are named
    # m[run-id]XXX.hifi_reads.[barcode].bam; thus to find the [barcode] we
    # can split on '.' and then the second to last element [-2].
    files = {f.split(".")[-2]: f for f in glob(f"{run_folder}/*/hifi_reads/*.bam")}
    # now, twisted files are within each smartcells, in a new folder - name can change but includes
    # the project_name, then call-lima/execution/[twisted-barcode] and
    # m[run-id]XXX.hifi_reads.[twisted-barcode].bam
    twisted_files = {
        f.split("/")[-2]: f
        for f in glob(f"{run_folder}/*/*/call-lima/execution/*/*.bam")
    }

    makedirs(outdir, exist_ok=True)

    commands, missing_files = [], []
    for _, row in df.iterrows():
        bc = row["barcode"]
        sn = row["sample_name"]
        pn = row["project_name"]
        lane = row["lane"]
        tai = row["twist_adaptor_id"]

        if tai is None or tai is np.nan:
            if bc not in files:
                missing_files.append(bc)
                continue
            ifile = files[bc]
        else:
            if tai not in twisted_files:
                missing_files.append(tai)
                continue
            ifile = twisted_files[tai]

        od = f"{outdir}/{pn}"

        makedirs(od, exist_ok=True)
        fn = f"{od}/{sn}_S000_L00{lane}_R1_001"
        cmd = (
            f"bam2fastq -j {threads} -o {fn} -c 9 "
            f"{ifile}; "
            f"fqtools count {fn}.fastq.gz > "
            f"{fn}.counts.txt"
        )
        commands.append(cmd)

    if missing_files:
        raise ValueError(f"{run_folder} is missing barcodes: {missing_files}")

    for cmd in commands:
        print(cmd)
