import glob
import gzip
import os
from sequence_processing_pipeline.util import (iter_paired_files,
                                               determine_orientation)


def split_similar_size_bins(data_location_path, max_file_list_size_in_gb,
                            batch_prefix):
    '''Partitions input fastqs to coarse bins

    :param data_location_path: Path to the ConvertJob directory.
    :param max_file_list_size_in_gb: Upper threshold for file-size.
    :param batch_prefix: Path + file-name prefix for output-files.
    :return: The number of output-files created, size of largest bin.
    '''
    # to prevent issues w/filenames like the ones below from being mistaken
    # for R1 or R2 files, use determine_orientation().
    # LS_8_22_2014_R2_SRE_S2_L007_I1_001.fastq.gz
    # LS_8_22_2014_R1_SRE_S3_L007_I1_001.fastq.gz

    # since the names of all fastq files are being scanned for orientation,
    # collect all of them instead of mistakenly pre-filtering some files.
    fastq_paths = glob.glob(data_location_path + '/*/*.fastq.gz')
    fastq_paths = [x for x in fastq_paths
                   if determine_orientation(x) in ['R1', 'R2']]

    # convert from GB and halve as we sum R1
    max_size = (int(max_file_list_size_in_gb) * (2 ** 30) / 2)

    split_offset = 0

    # ensure we are max-sized to start.
    current_size = max_size * 10
    fp = None

    bucket_size = 0
    max_bucket_size = 0

    for a, b in iter_paired_files(fastq_paths):
        r1_size = os.stat(a).st_size
        r2_size = os.stat(b).st_size

        output_base = os.path.dirname(a).split('/')[-1]
        if current_size + r1_size > max_size:
            # bucket is full.
            if bucket_size > max_bucket_size:
                max_bucket_size = bucket_size

            # reset bucket_size.
            bucket_size = r1_size + r2_size

            if fp is not None:
                fp.close()

            split_offset += 1
            current_size = r1_size
            fp = open(batch_prefix + '-%d' % split_offset, 'w')
        else:
            # add to bucket_size
            bucket_size += r1_size + r2_size
            current_size += r1_size

        fp.write("%s\t%s\t%s\n" % (a, b, output_base))

    if fp is not None:
        fp.close()

    if split_offset == 0:
        raise ValueError("No splits made")

    return split_offset, max_bucket_size


def demux_cmd(id_map_fp, fp_fp, out_d, task, maxtask):
    with open(id_map_fp, 'r') as f:
        id_map = f.readlines()
        id_map = [line.strip().split('\t') for line in id_map]

    # fp needs to be an open file handle.
    # ensure task and maxtask are proper ints when coming from cmd-line.
    with open(fp_fp, 'r') as fp:
        demux(id_map, fp, out_d, int(task), int(maxtask))


def demux(id_map, fp, out_d, task, maxtask):
    """Split infile data based in provided map"""
    delimiter = '::MUX::'
    mode = 'wt'
    ext = '.fastq.gz'
    sep = '/'
    rec = '@'

    openfps = {}

    for offset, (idx, r1, r2, outbase) in enumerate(id_map):
        if offset % maxtask == task:
            idx = rec + idx

            # setup output locations
            outdir = out_d + sep + outbase
            fullname_r1 = outdir + sep + r1 + ext
            fullname_r2 = outdir + sep + r2 + ext

            os.makedirs(outdir, exist_ok=True)
            current_fp_r1 = gzip.open(fullname_r1, mode)
            current_fp_r2 = gzip.open(fullname_r2, mode)
            current_fp = {'1': current_fp_r1, '2': current_fp_r2}
            openfps[idx] = current_fp

    # setup a parser
    seq_id = iter(fp)
    seq = iter(fp)
    dumb = iter(fp)
    qual = iter(fp)

    for i, s, d, q in zip(seq_id, seq, dumb, qual):
        # '@1', 'LH00444:84:227CNHLT4:7:1101:41955:2443/1'
        # '@1', 'LH00444:84:227CNHLT4:7:1101:41955:2443/1 BX:Z:TATGACACATGCGGCCCT' # noqa
        # '@baz/1

        # NB: from 6d794a37-12cd-4f8e-95d6-72a4b8a1ec1c's only-adapter-filtered results: # noqa
        # @A00953:244:HYHYWDSXY:3:1101:14082:3740 1:N:0:CCGTAAGA+TCTAACGC

        fname_encoded, sid = i.split(delimiter, 1)

        if fname_encoded not in openfps:
            continue

        current_fp = openfps[fname_encoded]

        # remove '\n' from sid and split on all whitespace.
        tmp = sid.strip().split()

        if len(tmp) == 1:
            # sequence id line contains no optional metadata.
            # don't change sid.
            # -1 is \n
            orientation = sid[-2]
            sid = rec + sid
        elif len(tmp) == 2:
            sid = tmp[0]
            metadata = tmp[1]
            # no '\n'
            orientation = sid[-1]
            # hexdump confirms separator is ' ', not '\t'
            sid = rec + sid + ' ' + metadata + '\n'
        else:
            raise ValueError(f"'{sid}' is not a recognized form")

        current_fp[orientation].write(sid)
        current_fp[orientation].write(s)
        current_fp[orientation].write(d)
        current_fp[orientation].write(q)

    for d in openfps.values():
        for f in d.values():
            f.close()
