# Why
# 1) cloudspades requires the index reads be inline in the record header
# 2) Ariadne requires the data are sorted by the barcodes
#
# Inlining is easy. Sorting is complex as the amount of data is large, and
# the ordering stems is determined external to the data being sorted. To
# determine order, all barcodes must be read in to gather the complete
# barcode <-> record association; if only partial data is read then
# associations to barcodes may be missed, and we cannot perform an insertion
# sort efficiently as we're writing to disk. Once we know an order for the
# records, we (currently) read in the entirety of the subsequent data (R1 then
# R2), reorder, and write. Performing this in blocks to minimize memory may be
# possible, but we have to assume access is random as a grouping barcode
# may be with any record along the file.
#
# A variety of approaches were considered, including:
# - indexing portions in a hashtable, reading inputs multiple times, and
#   writing in blocks. This was tested in both rust and python. The amount of
#   memory was large, and keeping it under control would be many many many
#   passes over data on disk or in memory
# - using pandas to do the grouping, which possibly avoids the memory burden
#   of a hashmap. it didn't
# - using mmap files. No go, these are large and we have to walk over them
#   a lot.
#
# Parsing this stuff adds a lot of overhead in Python. It will add some, if not
# a lot, in rust as well -- our test data had 65M sequences. So the current
# approach operates in the raw file data itself, using regex's to parse
# individual records. We use numpy for sorting and getting record orders.
# This is memory expensive but so far much less than the other approaches tried
# and it does not require multiple passes over files. We bottleneck on write
# IO, so to mitigate that, we are using a parallel gzip (pgzip), which still
# bottlenecks but gets better throughput.
#
# There probably are smarter ways to do this to reduce the memory burden.
# Right now, it's O(N) where N is the number of records. We load R1 and R2
# separately though so we at least halve the memory use. As for doing it
# faster, at the moment we appear to saturate time on gzip. Easiest solution
# would be to increase the number of threads, but then again, this process
# is expected to run in an array, and filesystem can only take so much.
#
# In addition to the inline tests, md5 checks to verify all record IDs are
# present in both R1 / R2, and relative to original input. Spot checks on
# an arbitrary set of records were performed on R1 / R2 to verify no apparent
# unusual modification. And spot checks were performed to verify that correct
# barcodes are incorporating as expected in output.
#
# author: Daniel McDonald (d3mcdonald@eng.ucsd.edu)
import numpy as np
import click
import re
import io
import pgzip
import gzip


RECORD = re.compile(rb'@\S+\n[ATGCN]+\n\+\n\S+\n')
BARCODE = re.compile(rb'@\S+\n([ATGCN]+)\n\+\n\S+\n')


def gather_order(i1_in_fp):
    """Determine record order

    This is a fancy way of saying: get all the barcodes, and sort them.

    We return the order of the sorted records, the unique barcodes,
    and the bounds for what barcode associated with what record
    """
    # determine barcode length
    _ = i1_in_fp.readline()
    b = i1_in_fp.readline()
    rec_len = len(b.strip())
    i1_in_fp.seek(0)

    # we need larger data in memory later anyway...
    i1 = i1_in_fp.read()
    start = 0
    end = len(i1)

    # get the number of records. we completely assume non-multiline fastq here
    newlines = i1.count(b'\n')
    assert newlines % 4 == 0
    barcodes = np.empty(newlines // 4, dtype='|S%d' % rec_len)

    # walk all index records
    # grab each barcode
    idx = 0
    while start < end:
        barcode_result = BARCODE.search(i1, pos=start)
        barcode = barcode_result.groups()[0]
        assert len(barcode) == rec_len  # get angry if the barcode is weird

        barcodes[idx] = barcode
        idx += 1
        start = barcode_result.end()

    # we no longer need the raw data so let's toss it
    del i1

    # determine the record order of a lexicographic sort
    # gather the unique barcodes so we can use them later, and the bounding
    # points in the sorted set
    record_order = barcodes.argsort()
    barcodes = barcodes[record_order]
    unique_barcodes, barcode_bounds = np.unique(barcodes, return_index=True)

    return record_order, unique_barcodes, barcode_bounds


def test_gather_order():
    i1data = [b'@foo', b'ATGC', b'+', b'!!!!',
              b'@bar', b'TTGG', b'+', b'!!!!',
              b'@baz', b'ATGC', b'+', b'!!!!',
              b'@oof', b'TTTT', b'+', b'!!!!',
              b'@rab', b'TTGG', b'+', b'!!!!',
              b'@zab', b'TTTT', b'+', b'!!!!',
              b'@ofo', b'TTTT', b'+', b'!!!!', b'']

    i1 = io.BytesIO(b'\n'.join(i1data))
    order, unique, bounds = gather_order(i1)

    exp_order = np.array([0, 2, 1, 4, 3, 5, 6])
    exp_unique = np.array([b'ATGC', b'TTGG', b'TTTT'])
    exp_bounds = np.array([0, 2, 4])

    assert (order == exp_order).all()
    assert (unique == exp_unique).all()
    assert (bounds == exp_bounds).all()


def troll_and_write(order, unique, bounds, in_, out_):
    """Walk over the raw data, spit out barcode amended records in order

    - read all data
    - get index boundaries for each record
    - pull out each record in order according to the barcode data
    - associate the barcode
    - write
    """

    data = in_.read()
    boundaries = np.empty([order.size, 2], dtype=np.uint64)

    stop = 0
    for idx in range(order.size):
        rec = RECORD.search(data, pos=stop)
        start, stop = rec.span()
        boundaries[idx] = np.array([start, stop], dtype=np.uint64)

    current_barcode_idx = 0
    current_barcode = unique[current_barcode_idx]
    current_barcode_bound_end = bounds[current_barcode_idx + 1]

    for order_idx, record_idx in enumerate(order):
        if order_idx >= current_barcode_bound_end:
            current_barcode_idx += 1

            if current_barcode_idx >= bounds.size:
                raise ValueError("should not happen?")
            current_barcode = unique[current_barcode_idx]

            if current_barcode_idx + 1 >= bounds.size:
                # run to the end
                current_barcode_bound_end = order.size
            else:
                current_barcode_bound_end = bounds[current_barcode_idx + 1]

        start, stop = boundaries[record_idx]
        record = data[start:stop]

        # in a one-off, these might pass by chance. It would be real weird
        # for them to always pass for all records in a large file.
        # n.b., b'foo'[0] is int, because yay, so we use a slice to maintain
        # a human readable character to test against as most mortals haven't
        # memorized the ascii table
        assert record[:1] == b'@'
        assert record[-1:] == b'\n'

        with_barcode = insert_barcode(record, current_barcode)
        out_.write(with_barcode)


def test_troll_and_write():
    i1data = [b'@foo', b'ATGC', b'+', b'!!!!',
              b'@bar', b'TTGG', b'+', b'!!!!',
              b'@baz', b'ATGC', b'+', b'!!!!',
              b'@oof', b'TTTT', b'+', b'!!!!',
              b'@rab', b'TTGG', b'+', b'!!!!',
              b'@zab', b'TTTT', b'+', b'!!!!',
              b'@ofo', b'TTTT', b'+', b'!!!!', b'']

    i1 = io.BytesIO(b'\n'.join(i1data))
    order, unique, bounds = gather_order(i1)

    # we assume records are in the same order, as that has previously been
    # observed w/ tellread and is the normal expectation
    r1data = [b'@foo', b'AATGC', b'+', b'!!!!!',
              b'@bar', b'ATTGG', b'+', b'!!!!!',
              b'@baz', b'AATGC', b'+', b'!!!!!',
              b'@oof', b'ATTTT', b'+', b'!!!!!',
              b'@rab', b'ATTGG', b'+', b'!!!!!',
              b'@zab', b'ATTTT', b'+', b'!!!!!',
              b'@ofo', b'ATTTT', b'+', b'!!!!!', b'']
    r1 = io.BytesIO(b'\n'.join(r1data))
    r1out = io.BytesIO()
    troll_and_write(order, unique, bounds, r1, r1out)
    r1out.seek(0)

    r1exp = [b'@foo BX:Z:ATGC-1', b'AATGC', b'+', b'!!!!!',
             b'@baz BX:Z:ATGC-1', b'AATGC', b'+', b'!!!!!',
             b'@bar BX:Z:TTGG-1', b'ATTGG', b'+', b'!!!!!',
             b'@rab BX:Z:TTGG-1', b'ATTGG', b'+', b'!!!!!',
             b'@oof BX:Z:TTTT-1', b'ATTTT', b'+', b'!!!!!',
             b'@zab BX:Z:TTTT-1', b'ATTTT', b'+', b'!!!!!',
             b'@ofo BX:Z:TTTT-1', b'ATTTT', b'+', b'!!!!!',
             b'']
    r1exp = b'\n'.join(r1exp)
    assert r1exp == r1out.read()


def create_tag(t):
    return b'BX:Z:%s-1' % t


def create_tag_no_suffix(t):
    return b'BX:Z:%s' % t


def insert_barcode(record, barcode):
    """Get the current ID, smash the needed tag in"""
    # @foo\nATGC\n+\n!!!!\n
    id_, remainder = record.split(b'\n', 1)
    tag = create_tag(barcode)
    return b'%s %s\n%s' % (id_, tag, remainder)


def readfq(fp):
    if fp.mode == 'rb':
        strip = bytes.strip
    else:
        strip = str.strip

    id_ = iter(fp)
    seq = iter(fp)
    dumb = iter(fp)
    qual = iter(fp)
    for rec in zip(id_, seq, dumb, qual):
        yield list(map(strip, rec))


def writefq(rec, out):
    for item in rec:
        out.write(item)
        out.write(b'\n')


@click.group()
def cli():
    pass


@cli.command()
def tests():
    test_gather_order()
    test_troll_and_write()


@cli.command()
@click.option('--r1-in', type=click.Path(exists=True), required=True)
@click.option('--r2-in', type=click.Path(exists=True), required=True)
@click.option('--i1-in', type=click.Path(exists=True), required=True)
@click.option('--r1-out', type=click.Path(exists=False), required=True)
@click.option('--r2-out', type=click.Path(exists=False), required=True)
@click.option('--threads', type=int, required=False, default=1)
@click.option('--no-sort', is_flag=True, default=False)
def integrate(r1_in, r2_in, i1_in, r1_out, r2_out, threads, no_sort):
    r1_in_fp = open(r1_in, 'rb')
    r2_in_fp = open(r2_in, 'rb')
    i1_in_fp = open(i1_in, 'rb')

    if no_sort:
        r1_out_fp = gzip.open(r1_out, mode='wb')
        r2_out_fp = gzip.open(r2_out, mode='wb')

        r1_sniff = r1_in_fp.readline().strip()
        r2_sniff = r2_in_fp.readline().strip()
        r1_in_fp.seek(0)
        r2_in_fp.seek(0)

        # outputs from tellread don't seem to have orientation information
        # some downstream programs hate this, so let's add if needed.
        if r1_sniff.endswith(b'/1'):
            if not r2_sniff.endswith(b'/2'):
                raise ValueError('unexpected endings: '
                                 f'{r1_sniff.decode("utf-8")} '
                                 f'{r2_sniff.decode("utf-8")}')
            orient_r1 = ''
            orient_r2 = ''
        else:
            assert b'/1' not in r1_sniff

            orient_r1 = b'/1'
            orient_r2 = b'/2'

        for (r1, r2, i1) in zip(*map(readfq, [r1_in_fp, r2_in_fp, i1_in_fp])):
            assert r1[0] == r2[0]
            assert r1[0] == i1[0]

            tag = create_tag_no_suffix(i1[1])
            r1[0] = b"%s%s %s" % (r1[0], orient_r1, tag)
            r2[0] = b"%s%s %s" % (r2[0], orient_r2, tag)
            writefq(r1, r1_out_fp)
            writefq(r2, r2_out_fp)
        r1_out_fp.close()
        r2_out_fp.close()
    else:
        # 200MB is what they use in their readme...
        r1_out_fp = pgzip.open(r1_out, mode='wb', thread=threads,
                               blocksize=2*10**8)
        r2_out_fp = pgzip.open(r2_out, mode='wb', thread=threads,
                               blocksize=2*10**8)

        order, unique, bounds = gather_order(i1_in_fp)

        for in_, out_ in zip([r1_in_fp, r2_in_fp], [r1_out_fp, r2_out_fp]):
            troll_and_write(order, unique, bounds, in_, out_)
            in_.close()
            out_.close()


if __name__ == '__main__':
    cli()
