#!/usr/bin/env python
import argparse
import collections
import logging
import os
import itertools
import json
import fileinput
import numpy as np
import pandas as pd
from datetime import datetime
import cytoolz as toolz
import gzip


def parse_slice(str_slice):
    """Parse a string representation of a python slice."""
    parts = str_slice.split(':')
    if len(parts) < 2 or len(parts) > 3:
        raise argparse.ArgumentTypeError("Invalid slice: %s" % str_slice)
    parts_type = []
    for part in parts:
        if not part:
            parts_type.append(None)
            continue
        try:
            parts_type.append(int(part))
        except ValueError:
            raise argparse.ArgumentTypeError("Invalid slice: %s" % str_slice)
    return slice(*parts_type)


def gen_argparse():
    parser = argparse.ArgumentParser(
        """Count shRNA seq reads."""
    )
    parser.add_argument("-l", "--library", type=str, required=True,
                        help="The library of snRNAs. (tsv)")
    parser.add_argument("-b", "--barcodes", type=str, required=True,
                        help="List of barcodes for demultiplexing")
    parser.add_argument("--barcode-range", default="-3:", type=parse_slice,
                        help="The slice of the reads that contain the barcode")
    parser.add_argument("--lib-range", default=":", type=parse_slice,
                        help="Only use a subsequence of the library sequences."
                        " You can specify the range as python slice")
    parser.add_argument("--seq-range", default=":", type=parse_slice,
                        help="Only use a subsequence of the reads. "
                        "You can specify the range as python slice.")
    parser.add_argument("--write-split", default=None, type=str,
                        help="Write all reads into files inside this "
                        "directory according to the barcode they contain")
    parser.add_argument("--stats", default=None, type=str,
                        help="Write statistics to json file.")
    parser.add_argument("input", nargs='+',
                        help="Fastq file(s) containing the reads.")
    parser.add_argument("output", help="Name of the output excel file.")
    return parser


def find_multiple(seq):
    """Return a list of elements that are not unique."""
    counter = collections.Counter(seq)
    return [key for key, count in counter.items() if count > 1]


def read_library(path, library_slice):
    library = pd.read_csv(path, sep='\t', names=["seq"], header=None)
    library = library['seq'].str.upper()

    if not library.index.is_unique:
        raise ValueError("Index of library sequences is not unique.")
    for seq in find_multiple(library.values):
        logging.warn("Library sequence is not unique: %s" % seq)
    library = library.drop_duplicates()

    library_trim = library.str[library_slice]
    #if find_multiple(library_trim.values):
    #    raise ValueError("Library is not unique after trimming.")

    return library.astype(bytes), library_trim.astype(bytes)


def read_barcodes(path):
    barcodes = pd.read_csv(path, sep='\t', names=["barcode"])
    barcodes = barcodes['barcode']
    barcodes = barcodes.str.upper()
    barcodes = barcodes.astype(bytes)
    if len(barcodes.unique()) != len(barcodes):
        raise ValueError("Barcodes are not unique.")
    if not barcodes.index.is_unique:
        raise ValueError("Barcode index is not unique.")

    return barcodes


class SplitWriter:
    """Sort reads into a set of fastq files."""
    def __init__(self, basedir, name_template):
        self._template = name_template
        self._basedir = basedir
        self._files = {}
        try:
            os.mkdir(basedir)
        except FileExistsError:
            pass

    def write(self, read, *args, **kwargs):
        fileobj = self._get_file(*args, **kwargs)
        for line in read:
            fileobj.write(line)

    def _get_file(self, *args, **kwargs):
        name = self._template.format(*args, **kwargs)
        if name in self._files:
            return self._files[name]
        path = os.path.join(self._basedir, name)
        self._files[name] = open(path, 'xb')
        return self._files[name]

    def __enter__(self):
        return self

    def __exit__(self, *args):
        for fileobj in self._files.values():
            fileobj.close()


def count_reads(reads, library, barcodes, loop_rc, adapters,
                oligo_offset, barcode_offset, adapter_offsets,
                expect_loop_idx, loop_idx_precision=2, split_writer=None,
                barcode_len=3, oligo_len=19):
    """ Count the number of reads that contain oligo sequences.

    Parameters
    ----------
    reads:
        Iterator of 4-tuples (the 4 lines of the fastq)
    library: pd.
        The oligo sequences.
    barcodes: pd.Series
        Barcodes for demultiplexing.
    loop_rc: upper case str
        The loop sequence should be in every read. All positions in the
        read are seen relative to this loop sequence.
    adapters: list of upper case str
        A list of possible adapters.
    *_offset: int
        Offsets for the sequences relative to the loop location.
    """

    assert len(adapters) == len(adapter_offsets)

    library_lookup = dict((v, i) for i, v in enumerate(library.values))
    barcode_lookup = dict((v, i) for i, v in enumerate(barcodes.values))

    count_data = np.zeros((len(library), len(barcodes)), dtype=int)
    counts = pd.DataFrame(count_data,
                          index=library.index,
                          columns=barcodes.index)

    count_vals = counts.values
    barcode_vals = list(barcodes.values)

    num_reads = 0
    adapter_counts = [0 for i in adapters]
    pos_dist = collections.Counter()
    oligo_counts = collections.Counter()
    no_loop = 0
    no_adapter = 0
    no_barcode = 0
    no_oligo = 0
    bad_loop_idx = 0

    for read in reads:
        num_reads += 1
        read_str = read[1].strip().upper()
        loop_idx = read_str.find(loop_rc)

        if loop_idx == -1:
            no_loop += 1
            continue

        pos_dist[loop_idx] += 1
        if abs(loop_idx - expect_loop_idx) > loop_idx_precision:
            bad_loop_idx += 1
            continue

        for i, (adapter, offset) in enumerate(zip(adapters, adapter_offsets)):
            if read_str[loop_idx + offset:].startswith(adapter):
                adapter_counts[i] += 1
                break
        else:
            no_adapter += 1
            continue

        barcode_idx = loop_idx + barcode_offset
        barcode = read_str[barcode_idx:barcode_idx+barcode_len]
        barcode_key = barcode_lookup.get(barcode, None)

        if barcode_key is None:
            no_barcode += 1
            if split_writer is not None:
                split_writer.write(read, barcode=barcode.decode())
            continue

        oligo_idx = loop_idx + oligo_offset
        oligo = read_str[oligo_idx:oligo_idx+oligo_len]
        oligo_key = library_lookup.get(oligo, None)

        if oligo_key is None:
            oligo_counts[oligo] += 1
            no_oligo += 1
            if split_writer is not None:
                split_writer.write(read, barcode="unknown_oligo")
            continue

        count_vals[oligo_key, barcode_key] += 1

    stats = {
        'num_reads': num_reads,
        'no_loop': no_loop,
        'loop__no_adapter': no_adapter,
        'loop__adapter__no_barcode': no_barcode,
        'loop__adapter__barcode__no_oligo': no_oligo,
        'bad_loop_idx': bad_loop_idx,
        'loop_pos_dist': pos_dist,
        'adapter_counts': adapter_counts,
    }

    return counts, stats, oligo_counts


def lines(name):
    if name == '-':
        yield from sys.stdin.buffer
    if name.endswith('.gz'):
        with gzip.open(name, 'rb') as file:
            yield from file
    with open(name, 'rb') as file:
        yield from file


def read_fastq(filename, slice=None):
    res = toolz.partition(4, lines(filename))
    if slice is not None:
        return itertools.islice(res, slice)
    return res


def main():
    args = gen_argparse().parse_args()

    library, library_trim = read_library(args.library, args.lib_range)
    barcodes = read_barcodes(args.barcodes)

    def lines(name):
        if name == '-':
            yield from sys.stdin.buffer
        if name.endswith('.gz'):
            with gzip.open(name, 'rb') as file:
                yield from file
        with open(name, 'rb') as file:
            yield from file

    reads = itertools.chain.from_iterable(lines(name) for name in args.input)
    reads = toolz.partition(4, reads)
    """reads = itertools.islice(reads, 10000000)"""

    if args.write_split:
        template = "reads_{barcode}_{source}.fastq"
        with SplitWriter(args.write_split, template) as writer:
            counts, stats = count_reads(
                reads, library_trim, barcodes, args.barcode_range,
                args.seq_range, writer
            )
    else:
        counts, stats = count_reads(
            reads, library_trim, barcodes, args.barcode_range,
            args.seq_range
        )

    counts.to_excel(args.output)

    counts.index.name = "gene"
    counts.columns.name = "barcode"

    counts = counts.unstack()
    counts.name = "count"

    groups = counts.reset_index().groupby("barcode")
    by_barcode = dict(
        (
            key,
            val.sort_values(by=["count", "gene"], ascending=False)
                    .reset_index()[["gene", "count"]]
        )
        for key, val in groups
    )
    counts_sorted = pd.concat(by_barcode, axis=1)
    counts_sorted.to_excel("counts_sorted.xlsx")

    if args.stats:
        stats['date'] = datetime.now().isoformat()
        with open(args.stats, 'w') as fileobj:
            json.dump(stats, fileobj, indent=4)


if __name__ == '__main__':
    main()
