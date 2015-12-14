#!/usr/bin/env python
"""
Count multiplexed shRNA reads.


"""
import argparse
import collections
import logging
import os
import itertools
import json
import skbio
import numpy as np
import pandas as pd


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
    library = pd.read_csv(path, sep='\t', names=["seq"])
    library = library['seq']
    library_trim = library.str[library_slice]
    for seq in find_multiple(library.values):
        logging.warn("Library sequence is not unique: %s" % seq)
    for seq in find_multiple(library_trim.values):
        logging.warn("Trimmed library sequence is not unique: %s" % seq)

    return library, library_trim


def read_barcodes(path):
    barcodes = pd.read_csv(path, sep='\t', names=["barcode"])
    barcodes = barcodes['barcode']
    barcodes = barcodes.str.upper()
    if len(barcodes.unique()) != len(barcodes):
        raise ValueError("Barcodes are not unique.")

    return barcodes


class SplitWriter:
    """Sort reads into a set of fastq files."""
    def __init__(self, basedir, name_template):
        self._template = name_template
        self._basedir = basedir
        self._files = {}
        os.mkdir(basedir)

    def write(self, read, *args, **kwargs):
        fileobj = self._get_file(*args, **kwargs)
        read.write(fileobj, format='fastq', variant="sanger")

    def _get_file(self, *args, **kwargs):
        name = self._template.format(*args, **kwargs)
        if name in self._files:
            return self._files[name]
        path = os.path.join(self._basedir, name)
        self._files[name] = open(path, 'x')
        return self._files[name]

    def __enter__(self):
        return self

    def __exit__(self, *args):
        for fileobj in self._files.values():
            fileobj.close()


def count_reads(reads, library, barcodes, barcode_range, seq_range,
                split_writer=None):
    library_lookup = dict((v, i) for i, v in enumerate(library.values))
    barcode_lookup = dict((v, i) for i, v in enumerate(barcodes.values))

    count_data = np.zeros((len(library), len(barcodes)), dtype=int)
    counts = pd.DataFrame(count_data,
                          index=library.index,
                          columns=barcodes.values)

    count_no_barcode = 0
    count_no_source = 0

    for read in reads:
        barcode_idx = barcodes.get(str(read[barcode_range]), None)
        source_idx = library.get(str(read[seq_range]), None)

        if source_idx is not None and barcode_idx is not None:
            counts.iat[source_idx, barcode_idx] += 1

        if barcode_idx is None:
            count_no_barcode += 1
            barcode = "no_valid_barcode"
        else:
            barcode = barcodes.iloc[barcode_idx]

        if source_idx is None:
            count_no_source += 1
            source = "unknown"
        else:
            source = ""

        if split_writer is not None:
            split_writer.write(read, barcode=barcode, source=source)

    stats = {
        'count_no_barcode': count_no_barcode,
        'count_not_in_lib': count_no_source,
    }

    return counts, stats


def main():
    args = gen_argparse().parse_args()

    library, library_trim = read_library(args.library, args.lib_range)
    barcodes = read_barcodes(args.barcodes)
    reads = [skbio.io.read(input, format='fastq', variant='sanger')
             for input in args.input]
    reads = itertools.chain(*reads)

    if args.write_split:
        template = "reads_{barcode}_{source}.fastq"
        with SplitWriter(args.write_split, template) as writer:
            counts, stats = count_reads(
                reads, library, barcodes, args.barcode_range,
                args.seq_range, writer
            )
    else:
        counts, stats = count_reads(
            reads, library, barcodes, args.barcode_range, args.seq_range
        )

    counts.to_excel(args.output)
    if args.stats:
        with open(args.stats, 'w') as fileobj:
            json.dump(stats, fileobj, indent=4)


if __name__ == '__main__':
    main()
