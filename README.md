# A small script for counting shRNA expression

Given a set of reads from shRNA sequencing, we count how often each
reference shRNA sequence occurs at the expected position in the
reads.

Each read should contain a barcode at a specified position. We divide
the reads according to those barcodes.

# Install

This library is pip-installable, or can be run from the source
directory.

# Formats

The reads should be in a set of `fastq`-files, that can optionally
be gzip compressed.

The shRNA-library and the barcodes should each be stored in a simple
headerless tsv with the two columns `id` and `sequence`.

We do not require library sequences to be unique. If a sequence is not
unique we print a warning and ignore all but one of the ids.

For usage see `./rnacount.py -h`.
