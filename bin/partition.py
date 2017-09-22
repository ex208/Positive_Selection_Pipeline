#!/usr/bin/env python
#
# partition.py
# # (c) The James Hutton Institute 2016
# Author: Eirini Xemantilotou

# pylint: disable=C0103

"""
Use Biopython to get the number of the nucleotide sequences for the
alignment, which can be uses at partition file for raxml
"""
import sys
from Bio import SeqIO


# Write out the file for the partition

input_file = sys.argv[1]
result_file = sys.argv[2]

# open and write the file with the right nucleotide number
with open(result_file, "w") as f:
    back_translations_nucleotides = list(SeqIO.parse(input_file, "clustal"))
    nucleotides = len(back_translations_nucleotides[0])
# Python strings use escape characters
# most commonly used are \t for tab and \n for new line. And \\ for slash.
# Using r"my string" will treat \ as plain \, but also \n as those two
# symbols
# The little r at the front is for raw string (do not parse the escape
# characters)
# Using "my\\string" will treat \\ as \
# Using second form because we want \n to be a new line
    f.write("DNA,p1=1-%i\\3, 2-%i\\3\nDNA,p2=3-%i\\3" %
            (nucleotides, nucleotides, nucleotides))
