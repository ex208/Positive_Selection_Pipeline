#!/usr/bin/env python
# pylint: disable=C0103
# clustal_to_rphylip.py
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard
"""
This python script uses Biopython to convert MSA format from clustal to
relaxed phylip-relaxed
(this preserves filename lengths)
"""
import os
import sys
from Bio import AlignIO

infname = sys.argv[1]
infstem = os.path.splitext(infname)[0]

print("Parsing %s CLUSTAL file" % infname)
data = list(AlignIO.parse(infname, "clustal"))

outfname = infstem + '.rphylip'
print("Writing %s relaxed PHYLIP output" % outfname)
AlignIO.write(data, outfname, 'phylip-relaxed')

sys.exit(0)
