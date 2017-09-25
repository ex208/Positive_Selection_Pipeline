#!/usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=C0103
#automation.py
#
# (c) The James Hutton Institute 2016
#Author: Eirini Xemantilotou

# This script does the following:
# [1] Retrieve UniProt Dickeya entries with a cross-reference to CAZy
#       Retrieves data in tabular format:
#       accession  CAZy_xrefs
#       Example:
#       P0C1A7	PL9;
#       P07103	CBM5;GH5;
# [2] Create a dictionary mapping CAZy families to accessions
#       e.g.
#       {GH3: [P11073, D2BXL2, ...], ...}
#       We use a default dict for this. As `caz` is just a string we use
#       io.StringIO to handle it like a file whose lines we iterate over.
#       StringIO reads and writes strings as files.
#       As the first line contains the column headers, we ignore it.
#       Each line is split into accession and a list of CAZy families; from
#       this we construct the dictionary.
# [3] For each CAZy family, retrieve the mapped UniProtKB accessions in fasta
#       format and write to a file with name <CAZy_family>.fasta. My simple
#       approach here retrieves each batch of sequences from uniprot.org. One
#       could also download all Dickeya sequences in one go and then get them
#       from that file. The bioservices package also does not provide the
#       batch download facility that the UniProt API has which could
#       be used here.
#
"""
Python script for getting sequences in FASTA format out from Uniprot, where
the taxonomy filter is `Dickeya` and the entry has a CAZy database cross-
reference.

Sequences are written to FASTA files named by the corresponding CAZy family
"""
import os
import shutil
from io import StringIO
from collections import defaultdict
import bioservices


os.chdir("../data")
os.mkdir("./cazy_dickeya")

u = bioservices.UniProt()

caz = u.search('database:(type:cazy) taxonomy:dickeya',
               frmt='tab',
               columns='id, database(CAZy)')

caz_map = defaultdict(list)
caz_file = StringIO(caz)
caz_iter = iter(caz_file)
next(caz_iter)  # this takes the iterator to line two
for line in caz_iter:
    acc, xref = line.split()
    fams = xref.split(';')[:-1]
    for fam in fams:
        caz_map[fam].append(acc)

for fam, members in caz_map.items():
    print('Working on CAZy family: {}'.format(fam))
    q = ' OR '.join(members)
    fa = u.search(q, frmt='fasta')
    fname = '{}.fasta'.format(fam)
    destination = "./cazy_dickeya"
    fasta = os.path.join(destination, fname)
    with open(fasta, 'w', encoding='ascii') as outfile:
        outfile.write(fa)
