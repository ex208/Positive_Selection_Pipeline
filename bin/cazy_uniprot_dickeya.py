#!/usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=C0103
"""cazy_uniprot_dickeya.py

(c) The James Hutton Institute 2016
Author: Eirini Xemantilotou

This script does the following:
[1] Retrieve UniProt Dickeya entries with a cross-reference to CAZy
       Retrieves data in tabular format:
       accession  CAZy_xrefs
       Example:
       P0C1A7	PL9;
       P07103	CBM5;GH5;
[2] Create a dictionary mapping CAZy families to accessions
       e.g.
       {GH3: [P11073, D2BXL2, ...], ...}
       We use a default dict for this. As `caz` is just a string we use
       io.StringIO to handle it like a file whose lines we iterate over.
       StringIO reads and writes strings as files.
       As the first line contains the column headers, we ignore it.
       Each line is split into accession and a list of CAZy families; from
       this we construct the dictionary.
[3] For each CAZy family, retrieve the mapped UniProtKB accessions in fasta
       format and write to a file with name <CAZy_family>.fasta. My simple
       approach here retrieves each batch of sequences from uniprot.org. One
       could also download all Dickeya sequences in one go and then get them
       from that file. The bioservices package also does not provide the
       batch download facility that the UniProt API has which could
       be used here.

Python script for getting sequences in FASTA format out from Uniprot, where
the taxonomy filter is `Dickeya` and the entry has a CAZy database cross-
reference.

Sequences are written to FASTA files named by the corresponding CAZy family
"""

import logging
import os
import shutil
import sys

from collections import defaultdict

import bioservices

# Set up logger
logging.basicConfig(filename='examples.log', level=logging.DEBUG)
logger = logging.getLogger("cazy_uniprot_dickeya.py")

# create handlers
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.INFO)

# create a logging format
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)

# add the handlers to the logger
logger.addHandler(handler)


# Define and create output directory
outdir = os.path.join('data', 'cazy_dickeya')
os.makedirs(outdir, exist_ok=True)
logger.info("Created %s", outdir)

# Make query at UniProt
logger.info("Running query at UniProt")
logger.debug("This is debug level")
uhandle = bioservices.UniProt()
caz = uhandle.search('database:(type:cazy) taxonomy:dickeya',
                     frmt='tab',
                     columns='id, database(CAZy)')

# Generate list of (accession, CAZy families) tuples from returned result
entries = [line.strip() for line in caz.split('\n')[1:] \
           if len(line.strip())]
logger.info("%d entries returned", len(entries))

# Parse returned entries into family: [accession1, accession2] dictionary
caz_map = defaultdict(list)  # key: CAZy family; value: list of accessions
for entry in entries:
    acc, xref = entry.split()
    fams = xref.split(';')[:-1]
    for fam in fams:
        caz_map[fam].append(acc)
logger.info("Grouped entries into %d families", len(caz_map))

# Obtain accession sequences from UniProt
for fam, members in caz_map.items():
    logger.info('Working on CAZy family: %s', fam)
    sequences = []  # Holds list of FASTA sequences for family accessions
    for accession in members:
        logger.debug("Getting sequence for %s", accession)
        try:
            sequences.append(uhandle.search(accession, frmt='fasta'))
        except:
            logger.warning("Could not get sequence for %s - skipping",
                            accession)
    outfname = os.path.join(outdir, '%s.fasta' % fam)
    with open(outfname, 'w') as outfile:
        outfile.write('\n'.join(sequences))

logger.info("Complete")
