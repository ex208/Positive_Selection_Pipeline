#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# back_translations.py
# (c) The James Hutton Institute 2016
# Author: Eirini Xemantilotou
# pylint: disable=C0103
"""
It generates the back translations for given MSA looking for the sequence
in the dickeya_cds_nt.fasta file
The fasta file is fixed to be the 29 Dickeya nucleotide genomes
"""
import os
import sys
from Bio import SeqIO

# Input fasta file from which we can retrieve info
fasta_file = sys.argv[1]


# Wanted file must be the rbbh fasta file we want to obtain the back
# translation for
# Input sequence IDs of interest
wanted_file = sys.argv[2]
infstem = os.path.splitext(wanted_file)[0]

# Output fasta file
result_file = infstem + '.back_translations.fasta'

# Lists for storing dna and protein sequences
dna_seqlist = list(SeqIO.parse(fasta_file, 'fasta'))
protein_fastalist = list(SeqIO.parse(wanted_file, 'fasta'))

#print(len(dna_seqlist))
#print(len(protein_fastalist))

# Store in a list those ids from the protein rbbh fasta file
# and search for them in the nucleitide fasta file
# Get the corresponding DNA sequence and write it out at the result file
seqs = []

with open(result_file, "w") as f:
    for protein in protein_fastalist:
        for dna in dna_seqlist:
            if dna.id == protein.id:
                seqs.append(dna)
                SeqIO.write([dna], f, "fasta")


# Remove the stop codon as it is included as * in the protein
# sequence and therefore its existance geenrates issus as
# protein.seq != dna.seq
back_tra_file = infstem + '_backtrans_nostops.fasta'
output_handle = open(back_tra_file, "w")
dna_seq_back_translations = list(SeqIO.parse(result_file, 'fasta'))
no_stops = (s[:-3] for s in dna_seq_back_translations)
dna_no_stops = SeqIO.write(no_stops, output_handle, 'fasta')
sys.exit(0)

