#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#pylint: disable=C0103
#RBBH_not_det_uniprot.py
#
# (c) The James Hutton Institute 2017
#Author: Eirini Xemantilotou
"""
Python script to automate the execution of the extracy rbbh python
script for seqs/CAZy families not detected by Uniprot for Dickeya
"""
import subprocess
import os

data_dir = "../data"
with open(os.path.join(data_dir, "locus_Dickeya_not_det.txt"), 'r') as fh:
    for line in fh:
        line = line.strip()
        cmd = "python3 ./extract_rbbh.py --db dickeya.db --seqfile dickeya_cds_aa.fasta --locus_tag %s  -v -l %s.log" %(line, line)
        subprocess.call(cmd, shell=True)
