#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#pylint: disable=C0103
#RBBH_not_det_uniprot.py
#
# (c) The James Hutton Institute 2017
#Author: Eirini Xemantilotou

"""
Python script to merge text files with locus tags for CAZ fams
were not detected by Uniprot
"""

import os


data_dir = "../data"
with open('../data/locus_Dickeya_not_det.txt', 'w') as outfile:
    for locus_txt in os.listdir(data_dir):
        if locus_txt.startswith("PL") or locus_txt.startswith("CE"):
            locus_list = locus_txt.split()
            for locus in locus_list:
                locus_dir = os.path.join(data_dir, locus)
                with open(locus_dir) as infile:
                    lines = infile.read()
                    outfile.write(lines)
f = open('../data/locus_Dickeya_not_det.txt',"r+")
d = f.readlines()
f.seek(0)
for i in d:
    if i != "A4U42_16960" + "\n":
        f.write(i)
f.truncate()
f.close()
