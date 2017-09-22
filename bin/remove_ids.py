#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# remove_ids.py
# (c) The James Hutton Institute 2016
# Author: Eirini Xemantilotou
# pylint: disable=C0103

"""
Python script to remove seqs from rbbh which are not closely related
"""

import os
import shutil
from Bio import SeqIO


def remove_ids(fasta_file, ids_remove, rbbh_dir):
    """
    Function which writes out a new fasta file excluding presenectes
    sequences which do not match
    fasta_file: the file we are accessing the seqs
    ids_remove: list of the identifiers for the seqs we want to remove
    rbbh_dir: the directory we have stored teh fasta_file
    Notice: the script works only for the CAZy families which
    has NOT split into subgroups
    """
    seqs = []
    wanted_file = os.path.join(rbbh_dir, fasta_file)
    infstem = wanted_file.split("/")[3].split(".")[0]
    result_file = infstem + "_2.rbbh.fasta"
    sequences = list(SeqIO.parse(wanted_file, 'fasta'))
    with open(result_file, "w") as f:
        for seqid in sequences:
            if seqid.id not in ids_remove:
                seqs.append(seqs)
                SeqIO.write([seqid], f, "fasta")
    shutil.move(os.path.join(".", result_file), rbbh_dir)


# CBM32
ids_remove = ["MK16_00577", "MK10_00583", "GBBC2040_00584", "IPO_2222_00582",
              "NCPPB_2511_03722", "Dd703_3437", "DW_0440_04032", "CSL_RW240_04053"]
remove_ids("CBM32.rbbh.fasta", ids_remove, "../data/cazy_rbbh")

#CBM63
ids_remove_1 = ["MK16_02185", "CSL_RW240_01936", "GBBC2039_02073", "NCPPB_2511_03183"]
remove_ids("CBM63.rbbh.fasta", ids_remove_1, "../data/cazy_rbbh")

#PL26
ids_remove_2 = ["NCPPB_2511_02550"]
remove_ids("PL26.rbbh.fasta", ids_remove_2, "../data/cazy_rbbh")

#GT41
ids_remove_3 = ["NCPPB_402_02760"]
remove_ids("GT41.rbbh.fasta", ids_remove_3, "../data/cazy_rbbh")


