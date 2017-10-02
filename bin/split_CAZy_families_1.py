#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#split_CAZy_families_1.py
# pylint: disable=C0103
# (c) The James Hutton Institute 2017
#Author: Eirini Xemantilotou
# This script has as purpose to split large CAZy families into subgroups to repeat
# CodeML analysis
"""
Script to split large CAZy fams into subgrous for
further analysis for positive selection
"""

import logging
import os
import sys
import re
import sys
import sqlite3
from io import StringIO
import shutil
import pandas as pd

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from bioservices import UniProt


# Set up logger
logging.basicConfig(filename='info.log', level=logging.DEBUG)
logger = logging.getLogger("split_CAZy_families_1.py")

# create handlers
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.INFO)

# create a logging format
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)

# add the handlers to the logger
logger.addHandler(handler)

# Create directories for each family teh analysis needs to be repeated 
fams = "../families/"+sys.argv[1]
os.makedirs(fams, exist_ok=True)

# Using UniProt BioServices to get data, rather than manual search
u = UniProt()
query_result = u.search('database:(type:cazy) AND taxonomy:dickeya',
                        frmt='tab',
                        columns="entry name, id, genes(OLN), protein names, database(CAZy), sequence")
result = pd.read_csv(StringIO(query_result), sep="\t")

# Get all proteins with a locus tag
locus = result[pd.notnull(result["Gene names  (ordered locus )"])]
locus = locus.reset_index()


def split_cazy_families(df):
    """Returns a DataFrame with one CAZy family per row.
    UniProt sequences will be duplicated, where they are members of more than one family
    """
    df = df.reset_index()
    for idx, entry in df.iterrows():
        cazy_families = [e for e in entry["Cross-reference (CAZy)"].strip().split(';') if e]
        df.set_value(idx, "Cross-reference (CAZy)", cazy_families[0])
        if len(cazy_families) > 1:
            for family in cazy_families[1:]:
                new_entry = entry.copy()
                new_entry["Cross-reference (CAZy)"] = family
                df = df.append(new_entry, ignore_index=True)
    return df
 
 # Getting all entries with loocus tags per CAZy family
locus_split = split_cazy_families(result[pd.notnull(result["Gene names  (ordered locus )"])])
# Turning the locus tags into a list
locus_tags_list = locus_split["Gene names  (ordered locus )"].tolist()

# Directory to store fasta files per CAZY family with locus tags and their seqs
locus_dir_fasta = '../data/locus/fasta'
os.makedirs(locus_dir_fasta, exist_ok=True)

for caz_family in locus_split["Cross-reference (CAZy)"].unique():
    locus_family_seq = list()
    for index, column in locus_split[locus_split["Cross-reference (CAZy)"] == caz_family].iterrows():
        locus_family_seq.append(SeqRecord(id=column["Gene names  (ordered locus )"],
                                          seq=Seq(column["Sequence"])))
    filename = os.path.join(locus_dir_fasta, "{0}_locus.fasta".format(caz_family))
    SeqIO.write(locus_family_seq, filename, "fasta")

# Directory with text files per family with teh names of the tags
locus_dir_txt = '../data/locus/txt'
os.makedirs(locus_dir_txt, exist_ok=True)

# No data on Uniprot for CAZy family CE8 so we need to provide the text file ourselves
shutil.copy("../data/CE8_locus.txt", locus_dir_txt)

# List of the locus tags per family and writ eout text files
locus_tags = []
for fasta in os.listdir(locus_dir_fasta):
    if fasta == sys.argv[1]+"_locus.fasta":
        fasta = os.path.join(locus_dir_fasta, fasta)
        with open(fasta) as handle_file:
            locus_ids = list(SeqIO.parse(handle_file, "fasta"))
            loc_tags = [loc.id for loc in locus_ids]
            locus_file_dir = os.path.join(locus_dir_txt, sys.argv[1]+"_locus.txt")
            with open(locus_file_dir, 'w') as file_handler:
                for locus_name in loc_tags:
                    locus_tags.append(locus_name)
                file_handler.write('\n'.join(loc_tags))


# Accessing dickyea database
database = "../data/dickeya.db"

# Function to split families into subgroups and then write out
# text files with the locus tags identifiers of each subgroup
# stored under the CAZy's family they belong to

def split_families():
    """
    Function to split big CAZy families into
    subgroups by connecting to the local database
    and grapping the rbbbh column
    Argument being CAZy family name
    """
    # connect to the sqlite database
    with sqlite3.connect(database) as conn:
        cursor = conn.cursor()
        # get the rbbh column
        sql = 'SELECT * FROM rbbh WHERE locus_tag_1=?'
        # List to store the locus tags seen in the RBBH
        egs = []
        # List to store the locus tags ids from the rbbh file
        locus_names = []
    locus_file_dir = os.path.join("../data/locus/txt", sys.argv[1]+"_locus.txt")
    with open(locus_file_dir) as handle:
        records = handle.readlines()
        tags = [r.strip() for r in records]
        for record in records:
            locus_names.append(record.strip())
            for locus_tag in locus_names:
                cursor.execute(sql, (locus_tag,))
                results = cursor.fetchall()
                for i in results:
                    if i[2] not in tags:
                        print("New sequence not in family!", i[2])
                eg = {i[2] for i in results}
                eg.add(locus_tag)
                if eg not in egs:
                    egs.append(eg)
                    for e in range(len(egs)):
                        groups = "group" + "_" + str(e)
                        groups_list = []
                        outfile = open("../families/"+sys.argv[1]+"/"+groups+'.txt', 'w')
                        #outfile = open(file+'.txt', 'w')
                        for data in egs[e]:
                            groups_list.append(data)
                        outfile.write('\n'.join(groups_list))

split_families()

# Write out those sequences based on the identifiers and
# extracting them from the corresponding rbbh the command
# line argument is the CAZy family we want to split into subgroups
# Assign the directory with the families
directory = "../families/"+sys.argv[1]
# Fasta as the input file to get the seqs
fasta_file = "../data/cazy_rbbh/"+sys.argv[1]+'.rbbh.fasta' # Input fasta file
# filter the text file
for text_files in os.listdir(directory):
    if text_files.endswith(".txt"):
        text_files = text_files.split(".")[0]
        # assign them as the wanted files
        # Input interesting sequence IDs, one per line
        wanted_file = "../families/"+sys.argv[1]+"/"+text_files+'.txt'
        # assign the result files with the sugroups
        # Output fasta file
        result_file = "../families/"+sys.argv[1]+"/"+text_files +'.fasta'
        # set for all wanted identifiers
        wanted = set()
        with open(wanted_file) as f:
            for line in f:
                line = line.strip()
                if line != "":
                    wanted.add(line)
        # open rbbh
        fasta_sequences = SeqIO.parse(open(fasta_file), 'fasta')
        # open result file
        with open(result_file, "w") as f:
            # if the id is in wanted then write it out from the fasta file
            # to result file named after the group it belongs to
            for seq in fasta_sequences:
                if seq.id in wanted:
                    SeqIO.write([seq], f, "fasta")

# Move text and fasta files in corresponding repos named as fasta and txt
fasta = "../families/"+sys.argv[1]+"/fasta"
os.makedirs(fasta)
family_dir = "../families/"+sys.argv[1]
for fasta_files in os.listdir(family_dir):
    if fasta_files.endswith(".fasta"):
        shutil.move(os.path.join(family_dir, fasta_files), fasta)

text = "../families/"+sys.argv[1]+"/txt"
os.makedirs(text)
for txt_files in os.listdir(family_dir):
    if txt_files.endswith(".txt"):
        shutil.move(os.path.join(family_dir, txt_files), text)

# Testing for duplicated seqs in the text files
sets = []
for text_files in os.listdir(text):
    if text_files.endswith(".txt"):
        text_files = os.path.join(text, text_files)
        with open(text_files) as f:
            lines = [line.strip() for line in f.readlines()]
            sets.append(set(lines))

for h in range(len(sets)):
    for j in range(h+1, len(sets)):
        duplicates = sets[h] & sets[j]
        if len(duplicates) > 0:
            print("%s group and %s group both contain %s" % (h, j, ', '.join(duplicates)))
            #logger.warning("%s group and %s group both contain %s" % (h, j, ', '.join(duplicates)
        else:
            print("No duplicates detected")
            logger.info("No duplicates detected")

# Testing for paralogs
for text_files in os.listdir(text):
    if text_files.endswith(".txt"):
        text_files = os.path.join(text, text_files)
        with open(text_files) as h:
            file = [line.strip() for line in h.readlines()]
            split_species = {re.split("_[0-9]*$", x)[0] for x in file}
            if len(file) != len(split_species):
                print("Paralogs found!")
                logger.warning("Paralogs detected!")
            else:
                print("No paralogs!")
                logger.info("No paralogs dtected!")

# Texting if teh fasta files has teh sequences from the text files
# Look at the fasta files 
for dirs in os.listdir(directory):
    if dirs != ".DS_Store":
        fasta = os.path.join(directory, dirs, "fasta")
        for fasta_files in os.listdir(fasta):
            if fasta_files.endswith(".fasta"):
                ins1 = fasta_files.split(".")[0]
                locus_tags = []
                fasta_files = os.path.join(fasta, fasta_files)
                with open(fasta_files, "r") as handle:
                    records = list(SeqIO.parse(handle, "fasta"))
                    tags = [r.id for r in records]
                    for record in records:
                        locus_tags.append(record.id)
            text = os.path.join(directory, dirs, "txt")
            for text_files in os.listdir(text):
                if text_files.endswith(".txt"):
                    ins2 = text_files.split(".")[0]
                    loc = []
                    text_files = os.path.join(text, text_files)
                    if ins1 == ins2:
                        for lines in open(text_files):
                            loc.append(lines.strip())
                        #print(text_files, len(loc), fasta_files, len(locus_tags) )
                        if len(loc) != len(locus_tags):
                            print(text_files) 
                            logger.warning("Text file and fasta file do not correspond to each other!")
                            logger.debug("Text file %s has %s locus and the %s has %s sequences" % (text_files, len(loc), fasta_files, len(locus_tags)))
logger.info("Complete")
