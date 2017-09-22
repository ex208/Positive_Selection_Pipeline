#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#split_CAZy_families.py
# pylint: disable=C0103
# (c) The James Hutton Institute 2017
#Author: Eirini Xemantilotou
# This script has as purpose to split large CAZy families into subgroups to repeat
# CodeML analysis
"""
Script to split large CAZy fams into subgrous for
further analysis for positive selection
"""
import os
import sys
import sqlite3
import shutil
from Bio import SeqIO


#Accessing dickyea database
database = "./dickeya.db"

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
        locus_tags = []
    with open("../data/cazy_rbbh/"+ sys.argv[1] +".rbbh.fasta", "r") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
        tags = [r.id for r in records]
        for record in records:
            locus_tags.append(record.id)
        # Call the database to obtain RBBH (egs) for each locus tag
        for locus_tag in locus_tags:
            #print("Locus tag", locus_tag)
            cursor.execute(sql, (locus_tag,))
            results = cursor.fetchall()
            for i in results:
                if i[2] not in tags:
                    print("New sequence not in family!", i[2])
            eg = {i[2] for i in results}
            #add the locus tag in the list too
            eg.add(locus_tag)
            #print(eg)
            # If the eg is not already in the list then append it
            if eg not in egs:
                print("New eg", locus_tag)
                #print(eg)
                egs.append(eg)
                # Create groups based on how many egs we have
                for e in range(len(egs)):
                    groups = "group" + "_" + str(e)
                    groups_list = []
                    outfile = open("../families/"+sys.argv[1]+"/"+groups+'.txt', 'w')
                    for data in egs[e]:
                        groups_list.append(data)
                    # write out the group list in the group text files seperated
                    # in new lines for easier processing later on
                    outfile.write('\n'.join(groups_list))
                    #print(groups)
                    #print(groups_list)
                    

# Call the function
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
