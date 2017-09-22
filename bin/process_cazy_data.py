#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# process_cazy_data.py
# (c) The James Hutton Institute 2016
# Author: Eirini Xemantilotou
# pylint: disable=C0103

"""
Python script to process data dowloaded by Uniprot by accessing CAZy through cross refernece
"""
import os
from io import StringIO
import shutil
import pandas as pd

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import bioservices
# Using BioServices, rather than a manual online search to get data
from bioservices import UniProt

# Using UniProt BioServices to get data, rather than manual search
u = UniProt()
query_result = u.search('database:(type:cazy) AND taxonomy:"Dickeya [204037]"',
                  frmt='tab',
                  columns="id, entry name, genes(OLN), protein names, database(CAZy), sequence")
result = pd.read_csv(StringIO(query_result), sep="\t")

# Get all proteins without a locus tag
nolocus = result[pd.isnull(result["Gene names  (ordered locus )"])]
nolocus = nolocus.reset_index()

# Get all proteins with a locus tag
locus = result[pd.notnull(result["Gene names  (ordered locus )"])]
locus = locus.reset_index()

# Write sequences with no locus tag to FASTA file
seqlist = list()
for idx, entry in result[pd.isnull(result["Gene names  (ordered locus )"])].iterrows():
    seqlist.append(SeqRecord(id=entry["Entry name"],
                             description=entry["Protein names"],
                             seq=Seq(entry["Sequence"])))
SeqIO.write(seqlist, "all_uniprot_nolocus_results.fasta", "fasta")
shutil.move("./all_uniprot_nolocus_results.fasta", "./data")

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

nolocus_split = split_cazy_families(result[pd.isnull(result["Gene names  (ordered locus )"])])
locus_split = split_cazy_families(result[pd.notnull(result["Gene names  (ordered locus )"])])

# Write nolocus seq for all families out in different files names $FAMILY_nolocus.fasta
outdir = 'nolocus_sequences'
os.makedirs(outdir, exist_ok=True)
for family in nolocus_split["Cross-reference (CAZy)"].unique():
    nolocus_family_seq = list()
    for idx, entry in nolocus_split[nolocus_split["Cross-reference (CAZy)"] == family].iterrows():
        nolocus_family_seq.append(SeqRecord(id=entry["Entry name"],
                                            description=entry["Protein names"],
                                            seq=Seq(entry["Sequence"])))
    filename = os.path.join(outdir, "{0}_nolocus.fasta".format(family))
    SeqIO.write(nolocus_family_seq, filename, "fasta")

# Write nolocus seq for all families out in different files names $FAMILY_nolocus.fasta
outdir = './locus_sequences'
os.makedirs(outdir, exist_ok=True)
for family in locus_split["Cross-reference (CAZy)"].unique():
    locus_family_seq = list()
    for idx, entry in locus_split[locus_split["Cross-reference (CAZy)"] == family].iterrows():
        locus_family_seq.append(SeqRecord(id=entry["Entry name"],
                                          description=entry["Protein names"],
                                          seq=Seq(entry["Sequence"])))
    filename = os.path.join(outdir, "{0}_nolocus.fasta".format(family))
    SeqIO.write(locus_family_seq, filename, "fasta")

shutil.move("./locus_sequences/", "../data")
shutil.move("./nolocus_sequences/", "../data")
#os.rmdir("./locus_sequences/")
#os.rmdir("./nolocus_sequences/")

# Getting the entries with two or mores values for CAZy and parsing into string
countdf = locus_split.groupby("Gene names  (ordered locus )").count()
result = locus_split.loc[locus_split["Gene names  (ordered locus )"].isin(countdf[countdf["Cross-reference (CAZy)"] > 1].index)]
resultlist = result.groupby("Gene names  (ordered locus )")["Cross-reference (CAZy)"].apply(list).tolist()

#  We turn into a list the dataframe for the nolocus tag
# so we can get the Uniprpot accession numbers of those seqs
no_locus_list = nolocus['Entry'].tolist()

# Write in a text file the uniprot accessions for sequences with no locus tags
with open('uniprot_no_locus.txt', 'w') as file_handler:
    for i in no_locus_list:
        file_handler.write("{}\n".format(i))
shutil.move("./uniprot_no_locus.txt", "../data")
#os.remove("./uniprot_no_locus.txt")

#  We turn into the list the dataframe for the locus tag
# so we can get the the locus tags for the RBBH analysis
locus_tags_list = locus_split["Gene names  (ordered locus )"].tolist()

# Write in a text file the  for sequences with  locus tags
with open('locus_Dickeya_CAZy.txt', 'w') as file_handler:
    for i in locus_tags_list:
        file_handler.write("{}\n".format(i))
shutil.move("./locus_Dickeya_CAZy.txt", "../data")
#os.remove("./locus_Dickeya_CAZy.txt")

# Run RBBH_CAZY python script which executes the extract_rbbh.py python script for all CAZy families
os.system('python3 ./RBBH_cazy.py')

# RMerge txt files with locus tags were not downloaded by accessing Uniprot
os.system('python3 ./merge_locus.py')

# Run RBBH_CAZY python script which executes the extract_rbbh.py python script for those CAZy families
# were not retrieved by uniprot
os.system('python3 ./RBBH_not_det_uniprot.py')

# Make directories based on the CAZy family the enzymes with known locus tags belong to
# Copying and moving the RBBH for the entries with more than one value
# making sure they exist in both families
path = "."
path_data = "../data"
cazy_locus_tags = os.path.join(path_data, "cazy_locus_tags")
os.makedirs(cazy_locus_tags)

mydict = {}
for x in range(len(locus_split)):
    currentid = locus_split.iloc[x, 3]
    currentvalue = locus_split.iloc[x, 5]
    mydict.setdefault(currentid, [])
    mydict[currentid].append(currentvalue)
for key, value in mydict.items():
    for i in os.listdir(path):
        if i.endswith("rbbh.fasta"):
            infstem = i.split('.', 1)[0]
            if key == infstem:
                if len(value) > 1:
                    value1 = value[0]
                    value2 = value[1]
                    string_value1 = ''.join(value1)
                    string_value2 = ''.join(value2)
                    dest_dir1 = os.path.join(cazy_locus_tags, string_value1)
                    dest_dir2 = os.path.join(cazy_locus_tags, string_value2)
                    dest_file1 = os.path.join(dest_dir1, i)
                    dest_file2 = os.path.join(dest_dir2, i)
                    if not os.path.exists(dest_dir1):
                        os.makedirs(dest_dir1)
                    if not os.path.exists(dest_dir2):
                        os.makedirs(dest_dir2)
                    shutil.copy(os.path.join(path, i), dest_dir1)
                    shutil.copy(os.path.join(path, i), dest_dir2)
                if len(value) == 1:
                    string_value = ''.join(value)
                    dest_dir = os.path.join(cazy_locus_tags, string_value)
                    dest_file = os.path.join(dest_dir, i)
                    if not os.path.exists(dest_file):
                        os.makedirs(dest_dir, exist_ok=True)
                    shutil.copy(os.path.join(path, i), dest_dir)
                
  
# sorting in dirs the rbbh for the seqs were not recognised by Uniprot 
# with cross reference CAZy 

for locus_txt in os.listdir(path_data):
    if locus_txt.startswith("PL") or locus_txt.startswith("CE"):
        with open(os.path.join(path_data,locus_txt)) as infile:
                    lines = infile.read()
                    locus_id = locus_txt.split("_")[0]
                    for rbbh in os.listdir("."):
                        if rbbh.endswith("rbbh.fasta"):
                            rbbh_id = rbbh.split(".")[0]
                            if rbbh_id in lines:
                                dest_dir3 = os.path.join(cazy_locus_tags,locus_id )
                                if not os.path.exists(dest_dir3):
                                    os.makedirs(dest_dir3)
                                shutil.copy(os.path.join(path, rbbh), dest_dir3)

for rbbh in os.listdir(path):
    if rbbh.endswith("rbbh.fasta") or rbbh.endswith(".log"):
        os.remove(rbbh)

# Make a directory to store all final rbbh files per CAZy family
final_rbbh = '../data/cazy_rbbh'
os.makedirs(final_rbbh, exist_ok=True)

cazy_locus = "../data/cazy_locus_tags/"
# We loop over the CAZy directories.
for family in os.listdir(cazy_locus):
# We create an empty list, we open the first file in the CAZy directory,
# we open a new fasta file to write the unique
# sequences, and we get all records and write them to the new file.
    seqlist = []
    cazy_fam = os.path.join(cazy_locus, family)
    with open(cazy_fam + "/" + os.listdir(cazy_fam)[0], "r") as f2:
        with open(family +".rbbh.fasta", "w") as output_handle:
            for records in SeqIO.parse(f2, "fasta"):
                output_handle.write(">" + records.id + "\n")
                output_handle.write(str(records.seq)+ "\n")
# We open the new fasta file in read mode and we add the seq ids
    with open(family +'.rbbh.fasta', "r") as output_handle:
        for records in SeqIO.parse(output_handle, "fasta"):
            seqlist.append(records.id)

# We  open the new fasta file in append mode and we loop over
# the fasta files within the directory, we add the sequences to the new fasta file if the record.id
# is not in the seqlist and once write it out then add the id
# to the seqlist with append
    with open(family +".rbbh.fasta", "a") as output_handle:
        for file in os.listdir(os.path.join(cazy_locus, family)):
            with open(cazy_fam + "/" + file, "r") as f1:
                for record in SeqIO.parse(f1, "fasta"):
                    if record.id not in seqlist:
                        output_handle.write(">" + record.id + "\n")
                        output_handle.write(str(record.seq)+ "\n")
                        seqlist.append(record.id)

for files in os.listdir(path):
    if files.endswith(".rbbh.fasta"):
        shutil.move(files, final_rbbh)
