#!/usr/bin/env python
# -*- coding: utf-8 -*-
# split_large_CAZy_fams.py
# 
# (c) The James Hutton Institute 2016
# Author: Eirini Xemantilotou
# pylint: disable=C0103


import os
import re
import sys
import shutil
from Bio import SeqIO
from Bio import Phylo

"""
functions to split a dendrogram into *disjoint* subtrees on two criteria:

1. do not allow a subtree to contain paralogues 
(two sequences from the same organism)
2. do not allow a subtree to have a depth greater than some value

If either criterion is violated,
the subtree in question is broken down into progressively smaller subtrees 
until a set of subtrees describing all original leaf nodes 
that do fit those criteria is obtained.
"""

# We split on paralogues (or lack of…) first, 
# then split *those* subtrees by subtree length, to give subtrees 
# that contain no paralogues *and* meet a minimum depth requiremend. 
# To do this, we can apply the recursive functions in turn:

os.makedirs(sys.argv[1])

species_pattern = re.compile("_[0-9]*$")

def has_paralogues(tree):
    """Returns True if the passed tree contains two or more paralogues.
    
    The identifiers for all leaves in the tree are processed with re
    to reduce them to only their organism identifier. The set() of these
    identifiers is nonredundant, so if it is the same length as the
    number of leaf nodes in the tree, there are no repeated organisms,
    and the tree contains no paralogues.
    """
    leaf_names = [re.split(species_pattern, e.name)[0] for e in tree.get_terminals()]
    return len(leaf_names) != len(set(leaf_names))

def is_top_level(clade, tree):
    """Returns True if clade is top level in tree"""
    return len(tree.get_path(clade)) == 1


def split_tree_noparalogues(tree, nplist=None):
    """Return a list of disjoint subtrees with no paralogues"""
    if nplist is None:  # Required to avoid state issues
        nplist = list()
    if not has_paralogues(tree):  # base case
        nplist.append(tree)
        return nplist
    else:  # recursive case
        for clade in [c for c in tree.find_clades()
                      if is_top_level(c, tree)]:
            subtree = Phylo.BaseTree.Tree(clade)
            split_tree_noparalogues(subtree, nplist)
    return nplist

def is_shorter_than(tree, length):
    return max([tree.distance(t) for t in tree.get_terminals()]) < length


def split_tree_depth(tree, branch_length=20, stlist=None):
    """Returns a list of subtrees with max len branch_length."""
    if stlist is None:  # Required to avoid problems with state
        stlist = list()
    if is_shorter_than(tree, branch_length):  # base case
        stlist.append(tree)
        return stlist
    else:  # recursive case
        for clade in [c for c in tree.find_clades() if is_top_level(c, tree)]:
            st = Phylo.BaseTree.Tree(clade)
            split_tree_depth(st, branch_length, stlist)
    return stlist
    
def split_tree(tree, branch_length=20):
    """Split tree into disjoint subtrees with maximum depth and
    containing no paralogues.
    """
    paralogue_splits = split_tree_noparalogues(tree)
    subtrees = []
    for ps in paralogue_splits:
        subtrees.extend(split_tree_depth(ps, branch_length))
    return subtrees


#for treefile in [f for f in os.listdir('.') if
                 #os.path.splitext(f)[-1] == '.new']:


treefile = sys.argv[1] +".new"
file_tree = treefile.split(".")[0]
tree = Phylo.read(treefile, format="newick")
print("Complete tree for {}:".format(treefile))
subtrees = split_tree(tree)
print("Subtrees ({})".format(len(subtrees)))
# We use Phylo to write out the new split trees in newich format
i = 1   
for st in subtrees:
    tdir = str(i)
    Phylo.write(st, treefile+tdir+'.nwk', 'newick')
    i +=1
        
    leaf_names = [e.name for e in st.get_terminals()]
    result_file = file_tree + "_"+tdir+ '.rbbh.fasta'
    with open("group_"+tdir+'.txt', 'w') as output:
        for name in leaf_names:
            output.write(name+ '\n')



directory = "../data"
# Fasta as the MSA file
fasta_file = "../data/cazy_rbbh/"+sys.argv[1]+".rbbh.fasta" # Input fasta file

        # assign them as the wanted files
for text_file in os.listdir("."):
    if text_file.endswith("txt"):
        file = text_file.split(".")[0]
        number = file.split("_")[1]
        wanted_file = "group_"+ number +".txt" # Input interesting sequence IDs, one per line
        # assign the result files with the sugroups 
        result_file = "group_" +number +'.fasta' # Output fasta file   
        # set for all wanted identifiers
        wanted = set()
        with open(wanted_file) as f:
            for line in f:
                line = line.strip()
                if line != "":
                    wanted.add(line)
                        
        # open msa 
        fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
        # open result file
        with open(result_file, "w") as f:
            # if the id is in wanted then write it out from the fasta file to result file named after the group ot belongs to
            for seq in fasta_sequences:
                if seq.id in wanted:
                    SeqIO.write([seq], f, "fasta")


# Move text and fasta files in corresponding repos named as fasta and txt
fasta = sys.argv[1]+"/fasta"
os.makedirs(fasta)
family_dir = sys.argv[1]
for fasta_files in os.listdir("."):
    if fasta_files.endswith(".fasta"):
       # shutil.move(os.path.join(family_dir, fasta_files), fasta)
        shutil.move(fasta_files, fasta)

text = sys.argv[1]+"/txt"
os.makedirs(text)
for txt_files in os.listdir("."):
    if txt_files.endswith(".txt"):
        shutil.move(txt_files, text)

dir_trees = "./trees"
os.makedirs(dir_trees, exist_ok=True)
for new_nwk in os.listdir("."):
    if new_nwk.endswith(".nwk"):
        shutil.move(new_nwk, dir_trees)
