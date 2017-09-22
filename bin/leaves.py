#!/usr/bin/env python
#
#*- coding: utf-8 -*-
#leaves.py
#
# (c) The James Hutton Institute 2017
#Author: Eirini Xemantilotou

"""
This python script includes three function named as:
get_leaves, external_nodes, cazy_leaves which assist with
getting labelled nodes, removing labelled leaves and creating
a dictionary of cazy family and the leaves under positive selection.
"""
import os
from collections import defaultdict
from ete3 import Tree
from ete3 import TreeStyle

#Function to get all descendents leaves for the nodes apppearing to be under positive selection
# We returt the leaves as a tuple. We are also applying a filter:
# '_'.join(leaf.name.split('_')[:2]in order to get the names for the species isolates
# and exclude the name of individual enzymes
# 1. load tree into ete Tree object
# 2. identify the #1 node
# 3. get the leaves beneath that node
# 4. return the leaves as an ordered tuple

def get_leaves(tree):
    """The get_leaves function returns
    the leaves of a given labelled tree stored in a tuple
    """
    tree = Tree(tree, format=1)
    midpoint = tree.get_midpoint_outgroup()
    tree.set_outgroup(midpoint)
    tree_style = TreeStyle()
    tree_style.show_leaf_name = False
    nodes = tree.search_nodes(name="#1")[0]
    leaves = sorted(['_'.join(leaf.name.split('_')[:2]) for leaf in nodes.iter_leaves()])
    leaves = tuple(leaves)
    return leaves

# Function to detect cazy families and tree directories marked at external nodes
def external_nodes(directory, dataframe):
    """
    external_nodes: returns a list of cazy families and tree directories which are
    marked at the external nodes.
    1. Directory: is the root directory which stores all CAZy families
    2. Dataframe: is the loaded csv with filtered qvals
    """
    # Cretae empty list for storing cazy families and tree directories for
    # labelled leaves
    cazy_fams = []
    tree_dirs = []
    # Get all CAZy families, subgroups and tree directoris
    resdict = defaultdict(list)
    for subd in [d for d in os.listdir(directory) if
                 os.path.isdir(os.path.join(directory, d))]:
        for enzymes in os.listdir(os.path.join(directory, subd)):
            if enzymes.startswith("g"):
                enzymes = enzymes.split(".")[0]
                subpath = os.path.join(directory, subd, enzymes, "positive_selection",
                                       "alternative")
                tdirs = dataframe[dataframe['CAZy_Family'] == subd+"_"+enzymes]['tdir']
                for tdir in tdirs:
                    try:
                        treefile = os.path.join(subpath, tdir, "tree.nw")
                        leaves = get_leaves(treefile)
                        resdict[subd].append(leaves)
    # Identify the directories with labelled leaves due to Index error
                    except IndexError:
    # Append to the empty lists the cazy fams and tree directories
                        cazy_fams.append(subd + "_" + enzymes)
                        tree_dirs.append(tdir)
     # In case Cazy families are not any further split into subgroups then have to state and
     # include those families too
            else:
                subpath = os.path.join(directory, subd, "positive_selection", "alternative")
                tdirs = dataframe[dataframe['CAZy_Family'] == subd]['tdir']
                for tdir in tdirs:
                    try:
                        treefile = os.path.join(subpath, tdir, "tree.nw")
                        leaves = get_leaves(treefile)
                        resdict[subd].append(leaves)
                    except IndexError:
                        cazy_fams.append(subd)
                        tree_dirs.append(tdir)
    return cazy_fams, tree_dirs

# Function to detect cazy families and tree directories marked at external nodes
def cazy_leaves(directory, dataframe):
    """
    cazy_leaves: returns a dictionary for the CAZy family members with the given dataframe
    which is the second argument and the leaves (which we can get by using
    the function get_leaves) cazy_fam:leaves. The first argument is the root directory in
    which we have stored all CAZy families with their results.
    """
    # Create empty dictionary to add cazy families and leaves
    resdict = defaultdict(list)
    # Get CAZy families, subgroups, tree dirs and tree files
    for subd in [d for d in os.listdir(directory) if
                 os.path.isdir(os.path.join(directory, d))]:
        for groups in os.listdir(os.path.join(directory, subd)):
            if groups.startswith("g"):
                groups = groups.split(".")[0]
                cazy_group = subd +"_" + groups
                subpath = os.path.join(directory, subd, groups, "positive_selection", "alternative")
                tdirs = dataframe[dataframe['CAZy_Family'] == subd+"_"+groups]['tdir']
                for tdir in tdirs:
                    treefile = os.path.join(subpath, tdir, "tree.nw")
                    # Get the labelled leaves using the get leaves function
                    leaves = get_leaves(treefile)
                    # Append cazy family and leaves to the dictionary
                    resdict[cazy_group].append(leaves)
# In case Cazy families are not any further split into subgroups then have to state and
# include those families too
            else:
                subpath = os.path.join(directory, subd, "positive_selection", "alternative")
                tdirs = dataframe[dataframe['CAZy_Family'] == subd]['tdir']
                for tdir in tdirs:
                    treefile = os.path.join(subpath, tdir, "tree.nw")
                    leaves = get_leaves(treefile)
                    resdict[subd].append(leaves)
    return resdict
