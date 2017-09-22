#!/usr/bin/env python3
# pylint: disable=C0103
# ete.py
# # (c) The James Hutton Institute 2016
# Author: Eirini Xemantilotou
#

"""
 This python script has as purpose to iterate over all unique nodes of a tree
 constructed by RAxML and add #1 so it can be used by
 codeML for positive selection analysis
"""
# IMPORTANT: Need to make sure that its tree file with marked branch
# is being stored in subdirectories with the appropriate control files

import sys
from ete3 import Tree

# The command line argument must be the RAxML tree
# This is the raxml tree generated previously
filename = sys.argv[1]

# Loads a tree structure from a newick string.
# The returned variable t is the root node for the tree.
t = Tree(filename)
u = Tree(filename, format=1)



# Function to modify a tree node label,
# print the tree code, and restore the node label
def modify_and_print_node(t, n):
    """
	Appends `#1` to the passed node's label, prints the tree
    in Newick format, then restores the original node label.

    - t  the Tree we are modifying
    - n  the TreeNode to modify
    """
    label = n.name
    n.name += '#1'
    outstr = t.write(format=1)
    if '#' in outstr:
        print(t.write(format=1))
        n.name = label


for node in t.traverse('postorder'):
    modify_and_print_node(t, node)

