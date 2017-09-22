#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# rename_dirs.py
# (c) The James Hutton Institute 2016
# Author: Eirini Xemantilotou
# pylint: disable=C0103
"""
function to rename the tree dirs
"""
import os

def rename_dirs(model, directory, fams):
    """
    rename_dirs function renames the treedirs so
    we can run the automation scrip without replacing
    the treedirs as they store info about mlc output from cluster
    1: is null or alternative model
    2: root directory CAZy families are stored in
    3: CAZy families we do not want to change
    the tree dir name(the ones split in subgroups)
    """
    for families in os.listdir(directory):
        if families not in fams:
            data_dir = os.path.join(directory, families, "positive_selection", model)
            for treedirs in os.listdir(data_dir):
                trees = os.path.join(data_dir, treedirs)
                os.rename(trees, trees + "_mlc")

# list of cazy families we dont want to apply the function on
group_fams = ["CBM48", "CBM50", "CE8", "GH1", "GH3", "GH13", "GH19",
              "GH23", "GH28", "GH33", "GH73", "GH103", "GH104",
              "GH105", "GT1", "GT2", "GT4", "GT9", "GT35", "GT51",
              "PL1", "PL3", "PL4", "PL9"]

# call teh function for teh null and alternative model
rename_dirs("null", "../families", group_fams)
rename_dirs("alternative", "../families", group_fams)
