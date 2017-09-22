#!/usr/bin/env python
# -*- coding: utf-8 -*-
#pylint: disable=C0103
#automation_split_CAZy_families.py
#
# (c) The James Hutton Institute 2017
#Author: Eirini Xemantilotou

"""
This script has as purpose to automate the
split_CAZy_families.py script
"""

import subprocess

# Arguments being the CAZy families we want to split into subgroups
arguments = ["CBM48", "CBM50", "CE8", "GH1", "GH3","GH13", "GH19", "GH23", "GH28", "GH33", "GH73", "GH103", "GH104", "GH105", "GT1", "GT2", "GT4","GT9", "GT35", "GT51", "PL1", "PL3", "PL4", "PL9"]
for argument in arguments:
    cmd = "python3 ./split_CAZy_families.py %s" %(argument)
    subprocess.call(cmd, shell=True)
