#!/usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=C0103
#automation.py
#
# (c) The James Hutton Institute 2016
#Author: Eirini Xemantilotou

"""
This script has as purpose to run the automation
shell script for all different CAZy families
"""

import os
import subprocess

groups_fams = ["CBM48", "CBM50", "CE8", "GH1", "GH3", "GH13", "GH19",
               "GH23", "GH28", "GH33", "GH73", "GH103", "GH104", "GH105",
               "GT1", "GT2", "GT4", "GT9", "GT35", "GT51", "PL1", "PL3", "PL4", "PL9"]

directory = "./families"
for families in os.listdir(directory):
    if families not in groups_fams:
        argument = families.split('.', 1)[0]
        cmd = "bash ps_automation.sh %s" %(argument)
        subprocess.call(cmd, shell=True)
