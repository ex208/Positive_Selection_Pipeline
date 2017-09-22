#!/usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=C0103
#automation_repeat_split_families.py
# (c) The James Hutton Institute 2017
#Author: Eirini Xemantilotou

'''
This script has as purpose to run the automation
shell script for CAZy families with big phylohgenetic trees
after splitting the families using the ./bin/split_CAZy_families.py s
script through the automation_split_CAZy_families.py
'''


import os
import subprocess
groups_fams = ["CBM48", "CBM50", "CE8", "GH1", "GH3", "GH13", "GH19", "GH23",
               "GH28", "GH33", "GH73", "GH103", "GH104", "GH105", "GT1", "GT2",
               "GT4", "GT9", "GT35", "GT51", "PL1", "PL3", "PL4", "PL9"]
directory = "./families"
for subd in os.listdir(directory):
    if subd in groups_fams:
        full_path = directory +"/" + subd
        final_path = os.path.join(full_path, "fasta")
        for fasta in os.listdir(final_path):
            argument_1 = subd
            argument_2 = fasta.split('.', 1)[0]
            cmd = "bash ps_automation_repeat_split_families.sh %s %s" %(argument_1, argument_2)
            subprocess.call(cmd, shell=True)
