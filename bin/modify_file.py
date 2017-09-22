#!/usr/bin/env python
# pylint: disable=C0103
# This python script modifies single lines in a textfile
# (c) The James Hutton Institute 2016
# Author: Eirini Xemantilotou
# modify_file.py
"""
Python script for modifying the content of a file by opening and then reading it
"""

import sys


infname = sys.argv[1]


def replace_line(file_name, line_num):
    """
    This function takes 2 arguments a file and a line within a file.
    It reads the file, and it manipulates a single line within the file.

    """
    with open(file_name, 'r') as fh:
        lines = fh.readlines()
        lines[line_num] = lines[line_num].strip('\n')+" I" + "\n"
        out = open(file_name, "w")
        out.writelines(lines)
        out.close()

replace_line(infname, 0)
