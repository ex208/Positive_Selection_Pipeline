#!/usr/bin/env python
# pylint: disable=C0103
# drzetree.py
# # (c) The James Hutton Institute 2016
# Author: Eirini Xemantilotou

"""
This script has as purpose to split tree files and annotate them by applying
drzetree.py "Directorize Tree", multitreefile is split off
into directories and a control file follows it with the tree name included.
"""
import sys
import os

# we use argments inside of hardcoding the the filenames in a script.
argquan = len(sys.argv)
if argquan != 3:
    print("This script requires two arguments: the name of the text file" +
          "and the template control file, whose MSA," +
          "tree and output will be replaced")
    sys.exit(2)

# functions are usually put on the top of a script
# (so all parts of the sript can use them if they want).
def fill_template(controlfile, dataset):
    """Replace $SEQFILE, $TREEFILE and $OUTFILE with data from data"""
    controlfile = controlfile.replace("$SEQFILE", dataset[0])
    controlfile = controlfile.replace("$TREEFILE", dataset[1])
    controlfile = controlfile.replace("$OUTFILE", dataset[2])
    return controlfile

# open combined tree file first
with open(sys.argv[1]) as f:
    fl = f.read()
# we are not interested in newlines because
# we are splitting on only semicolons.
flnn = fl.replace("\n", "")
treelst = flnn.split(";")
# get rid of final element, as ; marks the end and not the beginning of a tree.
del treelst[-1]
sz = len(treelst)

# We have a list called treelst with our individual trees
# though minus the semicolon.

data = ["backtrans_nostops.rphylip", "tree.nw", "control.mlc"]
data[1] = "tree.nw"


with open(sys.argv[2], 'r') as my_template_file:
    controlfile = my_template_file.read()


for i in range(sz):
    tdir = "tdir"+str(i+1)

    os.mkdir(tdir)
    with open(tdir+"/tree.nw", "w") as opf:
        opf.write(treelst[i]+";")
    with open(tdir+"/"+sys.argv[2], 'w') as f:
        f.write(fill_template(controlfile, data))
