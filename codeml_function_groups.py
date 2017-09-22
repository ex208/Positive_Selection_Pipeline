#!/usr/bin/env python
#*- coding: utf-8 -*-
#
# pylint: disable=C0103
#codeml_function_groups.py
#
# (c) The James Hutton Institute 2017
#Author: Eirini Xemantilotou

"""
This script is a function which allows us to obtain  data for the maximum
likelihood values of each cazy family and tree dir.
The function takes two arguments. The first argument corresponds to the directory
storing the CAZy families we are interested in and the second argument
corresponds to the model we are interested in (alternative or null).
The following function will only work for data for one CAZy family,
several groups fro each CAZy family and several tree dirs
"""

import os
from Bio.Phylo.PAML import codeml
import pandas as pd

def codeml_data_2(directory, model):
    """
    codeml_data function which takes as:
    1.first argument the directory storing the mlc files we want to extrcat
    the maximum likelihood
    2. second argument the model we are extracting the values from. Which can be either
    the null or alternative
    """
    columns = ['CAZy_Family', 'tdirs', model]
    data = pd.DataFrame(columns=columns)
    subdir = os.listdir(directory)
    for dirs in subdir:
        if dirs != ".DS_Store":
            path = os.path.join(directory, dirs)
            for group in os.listdir(path):
                if group.startswith("g"):
                    group = group.split(".")[0]
                    full_path = os.path.join(directory, dirs, group, "positive_selection", model)
                    for tdir in os.listdir(full_path):
                        if tdir != "tdir":
                            mlc = os.path.join(full_path, tdir, "control.mlc")
                            results = codeml.read(mlc)
                            lnL = (results.get("NSsites").get(2).get('lnL'))
                            data = data.append(pd.DataFrame({'CAZy_Family': dirs,
                                                             'sub_family': group,
                                                             'tdirs': tdir, model: lnL},
                                                            index=[0]), ignore_index=True)
    return data
    