#!/usr/bin/env python
#
#*- coding: utf-8 -*-
#ranking.py
# pylint: disable=C0103
# (c) The James Hutton Institute 2017
#Author: Eirini Xemantilotou


"""
This python script uses funcions within the **leaves.py** python script
in order to generate heatmap with internal labelled nodes which
are shown to be positive selected based on a pre fixed cut off on
q vals
"""

import subprocess
import pandas as pd
from leaves import external_nodes, cazy_leaves
import seaborn as sns

# Accessing the directory where all results are stored (to get CAZy families etc)
cazy_results = "../PhD_Project/PS_analysis/all_results"
# Assign result file with qvals included for all codeml results which will be loaded as dataframe
file_qvals = "./codeml_bind_data_qvals.csv"
qvals_dataframe = pd.DataFrame.from_csv(file_qvals, index_col=None)

# Call the function external nodes from leaves
remove_leaves = external_nodes(cazy_results, qvals_dataframe)

# Create a dataframe with the cazy families and tree directories
# of the labelled external nodes/leaves
nodes = [('CAZY_FAM', remove_leaves[0]),
         ('tree_dir', remove_leaves[1])]
external_nodes = pd.DataFrame.from_items(nodes)

# Add an identifier being the cazy family combined with
# the tree directory to both results and nodes dataframes
external_nodes["identifier"] = external_nodes["CAZY_FAM"] +"_"+ external_nodes["tree_dir"]
qvals_dataframe["identifier"] = qvals_dataframe['CAZy_Family'] +"_"+ qvals_dataframe["tdir"]
no_leaves_df = qvals_dataframe[~qvals_dataframe['identifier'].isin(external_nodes["identifier"])]
no_leaves_df.to_csv('no_external_nodes_ps.csv', mode='a', header=True)

# run r code to filter cazy fams based on qvals being (log(qvalue) < -20)
subprocess.call("Rscript filter_data.R ", shell=True)

# Assign the csv file which is filtered based on q values using the filter_data r code
final_ps = ("./final_ps.csv")
final_df = pd.DataFrame.from_csv(final_ps)

# Call the function cazy_leaves from leaves script
resdict = cazy_leaves(cazy_results, final_df)

# Get CAZy families
families = resdict.keys()

# get leaves
def flatten(l):
    """
	function to flatten leaves
	"""
    return [i for sublist in l for i in sublist]
leaves = set(flatten(resdict.values()))

# Create an empty dataframe with columns being the levaes and index the cazy famiilies.
# Fill it with zeros.
df = pd.DataFrame(0, columns=leaves, index=families)

#Loop over tge key:value in the dictionary to populate the data frame
for k, v in resdict.items():
    df.loc[k][v] += 1

# get the leaves as groups
for i in range(len(leaves)):
    key_alias = "group" + "_" + str(i)
# Change column names to group_1,2,3 etc for better readability
df.columns = ['group%s' % i for i in range(len(df.columns))]
df.columns = ['group%s' % i for i, c in enumerate(df.columns)]
df['Total'] = df.sum(axis=1)
df.loc['Total'] = df.sum()

# Save positive selected families as a csv file
df.to_csv('positeve_selected_qvals.csv')

# generate a heatmap
ax = sns.heatmap(df, annot=True, fmt="d", linewidths=.10)

# Save heatmpap as a figure and then print out as png
fig = ax.get_figure()
fig.savefig("heatmap_short_colnames.png")
