#!/usr/bin/env python
# -*- coding: utf-8 -*-
#pylint: disable=C0103
#codeml.py
#
# (c) The James Hutton Institute 2017
#Author: Eirini Xemantilotou
# This python script calls the codeml_data and codeml_data_2 functions
# and generate two csvs.
# Each csv file contain information about the null and alternative model's ML,
# the corresponding CAZy family and the tree directory.


import pandas as pd
from scipy import stats
# Import the codeml_data function from the codeml_function python script
from codeml_function import codeml_data
# Import the codeml_data function from the codeml_function_groups python script
from codeml_function_groups import codeml_data_2

groups_fams = ["CBM48", "CBM50", "CE8", "GH1", "GH3", + \
"GH13", "GH19", "GH23", "GH28", "GH33", "GH73", + \
 "GH103", "GH104", "GH105", "GT1", "GT2", "GT4", + \
  "GT9", "GT35", "GT51", "PL1", "PL3", "PL4", "PL9"]
directory = "./families"
for family in os.listdir(directory):
	if family not in groups_fams:
		# Call the function for the alternative model
		data_alt = codeml_data("./families", "alternative")
		# Call the function for the null model
		data_null = codeml_data("./families", "null")
		# Assign columns we want teh final dataframe to include
		cols_to_use = data_alt.columns - data_null.columns
		# Create a new dataframe which is the merge of the
		# null and alternative data frames adding the null column
		new_pd = pd.merge(data_null, data_alt[cols_to_use], left_index=True, right_index=True, how='outer')
		# Calculate the LRT and store values in a new column
		new_pd['LRT'] = (new_pd['alternative'] - new_pd['null']) *2
		# Calculate the p value
		new_pd['p-value'] = new_pd['LRT'].apply(stats.chisqprob, args=(1,))
		# Write out the table in a csv as codeml_results.csv
		new_pd.to_csv('codeml_results.csv', mode='a', header=False)
	else:
		# Call the function for the alternative model for split cazy families
		data_alt_split = codeml_data_2("./families", "alternative")
		# Call the function for the null model for split cazy families
		data_null_split = codeml_data_2("./families", "null")
		# Assign columns we want teh final dataframe to include
		cols_to_use = data_alt_split.columns - data_null_split.columns
		# Create a new dataframe which is the merge of the
		# null and alternative data frames adding the null column
		new_pd_2 = pd.merge(data_null_split, data_alt_split[cols_to_use], left_index=True,
					               		right_index=True, how='outer')
		# Calculate the LRT and store values in a new column
		new_pd_2['LRT'] = (new_pd_2['alternative'] - new_pd_2['null']) *2
		# Calculate the p value
		new_pd_2['p-value'] = new_pd_2['LRT'].apply(stats.chisqprob, args=(1,))
		# Write out the table in a csv as codeml_results.csv
		new_pd_2.to_csv('codeml_results_split.csv', mode='a', header=False)
