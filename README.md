
# README

## Positive selection Detection 

For packages, software, libraries required for the following code:
**pip install -r requirements.txt**

Note: In order to execute  the following the **dickeya.db**, **dickeya_cds_aa.fasta** and **dickeya_cds_nt.fasta** are required

#### 1. CAZy sequences for Dickeya
The first step is to obtain the input sequence data. Those data derive from CAZy. 
A python script named **cazy_uniprot_dickeya.py** (./bin) allows access to UniProt and returns data from CAZy based on a specific family which we can pre define.

The database has been specified to be CAZy and the taxonomy to be Dickeya.

By executing **python3 cazy_uniprot_dickeya.py** (navigating to the **./bin**) generates a new directory under the **./data** directory named as **./data/cazy_dickeya** and stores all CAZy sequences for Dickeya, separated per CAZy family and stored under the appropriate name. 

The script:

[1] Retrieve all Dickeya entries with a cross-reference to CAZy Retrieve in tabular format: accession  CAZy_xrefs
       
[2] Create a dictionary mapping CAZy families to accessions e.g. {GH3: [P11073, D2BXL2, ...], ...}
 We use a default dict for this. As `caz` is just a string we use io.StringIO to handle it like a file whose lines we iterate over. As the first line contains the column headers, we ignore it. Each line is split into accession and a list of CAZy families; from this we construct the dictionary.

 [3] For each CAZy family, retrieve the mapped UniProtKB accessions in fasta format and write to a file with name <CAZy_family>.fasta. My simple approach here retrieves each batch of sequences from uniprot.org. One could also download all Dickeya sequences in one go and then get them from that file.


### 2. Processing CAZy raw data
 Run **python3 process_cazy_data.py**

The following were performed in order to process the data retrieved by Uniprot:

1. Bioservices were used to get access to Dickeya data instead of downloading manually the tabular file from Uniprot. 
2. The next step was to pass the query as pandas data frame. The query is a string at the therefore StringIO was used  for reading it as a data frame.
3. Next we get all all proteins without a locus tag using the **isnull**
4. Get all proteins wit a locus tag using **notnull**
6. Next we write the sequences with no locus tag in a fasta file. 
For that we are creating an empty list initially and then we append the sequences and write them out in a fasta file called: **all_uniprot_nolocus_results.fasta**
7. Finally we write a function for splitting cazy families. We do that as some enzymes belong to more than one CAZy family and therefore we want them as two separate entries. 
9. Then with SeqIO we write out in a fasta file all sequences which belong to a single family. So we first write out the sequences with no locus based on CAZy family and then the sequences with known locus based on CAZy family again and store the sequence results under the : /Users/eirinixemantilotou/Documents/PhD/IBioIC/Project/cazy_locus_tags/results/locus_seqs or /no_locus_seqs
12. We use the groupby function of pandas data frame and we count the sequences by locus and then we make a dataframe (called result) which returns the values based on locus which are in our groupby dataframe and have one or more occurrence for CAZy family. We can also pass them in a list using apply(list).tolist()
13. We turn the data frame in a list for no locus tags and we write out the uniprot accession numbers for those sequences.
14. We turn into list the datafrmae for the seqs with known locus tags and then write in a text file their locus tags so it can be used as input for RBBH.
15. We run the **RBBH_cazy.py** 
16. The python script opens the text file with all locus tags and reads every line and then using subprocess invoke the extract_rbbh.py python script which generates all RBBH for the enzymes from CAZy with known locus tags.

However, Uniprot failed to detect and downlowland a couple of CAZy families with Dickeya entries. Those are the following CAZyy FAMILIES: **CE8, CE1, CE4, CE9, CE11, CE12, PL26**. For those CAZY families the sequences and locus tags were downloaded manually by accessing CAZy, The locus tags for the CAZy families were not detecting by accessing Uniprot were merged into a unique text files names as **locus_Dickeya_not_det.txt**. We merge the txt files by running **python3 merge_locus.py**
Then we obtain the RBBH for those locus tags by running the **RBBH_not_det_uniprot.py** python script. the script runs from within the **process_cazy_data.py** script.

17. Run **merge_locus.py** python script to merge all locus tags
18. Run **RBBH_not_det_uniprot.py** python script for those locus tags were not downloaded directly from accessing uniprot. 
19. We select all rbbh.fasta files from the working directory
20. We make a dictionary with keys being the locus tag and values the corresponding CAZy family from the locus_split data frame. Then we relate the locus tag to the filenames which include the locus name at the start of the filename (we get that by splitting the name of the file).
21. We make directories based on the value names and then we copy the RBBH files to the corresponding directory based on the values (which is the CAZy family the RBBH result came from)
22. We however have some entries (locus) which belong to more than one CAZy family. We can create a dictionary called mydict which has more than one values for a single key and stores as list so we can extract the values as 1, 2, 3 etc
23. We filter the RBBH files for the ones which have more than one value (CAZy family)
24. We extract the CAZy families for those ones as value1 and value2
25. Create dictionaries for those values if they do not exist and then we copy the RBBH files to the first destination which is the first CAZy family and then we copy to the second destination which is the second CAZy family they belong to

### 3. Keeping unique fasta sequences as RBBH ouptup grouped per CAZy family (process_cazy-data.py)

The next step was to to make sure that the files containing the RBBH output result fasta files do not contain any duplicate sequences. We can do that by looping over all fasta files in each CAZy directory first, then adding in a list all record.id for the first file and writing those sequences in a new file named after the CAZy family. Then we check all other fasta files stored within the directory for record.id which are not in the list already. If there are not in the list then we append them to the new fasta file and then append them to the list too. If they already exist in the list then we skip those sequences. 

The new fasta files containing unique sequences after RBBH analysis and named after the CAZy family they belong to were stored under the **data/cazy_rbbh** directory.


 

### 4. Splitting large diverse CAZy families into subgroups(**automation_split_CAZy_families.py** )

As some of those families have quite a few representatives enzymes, the next step  was to split those families based on the RBBH table from the dickeya database.
By running the **split_CAZy_families_1.py** we can use sys.argv[] to automate the process for all large CAZy families. Finally, a python script named as **automation_split_CAZy_families_1.py** was created to pass as arguments all individual CAZy families. 
The code run as follows: 

$ python3 automation_split_CAZy_families_1.py
The script generates the following : 
1. **./data/locus/txt**
2. **data/locus/fasta**
Those two directories containing the locus per CAZy family in text file and the corresponding seqs in fasta. 
3. Under each CAZy family which need to be split a txt directory with all equivalent groups for each family and fasta directory with tehe sequences. 

### 4. Amending CAZy families

Some CAZy families include sequences which are very diverse in comparison to the other sequences and therefore in order to obtain a more reliable output indicating positive selection we will have to exclude those sequences form the MSA. We do that for four  CAZy families by running the **remove_ids.py** python script. Run **rename_dies.py** and **rename_dirs_2.py**  to rename the directories and merge data. ### 5. Input data for Positive selection analysisa.  Generate Input for codeml
________
The next step is to obtain MSA and phylogenetic trees using RaxML in order to use those data as input for positive selection analysis. We generate those data by running the **automation_repeat_split_families.py** which runs the **ps_automation_repeat_split_families.sh** shell script and the **ps_automation.py** which runs the **ps_automation.sh** shell script.
Both shell scripts invoke several python scripts in order to obtain the back translations, MSA of protein and nucleotides  sequences , different formats of the MSA'S (FASTA, CLUSTAL, RPHYLIP). 
Split the nucleotide sequence into petition, generate phylogenetic trees using RaxML, splitting the trees into subtrees by alternatively assigning a different branch node each time, modify the rphylip MSA in order to be readable by RaxML and organise teach repo in that way that we get access to the phylogenetic tree, control file for CodemML and MSA. 
Finally, the generated input data were transferred otto the local cluster. 

b. Cluster 
_________Shell scripts were written for each CAZy family and subgroup in order to run codeml for all trees within each CAZy family/subgroup   for the alternative and null model. Finally, a codeml mlc output file was generated for each tree in within each CAZy family and stored under the corresponding tree directory storing information about the maximum likelihood 
### 6. Process for generating Codeml output and filtered data
Note that in order to analyse the data some of the code request ETE and therefore need to be installed in the directory in advance. 

1. The first step is to obtain a complete data set for the maximum likelihood values of all CAZY families and their subgroups for the null and alternative model. Two functions **codeml_data** and **codeml_data_2** were passed into 2 separate python scripts named as **codeml_function.py** and **codeml_function_groups.py**. We invoke those two scripts and corresponding functions from within the **codeml.py**. the script generates two cvs files named as **codeml_results.csv** and **codeml_results_split.csv**

For those families/treedirs which the codeml output was not complete then I run: **bash codeml.sh $CAZYFAMILY** or navigate within the directory the data is missing and simply run: codeml and the type of control file we are providing (depending if it the null or alternative model). 

2. The next step was to combine the two columns about  the CAZy family name and the group into one running the **combine_col.r** code and the output file is named as **codeml_results_split_COMB.csv**
4. The two files (codeml_results_split_COMB.csv and codeml_results.csv ) were bind into one by running the **bind.R** code #1 and names as **codeml_bind_data.csv**

5. Next step was to split the **codeml_bind_data.csv** into individual families in order to calculate the fdr and q values. I achieved that by running the **split.R** code. The cvs files were stored under the **split_fams** working directory. **Current location being PS directory** 
6. The fdr values were then calculated by running the **fdr.R** code which invokes the fdr_calc from within the **fdr_function.R** script. The cvs files including the fdr values were stored under the **fdr** directory. Current directory split_fams
7. The csv's were combined into one cvs named as **codeml_bind_data_fdr.csv** using the **bind.R** code #2
Current working directory **./PS/frd**. We move the **codeml_bind_data_fdr.csv** to PS directory. 
8. Then, the qvalues were obtained for each family by running the **qvals.R** script. the script invokes the **qval_calc** function from within the **qvalue_function.R** script. The input files are the csv's from the fdr directory
9. The csv files with the fdr and qvalues were combined into one by running the #3 from the **bind.R** code
and the and the file named as **codeml_bind_data_qvals.csv** Current directory **qvalue**. Move data in PS directory.
10. Finally, the **codeml_bind_data_qvals.csv** can be filtered based on either the q values of the fdr values. The **fdr.R** and **qvalue.R** include code for filtering the data. 
11. Before processing the **codeml_bind_data_qvals.csv**, we need to remove those families/trees/ branch sites which are marked externally. We do that through the **ranking.py** script. 
a. In order to get the labelled #1 branches we use ETE. The python script **leaves.py** does 3 things: Remove any external labelled nodes, get the labelled nodes and create a dictionary with CAZy_family : leaves in order to generate heatmaps etc of the positive selected cazy families and the nodes they seem to be under positive selection.
12. Firstly, I  remove all those entries which are labelled externally by using the **external_nodes** function from the **leaves.py**. Running the first 5 cells form within the **ranking.ipynb**
13. The code writes out to a cvs file named as **no_external_nodes_ps.csv** those families/trees which are labelled internally 
14. Next step is  to filter those families by using the **filter_data.r** R code and filtering based on filter(log(qvalue) < -20).
15. The **filtered_cazy_fams** cvs file was obtained by running the **filter_data.R** script based on **no_external_nodes_ps.csv**. The cvs contains all the CAzy families which are above the threshold I have preselected and also provides information for the amount of times the CAZy family has been observed. 
16. A cvs file named as **final_ps.csv** was generated by running the **filter_data.R** including all families with qvals etc with the preselected threshold
17. We make the final cvs a data frame (ranking.ipynb) and create a dictionary using the **cazy_leaves** function (from the **leaves.py** script) to get the CAZy families and the labelled leaves. 
18. We make the leaves into groups so we can have access and finally write out in a cvs a file named **positive_selected_qvals.csv** which shows the CAZy families and groups are above the preselected threshold
Finally, a heat map using seaboard was generated to represent those positive elected families and the png was named as **heatmap_full_analysis.png**
