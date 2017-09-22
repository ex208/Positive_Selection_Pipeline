#!/user/bin/sh
#
#ps_automation_repeat_split_families.sh
#
# (c) The James Hutton Institute 2017
#Author: Eirini Xemantilotou
#

```
The current shell script takes as input a locus tag and with the
assistance of several python scripts it generates
MSA and RaxML phylogenetic trees
```
# Make the folders with my tools (bin) and my fasta seqs more independent
SCRIPTS=./bin

# Under this folder I should keep the cds_aa.fasta file with all 29 Dickeya genomes and the Dickeya database
DATA=./data

# I also need a general directory  where I will keep all different subdirectories named after the locus tag name
RESULTS=./families

# Set the first commalnd line argument being the examining CAZy family
CAZYFAMILY=$1
# Set the second command line argument being the sugroup of teh CAZy family 
CAZYFAMILY_SUBDIR=$2

# Assign the working directories we will generate
CAZYFAMILY_DIR=./$RESULTS/$CAZYFAMILY
SUB_CAZYFAMILY_DIR=./$RESULTS/$CAZYFAMILY/$CAZYFAMILY_SUBDIR
RBBH_DIR=./$RESULTS/$CAZYFAMILY/$CAZYFAMILY_SUBDIR/rbbh
RAXML_DIR=./$RESULTS/$CAZYFAMILY/$CAZYFAMILY_SUBDIR/raxml
MSA_DIR=./$RESULTS/$CAZYFAMILY/$CAZYFAMILY_SUBDIR/msa
PS_DIR=./$RESULTS/$CAZYFAMILY/$CAZYFAMILY_SUBDIR/positive_selection
ALT_DIR=./$RESULTS/$CAZYFAMILY/$CAZYFAMILY_SUBDIR/ \
        positive_selection/alternative
NULL_DIR=./$RESULTS/$CAZYFAMILY/$CAZYFAMILY_SUBDIR/ \
        positive_selection/null
FASTA=./$RESULTS/$CAZYFAMILY/fasta

# Use a for loop to generate all required directories
OUTDIRS=($SCRIPTS $DATA $RESULTS $CAZYFAMILY_DIR \
        $MSA_DIR $RAXML_DIR $RBBH_DIR $PS_DIR \
        $ALT_DIR $NULL_DIR)
for directory in "${OUTDIRS[@]}"; do
    mkdir $directory
done

# Assign as prefix the subgroup
OUTPREFIX=$CAZYFAMILY_SUBDIR

# This command has as purpose to generate a MSA in standard format
echo "generate the MSA for the rbbh of $CAZYFAMILY utilising T-COFFEE"

# MSA for the protein sequence-standard output using t-coffee
t_coffee $FASTA/$OUTPREFIX.fasta
mv ./$OUTPREFIX.aln ./$OUTPREFIX.dnd ./$OUTPREFIX.html \
   ./$OUTPREFIX.fasta.aln  $MSA_DIR

# Set a variable for the output of the back_translation.py script
OUTPREFIX_2=$CAZYFAMILY_SUBDIR.rbbh_backtrans_nostops

# Copy the rbbh fasta files into the RBBH correspoing directory
cp $FASTA/$OUTPREFIX.fasta $RBBH_DIR

# Run python script and obtained corresponding nucleotide sequence for protein sequence 
python $SCRIPTS/back_translations_2.py \
        $SCRIPTS/dickeya_cds_nt.fasta \
        $RBBH_DIR/$OUTPREFIX.fasta

# Move the geberated files to the rbbh directory
mv $CAZYFAMILY_DIR/$OUTPREFIX.back_translations.fasta $RBBH_DIR
mv $CAZYFAMILY_DIR/$OUTPREFIX_2.fasta $RBBH_DIR

# MSA for the nucleotide sequence - clustalw output
t_coffee -other_pg seq_reformat \
         -in $RBBH_DIR/$OUTPREFIX_2.fasta \
         -in2 $MSA_DIR/$OUTPREFIX.aln -action +thread_dna_on_prot_aln \
         -output clustalw > $OUTPREFIX_2.aln

# Reformat clustalw alignment - phylip format
t_coffee -other_pg seq_reformat \
         -in ./$OUTPREFIX_2.aln \
         -output phylip > $OUTPREFIX_2.phylip

# Reformat clustalw alignment - fasta format
t_coffee -other_pg seq_reformat \
         -in ./$OUTPREFIX_2.aln \
         -output fasta_aln > $OUTPREFIX_2.fasta

# Move generated MSA files in MSA directory from the working directory
mv ./$OUTPREFIX_2.aln ./$OUTPREFIX_2.phylip ./$OUTPREFIX_2.fasta $MSA_DIR

# Python script to generate clustal_to_rphylip for 'relaxed PHYLIP' format.
# Files name is based on locus tag plus the rphylip extension
# An I was added at the first line
python $SCRIPTS/clustal_to_rphylip.py \
       $MSA_DIR/$OUTPREFIX_2.aln

cp $MSA_DIR/$OUTPREFIX_2.rphylip $RAXML_DIR

# python script for partitiotion of the dna alignment and write it out as a text file 
python $SCRIPTS/partition.py \
       $MSA_DIR/$OUTPREFIX_2.aln \
       dna_12_3.partition_$CAZYFAMILY.$CAZYFAMILY_SUBDIR.txt

# Move the text file in the raxml directory
mv dna_12_3.partition_$CAZYFAMILY.$CAZYFAMILY_SUBDIR.txt $RAXML_DIR

# Generate raxml tree. Move to the appropriate directory first as several files will be generated
cd $RAXML_DIR
raxmlHPC -f a -m GTRGAMMA -p 12345 \
         -q dna_12_3.partition_$CAZYFAMILY.$CAZYFAMILY_SUBDIR.txt \
         -x 12345 -# 100 \
         -s $OUTPREFIX_2.rphylip \
         -n $CAZYFAMILY.$CAZYFAMILY_SUBDIR.tree

# Back to the root
cd ../../../..

# Assign tree as a variable
RAXML_TREE=RAxML_bestTree.$CAZYFAMILY.$CAZYFAMILY_SUBDIR.tree

# So long as your script produces tree files with marked branches (maybe all
# in their own subdirectories, with the corresponding control files), doing
# this will completely automate the process of generating dN/dS analyses.
# The tricky bit is using ETE to mark the branches.

python3 $SCRIPTS/ete.py\
        $RAXML_DIR/$RAXML_TREE > all_trees.nw

# Move text file in palm directory
mv ./all_trees.nw  $PS_DIR

#Here I need to introduce some code for splitting the file all_trees.nw into individual files called #tree1.nw, tree2.nw, tree3.nw, tree4.nw etc

# Now I have all possible tree files in a file called all_trees.nw I can split this
# file into subfiles and store them in individual directories called Tree1, Tree2, Tree3
# etc

#Next step is to direct the files into directories called after the files name but without the .nw.
#We first have to create the directories and then we have to move the files in the right directories


python $SCRIPTS/drzetree.py \
       $PS_DIR/all_trees.nw ./control_file.ctl

mv  tdir* $ALT_DIR

python $SCRIPTS/drzetree.py\
       $PS_DIR/all_trees.nw ./control.fixed.ctl

mv  tdir* $NULL_DIR

# Add an I in the end of the phylip MSA file as codeml requires it
python $SCRIPTS/modify_file.py $MSA_DIR/$OUTPREFIX_2.rphylip

mv $MSA_DIR/$OUTPREFIX_2.rphylip $MSA_DIR/backtrans_nostops.rphylip

echo $NULL_DIR/tdir* | xargs -n 1 cp $MSA_DIR/backtrans_nostops.rphylip

echo $ALT_DIR/tdir* | xargs -n 1 cp $MSA_DIR/backtrans_nostops.rphylip
