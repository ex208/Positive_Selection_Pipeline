#!/user/bin/sh
#
#automation.sh
#
# (c) The James Hutton Institute 2016
#Author: Eirini Xemantilotou
#
# The current shell script takes as input a locus tag and with the
# assistance of several python scripts it generates MSA and RaxML phylogenetic trees

# I put the special variable in order to pass the LOCUS TAG that I want to
# calculate as a cmd parameter. Database is the second variable for the RBBH command


# Make the folders with my tools (bin) and my fasta seqs more independent
SCRIPTS=./bin

# Under this folder I should keep the cds_aa.fasta file with all 29 Dickeya genomes and the Dickeya database
DATA=./data


# I also need a general directory  where I will keep all different subdirectories named after the locus tag name
RESULTS=./families


CAZYFAMILY=$1

CAZYFAMILY_DIR=./$RESULTS/$CAZYFAMILY
RBBH_DIR=./$RESULTS/$CAZYFAMILY/rbbh
RAXML_DIR=./$RESULTS/$CAZYFAMILY/raxml
MSA_DIR=./$RESULTS/$CAZYFAMILY/msa
PS_DIR=./$RESULTS/$CAZYFAMILY/positive_selection
ALT_DIR=./$RESULTS/$CAZYFAMILY/positive_selection/alternative
NULL_DIR=./$RESULTS/$CAZYFAMILY/positive_selection/null


# Use a for loop to generate all required directories
OUTDIRS=($SCRIPTS $DATA $RESULTS $CAZYFAMILY_DIR $TEMPLATES \
        $MSA_DIR $RAXML_DIR $RBBH_DIR $PS_DIR $ALT_DIR $NULL_DIR)
for directory in "${OUTDIRS[@]}"; do
    mkdir $directory
done

# Assign as prefix the subgroup
OUTPREFIX=$CAZYFAMILY.rbbh


#This command has as purpose to generate a MSA in standard format

echo "generate the MSA for the rbbh of $CAZYFAMILY utilising T-COFFEE"

# MSA for the protein sequence-standard output
t_coffee $DATA/cazy_rbbh/$OUTPREFIX.fasta


mv ./$OUTPREFIX.aln ./$OUTPREFIX.dnd \
   ./$OUTPREFIX.html  ./$OUTPREFIX.fasta.aln $MSA_DIR


# Next step is the generation of the back translation for which there is
# a py script. The output name is based on the locus tag and backtrans as its the corresponding nucleotide sequence
# to the protein sequence obtained with rbbh and no stops as the file is excluded for stop codons
# stop codons should not be in the file as MSA will not be generated for nucleotides as it will not match with protein MSA

# Set a variable for the output of the back_translation.py script

OUTPREFIX_2=$CAZYFAMILY.rbbh_backtrans_nostops

# Run python script and obtained corresponding nucleotide sequence
python $SCRIPTS/back_translations.py \
       $DATA/dickeya_cds_nt.fasta \
       $DATA/cazy_rbbh/$CAZYFAMILY.rbbh.fasta

cp $DATA/cazy_rbbh/$OUTPREFIX.fasta $RBBH_DIR
mv $DATA/cazy_rbbh/$OUTPREFIX.back_translations.fasta $RBBH_DIR
mv $DATA/cazy_rbbh/$OUTPREFIX_2.fasta $RBBH_DIR


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
# Add an I at the first line
python $SCRIPTS/clustal_to_rphylip.py \
       $MSA_DIR/$OUTPREFIX_2.aln

cp $MSA_DIR/$OUTPREFIX_2.rphylip $RAXML_DIR



# The input file for raxml must be a dna alignmennt sequence in relaxed
# phylip format or fasta
#The nice thing about rapid bootstrapping is that it allows you to do a
# complete analysis
# (ML search + Bootstrapping) in one single step by typing
#If called like this RAxML will do 100 rapid Bootstrap searches, 20 ML
# searches and #return the best tree

# This shell script includes a text file which splits the codon positions.
#We infer distinct model parameters jointly for all 1st and 2nd positions
#in the alignment and separately for the 3rd
# Introduce a small python script that writes a partition text file with
# the corresponding number of nucleotides

python $SCRIPTS/partition.py \
       $MSA_DIR/$OUTPREFIX_2.aln dna_12_3.partition_$CAZYFAMILY.txt

# Move the text file in the raxml directory
mv dna_12_3.partition_$CAZYFAMILY.txt $RAXML_DIR

# Generate raxml tree. Move to the appropriate directory first as several files will be generated

cd $RAXML_DIR

raxmlHPC -f a -m GTRGAMMA -p 12345 -q dna_12_3.partition_$CAZYFAMILY.txt \
         -x 12345 -# 100 \
         -s $OUTPREFIX_2.rphylip \
         -n $CAZYFAMILY.tree

# Back to the root
cd ../../..

# Assign tree as a variable
RAXML_TREE=RAxML_bestTree.$CAZYFAMILY.tree

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

#for i in $(ls *.nw); do mkdir ${i% .*}; done #I have done it through the python drzetree script

#for i in $(ls *.nw); do mv $1 ${i%.*}; done #I have done it through the python drztree script



#Next I need to make sure that I have a control template file in each directory and that the tree file
#is substituted every time with the corresponding
#I have done this again with the drzetree.py

python $SCRIPTS/drzetree.py\
       $PS_DIR/all_trees.nw ./control_file.ctl

mv  tdir* $ALT_DIR

python $SCRIPTS/drzetree.py\
       $PS_DIR/all_trees.nw ./control.fixed.ctl

mv  tdir* $NULL_DIR
# Add an I in the end of the phylip MSA file as codeml requires it.

python $SCRIPTS/modify_file.py $MSA_DIR/$OUTPREFIX_2.rphylip

mv $MSA_DIR/$OUTPREFIX_2.rphylip $MSA_DIR/backtrans_nostops.rphylip

echo $NULL_DIR/tdir* | xargs -n 1 cp $MSA_DIR/backtrans_nostops.rphylip

echo $ALT_DIR/tdir* | xargs -n 1 cp $MSA_DIR/backtrans_nostops.rphylip
