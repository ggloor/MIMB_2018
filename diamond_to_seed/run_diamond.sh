#!/usr/bin/env bash

#24Jun2016 - JM
#update: Aug 15, 2017
# 1) Use diamond to compare all read files to the SEED fig.peg database
# 2) Merge all the best hits into a counts table

#--------------------------------------------------------------------------------------------------
# Setup
# This script should be in your working directory. Within your working directory should be the .gz
#	sequence files.
#--------------
# Paths to change:

DATA="" #the path to your reads (in .fastq.gz format)

DIAMOND="/Volumes/bin/diamond"			#diamond install location
DB="/Volumes/data/SEED_database/subsys4.dmnd"	#Path to SEED database (default: on Agrajag)
BIN="/Volumes/bin/diamond-program" #This is where blast_to_counts.pl and merge_counts.pl are. Default for Agrajag

#output paths
OUT1="diamond_output"
OUT2="diamond_output/counts_tables" #This is a temporary directory that will be deleted
#--------------------------------------------------------------------------------------------------
# Output
# There will be 1 output directory:
#	diamond_output will contain the diamond blast output (in .daa compressed format AND
#		decompressed to BLAST m8 format). Only one hit is retained per sequence by default.
#		Additionally, a single tab-delimited counts table will be built containing all samples:
#		all_counts.txt
#--------------------------------------------------------------------------------------------------

#----------------
# Part 1 - Use DIAMOND to compare reads per sample to the SEED database
#----------------

#make output directory if it doesn't already exist
mkdir -p $OUT1
echo -e "WARNING: output will be overwritten if files already exist\n"

for f in $DATA/*.gz; do	# e.g. F12_S17_L004_R2_001.fastq.gz

# Split on . and get the first field
	B=`basename $f`
	NAME=`echo $B | cut -d "." -f1`

#	echo $f
#	echo $B
#	echo $NAME
#	exit

#e.g. $DIAMOND blastx -d $DB -q ../data/sequence_files/F12_S17_L004_R2_001.fastq.gz -a diamond_output/F12_S17_L004_R2_001 --salltitles -k 3

	$DIAMOND blastx -d $DB -q ${f} -a ${OUT1}/${NAME} --salltitles -k 1
		# --salltitles           print full subject titles in output files
		# --max-target-seqs (-k) maximum number of target sequences to report alignments for
		# by default, diamond will use all available cores
		# -k is the number of hits to report

# Convert diamond output to blast tab-delimited format
	$DIAMOND view -a ${OUT1}/${NAME}.daa -o ${OUT1}/${NAME}.m8

done

#----------------
# Part 2 - From DIAMOND output, make counts tables
#----------------

# from DIAMOND output, make counts tables
# This is a per-file basis
# This will take a WHILE because Perl has to read line by line
#		Would have been smarter to do multiple files at a time - next time!

mkdir -p ${OUT2}	# Temporary directory to hold the individual counts tables

for f in ${OUT1}/*.m8; do

# Split on - and get the first field
	B=`basename $f`
	NAME=`echo $B | cut -d "." -f1`

#	echo $f
#	echo $NAME
# 	echo $OUT2
#	exit

	$BIN/blast_to_counts.pl $f $NAME $OUT2

done

# Merge all the counts tables into one
# Remove the individual tables when complete

	$BIN/merge_counts.pl $OUT2 > all_counts.txt

#echo $OUT2;exit

minsize=100
actualsize=$(wc -c <"all_counts.txt")

#Check the file is written (greater than 100kb) before deleting intermediary files
	if [[ $actualsize -gt $minsize ]]; then
		rm ${OUT2}/*_counts.txt
		rmdir $OUT2
		echo "output to all_counts.txt"
	else
		echo "all_counts.txt is empty. Please check individual counts tables and upstream scripts"
	fi

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#Future:

## Add checks before continuing
# Are there .gz files in the working dir
# Does the database exist
# Does diamond exist
# Are the perl scripts there
# Do not run diamond if output already exists
