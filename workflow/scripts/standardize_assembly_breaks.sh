#!/bin/bash

#Fasta names and breaklength values are passed from command line
FASTA_IN=$1
FASTA_OUT=$2
BREAKLEN=$3
NUM_REPETITIONS=$((BREAKLEN + 1))
REPEATED_STRING=$(printf "%*s" "$NUM_REPETITIONS" "" | tr " " "N")

#Remove line breaks from fasta
seqkit fx2tab $FASTA_IN > fasta_table.txt

#Replace each string of N's of any size to a string on N's 1 longer than the nucmer break length
sed -E "s/N+/$REPEATED_STRING/g" fasta_table.txt > fasta_table.newbreaks.txt

#Reformat as fasta
seqkit tab2fx fasta_table.newbreaks.txt > $FASTA_OUT

#Remove temporary tables
rm fasta_table.txt
rm fasta_table.newbreaks.txt
