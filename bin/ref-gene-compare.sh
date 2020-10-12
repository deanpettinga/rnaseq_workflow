#! /bin/bash

# get annotation field
awk '{print $9}' ../../refs/HanXRQr2.0-SUNRISE-2.1/annotation/HanXRQr2.0-SUNRISE-2.1.gff3 |\
 # get only rows with 'gene' annotation
 grep 'gene:' | awk 'BEGIN {FS = ":"} ; {print $2}' |\
 # extract just the gene name
 awk 'BEGIN {FS = ";"} ; {print $1}' |\
 # uniquify and output to file
 uniq > ../analysis/ref-gene-compare/HanXRQr2.0-SUNRISE-2.1.genes.tsv
