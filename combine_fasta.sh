#!/bin/bash

dir=$1
cat $dir/*.fasta >allseqs.fasta


num=$(grep '>' allseqs.fasta|wc -l)

echo "$num sequences identified and combined"
head allseqs.fasta

echo "sequences saved as allseqs.fasta"
