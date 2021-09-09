#!/bin/bash

multi_fasta=$1

bname=$(basename $multi_fasta .fasta)
apps/muscle -in $multi_fasta -out ./${bname}.msa.fasta
apps/muscle -in ./${bname}.msa.fasta -out ./${bname}.refined.phylip -refine -phyi

echo "File saved as ${bname}.refined.phylip"

