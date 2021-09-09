#!/bin/bash



POP1=$1
POP2=$2
vcf=$3

pop1bname=$(basename $POP1 .txt)
pop2bname=$(basename $POP2 .txt)
vcftools --gzvcf $vcf --weir-fst-pop $POP1 --weir-fst-pop $POP2 --out ${pop1bname}_vs_${pop2bname}.fst.txt

