#!/bin/bash

#Identifying population structure in Plasmodium falciparum genome sequence usinf vcf formats

conda creat -n vcfkit_Pop_Genetics

conda imstall -c bioconda vcftools -y

source activate vcfkit_Pop_Genetics

find usr/georgeboateng-sarfo/P.falciparum -iname '*.vcf' -exec mv '{}' usr/georgeboateng-sarfo/P.falciparum/Analysis_Sela \;

cpus=8

#filtering vcf file

vcftools --gzvcf Variant_orig_biallelic.vcf.gz --max-missing 0.5 --mac 3 --minQ 30 --recode --recode-INFO-all --out Variant_orig_biallelic.vcf

#Applying min depth for a genotype call and a min mean depth

vcftools --vcf Variant_orig_biallelic.vcf --minDP 3 --recode --recode-INFO-all --out Variant_orig_biallelic.vcf.g5mac3dp3

#removing poorly sequenced reads 

vcftools --vcf Variant_orig_biallelic.g5mac3dp3.vcf --missing-indv

#see the reads here:

cat out.imiss

#createing a list of individuals with more than 50% missing data
mawk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv

#Now that we have a list of individuals to remove, we can feed that directly into VCFtools for filtering
vcftools --vcf Variant_orig_biallelic.g5mac3dp3.vcf --remove lowDP.indv --recode --recode-INFO-all --out Variant_orig_biallelic.g5mac3dplm.vcf

#Now restrict the data to highly called variant reads

#vcftools --vcf Variant_orig_biallelic.g5mac3dp3.vcf --max-missing 0.95 --maf 0.05 --recode --recode-INFO-all --out Variant_orig_biallelic_DP3g95maf05.vcf --min-meanDP 20 

# BED files and associated plink format files

#slipt files into Asymptomatic and Asymptomatic VCFs

bcftools view Variant_orig_biallelic.g5mac3dp3.vcf -S Asymp.txt > AsympNT.g5mac3dp3.recode.vcf

bcftools view Variant_orig_biallelic.g5mac3dp3.vcf -S Clinical.txt > ClinicalNT.g5mac3dp3lm.recode.vcf.recode.vcf


#for Asymptomatic carraiers do

plink --vcf AsympNT.g5mac3dp3.recode.vcf --make-bed --chr-set 14 no-xy --out AsympNT.g5mac3dp3.recode

plink --bfile AsympNT.g5mac3dp3.recode.vcf --pca --chr-set 14 no-xy --out AsympNT.g5mac3dp3.recode

#for Symptomatic samples do

plink --vcf ClinicalNT.g5mac3dp3lm.recode.vcf.recode.vcf --make-bed --chr-set 14 no-xy --out ClinicalNT.g5mac3dp3lm.recode.vcf.recode

plink --bfile ClinicalNT.g5mac3dp3lm.recode.vcf.recode.vcf --pca --chr-set 14 no-xy --out ClinicalNT.g5mac3dp3lm.recode.vcf.recode

#Nucleotide diversity
#To calculate π over 100bp windows of the genome

vcftools --vcf AsympNT.g5mac3dp3.recode.vcf --window-pi 100 --out AsympNT.g5mac3dp3.recode_100
vcftools --vcf ClinicalNT.g5mac3dp3lm.recode.vcf.recode.vcf --window-pi 100 --out ClinicalNT.g5mac3dp3lm.recode.vcf.recode.vcf_100


#Tajima’s D
vcftools --vcf AsympNT.g5mac3dp3.recode.vcf --TajimaD 100 --out vcftools --vcf AsympNT.g5mac3dp3.recode_100
vcftools --vcf ClinicalNT.g5mac3dp3lm.recode.vcf.recode.vcf --TajimaD 100 --out vcftools --vcf ClinicalNT.g5mac3dp3lm.recode.vcf.recode_100

#Estimating population divergence with Fst

#vcftools --vcf AsympNT.g5mac3dp3.recode.vcf --weir-fst-pop population1 --weir-fst-pop population2 --fst-window-size 100 --out pop1_vs_pop2_FST_100
#vcftools --vcf ClinicalNT.g5mac3dp3lm.recode.vcf.recode.vcf --weir-fst-pop population1 --weir-fst-pop population2 --fst-window-size 100 --out pop1_vs_pop2_FST_100

echo "used R for Estimating MOI with Fws from"

conda deactivate 
