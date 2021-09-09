library(SeqArray)
library(Rcpp)
library(moimix)
library(ggplot2)
library(qqman)
library(tidyverse)
library(ggplot2)
library(PopGenome)

setwd('/Users/georgeboateng-sarfo/Desktop/Bioinformatics/Sela/Asym_Parameters')
seqVCF2GDS('AsymptoFiltered.vcf', 'Asymptomatic.gds')
Asymto_gds <- seqOpen('Asymptomatic.gds')
seqSummary(Asymto_gds)
Asym_sampleID <- seqGetData(Asymto_gds, 'sample.id')
Asymto_coord <- getCoordinates(Asymto_gds)
head(Asymto_coord)

#Estimate BAF matrix
Asym_BAF <- bafMatrix(Asymto_gds)
str(Asym_BAF)
plot(Asym_BAF$baf_site)
Asym_BAF < 0.95

#Draw histogram for MIO 
#Estimating MOI with binommix
seqSetFilter(Asymto_gds, 
             variant.id = Asymto_coord$variant.id[Asymto_coord$chromosome != '0985-7-SK-1'])
Asympt_counts <- alleleCounts(Asymto_gds)
m1 <- binommix(Asympt_counts, sample.id = "0985-7-SK-1", k = 2)
summary(m1)
param.estimates <- getTheta(m1)
param.estimates

#Estimating MOI with Fws
Asymp_Samples <- getFws(Asymto_gds)
hist(Asymp_Samples, col = 'gray', xlim = c(0.5,1), ylim = c(0,40), breaks = 5, border = 'black', main="Fws Statistics of Asymptomatic Patients", xlab="Within Sample F Stasticts ", ylab="Number of Samples")
abline(v = 0.95, lwd = 3, col = "blue")
qplot(Asymp_Samples, col= "grey", binwidth=0.1, main="Fws Statistics of Asymptomatics", xlab="Within Sample F Stasticts ", ylab="Number of Sample")
geom_abline(xintercept= 0.95)

ggplot()
Asymp_Samples > 0.95
write.csv (Asymp_Samples)
summary(Asymp_Samples)
summary(Asymp_Samples < 0.95)

#Population structure Analysis for Clinical Samples
setwd('/Users/georgeboateng-sarfo/Desktop/Bioinformatics/Sela/Clinical_Parameters/')
seqVCF2GDS('ClinicalFiltered.vcf', 'ClinicalFiltered.gds')
Clinical_gds <- seqOpen('ClinicalFiltered.gds')
seqSummary(Clinical_gds)
Clinical_coord <- getCoordinates(Clinical_gds)
Clinical_BAF <- bafMatrix(Clinical_gds)
str(Clinical_BAF)
qplot(Clinical_BAF, "0640-C-1-SK-1")

#Draw histogram for Clinical MIO 
Clinical_Samples <- getFws(Clinical_gds)
hist(Clinical_Samples, col = 'grey', xlim = c(0.3,1), ylim = c(0,80), breaks = 5, border = 'black', main="Fws Statistics for Symptomatic patients", xlab="Within Sample F Stasticts", ylab="Number of Samples")
abline(v = 0.95, col = 'red', lwd = 3)
qplot(Clinical_Samples, col = 'grey', main="Fws MOI Statistics for Symptomatic patients", xlab="Within Sample F Stasticts", ylab="Number of Samples")
Clinical_Samples > 0.95 
summary(Clinical_Samples > 0.95)
summary(Clinical_Samples < 0.95)

#Plotting nucleotide diversity for Asymptomatics
Asympto_ND <- read.table("Asympto100.windowed.pi", header = T)
hist(Asympto_ND$BIN_START,Asympto_ND$PI, br=1, col = 'grey', xlim = c(0.0,0.008), ylim = c(0,20), breaks = 5, border = 'red', main = 'Nucleotide Diversity amongst Asymtomatic Samples', xlab = 'Estimates Of Nucleotide diversity', ylab = 'Diversity')
qplot(Asympto_ND$BIN_START, Asympto_ND$PI, col = 'grey', main = 'Nucleotide Diversity amongst Asymptomatics Samples', xlab = "Chromsome Position", ylab = "Nucleotide diversity (π)")
PI <- subset(Asympto_ND, CHROM == 'chr1')
plot(Asympto_ND$BIN_START,Asympto_ND$PI, xlab="Chromsome Position",ylab="Nucleotide diversity (π)")
boxplot(Asympto_ND$BIN_START, Asympto_ND$PI, ylab = "Diversity")
summary(Asympto_ND)

##Plotting the PCA results in R for Asymtopmatics
getwd()
setwd('/Users/georgeboateng-sarfo/Desktop/Bioinformatics/Sela/Asym_Parameters')
Asymtopca1 <- read.table('AsymptomaticPCA.eigenvec', sep = ' ', header = F)

# read in data
newpca <- read_table2("AsymptomaticPCA.eigenvec", col_names = FALSE)
Asymptneweigenval <- scan("AsymptomaticPCA.eigenval")
# sort out the pca data
# remove nuisance column
newpca <- newpca[,-1]
# set names
names(newpca)[1] <- "ind"
names(newpca)[2:ncol(newpca)] <- paste0("PC", 1:(ncol(newpca)-1))
newpca <- as.tibble(data.frame(newpca))
# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = Asymptneweigenval/sum(Asymptneweigenval)*100)
# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light() +


# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

# sort out the individual species and pops
# spp
spp <- rep(length(newpca$ind))
spp[grep("Plasmod", newpca$ind)] <- "Plasmodium"

# location
loc <- rep(length(newpca$ind))
loc[grep("pf", newpca$ind)] <- "Plasmodium"

# combine - if you want to plot each in different colours
spp_loc <- paste0(spp, "_", loc)
# plot pca
b <- ggplot(newpca, aes(PC1, PC2, col = spp, shape = loc)) 
b <- b + scale_colour_manual(values = c("red", "blue"))
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))



Asymtopca1
plot(data=newpca, PC1~PC2, main= "Principal Component Analysis of Asymptomatic Pateints")
summary(Asymtopca1)

##Plotting Tajima's D 
Asymp_Tajima<- read.table("AsymptTD_100bp.Tajima.D", header = T)
hist(Asymp_Tajima$TajimaD, col = 'grey', border = 'red', main = "Tajima's D amongst Asymtomatic Samples", xlab = 'Position', ylab = "Tajima's D")
qplot(Asymp_Tajima$TajimaD, col = 'grey', xlim = c(-0.5,3), ylim = c(0,10), border = 'red', main = "Tajima's D amongst Asymptomatic Samples", xlab = 'Position', ylab = "Tajima's D")
qplot(Asymp_Tajima$BIN_START,Asymp_Tajima$TajimaD,xlab="chromosome position",ylab="Tajima's D", main = "Tajima's D estimates for Asympotmatic Patients")
hist(Asymp_Tajima$BIN_START, Asymp_Tajima$TajimaD, xlab="chromosome position", ylab="Tajima's D")



#Analysis for Clinical Samples
setwd('/Users/georgeboateng-sarfo/Desktop/Bioinformatics/Sela/Clinical_Parameters/')
seqVCF2GDS('ClinicalFiltered.vcf', 'ClinicalFiltered.gds')
Clinical_gds <- seqOpen('ClinicalFiltered.gds')
seqSummary(Clinical_gds)
Clinical_coord <- getCoordinates(Clinical_gds)
Clinical_BAF <- bafMatrix(Clinical_gds)
str(Clinical_BAF)
qplot(Clinical_BAF, "0640-C-1-SK-1")
View(Clinical_Samples)
names(Clinical_Samples)

#Draw histogram for MIO 
Clinical_Samples <- getFws(Clinical_gds)
hist(Clinical_Samples, col = 'blue', xlim = c(0.4,1), ylim = c(0,80), breaks = 5, border = 'red', main="Fws Statistics for Symptomatic patients", xlab="Within Sample F Stasticts", ylab="Number of Sample")
summary(Clinical_Samples < 0.95)

#Plotting nucleotide diversity
Clinical<- read.table("Clinical100.windowed.pi", header = T)
hist(Clinical$PI,br=10)
pi.chr1 <- subset(Clinical_ND, CHROM == "chr1")
boxplot(pi.chr1)
plot(pi.chr1$BIN_START, pi.chr1$PI, xlab = "positon", ylab = "Diversity")
qplot(Clinical$BIN_START, Clinical$PI, col = 'grey', main = 'Nucleotide Diversity amongst Clinical Samples', xlab = 'Estimates Of Nucleotide diversity', ylab = 'Diversity')
boxplot(Clinical_ND$PI,ylab = "Diversity")
qplot(Clinical$BIN_START, Clinical$PI, xlab="Chromsome Position", ylab="Nucleotide diversity π")
box(Clinical,ylab ="Nucleotide diversity (π)")

barplot()
library(dplyr)
library(tidyverse)
df <- data.frame(Clinical_ND, header = T)
df
ggplot(df)


##Plotting the PCA results in R
getwd()
setwd('/Users/georgeboateng-sarfo/Desktop/Bioinformatics/Sela/Clinical_Parameters/')
Clinicalpca1 <- read.table('ClinicalPCAPara.eigenvec', sep = ' ', header = F)
Clinicalpca1
plot(data=newpca, PC1~PC2, main = "Principal Component Analysis of Clinical Pateints")
summary(Clinicalpca1)
newpca <- read_table2("ClinicalPCAPara.eigenvec", col_names = FALSE)
Clinicaleigenval <- scan("ClinicalPCAPara.eigenval")
# sort out the pca data
# remove nuisance column
newpca <- newpca[,-1]
# set names
names(newpca)[1] <- "ind"
names(newpca)[2:ncol(newpca)] <- paste0("PC", 1:(ncol(newpca)-1))
# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = Clinicaleigenval/sum(Clinicaleigenval)*100)
# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()


##Plotting Tajima's D 
ClinicalTajima<- read.table("ClinicalTajD100.Tajima.D", header = T)
hist(ClinicalTajima$TajimaD, col = 'grey', xlim = c(-0.5,3), ylim = c(0,30), border = 'red', main = "Tajima's D amongst Clinical Samples", xlab = 'Position', ylab = "Tajima's D")
qplot(ClinicalTajima$TajimaD, col = 'grey', xlim = c(-0.5,3), ylim = c(0,10), border = 'red', main = "Tajima's D amongst Clinical Samples", xlab = 'Position', ylab = "Tajima's D")

# Fws Analysis of every variants 
setwd('/Users/georgeboateng-sarfo/Desktop/Bioinformatics/Sela')
WG_vcf <- seqOpen('Variant_org_biallelic.gds')
seqSummary(WG_vcf)
WG_coord <- getCoordinates(WG_vcf)
WG_baf <- bafMatrix(WG_vcf)
str(WG_baf)
plot(WG_baf, "0640-C-1-SK-1")

#plotting nucleotide diversity for clinical Samples
ClinicND <- read.table("Clinical100.windowed.pi", header = T)
hist(ClinicND$PI, br=10)
boxplot(ClinicND$PI, ylab="Diversity")
qplot(ClinicND$BIN_START,ClinicND$PI, xlab="Chromosome Position", ylab='Nucleotide Diversity π', main = "Nucleotide Diversity amongst Symptomatic Patients")
ClinicalTD <- read.table("ClinicalTajD100.Tajima.D", header = T)
qplot(ClinicalTD$BIN_START,ClinicalTD$TajimaD, xlab = 'Chromsome Position', ylab = "Tajima's D", main = "Tajima's D estimates for Symptomatic Patients")


#Draw histogram for MIO 
WG_Samples <- getFws(WG_vcf)
WG_Samples < 0.95
hist(WG_Samples, col = 'blue', xlim = c(0.4,1), ylim = c(0,200), breaks = 5, border = 'red', main="Fws Statistics for Total population", xlab="Within Sample F Stasticts", ylab="Number of Sample")

#Plotting barplot for Asymptomatic and Symptomatic samples
AsymClinic <- cbind(Asymp_Samples, Clinical_Samples)
barplot(AsymClinic, beside = T)
barchart(AsymClinic)


sessionInfo()
