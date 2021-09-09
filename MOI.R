library(SeqArray)
library(Rcpp)
library(gdsfmt)
gds.fn <- seqExampleFileName("gds")
gds.fn
seqSummary(gds.fn)
# open a GDS file
genofile <- seqOpen(gds.fn)
# display the contents of the SeqArray file in a hierarchical structure
genofile
# display all contents of the GDS file
print(genofile, all=TRUE, attribute=TRUE)

library(moimix)
setwd('/Users/georgeboateng-sarfo/Desktop/Bioinformatics/Sela')
system.file()
my_gds <- seqOpen('Variant_org_biallelic.gds')
seqSummary(my_gds)
# save sample identifiers
sample.id <- seqGetData(my_gds, "sample.id")
library(flexmix)
library(lattice)
NwCoord <- getCoordinates(my_gds)
head(NwCoord)
# first generate matrix of read counts for ref and alt alleles
counts <- alleleCounts(my_gds) 
counts
#################
seqOpen()
print(my_gds, all=TRUE, attribute=TRUE)
my_gds
# take out sample id
head(samp.id <- seqGetData(my_gds, "sample.id"))
# get "allele"
head(seqGetData(my_gds, "allele"))
z <- seqApply(my_gds, "allele",
              FUN=function(x) length(unlist(strsplit(x,","))), as.is="integer") 
table(z)
# covariance variable with an initial value
s <- 0

seqApply(my_gds, "$dosage", function(x)
{
  p <- 0.5 * mean(x, na.rm=TRUE)    # allele frequency
  g <- (x - 2*p) / sqrt(p*(1-p))    # normalization
  g[is.na(g)] <- 0                  # missing values
  s <<- s + (g %o% g)               # update the cov matrix 's' in the parent environment
}, margin="by.variant", .progress=TRUE)
# scaled by the number of samples over the trace
s <- s * (nrow(s) / sum(diag(s)))
# eigen-decomposition
eig <- eigen(s)

# eigenvalues
head(eig$value)
# eigenvectors
plot(eig$vectors[,1], eig$vectors[,2], xlab="PC 1", ylab="PC 2")
# get genomic coordinates of all variants
library(moimix, lexmix, lattice)
coords <- getCoordinates(my_gds)
head(coords)
#filter variant.ids not on apicoplast
seqSetFilter(my_gds, 
             variant.id = coords$variant.id[coords$chromosome != "Pf3D7_API_v3"])
isolate_baf <- bafMatrix(my_gds)
class(isolate_baf)
str(isolate_baf)
plot(isolate_baf, "0478-C-1-SK-1", main="Different types of variants", xlab="P. falciparum Chromosomes", ylab="Fws")
fws_all <- getFws(my_gds)
hist(fws_all, col = 'black', xlim = c(0.4,1), ylim = c(0,200), breaks = 5, border = 'red', main="Different types of variants", xlab="Within Sample F Stasticts ", ylab="Number of Sample")
qplot(fws_all, xlim = c(0.4,1), ylim = c(0,60), main="Different types of variants", xlab="Within Sample F Stasticts", ylab="Sample")
fws_all<0.95
head(fws_all)
###Drawing Hist of samples??
plot(fws_all)

########################################
# see if our sample that we estimated is multiclonal according to fws
fws_all["0456-C-1-SK-1"] < 0.95
isolate_baf
fws_all

pca1 <- read.table('VariantsPCA.eigenvec' ,sep=" ",header=F)
plot(data=pca1, PCA1~PCA2)
pca1

library(tidyverse)
library(broom)
library(ggfortify)
my_gds
eig
eig %>% tidy()
