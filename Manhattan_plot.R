library(qqman)
install.packages('qqman')
setwd('/Users/georgeboateng-sarfo/Desktop/Bioinformatics/Sela/manhattan-plots-main')
NewMan <- read.table('Asymptomatic_SampleID_vs_Clinical_SamplesIDs.fst.txt.weir.fst', header = TRUE )
NewManSub <- NewMan [complete.cases(NewMan),]
SNP<-c((1:nrow(NewManSub)))
mydf<- data.frame(SNP,NewManSub)
summary (mydf >1)
#Determine how many SNPs are on each Chromosome
as.data.frame(table(NewManSub$CHR))
#Plotting manhattan plots
manhattan(mydf, chr = 'CHROM', bp='POS', p='WEIR_AND_COCKERHAM_FST', snp="SNP", logp=FALSE, ylab= 'Fst Values', col = c("blue"))
manhattan(SNP)
manhattan(mydf)
str(gwasResults)
manhattan(gwasResults)

Trail <- read.table('BAN_vs_GAM.fst.weir.fst', header = TRUE)
Trail2 <- Trail [complete.cases(Trail),]
TSNP <- c((1:nrow(Trail2)))
mydf2 <- data.frame(TSNP,Trail2)
manhattan(mydf2, chr = 'CHROM', bp='POS', p='WEIR_AND_COCKERHAM_FST', snp="TSNP", logp=FALSE, ylab= 'Weir and Cockerham Fst', col = c("blue"))
