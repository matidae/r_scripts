library("DESeq")

counts <- read.table("/home/matias/dge/test.cds.counts", header=FALSE, row.names=1)

conds <- c(rep("E", 1), rep("M",1))
cds <- newCountDataSet(counts, conds)
esf <- estimateSizeFactors(cds)
esd <- estimateDispersions(esf, method="blind", sharingMode="fit-only")
resEM <- nbinomTest(esd, "E", "M")
resEBcds <- nbinomTest(esd, "E", "M")
resEBgeno <- nbinomTest(esd, "E", "M")
plotMA(resEM)
plotMA(resEBgeno)


plot(resEM$log2FoldChange ~ log(resEM$baseMean), pch=19, cex=.5,col=ifelse(resEM$padj < .05, "red", "black" ))

plot(resEBgeno$log2FoldChange ~ log(resEBgeno$baseMean), pch=19, cex=.5,col=ifelse( resEBgeno$padj < .1, "red", "black" ))

plotMA(resEBgeno[resEBgeno$baseMean >1 ,])
plotMA(resEBgeno[order(resEBgeno$padj),])
