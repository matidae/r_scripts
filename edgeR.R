library("edgeR")

counts <- read.table("/home/matias/sl.y486.map/blat.data/out/sl.all.counts", row.names=1)
conds <- c(rep("E", 5), rep("M", 4), rep("B",4))
cds <- DGEList(counts=counts, group=conds)

#Calculando el tamano de las librerias
cds <- calcNormFactors(cds)
cds <- estimateCommonDisp(cds)

de.com <- exactTest(cds, pair=c("E", "B"))

#Comparacion de grupos -1 down reg, 0 not significant, 1 up reg, p.value es de FDR
com <- decideTestsDGE(de.com, p.value=0.05)
summary(com)
topTags(de.com)

#Plot de las distancias de las muestras en terminos de log2fold_change 
etd <- estimateTagwiseDisp(cds)
etc <- estimateCommonDisp(cds)
plotBCV(etd)
plotMDS(etd, cex=.7)

#MA plot
de.eb.tgw <- exactTest(cds, dispersion = "auto", pair = c("E", "B"))
de.em.tgw <- exactTest(cds, dispersion = "auto", pair = c("E", "M"))
de.mb.tgw <- exactTest(cds, dispersion = "auto", pair = c("M", "B"))

resultsTbl.tgw <- topTags(de.eb.tgw , n = nrow(de.em.tgw$table))$table
de.genes.tgw <- rownames(resultsTbl.tgw)[resultsTbl.tgw$FDR < 0.05]

plotSmear(cds, de.tags=de.genes.tgw, main="Tagwise", pair = c("E","B") ,cex = .1,
           xlab="Log Concentration", ylab="Log Fold-Change", ylim=c(-7,9) )

abline(h=c(-2,2) , col="dodgerblue" )
length(de.genes.tgw )


sink("sl.eb.edger")
resultsTbl.tgw
sink()
