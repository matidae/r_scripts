library("DESeq")

counts <- read.table("/home/matias/sl.y486.map/blat.data/out/sl.all.counts", header=FALSE, row.names=1)
conds <- c(rep("E", 5), rep("M", 4), rep("B",4))

#Crea el objeto CountDataSet a partir de una matriz o dataframe
cds <- newCountDataSet(counts, conds)

#Normalizacion. Estimacion del tamaño efectivo de la liberia.
esf <- estimateSizeFactors(cds)

#Estimacion de la varianza. La dispersion el cuadrado del coeficiente de variacion biologica.
#Ej si un gen varia un 20% entre 2 muestras, su dispersion en 0.04
#Primero estima los valores de dispersion p/c/gen. Luego ajusta una curva y por ultimo
#asigna a cada gen un valor de dispersion eligiendo entre la curva y el estimado
esd <- estimateDispersions(esf, sharingMode="fit-only", fitType="local")
#Se puede explorar esta estimacion de dispersion con fitInfo, fData y plotDispEsts
fitInfo(esd)
plotDispEsts(esd)
head(fData(esd))

#Inferencia de expresion diferencial, asumiendo distribucion binomial negativa
resEB <- nbinomTest(esd, "E", "B")
resEM <- nbinomTest(esd, "E", "M")
resMB <- nbinomTest(esd, "M", "B")
#Crea un dataframe con las cols.: id, baseMean (all), baseMeanA (sampleA), baseMeanB (sampleB),
#foldchange, log2FoldChange, pval, padj (pval ajustado para muestras multiples por Benjamini-Hochberg, para controlar FDR)

#Plot log2FoldChange~baseMean (significativos a 5% FDR)
plotMA(resEB,col=ifelse(resEB$padj<=0.05, "red", "black"))
#Volcano plot
plot(resMB$log2FoldChange,-log10(resMB$padj),cex=.1,col=ifelse( resMB$padj < .05, "red", "black" ), 
     xlim=c(-6,9), ylim=c(0,150))

#Histograma de p-values
hist(resMB$padj, breaks=100, col="skyblue", border="slateblue", main="", xlab="p-value adj")

#Filtrado de genes segun FDR
resEBSig = resEB[resEB$padj < 0.05,]
#Lista de 50 genes con mas significancia estadistica de expresion diferencial
head(resEBSig[order(resEBSig$pval), ]$id, 50)
#Most strongly up-regulated genes
head( resEBSig[order( -resEBSig$foldChange, -resEBSig$baseMean ), ]$id, 50)
#Most strongly down-regulated genes
head( resEBSig[order( resEBSig$foldChange, -resEBSig$baseMean ), ]$id, 50)

##################################################
a<-subset(resEB, resEB$pval<0.05)
write(t(as.matrix(a)), file="asd2", append=FALSE, ncolumns=8)
nrow(subset(resMB, resMB$pval<0.05))
nrow(subset(resMB, resEB$padj<0.05))

dfit0 <- fitNbinomGLMs(cds, counts~conds)

##Plot log2foldChange
cn=resEBSig$log2FoldChange
plot(cn, col= ifelse((cn >= 2 | cn<=-2), "red", ifelse(cn <= 1,"black", "darkgreen")), pch=19, cex=.5)
abline(2,0)
resEBSig$foldChange

resEBSig2 = resEBSig[resEBSig$log2FoldChange>2,]
n<-na.omit(resEB)
b<-as.matrix(cbind(n$id,n$foldChange))
write(b, file="alal2")

hist(resEB$log2FoldChange, breaks=100, col="skyblue", 
     border="slateblue", main="", xlab="p-value adj", xlim=c(-4,4))
mean(na.omit(resEB$log2FoldChange))
