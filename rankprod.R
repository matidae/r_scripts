library(RankProd)
counts <- read.table("/home/matias/sl.y486.map/blat.data/out/sl.all.counts", header=FALSE, row.names=1)
data <- as.matrix(counts)
gnames <- rownames(data)
dataEB<-data[,-6:-9]
dataEM<-data[,-10:-13]
dataMB<-data[,-1:-5]

condsEB <- c(rep(0, 5), rep(1,4))
condsEM <- c(rep(0, 5), rep(1,4))
condsMB <- c(rep(0, 4), rep(1,4))

RP.out <- RP(dataMB, condsMB, gene.names=gnames, num.perm=100, logged=TRUE, na.rm=FALSE,plot=FALSE, rand=123)

topGene(RP.out, gene.names=gnames, cutoff=100)
head(RP.out)

sink("sl.mb.rankprod")
topGene(RP.out, gene.names=gnames, cutoff=100)
sink()
