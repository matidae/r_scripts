library("baySeq")
data <- read.table("/home/matias/sl.y486.map/blat.data/out/sl.all.counts", row.names=1)
replicates <- c(rep("E", 5), rep("M", 4), rep("B",4))
cl <- NULL

groups <-list(NDE = c(rep("N",13)), DE = c(rep("E", 5), rep("M", 4), rep("B",4)))
data<-as.matrix(data)
CDEB <- new("countData", data=data[,-6:-9], replicates=c(rep("E", 5),rep("B",4)), groups=list(NDE = c(rep("N",9)), DE = c(rep("E", 5),rep("B",4))))
CDEM <- new("countData", data=data[,-10:-13], replicates=c(rep("E", 5),rep("M",4)), groups=list(NDE = c(rep("N",9)), DE = c(rep("E", 5),rep("M",4))))
CDMB <- new("countData", data=data[,-1:-5], replicates=c(rep("M", 4),rep("B",4)), groups=list(NDE = c(rep("N",8)), DE = c(rep("M", 4),rep("B",4))))

libsizes(CDEB) <- getLibsizes(CDEB)
libsizes(CDEM) <- getLibsizes(CDEM)
libsizes(CDMB) <- getLibsizes(CDMB)

plotMA.CD(CDEB, samplesA = "E", samplesB = "B", 
          col = c(rep("red", 100), rep("black", 900)), cex=.3, pch=19)

CDEB <- getPriors.NB(CDEB, samplesize = 1000, estimation = "QL", cl=cl)
CDEM <- getPriors.NB(CDEM, samplesize = 1000, estimation = "QL", cl=cl)
CDMB <- getPriors.NB(CDMB, samplesize = 1000, estimation = "QL", cl=cl)

CDEB <- getLikelihoods.NB(CDEB, pET ='BIC', cl = cl)
CDEM <- getLikelihoods.NB(CDEM, pET ='BIC', cl = cl)
CDMB <- getLikelihoods.NB(CDMB, pET ='BIC', cl = cl)

CDEB@estProps
topCounts(CDEB, group="DE", number=50000)

plotPosteriors(CDEB, group = "DE", samplesA="E", sampleB="B" ,col = c(rep("red", 100), rep("black", 900)), cex=.3, pch=19)

sink("sl.mb.bayseq")
topCounts(CDMB, group="DE", number=50000)
sink()
