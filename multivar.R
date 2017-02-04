library(Biostrings)
library(ggplot2)
library(rgl)
library(RColorBrewer)
library(vegan)
library(MASS)

#Parametros a definir: tamaño de la palabra (kmer) y archivo fasta de entrada
kmer = 3
file_name = '/home/matias/data/y486.split.5.norepe2.fa'

# Función que determina la frecuencia de kmers para cada uno de los elementos de 
# una colección de strings (i.e. archivo multifasta)
calc_freqs <- function(file_name, kmer) {
  #Carga archivo fasta
  fasta <- readDNAStringSet(file_name, "fasta")
  #Cuenta la frecuencia de kmers en la secuencia original y su reverso 
  #complementaria para evitar sesgo de hebra
  count_plus <- oligonucleotideFrequency(fasta, kmer)
  count_minus <- oligonucleotideFrequency(reverseComplement(fasta), kmer)
  #Suma de ambos conteos y cálculo de frecuencias como dataframe.
  kmer_freqs <- as.data.frame(cbind(prop.table(count_plus + count_minus, 1)))
  #Obtengo nombre de las secuencias y los agrego al dataframe
  seq_names <- names(fasta)
  row.names(kmer_freqs) <- seq_names
  #Retorna el conteo de kmers para cada secuencia y el nombre de secuencia asociado
  return(kmer_freqs)
}

#-------------------------------------------------------------------------------

#Cálculo de la frecuencia de k-nucleótidos (en este caso 3)
kmer_freqs <- calc_freqs(file_name, kmer)

#Mini-hack para quedarme con los nombres de las 32 vars diferentes
uniq3mers <- rownames(unique(t(kmer_freqs[4,])))
kmer_freqs_uniq <- as.data.frame(kmer_freqs[, uniq3mers])

#PCA con la función princomp, escalado y usando matriz de correlación
pca <- princomp(kmer_freqs_uniq, scale=TRUE, cor = TRUE)

#Screeplot para visualizar el aporte de cada componente a la varianza
pvar <-as.data.frame(pca$sdev[1:10]^2/sum(pca$sdev^2))
colnames(pvar)<-c("prop_var")

ggplot(pvar, aes(x=seq(1,10),y=prop_var)) + geom_bar(stat="identity") + 
  scale_x_continuous(breaks=seq(1,10)) +
  labs(y= "proporcion de varianza", x = "componentes principales")

#Plot que muestra la varianza acumulada
ggplot(cumsum(pvar), aes(x=seq(1,10),y=prop_var)) + geom_line() + geom_point() +
  scale_x_continuous(breaks=seq(1,10)) +
  scale_y_continuous(limits=seq(0,1), breaks=seq(0,1,0.2)) +
  labs(y= "proporcion de varianza acumulada", x = "componentes principales")

#Tabla con valores de varianza observada y suma acumulada
head(cbind(pvar, cumsum(pvar)), n=3)

#Cargo datos de los 2 indices y contenido GC
z<-read.table(file='/home/matias/data/y486.split.5.norepe2.allindex', header=FALSE)
zindex<-cbind(z$V4,z$V7,z$V8)
colnames(zindex) <- c("stopOE","GC","orfindex")

#Construyo dataframe para ggplot
pca_data <- data.frame(PC1 = pca$scores[,1], PC2 = pca$scores[,2],
                       PC3 = pca$scores[,3], z_val = zindex)

#Creo un dataframe para los loadings (aporte de cada variable original)
x1<-summary(pca)$loadings[,1]
y1<-summary(pca)$loadings[,2]
loadings <- as.data.frame(cbind(x1, y1))
colnames(loadings)<-c('x1','y1')

#Renombro los 2 codones stop que estan como su reverso complementario 
#para identificarlos mejor
rownames(loadings)[rownames(loadings) == "TCA"] <- "TGA"
rownames(loadings)[rownames(loadings) == "CTA"] <- "TAG"

#Extraigo los 3 codones stop para facilitar identificación
loadings_stop <- loadings[c("TGA", "TAG", "TAA"),]
loadings <- loadings[!rownames(loadings) %in% c("TGA", "TAG", "TAA"),]

#-------------------------------------------------------------------------------

#Plot de PCA utilizando los 3 primeros componentes principales
colFunc <- colorRamp(c("darkgreen","green","orange", "darkred"))
cols <- colFunc(pca_data$z_val.orfindex/max(pca_data$z_val.orfindex))
plot3d(x=pca_data$PC1, y=pca_data$PC2, z=pca_data$PC3, type="p", 
       col=rgb(cols, maxColorValue=255), size=2, box=F)

#Plot de PCA + vectores de aporte de variables originales, coloreando cada punto segun valor de:
##stopOE
ggplot(pca_data, aes(PC1, PC2)) +  geom_point(shape=19, size=.5, alpha=.5, 
                                              aes(colour=z_val.stopOE)) +
        scale_colour_gradientn(colours=c("darkgreen","green","orange", "darkred"), 
                                limits=c(0.3,0.8), na.value="grey") +
                                labs(colour="stopOE", x="PC1", y="PC2")  +
  geom_segment(aes(xend=x1*20, yend=y1*20, x=rep(0,29), y=rep(0,29)), data=loadings, 
               arrow = arrow(length = unit(0.1,"cm")),alpha=0.5, lwd=.5) +
  geom_text(alpha=.9, data=loadings, size=2.8, fontface="bold", 
            aes(x=x1*25, y=y1*23, label=rownames(loadings))) +
  geom_segment(aes(xend=x1*20, yend=y1*20, x=rep(0,3), y=rep(0,3)), data=loadings_stop, 
               arrow = arrow(length = unit(0.1,"cm")),alpha=0.7, lwd=1, color="red") +
  geom_label(alpha=.5, data=loadings_stop, 
             aes(x=x1*30, y=y1*30, label=rownames(loadings_stop), fontface="bold"), color="red") +
  scale_y_continuous(limits = c(-9,8))

##orfindex
ggplot(pca_data, aes(PC1, PC2)) +  geom_point(shape=19, size=.5, alpha=.5, 
                                              aes(colour=z_val.orfindex)) +
       scale_colour_gradientn(colours=c("darkgreen","green","orange", "darkred"), 
                         limits=c(0,3), na.value="darkred") +
                         labs(colour="orfindex", x="PC1", y="PC2")  +
  geom_segment(aes(xend=x1*20, yend=y1*20, x=rep(0,29), y=rep(0,29)), data=loadings, 
               arrow = arrow(length = unit(0.1,"cm")),alpha=0.5, lwd=.5) +
  geom_text(alpha=.9, data=loadings, size=2.8, fontface="bold", 
            aes(x=x1*25, y=y1*23, label=rownames(loadings))) +
  geom_segment(aes(xend=x1*20, yend=y1*20, x=rep(0,3), y=rep(0,3)), data=loadings_stop, 
               arrow = arrow(length = unit(0.1,"cm")),alpha=0.7, lwd=1, color="red") +
  geom_label(alpha=.5, data=loadings_stop, 
             aes(x=x1*30, y=y1*30, label=rownames(loadings_stop), fontface="bold"), color="red") +
  scale_y_continuous(limits = c(-9,8))

##GC content
ggplot(pca_data, aes(PC1, PC2)) +  geom_point(shape=19, size=.5, alpha=.5, aes(colour=z_val.GC)) +
       scale_colour_gradientn(colours=c("darkgreen","green","orange", "darkred"), 
                         limits=c(0.26,0.7), na.value="grey") +
                         labs(colour="GC", x="PC1", y="PC2")  +
  geom_segment(aes(xend=x1*20, yend=y1*20, x=rep(0,29), y=rep(0,29)), data=loadings, 
               arrow = arrow(length = unit(0.1,"cm")),alpha=0.5, lwd=.5) +
  geom_text(alpha=.9, data=loadings, size=2.8, fontface="bold", 
            aes(x=x1*25, y=y1*23, label=rownames(loadings))) +
  geom_segment(aes(xend=x1*20, yend=y1*20, x=rep(0,3), y=rep(0,3)), data=loadings_stop, 
               arrow = arrow(length = unit(0.1,"cm")),alpha=0.7, lwd=1, color="red") +
  geom_label(alpha=.5, data=loadings_stop, aes(x=x1*30, y=y1*30, label=rownames(loadings_stop), + 
                                                 fontface="bold"), color="red") +
  scale_y_continuous(limits = c(-9,8))

#-------------------------------------------------------------------------------

#CCA con matriz de frecuencias de trinucleotidos y matriz de indices y contenido gc
df_zindex <- as.data.frame(zindex)
zcca<-cca(kmer_freqs_uniq~., data=df_zindex)

print(zcca)
spenvcor(zcca)

cca_data <- as.data.frame(cbind(zcca$CCA$wa, df_zindex))
cca_amb <- as.data.frame(zcca$CCA$biplot)
ggplot(cca_data, aes(CCA1*-1, CCA2)) +  geom_point(shape=19, size=.5, alpha=.5, 
                                                   aes(colour=orfindex)) +
  scale_colour_gradientn(colours=c("darkgreen","green","orange", "darkred"), 
                         limits=c(0,3), na.value="grey") +
  labs(colour="orfindex", x="CCA1", y="CCA2") +
  geom_segment(aes(xend=CCA1*-3, yend=CCA2*3, x=rep(0,3), y=rep(0,3)), data=cca_amb, 
               arrow = arrow(length = unit(0.1,"cm")), color="black" , alpha=0.7, lwd=1.5) +
  geom_label(alpha=.5, data=cca_amb, aes( x=CCA1*-3.7, y=CCA2*3.4, label=rownames(cca_amb))) +
  scale_x_continuous(limits=c(-3.5,4)) +
  scale_y_continuous(limits=c(-6,3))

#-------------------------------------------------------------------------------

#K-means utilizando los valores del CCA, le pide que agrupe en 2 clusters
km_cca <- kmeans(cca_data[,1:2], 2)
df_cca_data <- as.data.frame(cbind(cca_data[,1]*-1, cca_data[,2], km_cca$cluster))
colnames(df_cca_data) <-c("CCA1", "CCA2", "cluster")
kcenters <- as.data.frame(km_cca$centers)

#Grafico de clusters encontrados  por k-means y centroides
ggplot(df_cca_data, aes(CCA1,CCA2)) +
  geom_point(shape=19, size=1, alpha=.5,data=df_cca_data[df_cca_data$cluster==1,], 
             colour="darkred") +
  geom_point(shape=19, size=1, alpha=.5,data=df_cca_data[df_cca_data$cluster==2,], 
             colour="darkgreen") +
  geom_point(shape=18, size=4, data=kcenters, aes (x=kcenters[1,1]*-1, y=kcenters[1,2]), 
             colour="red") +
  geom_point(shape=18, size=4, data=kcenters, aes (x=kcenters[2,1]*-1, y=kcenters[2,2]), 
             colour="green") 

#-------------------------------------------------------------------------------

#LDA
##LDA utilizando los indices
#Tomo 5% muestras al azar de los cada uno de los clusters definidos por k-means
o_list <- sample(names(km_cca$cluster[km_cca$cluster==1]),250)
a_list <- sample(names(km_cca$cluster[km_cca$cluster==2]),420)
o_pca = pca_data[o_list,]
a_pca = pca_data[a_list,]
gpo<-"O"
o_freqs<-cbind(o_pca,gpo)
gpo<-"A"
a_freqs<-cbind(a_pca,gpo)
oa_freqs<-as.data.frame(rbind(o_freqs,a_freqs))

#Realizo LDA
df <- subset(oa_freqs, select = c("z_val.stopOE", "z_val.GC", "z_val.orfindex", "gpo"))
zlda<-lda(df$gpo~z_val.stopOE+z_val.GC+z_val.orfindex , data=df, method="moment")

##LDA utilizando los valores del PCA
#Tomo 100 muestras al azar de los cada uno de los clusters definidos por k-means
o_list <- sample(names(km_cca$cluster[km_cca$cluster==1]),100)
a_list <- sample(names(km_cca$cluster[km_cca$cluster==2]),100)
o_pca = pca_data[o_list,]
a_pca = pca_data[a_list,]
gpo<-"O"
o_freqs<-cbind(o_pca,gpo)
gpo<-"A"
a_freqs<-cbind(a_pca,gpo)
oa_freqs<-as.data.frame(rbind(o_freqs,a_freqs))

#Realizo LDA
df <- subset(oa_freqs, select = c("PC1", "PC2", "gpo"))
zlda<-lda(df$gpo~. , data=df, method="moment")


#De acá en adelante es lo mismo para ambos LDA
#Utilizo los resultados del LDA para predecir el grupo del resto de los datos
zpredict <- predict(object=zlda, newdata= pca_data)
clasif <- as.numeric(zpredict$posterior[,"A"] > zpredict$posterior[,"O"])
prob <- 1-zpredict$posterior[,2]
ndf<-cbind(pca$scores[,1], pca$scores[,2],clasif, prob)
colnames(ndf) <- c("PC1", "PC2", "clasif", "prob")
ndf <- as.data.frame(ndf)
zpredict_df <- as.data.frame(zpredict)

#Histograma de los valores que toma la función discriminante en los dos grupos
ggplot(zpredict_df, aes(LD1)) + 
  geom_histogram(data=zpredict_df[o_list,], bins=20, fill="red", alpha=.5) +
  geom_histogram(data=zpredict_df[a_list,], bins=20, fill="green", alpha=.5)

#Gráfico de las observaciones del PCA con el color de la probabilidad 
#de pertenecer a cada cluster
#Resaltados, los puntos usados para el entrenamiento del LDA
ggplot(ndf, aes(PC1, PC2)) +  geom_point(shape=19, size=.7, alpha=.4, aes(colour=prob)) +
  scale_colour_gradientn(colours=c("green3", "yellow","orange","red2")) +
  geom_point(shape=19, size=.7, alpha=1,data=ndf[v_list,], colour="darkred") +
  geom_point(shape=19, size=.7, alpha=1,data=ndf[a_list,], colour="darkgreen")

