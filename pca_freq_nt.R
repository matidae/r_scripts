require(Biostrings)
require (ggplot2)

#Define kmer size and file name
kmer = 3
file_name = '/home/matias/projects/r_scripts/mt1.all.fasta'

#Function that determines the frequency of kmer in a collection of strings (i.e. multifasta)
calc_freqs <- function(file_name, kmer) {
  #Load fasta file
  fasta <- readDNAStringSet(file_name, "fasta")
  #Count kmer and mirror sequences to avoid strand bias
  count_plus <- oligonucleotideFrequency(fasta, kmer)
  count_minus <- oligonucleotideFrequency(reverseComplement(fasta), kmer)
  #Calculate kmers' frequencies as dataframe.
  kmer_freqs <- as.data.frame(cbind(prop.table(count_plus + count_minus, 1)))
  #Get seq names
  seq_names <- names(fasta)
  #Add sequence names
  row.names(kmer_freqs) <- seq_names
  return(kmer_freqs)
}

kmer_freqs <- calc_freqs(file_name, kmer)
#Calculate PCA 
pca <- prcomp(kmer_freqs, scale = TRUE, center = TRUE)
#Load z values (i.e. sequence's properties)
z_value<-read.table(file='/home/matias/projects/r_scripts/gcindex', header=FALSE)$V2
#Create dataframe with the scores of the two first components
pca_data <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], z_val = z_value)
#Plot the data
ggplot(data = pca_data, mapping = aes(PC1, PC2)) + 
  geom_point(shape=19, size=1.5, aes(colour=z_val)) +
  scale_colour_gradientn(colours=c("darkgreen","green","orange", "darkred"), limits=c(min(z_value),max(z_value))) +
  labs(colour=" z val", x="PC1", y="PC2", title="PCA spp") +
  theme(plot.title=element_text(face="bold"))


