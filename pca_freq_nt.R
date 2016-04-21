require(Biostrings)

#Define kmer size and file name
kmer = 3
file_name = '/home/matias/projects/r_scripts/m.all.fasta'

#Function that determines the frequency of kmer in a collection of strings (aka multifasta)
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
#Create dataframe with the scores of the two first components
pca_data <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2])
#Plot the data
ggplot(pca_data, mapping = aes(PC1, PC2)) + geom_point()

