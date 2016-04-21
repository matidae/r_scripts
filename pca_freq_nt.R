require(Biostrings)

#Define kmer size and file name
kmer = 3
file_name = '/home/matias/projects/r_scripts/m.all.fasta'

#Function that determines the frequency of kmer in a given string
calc_freqs <- function(file_name, kmer) {
  #Load fasta file
  fasta <- readDNAStringSet(file_name, "fasta")
  #Mirror sequences to avoid strand bias
  count_plus <- oligonucleotideFrequency(orfs, kmer)
  count_minus <- oligonucleotideFrequency(reverseComplement(orfs), kmer)

