require(Biostrings)

kmer = 3
base = '/home/matias/assemblies/'
fa_list = c('file2.fa','file1.fa')

calc_freqs <- function(base, fa_list, kmer) {
  freq_all=NULL
  for (i in fa_list){
    orfs <- readDNAStringSet(paste(base,"/",i,sep=""), "fasta")
    count_plus <- oligonucleotideFrequency(orfs, kmer)
    count_minus <- oligonucleotideFrequency(reverseComplement(orfs), kmer)
    seq_names <- names(orfs)   
    if (!is.null(freq_all)) {      
      freq_tmp <- as.data.frame(cbind(prop.table(count_plus + count_minus, 1)))
      row.names(freq_tmp)<-seq_names
      freq_all <- rbind(freq_all,freq_tmp)
    } 
    else {
      freq_all <- as.data.frame(cbind(prop.table(count_plus + count_minus, 1)), grp=i)
      row.names(freq_all)<-seq_names
    }    
  }
  return(freq_all)
}
