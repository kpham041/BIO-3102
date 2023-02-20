library(ape)
library(stringr)
library(pegas)




dna1 = "TATGTCGCGTTAAATGGATCGTGATTGTTT"
dna2 = "TATGTCGGGGTAAATCGTTCGTGATTTTGT"
dna3 = "TATGTCGAGATAAATGGGTCGTGATTCTTT"
dna4 = "TATGTCGAGCTAAATTGGTCGTGATTTTGT"
dna5 = "TATGTCGCGGTAAATCGCTCGTGATTATGT"

list_cool = list(dna1, dna2, dna3, dna4, dna5)


find_difference <- function(sequences) {
  list_of_sequences <- lapply(sequences, str_split_1, pattern="")
  n <- length(list_of_sequences)
  seq_length <- length(list_of_sequences[[1]])
  N <- 0
  indices <- c()
  
  for (i in 1:seq_length) {
    base_counts <- table(sapply(list_of_sequences, function(seq) seq[i]))
    if (sum(base_counts > 1) > 1) {
      N <- N + 1
      indices <- c(indices, i)
    }
  }
  print(paste("Position of informative sites is: ", indices))
  return(indices)
}


find_difference(list_cool)







