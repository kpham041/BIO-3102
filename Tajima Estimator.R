dna1 ="CTATCGAGCTGTTTTGATCCCTCCTCCTCCAAAATTGATA"
dna2 ="CTATCGAGCTGTTTTGATCCCTCCTCCTCCAAAATTGATA"
dna3 ="CTATCGAGCTGTTTTGATCCCTCCTCCTCCAAAATTGATG"

list_cool = list(dna1, dna2, dna3)

########This function counts the pairwise N_polymorphic sites, return a vector with all the values!!



dna_difference_count <- function(sequences) {
  n <- length(sequences)
  count_vector <- numeric(n*(n-1)/2)
  index <- 1
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      seq1 <- sequences[[i]]
      seq2 <- sequences[[j]]
      count <- 0
      for (k in 1:length(seq1)) {
        if (seq1[k] != seq2[k]) {
          count <- count + 1
        }
      }
      count_vector[index] <- count
      index <- index + 1
    }
  }
  return(count_vector)
}

# Write a function for Tajima, take in DNA sequences as parameter

Tajima_estimator <- function(list_of_sequences){
  list_can_be <-lapply(list_of_sequences, str_split_1, pattern="")
  N_poly_pairwise <- dna_difference_count(list_can_be) 
  n <-length(list_can_be )
  L <- length(list_can_be [[1]])
  prob_pairwise <- N_poly_pairwise/L
  denominator<- n*(n-1)/2
  
  Theta_Tajima <- sum(prob_pairwise)/ denominator
  
  print(N_poly_pairwise)
  print(paste("This is ThetaT: ", Theta_Tajima))
  
}

## Testing the function out


Tajima_estimator(list_cool)

