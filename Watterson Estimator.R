
library(ape)
library(stringr)
library(pegas)


## This function will allow u to find N_poly 

find_difference <- function(sequences) {
  n <- length(sequences)
  seq_length <- length(sequences[[1]])
  N <- 0
  
  for (i in 1:seq_length) {
    base_counts <- table(sapply(sequences, function(seq) seq[i]))
    if (sum(base_counts != n) > 0) {
      N <- N + 1
    }
  }
  return(N)
}


find_difference(list_dna)


## If n is large, greater than 100, then we use a different formula

sum_of_reciprocals <- function(n) {
  sum <- 0
  for (i in 1:(n-1)) {
    sum <- sum + 1/i
  }
  return(sum)
}

### Easy formula for Watterson
Watterson_esimator_calculator <- function(n,L, N_poly){
  K_n = N_poly/L
  if (n <100){
    theta_watterson = K_n/sum_of_reciprocals(n)
  }else{
    theta_watterson= K_n/0.577+ln(n)
  }
  print(paste("This is N_poly: ", N_poly))
  print(paste("This is Kn: ", K_n))
  return(theta_watterson)
}

Watterson_esimator_calculator(4,20,1)




##### Lemmet try to upgrade thte thiug 


## NOw, this function should only have the 3 DNA sequences in it , the length shoudl be calcualted easy!!

dna1 = "CTATCGAGCTGTTTTGATCCCTCCTCCTCCAAAATTGATA"
dna2 = "CTATCGAGCTGTTTTGATCCCTCCTCCTCCAAAATTGATA"
dna3 = "CTATCGAGCTGTTTTGATCCCTCCTCCTCCAAAATTGATG"


list_dna = list(dna1, dna2, dna3)
list_dna = lapply(list_cool, str_split_1, pattern="")


find_difference <- function(sequences) {
  n <- length(sequences)
  seq_length <- length(sequences[[1]])
  N <- 0
  
  for (i in 1:seq_length) {
    base_counts <- table(sapply(sequences, function(seq) seq[i]))
    if (sum(base_counts != n) > 0) {
      N <- N + 1
    }
  }
  return(N)
}



## Advanced forumla for Waterson, where u just have to input DNA sequences as strings, and then make it a list

Watterson_esimator_calculator_2_upgraded<- function(list_of_sequences, u =1*10^(-9)){
  
  list_of_sequences <- lapply(list_of_sequences, str_split_1, pattern="")
  n <- length(list_of_sequences)
  L=length(list_of_sequences[[1]])
  N_poly= find_difference(list_of_sequences)
  K_n = N_poly/L
  
  if (n <100){
    theta_watterson = K_n/sum_of_reciprocals(n)
  }else{
    theta_watterson= K_n/0.577+ln(n)
  }
  
  print(paste("This the number of sequences:  ", n ))
  print(paste("This is the length of 1 sequence: ", L))
  print(paste("This is N_poly: ", N_poly))
  print(paste("This is Kn: ", K_n))
  print(paste("This is ThetaW: ", theta_watterson/(4*u)))
}





### Example of how to use this
# Step 1, type in the DNA that is given 

dna1 = "CTATCGAGCTGTTTTGATCCCTCCTCCTCCAAAATTGATA"
dna2 = "CTATCGAGCTGTTTTGATCCCTCCTCCTCCAAAATTGATA"
dna3 = "CTATCGAGCTGTTTTGATCCCTCCTCCTCCAAAATTGATG"


list_dna = list(dna1, dna2, dna3) # add dna4 or dna5 or.. if requried
## Step 2, use the functtion lol. The function would take 2 parameters, list_dna, and u = mutation rate. If u is not typed in, it is assumed to be 1*10^-9

Watterson_esimator_calculator_2_upgraded(list_dna)






