
###### needed package : "ape", "Biostrings", "tidyverse", "phangorn" 


library(tidyverse)
library(phangorn)
library(ape)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")

library(Biostrings)






######## Function started

K80_matrix_and_P_and_Q <- function(seq1, seq2, seq3) {
  
  # Convert input sequences to DNAStringSet object
  dna <- DNAStringSet(c(seq1, seq2, seq3))
  
  # Convert DNAStringSet object to DNAbin object
  x <- as.DNAbin(dna, "DNAStringset")
  
  # Calculate pairwise K80 distances using dist.dna function from ape package
  k80_matrix <- dist.dna(x, model = "K80", variance = FALSE,
                       gamma = FALSE, pairwise.deletion = FALSE,
                       base.freq = NULL, as.matrix = FALSE)
  
  # Function that Calculate P and Q
  titv<-function(dat){
    mat<-as.matrix(dat)
    res<-matrix(NA, ncol=dim(mat)[1], nrow=dim(mat)[1], dimnames=list(x=names(dat), y=names(dat)))
    for(i in 1:dim(mat)[1]){
      for(j in 1:dim(mat)[1]){
        vec<-as.numeric(mat[i,])+as.numeric(mat[j,])-8
        res[i,j]<-length(grep("200|56",vec)) #Transitions
        res[j,i]<-length(grep("152|168|88|104",vec)) #Transversions
      }
    }
    res
  }
  
  ti<-titv(x)
  tv <- t(ti)
  
  transition <-ti[lower.tri(ti)] #number of transition
  transversion <-tv[lower.tri(tv)] #Number of transversions
  
  print(paste("P for S1_2 then S1_3 then then S2_3 is:", transition/nchar(seq1)))
  print(paste("Q for S1_2 then S1_3 then S2_3 is:", transversion/nchar(seq1)))
  tree <-upgma(lord)
  plot(tree)
  return(k80_matrix)
}

#### Test function 

seq1 <- "TATATGGCGTTAAAGGGATAGAGGTTGTTT"
seq2 <- "TATGCCGGGGTAGCTCGTTGCTCATTTTGT"
seq3 <- "TATGTAGAGATAATAGGGTCTAGTTTCTTT"


K80_matrix_and_P_and_Q(seq1, seq2, seq3)







