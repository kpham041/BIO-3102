
library(ape) # will be used to calculate distance
library(stringr) ## use to create the x matrix ACGT, with A,C,G,T split out 

dna1 <- "GAGACTACGACTAGAGCTAGACGGTACAC"
dna2 <- "GAGGCCACTACCCGAGTTGGACAGAACAC"

x=rbind(str_split_1("GAGACTACGACTAGAGCTAGACGGTACAC",""),str_split_1("GAGGCCACTACCCGAGTTGGACAGAACAC",""))
x <-as.DNAbin(x)

dist.dna(x, model = "TN93", variance = FALSE,
         gamma = FALSE, pairwise.deletion = FALSE,
         base.freq = NULL, as.matrix = FALSE)

## We can write a function

## THis is the function that will be used to calculate the pairwise distance

## This work for 2 sequences, but we could add more


michael_pairwise_distance <- function (seq1, seq2, type, t=1){
  x=as.DNAbin(rbind(str_split_1(seq1,""),str_split_1(seq2,"")))
  y=dist.dna(x, model = type, variance = FALSE,
           gamma = FALSE, pairwise.deletion = FALSE,
           base.freq = NULL, as.matrix = FALSE)
  r= y/(2*t)
  
  print(paste("This is the substituion rate base on ", type, ":", "r=", r))
  print(paste("This is the distance base on ", type, ":", "D=", y))
  
  
}
 
## Example in powerpoints, slide 33
dna1 <- "GAGACTACGACTAGAGCTAGACGGTACAC"
dna2 <- "GAGGCCACTACCCGAGTTGGACAGAACAC"

michael_pairwise_distance(dna1, dna2, "JC69", t=2.3)

michael_pairwise_distance(dna1, dna2, "K80")

michael_pairwise_distance(dna1, dna2, "TN93")





### What's next ? 


## WE need to calculate P, Q ## NAh bro, I would just do this on paper


## Substituion rate based on each model, I need to find the sub

## Bro, I am gonna fall alsleep so may times today dude!! I am gonn aeraking die of sleeping dude!!! what the hell bro!! The firs tthing u wak up up sis that you have to check your honour thesisi stuff  udd



### SO we have Matrix M 




# library(ape)
# 
# # calculate the Jukes-Cantor rate matrix
# jc69.rate <- function(pi, k) {
#   pi <- pi[pi != 0]
#   pi <- pi / sum(pi)
#   rate <- -k * (diag(3) - matrix(pi, 3, 3))
#   return(rate)
# }
# 
# # example: k = 1
# k <- 1
# pi <- c(1,0,0,0)
# rate <- jc69.rate(pi, k)
# rate

#dsfdf This is very good!! The



#########


### N POLY IS THE NUMB

