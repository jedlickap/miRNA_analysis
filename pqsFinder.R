#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("You should provide two arguments: 1. input fasta file and 2. output file name.", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "pqsOut.txt"
}

library(pqsfinder)
multFa <- readDNAStringSet(args[1])
pqsList <- lapply(multFa, pqsfinder)
capture.output(pqsList, file=args[2])
