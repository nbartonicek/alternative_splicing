library(edgeR)
library(Rsubread)
library(org.Mm.eg.db)
#first build the index
buildindex(basename="mm39",reference="../annotation/Mus_musculus.GRCm39.dna.toplevel.fa.gz")

