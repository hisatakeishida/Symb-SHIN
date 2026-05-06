print("**********************************************************************************")
print("Loading packages")
print("**********************************************************************************")

Rlibdir <- "/scratch/project_mnt/S0026/ishida/0_scripts/R/Rlibs"
# dir.create(Rlibdir, recursive = TRUE, showWarnings = FALSE)
.libPaths(Rlibdir)

options(repos = c(CRAN = "https://cloud.r-project.org"))

# install.packages("phangorn")
# remotes::install_github("liamrevell/phytools")

library(readr)
library(ape)
library(phangorn)
library(tidyverse)
library(phytools)
?phytools

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 3) {
  stop("Usage: Rscript 01_avg_tree.R <ref_tree_file> <trees_dir> <output_file>")
}

ref_tree_file <- args[1]
trees_dir     <- args[2]
output_file   <- args[3]

# Read reference tree
ref_tree <- read.tree(ref_tree_file)
ref_tree <- makeNodeLabel(ref_tree, prefix = "")

# Read all trees from directory
list_files <- list.files(trees_dir, full.names = TRUE)
list_files

trees <- list()
for (i in list_files){
  i <- read.tree(i)
  trees[[length(trees)+1]] <- i
}
trees

class(trees)<-"multiPhylo"

rf.tree<-averageTree(trees,method="symmetric.difference")
# qpd.tree<-averageTree(trees,method="quadratic.path.difference")

plot(rf.tree)
nodelabels()
tiplabels()

A <- makeNodeLabel(rf.tree, prefix="")
A$node.label <- prop.clades(rf.tree, trees)
list = prop.clades(rf.tree,trees)
list

sum(list)/length(list)
plot(A,show.node.label=TRUE)

write.tree(A, file=output_file)
