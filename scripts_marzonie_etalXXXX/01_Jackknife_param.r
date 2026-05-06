print("**********************************************************************************")
print("Loading packages")
print("**********************************************************************************")

Rlibdir <- "/R/Rlibs"
# dir.create(Rlibdir, recursive = TRUE, showWarnings = FALSE)
.libPaths(Rlibdir)

options(repos = c(CRAN = "https://cloud.r-project.org"))

# install.packages("phytools")
# install.packages("Hmisc")
# install.packages("phytools")

library(ape)
library(phytools)
library(Hmisc)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 3) {
  stop("Usage: Rscript 01_Jackknife_param.r <ref_tree_file> <trees_dir> <output_file>")
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

plot(ref_tree)
nodelabels()
tiplabels()

A <- makeNodeLabel(ref_tree, prefix="")
A$node.label <- prop.clades(ref_tree, trees)
list = prop.clades(ref_tree,trees)
list

sum(list)/length(list)
plot(A,show.node.label=TRUE)

write.tree(A, file=output_file)
