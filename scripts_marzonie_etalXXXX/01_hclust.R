print("**********************************************************************************")
print("Loading packages")
print("**********************************************************************************")

Rlibdir <- "/scratch/project_mnt/S0026/ishida/0_scripts/R/Rlibs"
# dir.create(Rlibdir, recursive = TRUE, showWarnings = FALSE)
.libPaths(Rlibdir)

options(repos = c(CRAN = "https://cloud.r-project.org"))

library(readr)
library(ape)
library(phangorn)
library(tidyverse)
library(phytools)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    stop("Usage: Rscript 01_hclust.R <distance_matrix_file> <output_tree_file>")
}

print("**********************************************************************************")
print("Loading data and generate nj tree")
print("**********************************************************************************")

input_file <- args[1]
output_file <- args[2]

# Read distance matrix
d2s_table <- read_table(input_file, col_names = F) %>%  column_to_rownames("X1") %>% as.matrix()
colnames(d2s_table) <-rownames(d2s_table)
d2s_dist <- as.dist(d2s_table)

d2s_dist_hclust <- unroot(as.phylo(hclust(d2s_dist, method = "ward.D")))
# d2s_dist_hclust <- as.phylo(hclust(d2s_dist, method = "ward.D"))

# Write tree in Newick format
ape::write.tree(d2s_dist_hclust, file=output_file)

cat(paste("Hclust tree saved to", output_file, "\n"))



