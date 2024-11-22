### final analysis; with mangabey
# EVS 10/2024
library(ape)
library(phyloseq)
library(tidyverse)
library(microViz)

# fungi
load("private/hominid_phyloITS.RData")
source("helpers/fx_myDend.R")

## ---- Bray Curtis ----

d <- myDend(psf, distance = "bray", taxlevel = "Genus")
write.tree(d, file = "data/fung_bray_genus_dendrogram.newick")

# read in tree rooted in Figtree
fung <- read.tree("data/figtree_fung_bray_genus_dendrogram.newick")

# host
hr <- read.tree("data/rooted_hominid.newick")
# rotate
hr.rot <- rotateConstr(hr, constraint = c("Mangabey", "Gorilla_gorilla", "Gorilla_beringei", "Homo_sapien", "P_troglodytes_schweinfurthii"))
# write
write.tree(hr.rot, file = "data/hominid_rotated.newick")

# set path to TreeCMP
TC <- "/bin/TreeCmp_v2.0-b76/bin/"


# run TreeCMP on the random trees versus the host tree 
system(paste0("java -jar ", TC, "/treeCmp.jar -r ", 
              "data/hominid_rotated.newick -i data/random_hominidtrees.newick -d mc rc ns tt mp mt co rcw nsw gdr cow -N -o updated_methods/wmangabey/hominid_v_random_rotated.txt"))


# run on fungal dendrogram
system(paste0("java -jar ", TC, "/treeCmp.jar -r data/hominid_rotated.newick -i data/figtree_fung_bray_genus_dendrogram.newick -d mc rc ns tt mp mt co rcw nsw gdr cow -N -o data/obs_fung_bray_rotated.txt"))

# get results
res <- read.table("data/obs_fung_bray_rotated.txt", sep = "\t", header = TRUE)

# calculate p values
source("helpers/fx_myPval.R")
myPval2(obs = "data/obs_fung_bray_rotated.txt", stoch = "data/hominid_v_random_rotated.txt", normalized = TRUE)

## ---- Jaccard ----


j <- myDend(psf, distance = "jaccard", taxlevel = "Genus")
write.tree(j, file = "data/fung_jaccard_genus_dendrogram.newick")


# run on fungal dendrogram
system(paste0("java -jar ", TC, "/treeCmp.jar -r data/hominid_rotated.newick -i data/figtree_fung_jaccard_genus_dendrogram.newick -d mc rc ns tt mp mt co rcw nsw gdr cow -N -o data/obs_fung_jaccard_rotated.txt"))

# get results
res <- read.table("data/obs_fung_jaccard_rotated.txt", sep = "\t", header = TRUE)

# calculate p values
source("helpers/fx_myPval.R")
myPval2(obs = "data/obs_fung_jaccard_rotated.txt", stoch = "data/hominid_v_random_rotated.txt", normalized = TRUE)

