# Sequence divergence analysis and figure - final
# EVA 11/2024

library(multcompView)
library(ape)
library(tidyverse)
library(usedist)
library(ggpubr)
library(ggtext)
library(cowplot)
library(vegan)
#library(phyloseq)
#library(microViz)
library(cowplot)

# get sequence divergence function

source("helpers/fx_MyDiv.R")

# get parafit/PACo results for significant fungi
load("data/all_genera_list.RData")

# get all genera and OTUs
allgen <- read.table("data/all_pacoparafit_results_wr2.txt", sep = "\t", header = TRUE)

## ---- divergence: raw nucleotide differences ----


# loop through and get divergence estimates in a big dataframe
myotus <- unique(allgen$otu)
alldf <- data.frame()
for(i in 1:length(myotus)) {
  
  # get genus name
  sub <- allgen %>% filter(otu == myotus[i])
  genname <- unique(sub$labname)
  # loop
  out <- myDiv(alignment = read.dna(paste0("data/unique_clipkit_derep_haplotypes_", genname, ".fasta"), format = "fasta"),
               OTU = myotus[i],
               dist.type = "N")
  
  # add to dataframe
  alldf <- rbind(alldf, out[[1]] %>% mutate(otu = myotus[i], Genus = genname))
}

# get summary stats for total sequence divergence across all OTUs
alldf %>% group_by(otu, Genus) %>% get_summary_stats(Distance, type = "common") %>% arrange(desc(mean)) %>% print(n=45)

## ---- divergence: proportions ----


propdf <- data.frame()
for(i in 1:length(myotus)) {
  
  # get genus name
  sub <- allgen %>% filter(otu == myotus[i])
  genname <- unique(sub$labname)
  # loop
  out <- myDiv(alignment = read.dna(paste0("data/unique_clipkit_derep_haplotypes_", genname, ".fasta"), format = "fasta"),
               OTU = myotus[i],
               dist.type = "raw")
  
  # add to dataframe
  propdf <- rbind(propdf, out[[1]] %>% mutate(otu = myotus[i], Genus = genname))
}

## ---- use proportional divergence to calibrate clock ----


# start with averages and use confidence intervals 
avgdf2 <- propdf %>% 
  group_by(Label, otu, Genus) %>% 
  get_summary_stats(Distance, type = "mean_sd") %>% 
  left_join(allgen, by = "otu") %>% drop_na(pacor2)

# using the average between HS and PT to 'calibrate' for each OTU
hs.pt1 <- avgdf2 %>% 
  ungroup() %>% 
  filter(Label %in% c("Between H. sapiens and P. troglodytes")) %>% 
  # assume that host speciation is 6MYA
  mutate(hspt.div = 6) %>% 
  mutate(hspt.clock = mean / hspt.div) %>% # mutation rate every 1MY 
  dplyr::select(otu, Genus, hspt.clock)
summary(hs.pt1$hspt.clock) # mean 0.2% every MYA

# add this back to the dataset and use it to calculate the other 'clocks'
clockdf2 <- propdf %>% 
  ungroup() %>% 
  left_join(allgen, by = "otu") %>% drop_na(pacor2) %>% 
  filter(issig == TRUE) %>% 
  filter(!Label %in% c("Between H. sapiens and P. troglodytes")) %>% 
  left_join(hs.pt1) %>% 
  # calculate the divergence time of other speciations
  mutate(divtime = Distance / hspt.clock)

# get confidence intervals around the estimates
clockest <- clockdf2 %>% 
  group_by(otu, Genus, Species, Label) %>% 
  get_summary_stats(divtime, type = "full") %>% 
  mutate(cihi = mean + ci,
         cilo = mean - ci)

### wrangle and save for summary statistics
ssum <- clockest %>% 
  # make column of pretty-fied names that match codiv fig
  mutate(taxname = str_replace(Species, " Genus", " spp."),
         taxname = str_replace(taxname, " Order", " spp."),
         taxname = str_replace(taxname, "_", " ")) %>% 
  mutate(taxname = if_else( otu == "OTU_456", "Xylaria spp. (1)", taxname),
         taxname = if_else(otu == "OTU_384", "Xylaria spp. (2)", taxname))  %>% 
  dplyr::select(Label, otu, taxname, mean, ci)


### ---- save image ----

# to avoid re-running the loop, save as RData image
#save.image(file = "data/image_sequencedivergence.RData") ### not saved; big
