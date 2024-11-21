## VSEARCH to phyloseq
# EVS 2/2024

library(tidyverse)
library(phyloseq)
library(microViz)
library(ggpubr)
library(vegan)
library(EcolUtils)
library(Biostrings)

#### ----- make phyloseq object ----

# get OTU table
otab <- read.table("data/hominid_otutab.txt", sep = "\t", header = TRUE, comment.char = "")
names(otab) <- str_remove(names(otab), "^X")
names(otab) <- str_remove(names(otab), "_S(\\d){1,4}")
otab <- otab %>% column_to_rownames(var = ".OTU.ID")

# get sintax and wrangle into submission
sin <- read.table("data/hominid_sintax50.txt", sep = "\t", header = FALSE, comment.char = "", na.strings = "") %>% dplyr::select(V1, V4) %>% separate_wider_delim(cols = V4, names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), delim = ",", too_few = "align_start") %>% mutate(OTU_ID = str_extract(V1, "OTU_(\\d){1,4}")) %>% mutate(across(!OTU_ID, ~str_remove(.x, "[:alpha:]:"))) %>% column_to_rownames(var = "OTU_ID") %>% dplyr::select(-V1)

# get sample metadata
##### not publicly shared
meta <- readxl::read_xlsx("private/hominid_metadataITS.xlsx")
### add two blank files to the metadata file
bls <- data.frame(SampleID = c("Blank1", "Blank2"),
                  FileID = c("Blank1_009_A09_011", "Blank2_009_A10_011"),
                  SampleName = "Blanks",
                  Group = "Blanks",
                  Captivity_Status = NA,
                  Wonky = "",
                  Dataset = "Blanks")
meta <- meta %>% rbind(bls)

# get reference sequences; big to share on Github
#seqs <- readDNAStringSet("data/hominid_otus.fasta")
#names(seqs) <- str_extract(names(seqs), "OTU_(\\d){1,100000}")

#### make phyloseq
ps <- phyloseq(
  otu_table(otab, taxa_are_rows = TRUE),
  tax_table(sin %>% as.matrix()),
  sample_data(meta %>% filter(!str_detect(FileID, "Unknown")) %>% column_to_rownames(var = "FileID"))
)
ps@refseq <- seqs

## wrangle names
ps <- ps %>% 
  ps_mutate(Species = case_when(
    str_detect(Group, "BaAka") ~ "Human_BaAka",
    str_detect(Group, "Bantu") ~ "Human_Bantu",
    str_detect(Group, "Lowland") ~ "Lowland Gorilla",
    str_detect(Group, "Mountain") ~ "Mountain Gorilla",
    str_detect(Group, "Chimp") ~ "Chimp",
    str_detect(Group, "Blanks") ~ "Blanks",
    str_detect(Group, "Mangabey") ~ "Mangabey"
  ),
  SpeciesCaptive = paste(Captivity_Status, Species, sep = "_")) %>% 
  ps_mutate(SpeciesCaptive = if_else(str_detect(SpeciesCaptive, "Human"), Species, SpeciesCaptive))

### ---- filter and rarefy ----

## are there any not assigned at Kingdom or phylum?
ps1 <- ps %>% microViz::tax_select("Fungi", "Kingdom") 

## remove samples found in blanks
pbl <- ps1 %>% ps_filter(Group %in% "Blanks")
bl <- taxa_names(pbl)
psdecon <- ps1 %>% tax_select(tax_list = bl, deselect = TRUE)

# how many unassigned at phylum level?
nas <- psdecon %>% tax_fix() %>% tax_select("Fungi Kingdom", "Phylum") # 186

### rarefy and filter, then resolve NAs
## rarefy
otab1 <- psdecon %>% otu_get() %>% data.frame()
rarecurve(otab1, step = 1000)
rtab <- otab1[rowSums(otab1) > 1500,]
rare <- rrarefy.perm(rtab, sample = 1500, n = 100, round.out = TRUE)
rarecurve(rare, step = 100)
# re-phyloseq
psrare <- phyloseq(otu_table(rare, taxa_are_rows = FALSE), tax_table(psdecon), sample_data(psdecon), refseq(psdecon))
# remove empty samples and taxa
psrare <- psrare %>% tax_filter(min_total_abundance = 1)
# filter
psfilt <- psrare %>% tax_filter(min_prevalence = 3, min_total_abundance = 1)
psfilt <- subset_samples(psfilt, sample_sums(psfilt) > 0)
sample_names(psfilt) <- psfilt@sam_data$SampleName
sample_names(psfilt) <- str_replace_all(sample_names(psfilt), "-", "_")

## ---- resolve NA phyla ----

# get NAs
nas <- psfilt %>% tax_fix() %>% tax_select("Fungi Kingdom", "Phylum") # 25
taxa_names(nas)

#### FROM HERE: manually resolve NA-assigned phylum using BLAST search against the Unite database and the NCBI nt database

# from this; taxa to remove
torem <- c("OTU_1077", "OTU_1409", "OTU_1438", "OTU_1472", "OTU_148", "OTU_1548", "OTU_1615", "OTU_1811", "OTU_2061", "OTU_2372", "OTU_2473", "OTU_2556", "OTU_2828", "OTU_3018", "OTU_3084", "OTU_3112", "OTU_344", "OTU_651", "OTU_823", "OTU_964", "OTU_966", "OTU_990")

# create dataframe of OTUs to replace
torep <- taxa_names(nas)[!taxa_names(nas) %in% torem]
newnames <- data.frame(OTU = c("OTU_1190", "OTU_13", "OTU_131", "OTU_1647", "OTU_2040", "OTU_2141", "OTU_218", "OTU_384", "OTU_439", "OTU_58", "OTU_663", "OTU_693"),
                       Tax = c("k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Agaricales;f__Psathyrellaceae;g__Coprinellus",
                               "k__Fungi;p__Basidiomycota;c__Malasseziomycetes;o__Malasseziales;f__Malasseziaceae;g__Malassezia", 
                               "k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Trichosphaeriales;f__Trichosphaeriaceae;g__Nigrospora", 
                               "k__Fungi;p__Ascomycota;c__Sordariomycetes",
                               "k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Hypocreales",
                               "k__Fungi;p__Ascomycota;c__Eurotiomycetes;o__Eurotiales;f__Aspergillaceae;g__Aspergillus", 
                               "k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Tremellales",
                               "k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Xylariales;f__Xylariaceae;g__Xylaria",
                               "k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Hypocreales;f__Nectriaceae;g__Fusicolla",
                               "k__Fungi;p__Ascomycota;c__Dothideomycetes;o__Capnodiales",
                               "k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Trechisporales;f__Hydnodontaceae;g__Trechispora",
                               "k__Fungi;p__Ascomycota;c__Eurotiomycetes;o__Eurotiales" )) %>% 
  separate_wider_delim(cols = Tax, delim = ";", too_few = "align_start",
                       names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), cols_remove = TRUE) %>% mutate(Species = NA) %>% 
  mutate(across(!OTU, ~str_remove(.x, "[:alpha:]__")))


## make new phylo
psname <- psfilt %>% tax_select(tax_list = torem, deselect = TRUE) 
newtt <- psname %>% tt_get() %>% as.data.frame() %>% rownames_to_column(var = "OTU") %>% filter(!OTU %in% torep) %>% 
  full_join(newnames) %>% column_to_rownames(var = "OTU") %>% as.matrix()
tax_table(psname) <- newtt

# rename and save
psf <- psname
save(psf, file = "private/hominid_phyloITS.RData") # not shared publicly because it has the accompanying metadata in it
