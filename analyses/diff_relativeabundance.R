## differential relative abundance between humans and non-human primate hominids
# EVS 10/2024


library(microViz)
library(phyloseq)
library(tidyverse)

# get data
load("/private/hominid_ITSphylo.RData")

# aggregate at Genus, remove Mangabey, and prettyfy names
psfilt <- psf %>% 
  tax_fix() %>% 
  ps_filter(!Group == "Mangabey") %>% 
  tax_glom("Genus") 
# rename to scientific name so tip labels match
psfilt <- psfilt %>% 
  ps_mutate(SciName = case_when(
    SpeciesCaptive %in% "Wild_Chimp" ~ "P_troglodytes_schweinfurthii",
    SpeciesCaptive %in% "Wild_Lowland Gorilla" ~ "Gorilla_gorilla",
    SpeciesCaptive %in% "Wild_Mountain Gorilla" ~ "Gorilla_beringei",
    SpeciesCaptive %in% "Human_Bantu" ~ "Homo_sapien",
    SpeciesCaptive %in% "Human_BaAka" ~ "Homo_sapien"
  )) %>% 
  ps_mutate(ItalName = case_when(
    SciName %in% "Gorilla_beringei" ~ "*G. beringei*",
    SciName %in% "Gorilla_gorilla" ~ "*G. gorilla*",
    SciName %in% "Homo_sapien" ~ "*H. sapien*",
    SciName %in% "P_troglodytes_schweinfurthii" ~ "*P. troglodyte*"
  ))
taxa_names(psfilt) <- psfilt@tax_table[,'Genus']

## ---- chi squared test for number of unique genera ----

# convert to presence-absence
pa <- psfilt %>% tax_transform("pa") %>% ps_melt()

# get total numbers for each species
hum <- pa %>% group_by(OTU, SciName) %>% summarize(tot = sum(Abundance)) %>% 
  pivot_wider(names_from = SciName, values_from = tot) %>% 
  mutate(ape.tot = sum(Gorilla_beringei, Gorilla_gorilla, P_troglodytes_schweinfurthii)) %>% 
  mutate(onlyhuman = if_else(Homo_sapien > 1 & ape.tot == 0, 1, 0)) %>% 
  filter(onlyhuman == 1)

# get numbers for a contingency table
totgen <- length(unique(hum$OTU))
# how many total shared taxa?
totgen - 36-23-46-0

# make table (by hand... icky)
ctab <- matrix(c(36, 23, 46, 0, 196, 196, 196, 196), byrow = TRUE, nrow = 2)
rownames(ctab) <- c("Unique", "Non-Unique")
colnames(ctab) <- c("human", "chimp", "lgorilla", "mgorilla")
chisq.test(ctab) # significant

# do pairwise tests
p <- pairwise.prop.test(t(ctab))

## ---- compare across core prevalent taxa ----

# convert to presence-absnece
pa <- psfilt %>% tax_transform("pa") %>% ps_melt()

# calculate prevalence
prev <- pa %>% group_by(OTU, SciName) %>% summarize(tot = sum(Abundance)) %>% ungroup() %>% 
  mutate(totprev = case_when(
    SciName == "Gorilla_beringei" ~ tot / stot$n[stot$SciName == "Gorilla_beringei"],
    SciName == "Gorilla_gorilla" ~ tot / stot$n[stot$SciName == "Gorilla_gorilla"],
    SciName == "P_troglodytes_schweinfurthii" ~ tot / stot$n[stot$SciName == "P_troglodytes_schweinfurthii"],
    SciName == "Homo_sapien" ~ tot / stot$n[stot$SciName == "Homo_sapien"]
  )) %>% dplyr::select(-tot) %>% 
  filter(totprev > 0.1) %>% group_by(OTU) %>% count() %>% filter(n == 4)

# get a list of the prevalent genera
prevgenera <- prev$OTU

## compare across these prevalent genera
prevps <- psfilt %>% tax_select(prevgenera, "Genus", strict_matches = TRUE) %>% 
  # remove Sacch order (not really a genus)
  tax_select("Saccharomycetales Order", ranks_searched = "Genus", strict_matches = TRUE, deselect = TRUE )

# make three separate tests for each host species
hc1 <- taxatree_models2stats(
  prevps %>% ps_filter(SciName %in% c("Homo_sapien", "P_troglodytes_schweinfurthii")) %>% 
    tax_fix() %>% 
    tax_transform("compositional", rank = "Genus") %>% 
    #tax_filter(min_prevalence = 0.1, use_counts = TRUE) %>% 
    taxatree_models(
      type = lm,
      trans = "log2", trans_args = list(zero_replace = "halfmin"),
      ranks = "Genus",
      variables = "ishuman"
    )) %>% taxatree_stats_get()
hg1 <- taxatree_models2stats(
  prevps %>% ps_filter(SciName %in% c("Homo_sapien", "Gorilla_gorilla")) %>% 
    tax_fix() %>% 
    tax_transform("compositional", rank = "Genus") %>% 
    #tax_filter(min_prevalence = 0.1, use_counts = TRUE) %>% 
    taxatree_models(
      type = lm,
      trans = "log2", trans_args = list(zero_replace = "halfmin"),
      ranks = "Genus",
      variables = "ishuman"
    )) %>% taxatree_stats_get()  
hm1 <- taxatree_models2stats(
  prevps %>% ps_filter(SciName %in% c("Homo_sapien", "Gorilla_beringei")) %>% 
    tax_fix() %>% 
    tax_transform("compositional", rank = "Genus") %>% 
    #tax_filter(min_prevalence = 0.1, use_counts = TRUE) %>% 
    taxatree_models(
      type = lm,
      trans = "log2", trans_args = list(zero_replace = "halfmin"),
      ranks = "Genus",
      variables = "ishuman"
    )) %>% taxatree_stats_get()

# combine
all1 <- hc1 %>% mutate(comp = "P_troglodytes") %>% 
  rbind(hg1 %>% mutate(comp = "G_gorilla")) %>% 
  rbind(hm1 %>% mutate(comp = "G_beringei"))

# adjust for multiple comparisons
all1 <- all1 %>% mutate(padj = p.adjust(p.value, method = "fdr"))%>% 
  mutate(combo = paste(taxon, comp))

# get sigs
sig1 <- all1 %>% filter(padj < 0.05)
# get sig combinations
sigs <- sig1 %>% 
  dplyr::select(OTU = taxon, SciName = comp) %>% 
  mutate(combo = paste(OTU, SciName))

# calculate log2 change
lg <- prevps %>% 
  tax_transform("compositional") %>% 
  tax_transform("log2", zero_replace = "halfmin", chain = TRUE) %>% 
  ps_melt()

# get the means
lgm <- lg %>% group_by(OTU, SciName) %>% summarize(avg = mean(Abundance)) %>% ungroup() %>% 
  pivot_wider(names_from = SciName, values_from = avg) %>% 
  # to calculate log2 fold change, want log2(A) - log2(B)
  mutate(GB = Homo_sapien - Gorilla_beringei,
         GG = Homo_sapien - Gorilla_gorilla,
         PT =  Homo_sapien - P_troglodytes_schweinfurthii) %>% 
  dplyr::select(OTU, G_beringei = GB, G_gorilla = GG, P_troglodytes = PT) %>% 
  pivot_longer(!OTU, names_to = "SciName", values_to = "fc") %>% 
  mutate(combo = paste(OTU, SciName)) %>% 
  left_join(all1, by = "combo") %>% 
  mutate(issig = if_else(padj < 0.05, TRUE, FALSE)) %>% 
  mutate(ItalName = case_when(
    SciName %in% "G_gorilla" ~ "*G. gorilla*",
    SciName %in% "G_beringei" ~ "*G. beringei*",
    SciName %in% "P_troglodytes" ~ "*P. troglodyte*"
  )) %>% 
  mutate(SciName = factor(SciName, ordered = TRUE, levels = c("P_troglodytes", "G_gorilla", "G_beringei")))

# write out for summary statistics (supplementary table)
sums <- lgm %>% 
  dplyr::select(SciName, Genus = OTU, fc, estimate, std.error, statistic, padj) %>% 
  arrange(SciName, Genus)

### ---- compare across the less prevalent genera ----


## remove the 10 core genera and compare everything else
noprevps <- psfilt %>% tax_select(prevgenera, "Genus", strict_matches = TRUE, deselect = TRUE) 

# make three separate comparisons for each host species
hc2 <- taxatree_models2stats(
  noprevps %>% ps_filter(SciName %in% c("Homo_sapien", "P_troglodytes_schweinfurthii")) %>% 
    tax_fix() %>% 
    tax_transform("compositional", rank = "Genus") %>% 
    #tax_filter(min_prevalence = 0.1, use_counts = TRUE) %>% 
    taxatree_models(
      type = lm,
      trans = "log2", trans_args = list(zero_replace = "halfmin"),
      ranks = "Genus",
      variables = "ishuman"
    )) %>% taxatree_stats_get()
hg2 <- taxatree_models2stats(
  noprevps %>% ps_filter(SciName %in% c("Homo_sapien", "Gorilla_gorilla")) %>% 
    tax_fix() %>% 
    tax_transform("compositional", rank = "Genus") %>% 
    #tax_filter(min_prevalence = 0.1, use_counts = TRUE) %>% 
    taxatree_models(
      type = lm,
      trans = "log2", trans_args = list(zero_replace = "halfmin"),
      ranks = "Genus",
      variables = "ishuman"
    )) %>% taxatree_stats_get()  
hm2 <- taxatree_models2stats(
  noprevps %>% ps_filter(SciName %in% c("Homo_sapien", "Gorilla_beringei")) %>% 
    tax_fix() %>% 
    tax_transform("compositional", rank = "Genus") %>% 
    #tax_filter(min_prevalence = 0.1, use_counts = TRUE) %>% 
    taxatree_models(
      type = lm,
      trans = "log2", trans_args = list(zero_replace = "halfmin"),
      ranks = "Genus",
      variables = "ishuman"
    )) %>% taxatree_stats_get()

# combine
all2 <- hc2 %>% mutate(comp = "P_troglodytes") %>% 
  rbind(hg2 %>% mutate(comp = "G_gorilla")) %>% 
  rbind(hm2 %>% mutate(comp = "G_beringei"))

# adjust for multiple comparisons
all2 <- all2 %>% mutate(padj = p.adjust(p.value, method = "fdr"))%>% 
  mutate(combo = paste(taxon, comp))

# get sigs
sig2 <- all2 %>% filter(padj < 0.05)
# get sig combinations
sigs2 <- sig2 %>% 
  dplyr::select(OTU = taxon, SciName = comp) %>% 
  mutate(combo = paste(OTU, SciName))

sigs2 %>% group_by(OTU) %>% count() %>% arrange(desc(n))


# calculate log2 change
lg2 <- noprevps %>% 
  tax_transform("compositional") %>% 
  tax_transform("log2", zero_replace = "halfmin", chain = TRUE) %>% 
  ps_melt()


#### get means of all taxa for summary stats
sums2 <- lg2 %>% group_by(OTU, SciName) %>% summarize(avg = mean(Abundance)) %>% ungroup() %>% 
  pivot_wider(names_from = SciName, values_from = avg) %>% 
  # to calculate log2 fold change, want log2(A) - log2(B)
  mutate(GB = Homo_sapien - Gorilla_beringei,
         GG = Homo_sapien - Gorilla_gorilla,
         PT =  Homo_sapien - P_troglodytes_schweinfurthii) %>% 
  dplyr::select(OTU, G_beringei = GB, G_gorilla = GG, P_troglodytes = PT) %>% 
  pivot_longer(!OTU, names_to = "SciName", values_to = "fc") %>% 
  mutate(combo = paste(OTU, SciName)) %>% 
  left_join(all2, by = "combo") %>% 
  mutate(issig = if_else(padj < 0.05, TRUE, FALSE)) %>% 
  dplyr::select(SciName, Genus = OTU, fc, estimate, std.error, statistic, padj) %>% 
  arrange(SciName, Genus)

