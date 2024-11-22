## pairwise B-C distances and variance in evenness
# EVS 12/2023-11/2024


library(vegan)
library(phyloseq)
library(microViz)
library(tidyverse)
library(ggpubr)
# to get pairwise distances into a nice dataframe:
#devtools::install_github("kylebittinger/usedist")
library(usedist)
library(rstatix)
library(pairwiseAdonis)
library(car)

# get phyloseq
load("/private/hominid_phyloITS.RData")

# aggregate at Genus and rename
psfilt <- psf %>% 
  ps_filter(!Group == "Mangabey") %>% 
  tax_fix() %>% 
  tax_glom("Genus")
# rename to scientific name so tip labels match
psfilt <- psfilt %>% 
  ps_mutate(SciName = case_when(
    SpeciesCaptive %in% "Wild_Chimp" ~ "P. troglodytes",
    SpeciesCaptive %in% "Wild_Lowland Gorilla" ~ "G. gorilla",
    SpeciesCaptive %in% "Wild_Mountain Gorilla" ~ "G. beringei",
    SpeciesCaptive %in% "Human_Bantu" ~ "H. sapiens",
    SpeciesCaptive %in% "Human_BaAka" ~ "H. sapiens"
  ))


# calculate bray
uf <- psfilt %>% dist_calc("bray") %>% microViz::dist_get()

# use 'usedist' package to convert to nice dataframe
ufdf <- dist_groups(d = uf, g = psfilt@sam_data$SciName)

### ---- permanova for within and between groups ----

# get sample data
sampdf <- psfilt %>% samdat_tbl()
# compare all pairwise
pairwise.adonis2(uf ~ SciName, data = sampdf) # all
# double check global statistic
mod <- adonis2(uf ~ SciName, data = sampdf)
mod

# get beta dispersion
bd <- betadisper(uf, group = sampdf$SciName)
anova(bd)
TukeyHSD(bd, which = "group")

### ---- anova ----

## pretty-fy names
forplot <- ufdf %>% filter(!str_detect(Label, "Within")) %>% 
  mutate(labplot = case_when(
    Label == "Between G. gorilla and P. troglodytes" ~ "Between *P. troglodytes* and *G. gorilla*",
    Label == "Between H. sapiens and P. troglodytes" ~ "Between *H. sapiens* and *P. troglodytes*",
    Label == "Between G. beringei and P. troglodytes" ~ "Between *P. troglodytes* and *G. beringei*",
    Label == "Between G. gorilla and H. sapiens" ~ "Between *H. sapiens* and *G. gorilla*",
    Label == "Between G. beringei and G. gorilla" ~ "Between *G. gorilla* and *G. beringei*",
    Label == "Between G. beringei and H. sapiens" ~ "Between *H. sapiens* and *G. beringei*"
  ))

# test with ANOVA
mod <- aov(Distance ~ labplot, data = forplot)
t <- TukeyHSD(mod)

## ---- Shannons div and Levene test ----

# get adiv
adiv <- estimate_richness(psfilt, measures = c("Shannon")) %>% 
  rownames_to_column(var = ".sample_name") %>% 
  full_join(psfilt %>% samdat_tbl())

# get global statistic for Levene test
lt <- leveneTest(Shannon ~ SciName, data = adiv, center = median) # significant
TukeyHSD(lt)

# calculate medians (same as Levene test, just done by hand)
meds <- aggregate(Shannon ~ SciName, data = adiv, FUN = mean)
colnames(meds)[2] <- "sh_mean"
adiv1 <- adiv %>% full_join(meds)
adiv1$res <- abs(adiv1$Shannon - adiv1$sh_mean)

# double check identical test statistic
res.aov <- aov(res ~ SciName, data = adiv1) # same p, F, and DF as leveneTest
TukeyHSD(res.aov) 