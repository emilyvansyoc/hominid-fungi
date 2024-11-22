### figure 3


library(tidyverse)
library(ape)
library(vegan)
library(ggtree)
library(ggpubr)
library(ggtext)
library(scales)
library(treeio)


#### ---- panel A: frankenstein a random tree for schematic of codiv ----

source("R/colors.R")

# get a tree
tr <- rtopology(n=24, rooted = TRUE)
plot(tr)
tr1 <- keep.tip(tr, tip = c("t4", "t3", "t24", "t11", "t10", "t8", "t21", "t15", "t2", "t9", "t1", "t18"))
plot(tr1)
tr1$tip.label <- c("*G. beringei*", "*G. beringei*", "*P. troglodyte*", "*P. troglodyte*", "*G. gorilla*", "*G. beringei*", "*G. gorilla*", "*H. sapien*", "*G. gorilla*", "*G. beringei*", "*G. gorilla*", "*H. sapien*")

# need unique tip labels
tr1$tip.label <- paste0("tip_", seq(1:12))

plot(tr1)

## get phylopic of a yeast
flabdf <- data.frame(
  row.names = tr1$tip.label,
  name = tr1$tip.label,
  uid = c("35db2572-a27d-4f83-9b45-11195d6fd5af"),
  forcolor = c("*G. beringei*", "*G. beringei*", "*P. troglodyte*", "*P. troglodyte*", "*G. gorilla*", "*G. beringei*", "*G. gorilla*", "*H. sapien*", "*G. gorilla*", "*G. beringei*", "*G. gorilla*", "*H. sapien*")
)

### plot tree-style aka topo cong pic style
tr1 %>% ggtree(size = 6, branch.length = "none", ladderize = FALSE) %<+% 
  flabdf + 
  geom_tiplab(aes(image = uid, color = forcolor), geom = "phylopic", size = 0.08, offset = -0.2) +
  scale_x_reverse(limits = c(11, 0)) + # artficially smoosh this so it fits better
  # manual colors
  scale_color_manual(values = italic.cols) +
  theme(legend.position = "none",
        plot.margin = margin(l=0)) 


#### plot circle tree
### takes some tweaking to get the correct aspect ratios of our little yeast guys
df <- flabdf
df$svg <- lapply(df$uid, get_phylopic)

# make plot
ggtree(tr1, size = 6, layout = "circular", branch.length = "none") %<+% df +
  geom_phylopic(aes(img = svg, color = forcolor), size = 0.5) +
  scale_color_manual(values = italic.cols) +
  theme(#legend.text = element_markdown(size = 18),
    #legend.title = element_blank(),
    legend.position = "none",
    plot.margin = margin(t=0, r=0, b=0, l=0, unit = "pt"),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA))

## ---- panels B-G and supplementary figure: neighbor joining trees ----

# get trees
load("data/all_otu_subtrees.RData")

# get significant PACo and parafit results
load("data/final_sigs_otus_parafitpaco.RData")

# get colors
source("helpers/colors.R")

### define function to make trees

myPlot <- function(mytree) {
  
  # break down hominid representation
  tipdf <- data.frame(tips = mytree$tip.label) %>% 
    mutate(tip1 = str_remove(tips, "_(\\d){1,4}$")) %>% 
    mutate(hominid = str_extract(tip1, "[:alpha:]{1,12}$")) %>% 
    column_to_rownames(var = "tips")
  tipdf %>% group_by(hominid) %>% count()
  
  # add sci name to tipdf
  tipdf <- tipdf %>% 
    mutate(italname = case_when(
      hominid == "Chimp" ~ "*P. troglodyte*",
      hominid == "Human" ~ "*H. sapien*",
      hominid == "LowGorilla" ~ "*G. gorilla*",
      hominid == "MountGorilla" ~ "*G. beringei*"
    )) 
  
  # simplify tip labels
  plottree <- mytree
  plottree$tip.label <- tipdf$italname
  td <- data.frame(rownames = plottree$tip.label,
                   hominid = plottree$tip.label)
  
  # plot circular
  tplot <- plottree %>% ggtree(ladderize = FALSE, layout = "fan", branch.length = "none") %<+% td + 
    geom_tippoint(aes(color = hominid), size = 3) +
    scale_color_manual(values = italic.cols) +
    theme(#legend.text = element_markdown(size = 18),
      #legend.title = element_blank(),
      legend.position = "none",
      plot.margin = margin(t=0, r=0, b=0, l=0, unit = "pt"),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA))
  
  # return plot
  return(tplot)
  
}

### make trees
# do by hand to get ordering right for now
### 11/13/24- update to keep only top 6 trees, move rest to supplementary
bigsig <- sigps[1:6,]
sigtrees <- alltrees[names(alltrees) %in% bigsig$otu]
sigtrees <- list(alltrees$OTU_169, alltrees$OTU_456, alltrees$OTU_329, alltrees$OTU_5023, alltrees$OTU_2339, alltrees$OTU_86)
names(sigtrees) <- bigsig$otu

# get "little sigs"
lsig <- sigps[7:11,]
lsigtrees <- alltrees[names(alltrees) %in% lsig$otu]
lsigtrees <- list(alltrees$OTU_559, alltrees$OTU_384, alltrees$OTU_282, alltrees$OTU_89, alltrees$OTU_12)

# plot
plist <- lapply(sigtrees, myPlot)
names(plist) <- names(sigtrees)

# save
ggarrange(plotlist = plist, ncol = 3, nrow = 2)

## ---- panel H: fungal tree with cophy taxa ----

### takes some Frankenstein-ization to figure out which nodes to highlight (node labels in ape... ugh)

# get the tree pulled from another paper
ftree <- read.tree("/li2021_fungaltree/1672taxa_290genes_bb_1.treefile")


# get tip labels
ftips <- data.frame(tiplab = ftree$tip.label) %>% 
  mutate(genus = sapply(str_split(tiplab, "_"), `[`, 1))


# get the genera present in our dataset
load("data/all_genera_list.RData")

# subset
fsub <- ftips %>% 
  filter(genus %in% gen)

## individual node labels were added with prefix "mynode"
tr <- read.newick("updated_methods/codiv/li2021_fungaltree/li2021_tree_mynodes.newick") %>% as.treedata()

# make tree data 
td <- as_tibble(tr)

#### get Pleosporales nodes
# BLAST suggest these may belong to the Thyridariaceae family
tab <- readxl::read_xlsx(path = "/li2021_fungaltree/1-s2.0-S0960982221001391-mmc3.xlsx", sheet = "B")
# get Pleosporale
pleo <- tab %>% filter(order == "Pleosporales")
pleosub <- pleo %>% filter(tip_id %in% ftree$tip_id)
pleotree <- keep.tip(ftree, tip = pleo$tip_id)

# make dataframe of node labels by hand
# add labels to nodes
mynodes <- data.frame(
  genus = c("Aureobasidium", "Saturnispora", "Malassezia", 
            "Talaromyces", #"Geotrichum", 
            "Xylaria/Nigrospora", "Pleosporales", "Cladosporium",
            # 11/13 add
            "Lasiodiplodia"),
  mynode = c("node475", # aureobasidium
             "node247", # saturnispora
             "node414", # malassezia
             "node589", # talaro#"node327", 
             "node159", # xylaria/nigrospora
             # change pleosporales node to fit
             #"node449", 
             "node455", # pleo
             "node467", # cladosporium
             #"node448" falls in a weird spot on the tree - shift slightly
             "node463" # lasioplodia
  )) %>%  
  left_join(td, by = c("mynode" = "label"))

# get the default colors
library(scales)
cols <- hue_pal()(8)
names(cols) <- mynodes$genus


### build plot
tr %>% ggtree() +
  geom_hilight(data = mynodes, aes(node = node, fill = genus),extend = 3.4, alpha = 0.9) +#,  type = "gradient", gradient.direction = "tr") +
  xlim(0, 3.4) +
  ggpubr::rotate() +
  theme(legend.position = "none") +
  scale_fill_manual(values = cols)
# save

##### 11/19; only the fungi shown in figure 3
mynodes1 <- mynodes %>% filter(genus %in% c("Aureobasidium", "Xylaria/Nigrospora", "Talaromyces", "Pleosporales", "Saturnispora"))
tr %>% ggtree() +
  geom_hilight(data = mynodes1, aes(node = node, fill = genus),extend = 3.4, alpha = 0.9) +#,  type = "gradient", gradient.direction = "tr") +
  xlim(0, 3.4) +
  ggpubr::rotate() +
  theme(legend.position = "none") +
  scale_fill_manual(values = cols)





