### figure 1

library(ggpubr)
library(ape)
library(tidyverse)
library(ggtree)
library(gridExtra)
library(rphylopic)
library(ggimage)
library(usedist)
library(cowplot)
library(rstatix)
library(ggtext) # for element markdown
library(multcompView) # for letters in boxplot
library(ggh4x) # for colored strip backgrounds



# get colors
source("helpers/colors.R")

## ---- panel A: topological congruency ----

## load hominid tree
hom <- read.tree("data/hominid_rotated.newick")
host <- hom

### read in fungal tree
figtree <- read.tree("data/figtree_fung_bray_genus_dendrogram.newick")

### add mangabey to colors
mang <- "#000000"
names(mang) <- "Mangabey"
manga.cols <- c(unders.cols, mang)

# get images
h.img <- get_phylopic(uuid = "036b96de-4bca-408e-adb6-2154fcd724ef", preview = TRUE)
g.img <- get_phylopic(uuid = "142e0571-3b5f-443d-a887-b572a224ea22", preview = TRUE)
c.img <- get_phylopic(uuid = "7133ab33-cc79-4d7c-9656-48717359abb4", preview = TRUE)
# 5/8/2024; get phylopic for fungi
f.img <- get_phylopic(uuid = "aaecd181-feb8-4203-8c64-f46384257e59", preview = TRUE) # hyphae
f.img1 <- get_phylopic(uuid = "e602729e-044f-4d2b-bab2-64a87a0b48c7", preview = TRUE) # yeast
# mangabey
m.img <- get_phylopic(uuid = "eccbb404-c99f-41f9-8785-01a7f57f1269", preview = TRUE)

#  build dataframe for fungi 
flabdf <- data.frame(
  row.names = host$tip.label,
  name = host$tip.label,
  hyphae = c("aaecd181-feb8-4203-8c64-f46384257e59"),
  yeast = c("e602729e-044f-4d2b-bab2-64a87a0b48c7"),
  forcolor = host$tip.label
)

# build dataframe for host
labdf <- data.frame(
  row.names = host$tip.label,
  name = host$tip.label,
  # uid = c("142e0571-3b5f-443d-a887-b572a224ea22", "142e0571-3b5f-443d-a887-b572a224ea22", ## GORILLA IS TWICE
  #     "036b96de-4bca-408e-adb6-2154fcd724ef", "7133ab33-cc79-4d7c-9656-48717359abb4"),
  forcolor = host$tip.label
) %>% 
  mutate(uid = case_when(
    name %in% "Gorilla_beringei" ~ "142e0571-3b5f-443d-a887-b572a224ea22",
    name %in% "Gorilla_gorilla" ~ "142e0571-3b5f-443d-a887-b572a224ea22",
    name %in% "Homo_sapien" ~ "036b96de-4bca-408e-adb6-2154fcd724ef",
    name %in% "P_troglodytes_schweinfurthii" ~ "7133ab33-cc79-4d7c-9656-48717359abb4",
    name %in% "Mangabey" ~ "eccbb404-c99f-41f9-8785-01a7f57f1269"
  ))


### make host plot
plot_host <- ggtree(host, size = 4, color = "black"#, branch.length = "none"
) %<+% labdf + ## the weird symbol is ggtree "attacher"
  geom_tiplab(aes(image = uid, color = forcolor), geom = "phylopic", 
              offset = 0.1,
              size = c(0.09, # human 
                       0.14, 0.14, # gorillas
                       0.11, # chimp
                       0.16
              )) +
  # show root
  geom_rootedge(rootedge = 0.05, size = 4) +
  xlim(-0.05,0.7) +
  scale_color_manual(values = manga.cols) +
  theme(legend.position = "none",
        plot.margin = margin(r=0))

# make fungi plot
plot_fungi <- figtree %>% ggtree(size = 3, color = "black", branch.length = "none"
)%<+% #+ 
  #geom_tippoint(aes(color = label), size = 15) +
  flabdf + 
  geom_tiplab(aes(image = hyphae, color = forcolor), geom = "phylopic",
              offset = -0.7, size = 0.16, alpha = 0.7)  +
  geom_tiplab(aes(image = yeast, color = forcolor), geom = "phylopic",
              offset = -0.7, size = 0.16) +
  
  geom_rootedge(rootedge = 0.5, size = 3)+
  scale_x_reverse(limits = c(6, -1)) +
  # manual colors
  scale_color_manual(values = manga.cols) +
  theme(legend.position = "none",
        plot.margin = margin(l=0),
        panel.margin = margin(l=0))

## arrange together
fin.plot <- ggarrange(plot_host, plot_fungi, ncol = 2) +
  theme(plot.background = element_rect(color = NA, fill = NA),
        plot.margin = margin(l=0, r=0))

## ---- panel B: pairwise comparisons ----

## get objects from analysis script
source("analyses/pairwise_braycurtis.R")

# remove duplicates in "within" from making longer
dat <- ufdf %>% 
  filter(str_detect(Label, "Within")) %>% 
  mutate(newLab = Group1,
         is.within = "Within") %>% 
  dplyr::select(Item1, Item2, is.within, Distance, newLab) %>% 
  rbind(ufdf %>% 
          filter(str_detect(Label, "Between")) %>% 
          mutate(is.within = "Between") %>% 
          dplyr::select(Item1, Item2, Group1, Group2, is.within, Distance) %>% 
          pivot_longer(cols = c(Group1, Group2), names_to = "oldGroup", values_to = "newLab") %>% dplyr::select(-oldGroup)) %>% 
  # make colored panels
  mutate(fulllab = interaction(is.within, newLab, sep = "-", lex.order = TRUE)) %>% 
  ### 5/23; order by phylogeny not median
  mutate(newLab = factor(newLab, ordered = TRUE, levels = c("P. troglodytes", "H. sapiens", "G. gorilla", "G. beringei")))

# make color vector
sci.cols <- italic.cols
names(sci.cols) <- c("Within-G. beringei", "Within-G. gorilla", "Within-H. sapiens", "Within-P. troglodytes")

# add "between" colors (all white)
f.cols <- c(sci.cols, rep("white", 4))
names(f.cols)[5:8] <- c("Between-H. sapiens", "Between-P. troglodytes", "Between-G. gorilla", "Between-G. beringei")

# make plot
p1 <- ggplot(data = dat1, aes(x = is.within, y = Distance)) +
  geom_jitter(
    width = 0.2, color = "darkgrey", alpha = 0.3) +
  geom_boxplot(#aes(fill = is.within),
    aes(fill = fulllab),
    size = 1.5, outlier.shape = NA, notch = FALSE, alpha = 0.7) +
  coord_flip() +
  ggh4x::facet_wrap2(~newLab, strip = strip_themed(background_y = elem_list_rect(fill = c("#FDE725FF", "#35B779FF", "#31688EFF", "#440154FF"))), 
                     ncol = 1, strip.position = "left") +
  labs(y = "Bray-Curtis pairwise distances", x = "") +
  stat_pvalue_manual(data = stats, label = "p.signif", hide.ns = TRUE, tip.length = 0, size = 10, coord.flip = TRUE) + # brackets
  scale_fill_manual(values = f.cols) +
  #ylim(c(0, 1.15)) +
  scale_x_discrete(position = "top") +
  theme_pubr(base_size = 22) +
  theme(strip.background = element_rect(fill = c("#FDE725FF", "#35B779FF", "#31688EFF", "#440154FF", alpha = 0.9), color = NA),
        strip.text = element_text(size = 22, face = "italic", color = "white"),
        axis.text.y = element_text(size = 22),
        legend.position = "none",
        panel.spacing = unit(0.05, "cm")) 

## ---- panel C: between species ----

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

# make comparisons
mod <- aov(Distance ~ labplot, data = forplot)
t <- TukeyHSD(mod)

# use multcompview to get letters easily
cld <- multcompLetters4(mod, t, reversed = TRUE)
cldf <- as.data.frame.list(cld$labplot) %>% 
  rownames_to_column(var = "labplot") %>% 
  dplyr::select(labplot, Letters) %>% 
  # get positions
  full_join(forplot %>% group_by(labplot) %>% get_summary_stats(Distance, type = "common") %>% dplyr::select(median, labplot)) 

# make fancy colored text
labs <- c("Between <b style='color:#FDE725FF'>*P. troglodytes*</b> and <b style='color:#31688EFF'>*G. gorilla*</b>",
          "Between <b style='color:#FDE725FF'>*P. troglodytes*</b> and <b style='color:#440154FF'>*G. beringei*</b>",
          "Between <b style='color:#31688EFF'>*G. gorilla*</b> and <b style='color:#440154FF'>*G. beringei*</b>",
          "Between <b style='color:#35B779FF'>*H. sapiens*</b> and <b style='color:#FDE725FF'>*P. troglodytes*</b>",
          "Between <b style='color:#35B779FF'>*H. sapiens*</b> and <b style='color:#31688EFF'>*G. gorilla*</b>",
          "Between <b style='color:#35B779FF'>*H. sapiens*</b> and <b style='color:#440154FF'>*G. beringei*</b>")

# do some wrangling
fp <- forplot %>% 
  mutate(fulllab = case_when(
    labplot == "Between *P. troglodytes* and *G. gorilla*" ~ labs[1],
    labplot == "Between *P. troglodytes* and *G. beringei*" ~ labs[2],
    labplot == "Between *G. gorilla* and *G. beringei*" ~ labs[3],
    labplot == "Between *H. sapiens* and *P. troglodytes*" ~ labs[4],
    labplot == "Between *H. sapiens* and *G. gorilla*" ~ labs[5],
    labplot == "Between *H. sapiens* and *G. beringei*" ~ labs[6]
  ))

# more wrangling
cldf1 <- cldf %>% left_join(fp %>% dplyr::select(fulllab, labplot) %>% distinct()) 

# make the plot 
p2 <- ggplot(data = fp, aes(x = fct_reorder(fulllab, desc(Distance)), y = Distance)) +
  
  geom_point(position = position_jitter(width = 0.2, seed = 123), alpha = 0.5, color = "darkgrey") +
  geom_boxplot(size = 1.5, outlier.shape = NA, notch = FALSE, fill = "lightcyan4", alpha = 0.7) +
  geom_text(data = cldf1, aes(x = fulllab, label = Letters, y = 0.05),
            size = 10) +
  #ylim(-0.1, 1.1) +
  scale_y_continuous(limits = c(-0.1, 1.1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  coord_flip() +
  theme_pubr(base_size = 24) +
  labs(x = "", y = "Bray-Curtis pairwise distances") +
  theme(axis.text.y = element_markdown())

## ---- add B and C together (A added in Illustrator for aspect ratios) ----

ggarrange(p1, p2, ncol = 2, widths = c(0.8, 1), labels = c("A.", "B."), font.label = list(size = 26),
          vjust = 0.2) +
  theme(plot.margin = margin(t=1, unit = "cm"),
        plot.background = element_rect(color = "white", fill = "white"))
