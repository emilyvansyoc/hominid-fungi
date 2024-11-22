# create figure 2

library(microViz)
library(phyloseq)
library(tidyverse)
library(ggpubr)
library(ggtext)
library(ggvenn)
library(rphylopic)

# get colors
source("helpers/colors.R")

# get data objects
source("analyses/diff_relativeabundance.R")

## ---- panel A: venn diagram ----
# some text labels are added in Illustrator for convenience

## wrangle sum of unique genera for venn diagram
# make nested list
tot <- pa %>% group_by(OTU, SciName) %>% summarize(tot = sum(Abundance)) 
mylist <- list('*H. sapiens*' = tot$OTU[tot$SciName == "Homo_sapien" & tot$tot > 0],
               '*P. troglodytes*' = tot$OTU[tot$SciName == "P_troglodytes_schweinfurthii" & tot$tot > 0],
               '*G. gorilla*' = tot$OTU[tot$SciName == "Gorilla_gorilla" & tot$tot > 0],
               '*G. beringei*' = tot$OTU[tot$SciName == "Gorilla_beringei" & tot$tot > 0])

# make dataframe input
tf <- tot %>% 
  mutate(ItalName = case_when(
    SciName %in% "Gorilla_beringei" ~ "*G. beringei*",
    SciName %in% "Gorilla_gorilla" ~ "*G. gorilla*",
    SciName %in% "Homo_sapien" ~ "*H. sapien*",
    SciName %in% "P_troglodytes_schweinfurthii" ~ "*P. troglodyte*"
  )) %>% 
  mutate(tf = if_else(tot > 0, TRUE, FALSE)) %>% 
  dplyr::select(value = OTU, ItalName, tf) %>%
  pivot_wider(names_from = ItalName, values_from = tf)

## build plot 
p2 <- ggplot(tf) +
  geom_venn(aes(A =  `*G. beringei*`, B = `*G. gorilla*`, C = `*H. sapien*`, D = `*P. troglodyte*`),
            show_stats = "c",
            show_set_totals = "c",
            fill_color = c("#440154FF", "#31688EFF",
                           "#35B779FF", "#FDE725FF"),
            stroke_color = NA,
            fill_alpha = 0.6,
            #set_names= c("","","",""),
            set_name_color = "white") + ### hack to "remove" labels (?!)
  coord_fixed() +
  # add phylopics for the plot
  add_phylopic(uuid = "142e0571-3b5f-443d-a887-b572a224ea22", x = -1.5, y = -1.4, ysize = 0.6, fill = "#440154FF") +
  add_phylopic(uuid = "7133ab33-cc79-4d7c-9656-48717359abb4", x = 1.55, y = -1.4, ysize = 0.6, fill = "#FDE725FF") +
  add_phylopic(uuid = "142e0571-3b5f-443d-a887-b572a224ea22", x = -0.8, y = 1.4, ysize = 0.6, fill = "#31688EFF") +
  add_phylopic(uuid = "036b96de-4bca-408e-adb6-2154fcd724ef", x = 0.8, y = 1.5, ysize = 0.8, fill = "#35B779FF") +
  ## add legend by hand
  
  add_phylopic(uuid = "7133ab33-cc79-4d7c-9656-48717359abb4", x = -2.2, y = 3, ysize = 0.3, fill = "#FDE725FF") +
  add_phylopic(uuid = "036b96de-4bca-408e-adb6-2154fcd724ef", x = -2.2, y = 2.5, ysize = 0.5, fill = "#35B779FF") +
  add_phylopic(uuid = "142e0571-3b5f-443d-a887-b572a224ea22", x = -2.2, y = 2, ysize = 0.3, fill = "#440154FF") +
  add_phylopic(uuid = "142e0571-3b5f-443d-a887-b572a224ea22", x = -2.2, y = 1.5, ysize = 0.3, fill = "#31688EFF") +
  geom_text(x = -2, y = 3, label = "P. troglodytes", fontface = "italic", size = 5, hjust = "left", check_overlap = TRUE) +
  geom_text(x = -2, y = 2.5, label = "H. sapiens", fontface = "italic", size = 5,hjust = "left", check_overlap = TRUE) +
  geom_text(x = -2, y = 2, label = "G. beringei", fontface = "italic", size = 5, hjust = "left",check_overlap = TRUE) +
  geom_text(x = -2, y = 1.5, label = "G. gorilla", fontface = "italic", size = 5, hjust = "left",check_overlap = TRUE) +
  theme_void() +
  theme(plot.margin = margin(t=-2, unit = "cm"))

## ---- panel B: phylopic dot plot ----

# build plot
p1 <- ggplot(lgm, aes(x = ItalName, y = fc, fill = issigcol)) +
  facet_wrap(~gen.order, scales = "fixed", strip.position = "left", ncol = 1,
             dir = "v") +
  geom_hline(yintercept = 0) +
  geom_segment(aes(y = 0, yend = fc, color = issigcol), position = position_dodge2(0.9), linetype = "twodash", linewidth = 0.6) +
  geom_phylopic(aes(uuid = uuid), size = 1.8, position = position_dodge(0.9)) +
  scale_fill_manual(values = mycols) +
  scale_color_manual(values = mycols, guide = "none") +
  theme_pubr()+
  labs(y = "Log2 fold change") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA, color = NA),
        strip.placement = "inside",
        strip.text.y.left = element_text(angle = 0, face = "bold.italic", hjust = 1),
        text = element_text(size = 18)) +
  coord_flip()

# ---- add together and save ----

# 
ggarrange(p2,p1,  widths = c(0.7, 1), labels = c("A.", "B.")) +
  theme(plot.background =  element_rect(fill = "white", color = "white"))
