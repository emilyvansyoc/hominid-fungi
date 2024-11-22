### sequence divergence function
# EVS 11/2024

library(tidyverse)
library(vegan)


myDiv <- function(alignment, # fasta file of genus level alignment
                  OTU, # character of OTU to pull
                  dist.type # from dist.dna 
) {
  
  # print progress message
  cat("\n ##### Working on sequence divergence for ", OTU, " ######## \n")
  
  # subset alignment to get the OTU of interest
  sub <- alignment[str_detect(labels(alignment), paste0(OTU, ";")), ]
  # print synopsis to screen
  print(sub)
  
  # get the divergence
  div <- dist.dna(sub, model = dist.type)
  
  # data frame-ize
  groupdf <- data.frame(ids = labels(div)) %>% 
    mutate(hominid = case_when(
      str_detect(ids, "MountGorilla") ~ "G. beringei",
      str_detect(ids, "LowGorilla") ~ "G. gorilla",
      str_detect(ids, "Human") ~ "H. sapiens",
      str_detect(ids, "Chimp") ~ "P. troglodytes"
    ))
  disdf <- dist_groups(div, g = groupdf$hominid)
  
  # print the number of exact sequence variants for each hominid
  cat("\n Number of sequence variants for each hominid: \n")
  print(groupdf %>% group_by(hominid) %>% count())
  
  # set up for plot
  forplot <- disdf %>% 
    filter(!str_detect(Label, "Within")) %>% 
    mutate(labplot = case_when(
      Label == "Between G. gorilla and P. troglodytes" ~ "Between *P. troglodytes* and *G. gorilla*",
      Label == "Between H. sapiens and P. troglodytes" ~ "Between *H. sapiens* and *P. troglodytes*",
      Label == "Between G. beringei and P. troglodytes" ~ "Between *P. troglodytes* and *G. beringei*",
      Label == "Between G. gorilla and H. sapiens" ~ "Between *H. sapiens* and *G. gorilla*",
      Label == "Between G. beringei and G. gorilla" ~ "Between *G. gorilla* and *G. beringei*",
      Label == "Between G. beringei and H. sapiens" ~ "Between *H. sapiens* and *G. beringei*"
    )) 
  
  # make test
  mod <- aov(Distance ~ labplot, data = forplot)
  t <- TukeyHSD(mod)
  cld <- multcompLetters4(mod, t, reversed = TRUE)
  cldf <- as.data.frame.list(cld$labplot) %>% 
    rownames_to_column(var = "labplot") %>% 
    dplyr::select(labplot, Letters) %>% 
    # get positions
    full_join(forplot %>% group_by(labplot) %>% get_summary_stats(Distance, type = "common") %>% dplyr::select(mean, sd, labplot)) 
  # print stat test
  cat("\n ANOVA test: \n")
  print(cldf)
  
  # plot
  #myplot <- ggplot(data = forplot, 
  #                aes(x = fct_reorder(labplot, desc(Distance)), y = Distance)) +
  
  # geom_point(position = position_jitter(width = 0.2, seed = 123), alpha = 0.5, color = "darkgrey") +
  # geom_boxplot(size = 1.5, outlier.shape = NA, notch = FALSE, fill = "lightcyan4", alpha = 0.7) +
  # geom_text(data = cldf, aes(x = labplot, label = Letters, y = 0),
  # size = 10) +
  #ylim(-0.1, 1.1) +
  #scale_y_continuous(limits = c(-0.1, 1.1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  #coord_flip() +
  # theme_pubr() +
  #labs(x = "", y = "Nucleotide sequence divergence", title = OTU) +
  #theme(axis.text.y = element_markdown())
  
  # plot in the window
  #plot(myplot)
  
  # return stat test and plot
  #return(list(forplot, myplot))
  return(list(forplot))
  
  
  
}
