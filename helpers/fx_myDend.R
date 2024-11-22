## inputting a phyloseq object, do species aggregation and clustering

myDend <- function(physeq, distance, taxlevel # use NA if no glom
) {
  
  # glom 
  if(!is.na(taxlevel)) {
    psg <- physeq %>% 
      tax_fix() %>% 
      tax_glom(taxlevel)  %>% 
      ps_mutate(SciName = case_when(
        SpeciesCaptive %in% "Wild_Chimp" ~ "P_troglodytes_schweinfurthii",
        SpeciesCaptive %in% "Wild_Lowland Gorilla" ~ "Gorilla_gorilla",
        SpeciesCaptive %in% "Wild_Mountain Gorilla" ~ "Gorilla_beringei",
        SpeciesCaptive %in% "Human_Bantu" ~ "Homo_sapien",
        SpeciesCaptive %in% "Human_BaAka" ~ "Homo_sapien",
        SpeciesCaptive %in% "Wild_Mangabey" ~ "Mangabey"
      ))
    taxa_names(psg) <- psg@tax_table[,taxlevel]
  } else {
    psg <- physeq %>% #tax_fix()  %>% 
      ps_mutate(SciName = case_when(
        SpeciesCaptive %in% "Wild_Chimp" ~ "P_troglodytes_schweinfurthii",
        SpeciesCaptive %in% "Wild_Lowland Gorilla" ~ "Gorilla_gorilla",
        SpeciesCaptive %in% "Wild_Mountain Gorilla" ~ "Gorilla_beringei",
        SpeciesCaptive %in% "Human_Bantu" ~ "Homo_sapien",
        SpeciesCaptive %in% "Human_BaAka" ~ "Homo_sapien",
        SpeciesCaptive %in% "Wild_Mangabey" ~ "Mangabey"
      ))
  }
  
  # collapse
  df <- psg %>% ps_melt()
  dfc <- df %>% 
    # get median
    group_by(OTU, SciName) %>% 
    summarize(avg = round(mean(Abundance), 0),
              med = round(median(Abundance, 0)))
  
  # re-phyloseq
  psmed <- phyloseq(
    sample_data(data.frame(
      row.names = c("Gorilla_beringei", "Gorilla_gorilla", "Homo_sapien", "P_troglodytes_schweinfurthii", "Mangabey"),
      SciName = c("Gorilla_beringei", "Gorilla_gorilla", "Homo_sapien", "P_troglodytes_schweinfurthii", "Mangabey")
    )),
    otu_table(dfc %>% dplyr::select(-med) %>% 
                pivot_wider(names_from = OTU, values_from = avg) %>% 
                column_to_rownames(var = "SciName"), taxa_are_rows = FALSE),
    tax_table(df %>% 
                dplyr::select(OTU, Kingdom, Phylum, Class, Order, Family, Genus) %>% 
                distinct() %>% 
                column_to_rownames(var = "OTU") %>% 
                as.matrix())
  )
  
  # cluster and return dendrogram object
  if(distance == "jaccard") {
    clus <- psmed %>% dist_calc("jaccard", binary = TRUE) %>% microViz::dist_get() %>% hclust(method = "average") %>% as.phylo()
  } else {
    clus <- psmed %>% dist_calc(distance) %>% microViz::dist_get() %>% hclust(method = "average") %>% as.phylo()
  }
  
  # 
  return(clus)
  
  
  
}
