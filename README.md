# hominid-fungi
Code for analyses and figures to support the manuscript Van Syoc et al., currently under review. Work in progress and will be updated as revisions progress. Files without public sharing access are marked with the psuedo-pathfile, "/private/"

## Folders:  

### amplicon-processing
Steps to process raw sequencing data into a phyloseq object. Assuming SRA data has already been downloaded. An accompanying metadata file was obtained from collaborators.  

1. `ITS_qualitycontrol_hominid.R`: inputs raw sequencing files, then removes N's, primer sequences, and does basic quality filtering using the dada2 platform.  
2. `ITS_vsearch.sh`: inputs a directory of filtered reads, then runs the VSEARCH pipeline to merge and dereplicate reads, then clusters OTUs and detects de novo chimeras.  
5. `vsearch-to-phyloseq_rarefy_filter.R`: puts together the VSEARCH output into phyloseq and runs a rarefaction and filtering pipeline; output is the finalized phyloseq object for downstream analyses.

### analyses  
Scripts used for analyses described in the paper, including phylosymbiosis, cophylogeny, and nucleotide sequence divergence. 

1. `run_randomtrees.sh`: uses script downloaded from https://github.com/awbrooks19/phylosymbiosis (written in Python 2) to create random trees for statistical testing
2. `topological_congruency.R`: runs topological congruency tests using the framework implemented in Brooks et al 2016. 
3. `pairwise_braycurtis.R`: pairwise tests of Bray-Curtis distances between hominids to support phylosymbiosis.  
4. `diff_relativeabundance.R`: differential relative abundance tests between humans and non-human primate hominids.  
5. `cophylogeny.R`: applies the cophylogeny framework to the fungal genera. 
6. `sequence_divergence.R`: calculates nucleotide sequence divergence and creates a calibration to estimate hominid speciation times 

#### subfolder: "cophylogeny-scripts" 
Runs a bunch of tiny little scripts to build the cophylogeny framework. See `cophylogeny.R` in the 'analyses' folder that puts them all together in an R wrapper.  
5a. `run-vsearch-derep.sh`: dereplicate exact sequence variants by length and identifier (decreases alignment size)  
5b. `mafft-alignments.sh`: align sequences with MAFFT and trim with clipkit  
5c. `run-seqkit-rename.sh`: uses seqkit to assign unique identifiers to all sequences; required by downstream tree building softwares  
5d. `run-PACO.R`: R script that is submitted as an HPC job to run PACo in a parallel loop (PACo is very slow on bigger trees)  
5e. `run-parafit.R`: run Parafit on each OTU tree  
5f. `results_paco-parafit.R`: combines results from PACo and Parafit, adjusts for multiple comparisons  
5g. `determine_randomchance.R`: shuffles tip labels from the fungal OTU trees to create 100 random trees, then run both PACo and Parafit to determine the odds of observing cophylogeny given random chance  



### figures  
Scripts used to generate each of the main and supplementary figures. Most source the analysis scripts for data objects. Some figures are paneled in Adobe Illustrator or other schematic objects added in Adobe Illustrator; these are noted below. All other plot elements are generated in the R scripts.  

1. `fig1.R`: figure 1. Some text on the trees added in Illustrator.
2. `fig2.R`: figure 2. The green human icon and arrows are added in Illustrator.  
3. `fig3.R`: figure 3. The text on the schematic and fungal tree are added in Illustrator.  
4. `fig4.R`: figure 4.


### data
Data objects as described throughout the analyses scripts 

### helpers
Helper scripts with custom functions etc.  

1. `fx_myDend.R`: inputs phyloseq object and calculates a distance matrix and hierarchical clustering to return a tree object (dendrogram) 
2. `fx_myPval.R`: calculates a P value for topological congruency testing (formatted for R following framework of Brooks et al 2016).  
3. `colors.R`: sets colors for host species for figures  
4. `fx_myDiv.R`: calculates sequence divergence using dist.dna in ape
