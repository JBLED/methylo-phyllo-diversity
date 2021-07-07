# methylo-phyllo-diversity
FINE-SCALE ADAPTATIONS TO ENVIRONMENTAL VARIATION AND GROWTH STRATEGIES DRIVE PHYLLOSPHERE METHYLOBACTERIUM DIVERSITY.

https://github.com/JBLED/methylo-phyllo-diversity.git

# IF USING ANY INFORMATION RELATED TO THIS PROJECT, PLEASE READ AND CITE:

Leducq et al. (2021) Fine-scale adaptations to environmental variation and growth strategies drive phyllosphere Methylobacterium diversity.

PRE-PRINT VERSION: https://www.biorxiv.org/content/10.1101/2021.06.04.447128v4 and updated versions

# DESCRIPTION

# 1-Phylo-of-plant-ass-Methylo-diversity

- Script-Reference-Genome-Ecology.R
Methylobacterium diversity associated with plants, especially the phyllosphere: map the ecological origin of Methylobacteriaceae genomes publicly available in September 2020, including 153 Methylobacteria, 30 Microvirga and 2 Enterovirga on a consensus phylogenetic tree built from Methylobacteriaceae rpoB gene nucleotide sequences. Show the distribution of pairwise nucleotide similarity accross and within Methylobacteriaceae genus and Methylobacterium clades. 

# 2-Bacteria-community-timeline-16s

Bacterial phyllosphere diversity assessed by Illumina sequencing of the 16S RNA ribosomal gene using primers 799F-1115R targeting the V5-V6 region, performed on 46 phyllosphere samples from 13 trees from two forests sampled 3-4 times throughout the 2018 growth season.

- 1-Rarefaction-and-control-16S.R
ASV (Amplicon sequence variants) filtering and rarefaction. 

- 2-ASV-Taxonomy-assignation-16S-SILVA.R
ASV taxonomic assignation using SILVA v138 16s database. 

- 3-Meta-Analyses-16s.R
Community analysis (Diversity, PERMANOVA, PCA)

# 3-culturable-Methylo-diversity-timeline

- Script-Strains-2018.R
Community analysis on the isolable part of Methylobacterium diversity from 36 trees from two forests sampled 3-4 times throughout the 2018 growth season.

# 4-Methylobacterium-community-timeline-rpoB

Community analysis of Methylobacterium diversity assessed by barcoding from rpoB amplicon sequencing of 179 leaf samples from 53 trees in two forests sampled 3-4 times throughout the 2018 growth season.

- 1-Rarefaction-and-control-rpoB.R
ASV (Amplicon sequence variants) filtering and rarefaction (using a synthetic METH community as positive control)

- 2-ASV-Taxonomy-assignation-phylogeny-rpoB.R
ASV taxonomic assignation using a homemade rpoB database, curated by phylogeny. 

- 3-Filtering-Methylobacterium-ASV+scaling-tree.R
Filtering of Methylobacteriaceae ASV (Step 1), scaling of the Methylobacteriaceae ASV phylogeny (Step 2) on pairwise nucleotide similarity (Step 3) and comparison of Methylobacterium diversities evaluated by barcoding and isolation (Step 4). 

- 4-Meta-Analyses-rpoB.R
Community analysis (Diversity, PERMANOVA, PCA) and spatial/temporal autocorrelation analyses.

# 5-Growth-performance

Growth monitoring of 79 Methylobacterium isolates (sampled in 2018 in two forests) for four temperature treatments mimicking temperature variations during the growing season.

- 1-Correction_photos
Correction for background noise of pixel intensity of original plate pictures (three time points x temperature x 24 plates)

- 2-Measure_intensity
Determination of average intensity of bacterial spots (36 spots per petri).

- 3-Predict-Growth-curves-Meta-analysis.R
Correction of average pixel intensity for edge effects (competition), prediction of log normal growth curves, rate growth rate and maximum growth (yield) from tree time points and meta-analyses on growth rate and yield (PCA, ANOVA).
