library(tidyverse)
library(openxlsx)
library(knitr)
library(kableExtra)


genes_0y <- read.xlsx("../../metadata/0y/0y_genes.xlsx")
load(file="../../temp_saves/HierarchicalFunctionalAnnotation.Rdata")
load(file="../../output/PathwayResults/HierarchicalFunctionalAnnotation_stats.Rdata")
load(file="../../output/PathwayResults/DMPR2Gene_PC.Rdata")
KEGG_pathways_annotated <- read.xlsx("../../metadata/KEGG/KEGGPathwayAnnotation_Categorisation-Final.xlsx")

# Differentially mehtylated genes 0y and overlap with 2y
c(genes_0y %>% filter(!is.na(Gene.symbol)) %>% pull(Gene.symbol) %>% unique()) %in% c(DMPR2Gene_PC %>% filter(!is.na(gene_name)) %>% pull(gene_name) %>% unique()) %>% table()

GenesRemainedDM <- c(genes_0y %>% filter(!is.na(Gene.symbol)) %>% pull(Gene.symbol) %>% unique())[c(genes_0y %>% filter(!is.na(Gene.symbol)) %>% pull(Gene.symbol) %>% unique()) %in% c(DMPR2Gene_PC %>% filter(!is.na(gene_name)) %>% pull(gene_name) %>% unique())]

GenesRemainedDM

# Differentially methylated pathways 0y and overlap 2Y
KeggGenes2Hierarchy %>% 
  filter(!is.na(gene_abbr),
         gene_abbr %in% genes_0y$Gene.symbol,
         !parent_name %in% c('Brite Hierarchies', 'Not Included in Pathway or Brite')) %>% 
  select(parent_name, child_name, path, gene_abbr) %>% distinct() %>%
  group_by(path) %>% summarise(n_DMGenes = n()) %>%
  arrange(-n_DMGenes)%>%
  kbl() %>% kable_styling()

pathways_0y <- KeggGenes2Hierarchy %>% 
  filter(!is.na(gene_abbr),
         gene_abbr %in% genes_0y$Gene.symbol,
         !parent_name %in% c('Brite Hierarchies', 'Not Included in Pathway or Brite')) %>% 
  select(parent_name, child_name, path, gene_abbr) %>% distinct()

pathways_0y$path %in% c(FuncAnn_stats_path %>% filter(n_DMGenes>0) %>% pull(path)) %>% table()
pathways_0y$path[!pathways_0y$path %in% c(FuncAnn_stats_path %>% filter(n_DMGenes>0) %>% pull(path))]
PathwaysRemainedDM <- pathways_0y$path[pathways_0y$path %in% c(FuncAnn_stats_path %>% filter(n_DMGenes>0) %>% pull(path))]

# Classification distribution 0y
## All 0y pathways
KEGG_pathways_annotated %>% filter(path %in% pathways_0y$path) %>% 
  summarise(Neurocognitive_or_Growth = sum(Neurocognitive_or_Growth=="Yes", na.rm = T),
            ICU_stay = sum(ICU_stay=="Yes", na.rm = T),
            PreAdmission = sum(PreAdmission=="Yes", na.rm = T))

## 0y pathways that remained Differentially methylated
KEGG_pathways_annotated %>% filter(path %in% PathwaysRemainedDM) %>% 
  summarise(Neurocognitive_or_Growth = sum(Neurocognitive_or_Growth=="Yes", na.rm = T),
            ICU_stay = sum(ICU_stay=="Yes", na.rm = T),
            PreAdmission = sum(PreAdmission=="Yes", na.rm = T))
