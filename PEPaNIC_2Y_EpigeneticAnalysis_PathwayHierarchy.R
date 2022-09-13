##### Setup ##### 
## Load libraries
library(tidyverse)
library(rentrez)
library(jsonlite)
source("CpG2Gene_functions.R")

## DMP and DMR results
load("../../temp_saves/DMPResults/DMP_PtCtrl.Rdata")
load("../../output/DMRResults/DMR_PtCtrl.Rdata")

## Get CpG sites that remained after preprocessing
load("../../minfi_sets/Swabs_2yFU_interSampleQC_funnorm_Convertedmethylset.Rdata")
rownames_GRcset <- rownames(GRcset)
rm(GRcset)

## Get UCSC Refgene data from annotation file
UCSC_RefGene_Group <- data.frame(UCSC_RefGene_Accession = data.frame(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other) %>%
                                   separate_rows(UCSC_RefGene_Accession, sep=";") %>% 
                                   pull(UCSC_RefGene_Accession), 
                                 UCSC_RefGene_Group = data.frame(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other) %>%
                                   separate_rows(UCSC_RefGene_Group, sep=";") %>% 
                                   pull(UCSC_RefGene_Group))
##### DMP DMR Gene names ##### 
## Patient Control DMPs
DMP_PC <- CpG2Gene(DMP_PtCtrl$Sign_CpG)
DMP_PC_genes <- DMP_PC %>% pull(UCSC_RefGene_Name) %>% str_split(pattern = ";") %>% unlist() %>% unique()
DMP_PC_tidyAcession <- DMP_PC %>% select(CpG, UCSC_RefGene_Accession) %>% separate_rows(UCSC_RefGene_Accession, sep=";")
DMP_PC_Accession <- DMP_PC_tidyAcession %>% filter(UCSC_RefGene_Accession != "") %>% pull(UCSC_RefGene_Accession) %>% unique()

## Patient Control DMRs
CpG_DMRPC <- DMR2CpG(data.frame(DMR_PtCtrl$DMRRanges)[,1:3], rownames_GRcset)
DMR_PC <- CpG2Gene(CpG_DMRPC$CpG)
DMR_PC_genes <- DMR_PC %>% pull(UCSC_RefGene_Name) %>% str_split(pattern = ";") %>% unlist() %>% unique()
DMR_PC_tidyAcession <- DMR_PC %>% select(CpG, UCSC_RefGene_Accession) %>% separate_rows(UCSC_RefGene_Accession, sep=";")
DMR_PC_Accession <- DMR_PC_tidyAcession %>% filter(UCSC_RefGene_Accession != "") %>% pull(UCSC_RefGene_Accession) %>% unique()


# df_DAVID <- merge_genelists(list(DMP_PC = DMP_PC_genes, DMR_PC = DMR_PC_genes))
# write_tsv(df_DAVID, file = "../../output/PathwayResults/DMPR_DavidInput.tsv")


DMP_PC_tidyAcession$gene_ids <- get_gene_ids(DMP_PC_tidyAcession$UCSC_RefGene_Accession) %>% unlist
DMR_PC_tidyAcession$gene_ids <- get_gene_ids(DMR_PC_tidyAcession$UCSC_RefGene_Accession) %>% unlist

DMP_PC_entrezdf <- bind_rows(get_gene_info(DMP_PC_tidyAcession$gene_ids)) %>% 
  mutate(gene_link = paste0("https://www.ncbi.nlm.nih.gov/gene/", gene_id)) 
DMR_PC_entrezdf <- bind_rows(get_gene_info(DMR_PC_tidyAcession$gene_ids)) %>% 
  mutate(gene_link = paste0("https://www.ncbi.nlm.nih.gov/gene/", gene_id)) 

DMP_PC_entrezdf %>% filter(gene_name =="Not found")
DMR_PC_entrezdf %>% filter(gene_name =="Not found")

DMP2Gene_PC <- left_join(DMP_PC_tidyAcession, DMP_PC_entrezdf, by=c("gene_ids"="gene_id"))
DMR2Gene_PC <- left_join(DMR_PC_tidyAcession, DMR_PC_entrezdf, by=c("gene_ids"="gene_id"))
DMPR2Gene_PC <- bind_rows(DMP2Gene_PC %>% mutate(origin = "DMP"), DMR2Gene_PC %>% mutate(origin = "DMR"))


# Get CpG info from DMPs and DMRs
CpGAnnot_ID <- as.data.frame(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other) %>%rownames_to_column("CpG") %>%
  select(CpG, UCSC_RefGene_Accession) %>%
  separate_rows(UCSC_RefGene_Accession, sep = ";")
CpGAnnot_group <- as.data.frame(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other) %>%
  select(UCSC_RefGene_Group) %>%
  separate_rows(UCSC_RefGene_Group, sep = ";")
CpGAnnot <- bind_cols(CpGAnnot_ID, CpGAnnot_group) %>% distinct()

CpG_DMP <- left_join(DMP_PC_tidyAcession %>% distinct(CpG, gene_ids, .keep_all = T),
                     DMP_PtCtrl$DMPresult %>% rownames_to_column("CpG"),
                     by="CpG") %>%
  left_join(CpGAnnot, by=c("CpG", "UCSC_RefGene_Accession")) %>%
  mutate(Hypomethylated_Patient = sign(t)<0,
         UCSC_RefGene_Group_short = case_when(UCSC_RefGene_Group %in% c("1stExon")~"1stExon",
                                              UCSC_RefGene_Group %in% c("3'UTR", "5'UTR")~"UTR",
                                              UCSC_RefGene_Group %in% c("Body", "ExonBnd")~"Body",
                                              UCSC_RefGene_Group %in% c("TSS1500", "TSS200")~"Promotor",
                                              TRUE ~ "Intergenic"))

CpG_DMR <- left_join(DMR_PC_tidyAcession %>% distinct(CpG, gene_ids, .keep_all = T), CpG_DMRPC %>% select(-c(UCSC_RefGene_Accession, UCSC_RefGene_Group, strand)), by="CpG") %>%
  left_join(DMR_PtCtrl$DMRRanges %>% as.data.frame() %>% rownames_to_column("DMR") %>% mutate(DMR = paste0("DMR", DMR)), by="DMR") %>%
  left_join(CpGAnnot, by=c("CpG", "UCSC_RefGene_Accession")) %>%
  mutate(UCSC_RefGene_Group_short = case_when(UCSC_RefGene_Group %in% c("1stExon")~"1stExon",
                                              UCSC_RefGene_Group %in% c("3'UTR", "5'UTR")~"UTR",
                                              UCSC_RefGene_Group %in% c("Body", "ExonBnd")~"Body",
                                              UCSC_RefGene_Group %in% c("TSS1500","TSS200")~"Promotor",
                                              TRUE ~ "Intergenic"))


# Save results
save(DMP2Gene_PC, DMR2Gene_PC, DMPR2Gene_PC, file = "../../output/PathwayResults/DMPR2Gene_PC.Rdata")
write.xlsx(CpG_DMP, "../../output/DMPResults/CpGinDMP_PC.xlsx")
write.xlsx(CpG_DMR %>%
             dplyr::select(-c(`strand`, overlapping.genes)) %>%
             dplyr::rename("DMR_start"=start,
                           "DMR_end"=end,
                           "DMR_width"=width,
                           "DMR_n_CpGs" = no.cpgs),
           "../../output/DMRResults/CpGinDMR_PC.xlsx")

##### Get Pathway info from KEGG database #####
kegg_brite_data <- readLines("http://rest.kegg.jp/get/br:hsa00001/json")
kegg_brite <- fromJSON(kegg_brite_data)
#kegg_brite <- fromJSON("http://rest.kegg.jp/get/br:hsa00001/json")

df_KEGG_hierarchy_all <- data.frame(parent = character(), child = character(), n_paths = numeric(), hsa_ids = character())
df_KEGG_pathwayGenes_all <- data.frame(path = character(), genes = character(), n_genes = integer())

for(l1 in 1:length(kegg_brite$children$children)){
  for(l2 in 1:length(kegg_brite$children$children[[l1]]$children)){
    df_KEGG_hierarchy_all_temp <- data.frame(parent = kegg_brite$children$name[[l1]],
                                             child = kegg_brite$children$children[[l1]]$name[[l2]],
                                             n_paths = length(kegg_brite$children$children[[l1]]$children[[l2]]$children),
                                             hsa_ids = kegg_brite$children$children[[l1]]$children[[l2]]$name %>% 
                                               lapply(function(i) str_split(i," ", n=2)[[1]][[1]]) %>%
                                               unlist() %>%
                                               paste0(collapse = " // "))
    for(l3 in 1:length(kegg_brite$children$children[[l1]]$children[[l2]]$children)){
      if(!is.null(kegg_brite$children$children[[l1]]$children[[l2]]$children[[l3]])){
        df_KEGG_pathwayGenes_all_temp <- data.frame(path = kegg_brite$children$children[[l1]]$children[[l2]]$name[[l3]],
                                                    genes = kegg_brite$children$children[[l1]]$children[[l2]]$children[[l3]]$name %>% paste0(collapse = " // "),
                                                    n_genes = length(kegg_brite$children$children[[l1]]$children[[l2]]$children[[l3]]$name))
        df_KEGG_pathwayGenes_all <- bind_rows(df_KEGG_pathwayGenes_all, df_KEGG_pathwayGenes_all_temp)
      }
    }
    df_KEGG_hierarchy_all <- bind_rows(df_KEGG_hierarchy_all, df_KEGG_hierarchy_all_temp)
  }
}

df_KEGG_pathwayGenes_all <- df_KEGG_pathwayGenes_all %>%
  separate(path, c("path_id", "path"), sep = " ", extra="merge")

df_KEGG_hierarchy_all <- df_KEGG_hierarchy_all %>% 
  group_by(parent)%>% 
  summarise(n_paths = sum(n_paths)) %>% 
  mutate(child = parent, parent = paste0(" ",kegg_brite$name)) %>% 
  select(parent, child, n_paths) %>%
  bind_rows(df_KEGG_hierarchy_all) %>%
  separate(parent, into = c("parent_ID", "parent_name"), sep = "\\s", extra = "merge", convert = T) %>% 
  separate(child, into = c("child_ID", "child_name"), sep = "\\s", extra = "merge", convert = T)

all_hsa_ids <- unique(df_KEGG_pathwayGenes_all$path_id)
df_KEGG_hierarchy <- data.frame(hsa_id = character(), name = character(), term = character(), description = character(), depth = numeric())

for(hsa in all_hsa_ids){
  print(paste0("http://rest.kegg.jp/get/pathway+hsa",hsa))
  df <- NULL
  #try(df <- read_table(paste0("http://rest.kegg.jp/get/pathway+hsa",hsa), progress = F, col_names = c("ENTRY", paste0("X", 1:1000))))
  try(df <- read_table(readLines(paste0("http://rest.kegg.jp/get/pathway+hsa",hsa)), progress = F, col_names = c("ENTRY", paste0("X", 1:1000))))
  if(is.null(df)) df <- data.frame(ENTRY = character())
  
  if(nrow(df %>% filter(ENTRY == "PATHWAY_MAP"))){
    name <- unlist(df %>% filter(ENTRY == "PATHWAY_MAP") %>% pull(2))
  }else{
    name <- NA
  }
  
  if(nrow(df %>% filter(ENTRY == "CLASS"))){
    hierarchy <- str_split(df %>% filter(ENTRY == "CLASS") %>% pull(2), pattern = "; ") %>% unlist()
  }else{
    if(!is.na(name)){
      hierarchy <- "Not Available"
    }else{
      hierarchy <- NA 
    }
  }
  
  if(nrow(df %>% filter(ENTRY == "DESCRIPTION"))){
    description <- unlist(df %>% filter(ENTRY == "DESCRIPTION")) %>% .[2:1001]
    description <- paste0(description[!is.na(description)],collapse=" ")
    
  }else{
    description <- NA
  }
  
  
  df_KEGG_hierarchy_temp <- data.frame(hsa_id = hsa, name = name, term = hierarchy, description = description, depth = 1:length(hierarchy))
  df_KEGG_hierarchy <- bind_rows(df_KEGG_hierarchy, df_KEGG_hierarchy_temp)
}

DMPR2KeggGenes <- df_KEGG_pathwayGenes_all %>%
  separate_rows(genes, sep = " // ") %>%
  separate(genes, c("Abbr", "gene_name_KEGG"), sep="; ", extra = "merge") %>%
  separate(Abbr, c("PMID_gene_ID", "gene_abbr"), sep = " ", extra="merge") %>%
  left_join(DMPR2Gene_PC, by=c("PMID_gene_ID"="gene_ids")) %>%
  left_join(UCSC_RefGene_Group, by="UCSC_RefGene_Accession")

KeggGenes2Hierarchy <- left_join(df_KEGG_hierarchy_all %>% separate_rows(hsa_ids, sep = " // "),
                                 DMPR2KeggGenes %>% filter(!is.na(path_id)),
                                 by=c("hsa_ids"="path_id")) %>% 
  left_join(df_KEGG_hierarchy %>%
              filter(!is.na(hsa_id)) %>%
              select(hsa_id, pathway_description = description, depth),
            by=c("hsa_ids"="hsa_id")) %>% 
  select(1:4, depth, hsa_ids, path,  n_paths, pathway_description,
         PMID_gene_ID, Accession_gene_ID = UCSC_RefGene_Accession,
         gene_abbr,
         gene_name, gene_fullname, gene_description, gene_link,
         CpG_location = UCSC_RefGene_Group, CpG_origin = origin)




##### Statististics from hierarchical distribution of differentially methylated pathways #####
FuncAnn_stats_path <- KeggGenes2Hierarchy %>% 
  filter(!parent_name %in% c('Brite Hierarchies', 'Not Included in Pathway or Brite', 'hsa00001')) %>%
  select(parent_ID, parent_name, child_ID, child_name, path, hsa_ids, gene_abbr, Accession_gene_ID) %>%
  distinct() %>% mutate(IsDMGene = !is.na(Accession_gene_ID)) %>%
  group_by(parent_ID, parent_name, child_ID, child_name, hsa_ids, path) %>%
  summarise(n_DMGenes = sum(IsDMGene),
            n_Genes = n(),
            p_DMGenes = round(n_DMGenes/n_Genes*100,3),
            wp_DMGenes = ifelse(n_DMGenes, round(log(n_DMGenes)*p_DMGenes, 3), -1)) %>%
  ungroup() %>% filter(n_Genes>1) %>%
  arrange(-wp_DMGenes)

FuncAnn_stats_subdomain <- FuncAnn_stats_path %>%
  group_by(parent_ID, parent_name, child_ID, child_name) %>%
  summarise( n_DMPaths = sum(n_DMGenes>0),
             n_Paths = n(),
             p_DMPaths = round(n_DMPaths/n_Paths*100, 3),
             wp_DMPaths = ifelse(n_DMPaths, round(log(n_DMPaths)*p_DMPaths, 3), -1),
             n_DMGenes = sum(n_DMGenes),
             n_Genes = sum(n_Genes),
             p_DMGenes = round(n_DMGenes/n_Genes*100, 3),
             wp_DMGenes = ifelse(n_DMGenes, round(log(n_DMGenes)*p_DMGenes, 3), -1)) %>%
  ungroup %>%
  arrange(-wp_DMPaths)


FuncAnn_stats_domain <- FuncAnn_stats_subdomain %>% group_by(parent_ID, parent_name) %>%
  summarise(n_DMPaths = sum(n_DMPaths),
            n_Paths = sum(n_Paths),
            p_DMPaths = round(n_DMPaths/n_Paths*100, 3),
            wp_DMPaths = ifelse(n_DMPaths, round(log(n_DMPaths)*p_DMPaths, 3), -1),
            n_DMGenes = sum(n_DMGenes),
            n_Genes = sum(n_Genes),
            p_DMGenes = round(n_DMGenes/n_Genes*100, 3),
            wp_DMGenes = ifelse(n_DMGenes, round(log(n_DMGenes)*p_DMGenes, 3), -1)) %>%
  arrange(-wp_DMPaths) %>%
  ungroup()


save(FuncAnn_stats_domain,
     FuncAnn_stats_subdomain,
     FuncAnn_stats_path,
     file="../../output/PathwayResults/HierarchicalFunctionalAnnotation_stats.Rdata")

save(KeggGenes2Hierarchy, CpG_DMP, CpG_DMR, file="../../temp_saves/HierarchicalFunctionalAnnotation.Rdata")

FuncAnn_stats_path %>%
  select(-c(1,3,10)) %>%
  mutate(link = paste0("https://www.kegg.jp/entry/pathway+hsa", hsa_ids)) %>%
  write.xlsx("../../output/PathwayResults/KEGGPathwayAnnotation_Categorisation.xlsx")