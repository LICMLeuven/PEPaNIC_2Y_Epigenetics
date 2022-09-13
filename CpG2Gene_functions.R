library(tidyverse)
library(minfi)

load("../../metadata/Variable_names/entrez_api_key.Rdata")
# entrez_api_key <- ""

annot_cols <- c('UCSC_RefGene_Name','UCSC_RefGene_Accession', 'UCSC_RefGene_Group')

# Get Gene info from CpG names
CpG2Gene <- function(CpG_selection){
  df <- as.data.frame(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other[CpG_selection, annot_cols]) %>%
    rownames_to_column("CpG") %>%
    separate_rows(UCSC_RefGene_Name, UCSC_RefGene_Group, sep=";") %>%
    distinct() %>%
    left_join(as.data.frame(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations) %>% rownames_to_column("CpG"), by="CpG")
  return(df)
}

# Get CpG names from DMRs
DMR2CpG <- function(DMR_Range, available_CpGs){
  CpG_positions <- as.data.frame(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations)[available_CpGs,] %>% rownames_to_column("CpG")
  lapply(1:nrow(DMR_Range), function(DMR){
    CpG_positions %>%
      filter(chr == DMR_Range[DMR,'seqnames'], pos >= DMR_Range[DMR,'start'], pos <= DMR_Range[DMR, 'end']) %>%
      mutate(DMR = paste0("DMR", DMR))
  }) %>% bind_rows() %>% 
    left_join(as.data.frame(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other[, annot_cols]) %>% rownames_to_column("CpG"), by="CpG")
}

# Merge multiple gene lists (practical for DAVID import)
merge_genelists <- function(genelist){
  maxLength <- max(sapply(genelist, length))
  df <- lapply(1:length(genelist), function(i){
    data.frame(gene_list = c(genelist[[i]], rep(NA,(maxLength-length(genelist[[i]])))))
  }) %>% bind_cols()
  colnames(df) <- names(genelist)
  return(df)
}

# Get gene IDs from NCBI gene database using RefGene Accession IDs (entrez)
get_gene_ids <- function(UCSC_RefGene_Accession){
  return(lapply(UCSC_RefGene_Accession, function(acc_id){
    #print(acc_id)
    if(acc_id=="") return("CpG not within gene")
    gene_id <- entrez_search( # Search for genes in NCBI gene database
      db = "gene",
      term=paste0(
        "Human[ORGN] AND alive[prop] AND ", # Only human genome
        acc_id),
      api_key = entrez_api_key)$ids
    if(length(gene_id)==1){
      return(gene_id)
    } else {
      if(length(gene_id)==0){
        return("No gene ID Found")
      } else {
        gene_id2 <- paste0(gene_id, ";")
        return(gene_id2)
      }
    }
  }))
}

# Get gene info from NCBI gene database using NCBI gene IDs 
get_gene_info <- function(gene_ids){
  list_entrez <- lapply(unique(gene_ids), function(gene_id){
    if(gene_id %in% c("CpG not within gene", "No gene ID Found")){
      return(data.frame(gene_id = NA,
                        gene_name = NA,
                        gene_fullname = NA,
                        gene_description = NA))
    }
    gene_summary <- NULL
    try(gene_summary <- entrez_summary("gene", gene_id, api_key = entrez_api_key))
    #Sys.sleep(1)
    
    if(!is.null(gene_summary)){
      return(data.frame(gene_id = gene_id,
                        gene_name = gene_summary$name,
                        gene_fullname = gene_summary$description,
                        gene_description = gene_summary$summary))
    } else {
      return(data.frame(gene_id = gene_id,
                        gene_name = "Not found",
                        gene_fullname = NA,
                        gene_description = NA))
    }
  })
}