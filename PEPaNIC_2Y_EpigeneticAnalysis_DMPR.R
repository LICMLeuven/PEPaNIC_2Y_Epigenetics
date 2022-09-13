  ##### Setup #####
set.seed(1234)
library(LICMEpigenetics)
library(magrittr)
library(dplyr)
library(openxlsx)



##### Preprocessing #####
Heatmap_vars_cat <- c("Center",	"Gender", "Ptn_Cntrl", "Early_Late", "Malignancy", "Syndrome", "Non_European", "Non_Caucasian", "Plate", "Well", "Chip")
Heatmap_vars_cont <- c("Age")

df_metadata <- read.xlsx("../metadata/metadata_swabs2yFU.xlsx")
df_plate <- read.xlsx("../metadata/swabinformation_swabs2yFU.xlsx", sheet = "Identification-swab-final") %>% filter(is.na(Opmerking))
df_swabs <- left_join(df_metadata, df_plate %>% select(Swabnr, Plate = Plate.DNA.methylation), by = "Swabnr")

## Local
memory.limit(size=56000) # in case of 64-bit system
defectSamples <- c("201868500039_R06C01_400989_SWAB_Patient_Pepanic10418_Late_2y", "201868500039_R07C01_400995_SWAB_Patient_Pepanic10413_Early_2y", "201868590200_R02C01_301343_SWAB_Patient_Pepanic534_Early_2y") # 2 of the 3 defect samples were rerun and are included in the final set
badSamples <- c("201868590162_R06C01_301945_SWAB_Patient_Pepanic577_Late_2y","202226320059_R02C01_400710_SWAB_Patient_Pepanic10429_Late_2y")

set <- Make_RGset("../idats/swabs_2yFU/", "Swabs_2yFU", save = TRUE) # Make RGset from idats
set$rgset <- set$rgset[,-which(colnames(set$rgset) %in% defectSamples)] # 

set_clean <- RemoveBadSamples(set, badSamplesNames = badSamples)
save(set_clean, file=getPath(set_clean, "set"))

## Cluster
timestamp()
load("../minfi_sets/Swabs_2yFU_interSampleQC_set.Rdata")

### Normalisation
timestamp()
set <- ProbeNormalisation(set_clean) %>% # functional normalisation
  ProbeExclusion() %>%
  ConvertSet(M=TRUE, save=FALSE) %>%
  BatchEffect_prc() %>%
  BatchEffect_Visualise(df_swabs, "Swabnr", "Chip", Heatmap_vars_cat, Heatmap_vars_cont)

GRcset <- ratioConvert(set$mset)
save(GRcset, file="../minfi_sets/Swabs_2yFU_interSampleQC_funnorm_Convertedmethylset.Rdata")


##### DMP and DMR analysis #####
###### Set variables ######

DMR_lambda <-  1000
cat_vars_PtCtrl <- c("Group", "Center", "Gender", "Malignancy", "Syndrome", "Non_European", "Non_Caucasian", "Language")
cont_vars_PtCtrl <- c("Age", paste0("PC", 1:30))

## Load Epigenetic data
load("../../minfi_sets/Swabs_2yFU_interSampleQC_funnorm_Convertedmethylset.Rdata")

## Load adjustment variables
load("../../temp_saves/2y_Patient_info/2yFU_AdjVar_FunctOutcomes.Rdata")

## Load PCA batch effect data
df_PCA <- read.csv("../../output/BatchEffect/Swabs_2yFU_interSampleQC_funnorm_prc_df_extended_ControlProbes.csv") %>%
  select(Patient_Number, paste0("PC", 1:30), filenames) %>%
  separate(Patient_Number, c(NA, "Patnr"), convert = T, sep = "panic")

## Join PCA data with adjustment variables
df <- left_join(df_PCA, df_2Y, "Patnr")

## Tests
### Test that all swabs have metadata
df_missing <- df %>% filter(is.na(Group)) %>% select(Patnr) %>% arrange(Patnr)
if(nrow(df_missing)) stop("Missing Patnr")

### Test that all IDs (Swabs vs Metadata) are aligned.
matching_ids <- data.frame(meth_Ids = colnames(GRcset) %>% lapply(function(i) strsplit(i,'_')[[1]][[6]]) %>% unlist(), df_IDs = df$Patnr) %>%
  separate(meth_Ids, c(NA, "meth_Ids"), sep = "Pepanic", convert = TRUE) %>%
  mutate(test = meth_Ids==df_IDs)
table(matching_ids$test)
if(!all(matching_ids$test)) stop("Patnr does not allign with the swabs")


###### Patient-Control ######
### DMP analysis
DMP_PtCtrl <- DMP_limma(GRcset, df, cat_vars_PtCtrl, cont_vars_PtCtrl)
save(DMP_PtCtrl, file="../../temp_saves/DMPResults/DMP_PtCtrl.Rdata")

### DMR analysis
DMR_PtCtrl <- DMR_dmrcate(GRcset, DMP_PtCtrl)
save(DMR_PtCtrl, file="../../output/DMRResults/DMR_PtCtrl.Rdata")

###### Early-Late PN ######
### Data prep 
cat_vars_EL <- c("RandomisationGroup", "Center", "Gender", "Malignancy", "Syndrome", "Non_European", "Non_Caucasian", "Language", "Diagnosis", "STRONGkids")
cont_vars_EL <- c("Age", "PIM3", "PeLOD", paste0("PC", 1:30))

df_EL <- df %>% filter(Group=="Patient")
filenames_EL <-  df_EL %>% pull(filenames)

#### Select significant CpG sites and exclude healthy controls
load("../../minfi_sets/Swabs_2yFU_interSampleQC_funnorm_Convertedmethylset.Rdata")
GRcset_EL_DMP <- GRcset[DMP_PtCtrl$Sign_CpG, filenames_EL]

#### Select CpG sites within DMRs
DMR_PtCtrl_Locations <- data.frame(DMR_PtCtrl$DMRRanges)[,1:3]
IlluminaPosition_data <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations %>% data.frame() %>% rownames_to_column("CpG")
DMR_PtCtrl_CpGs <- lapply(1:nrow(DMR_PtCtrl_Locations), function(DMR_i){
  IlluminaPosition_data %>% 
    filter(chr == DMR_PtCtrl_Locations[DMR_i, "seqnames"],
           pos >= DMR_PtCtrl_Locations[DMR_i, "start"],
           pos <= DMR_PtCtrl_Locations[DMR_i, "end"]) %>%
    mutate(DMR = DMR_i)
}) %>% bind_rows() %>% filter(CpG %in% rownames(GRcset))

GRcset_EL_DMR <- GRcset[DMR_PtCtrl_CpGs$CpG, filenames_EL]

### DMP Analysis Early vs Late PN
DMP_EL <- DMP_limma(GRcset_EL_DMP, df_EL, cat_vars_EL, cont_vars_EL)
save(DMP_EL, file="../../temp_saves/DMPResults/DMP_EL.Rdata")

### DMR Analysis Early vs Late PN
DMR_DMP_EL <- DMP_limma(GRcset_EL_DMR, df_EL, cat_vars_EL, cont_vars_EL)
try(DMR_EL <- DMR_dmrcate(GRcset_EL_DMR, DMR_DMP_EL))
save(DMR_EL, file="../../temp_saves/DMPResults/DMR_EL.Rdata")





