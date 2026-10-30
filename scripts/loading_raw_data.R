library(GEOquery)
library(limma)
library(DESeq2)
# library(AnnotationDbi)
# library(org.Hs.eg.db)

source("C:/Users/lrando11/Documents/projects/human_bpd_meta-analysis/Functions/opening_files_handling_duplicates.R")


# DNA microarray
dna_list <- list(GSE8586 = "GSE8586",
                 GSE32472 = "GSE32472",
                 GSE225881 = "GSE225881",
                 GSE188944 = "GSE188944")


# Map loops through both the elements (gse_id) and their list names (gse_name)
Map(function(gse_id, gse_name) {
  
  # Download and extract the expression set
  geo_obj <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]
  
  # Extract the expression matrix
  raw_mtx <- exprs(geo_obj)
  
  # Construct the precise file path
  output_path <- paste0("C:/Users/lrando11/Documents/projects/human_bpd_meta-analysis/Manuscript/final_raw_data/raw_", gse_name, ".rds")
  
  # Save the matrix
  saveRDS(raw_mtx, file = output_path)
  
}, dna_list, names(dna_list))




# RNA-seq
raw_gse106910 <- read.delim("C:/Users/lrando11/Documents/projects/human_bpd_meta-analysis/Manuscript/final_raw_data/GSE106910_combined_counts_geo.txt.gz")
saveRDS(raw_gse106910, "C:/Users/lrando11/Documents/projects/human_bpd_meta-analysis/Manuscript/final_raw_data/raw_GSE106910.rds")

raw_gse305315 <- read.csv("C:/Users/lrando11/Documents/projects/human_bpd_meta-analysis/Manuscript/final_raw_data/GSE305315_counts_all_samples.csv.gz")
saveRDS(raw_gse305315, "C:/Users/lrando11/Documents/projects/human_bpd_meta-analysis/Manuscript/final_raw_data/raw_GSE305315.rds")

raw_gse125873 <- read_files("C:/Users/lrando11/Documents/projects/human_bpd_meta-analysis/RNA_seq/GSE125873_Ryan")
saveRDS(raw_gse125873, "C:/Users/lrando11/Documents/projects/human_bpd_meta-analysis/Manuscript/final_raw_data/raw_GSE125873.rds")

raw_gse220135 <- read_files_gse220135("C:/Users/lrando11/Documents/projects/human_bpd_meta-analysis/RNA_seq/GSE220135_Bell/Day_1_14_28")
saveRDS(raw_gse220135, "C:/Users/lrando11/Documents/projects/human_bpd_meta-analysis/Manuscript/final_raw_data/raw_GSE220135.rds")
  
  