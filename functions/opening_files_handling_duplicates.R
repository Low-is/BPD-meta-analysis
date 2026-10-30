#### Open GSE220135 files (tar files) ####
read_files_gse220135 <- function(file_path) {
  # Select only gzipped text files
  files <- list.files(path = file_path, pattern = "*.txt.gz$", full.names = TRUE)
  
  exprs_list <- list()
  
  for (i in seq_along(files)) {
    file <- files[i]
    # Use read.delim with gzfile to read gzipped files
    rna_file <- read.delim(gzfile(file), header = FALSE)
    exprs_list[[i]] <- rna_file[, -1]  # exclude the first column
  }
  
  exprs_mtx <- do.call(cbind, exprs_list)
  
  # Extract the sample IDs from the file names
  geo <- sapply(basename(files), function(x) {
    strsplit(x, split = "_")[[1]][1]
  })
  
  # Extract the gene IDs from the first file
  genes <- read.delim(gzfile(files[1]), header = FALSE)[, 1]
  
  colnames(exprs_mtx) <- geo
  rownames(exprs_mtx) <- genes
  
  # Add the gene names as the first column
  exprs_mtx <- cbind(Genes = genes, exprs_mtx)
  
  return(exprs_mtx)
}

#### Open RNA-seq TAR files ####
read_files <- function(file_path) {
  files <- list.files(path = file_path, pattern = "GSM.*\\.txt\\.gz$", full.names = TRUE)
  exprs_list <- list()
  for(i in seq_along(files)) {
    file <- files[i]
    rna_file <- read.delim(gzfile(file), header = FALSE)[-c(1:2), 7]
    exprs_list[[i]] <- rna_file
  }
  exprs_mtx <- do.call(cbind, exprs_list)
  
  geo <- sapply(basename(files), function(x) {
    strsplit(x, split = "_")[[1]][1]
  })
  
  genes <- read.delim(gzfile(files[1]), header = FALSE)[-c(1:2), 1]
  colnames(exprs_mtx) <- geo
  
  exprs_mtx <- cbind(GeneID = genes, exprs_mtx)
  
  return(exprs_mtx)
}


#### Handling duplicated genes ####
handling_duplicates <- function(data_table) {
  
  if (is.list(data_table$Genes)) {
    data_table[, Genes := vapply(Genes, function(x) {
      # Check if the list is empty or NULL
      if (is.null(x) || length(x) == 0 || all(is.na(x))) {
        return(NA_character_)  # Return NA for empty or NULL elements
      } else if (is.list(x)) {
        return(as.character(unlist(x)))  # Unlist nested lists
      } else {
        return(as.character(x))  # Directly return non-list elements
      }
    }, FUN.VALUE = character(1))] # Ensures results is a character vector
  }
    
  data_table <- na.omit(data_table[(Genes != "" & !is.na(Genes)), ])
  
  if (any(duplicated(data_table$Genes))) {
    exprs_data.table <- data_table[, lapply(.SD, mean), by = Genes, .SDcols =! 'Genes']
  }
  
  exprs_mtx <- as.matrix(exprs_data.table[, -1])
  rownames(exprs_mtx) <- exprs_data.table$Genes
  
  return(exprs_mtx)
}
