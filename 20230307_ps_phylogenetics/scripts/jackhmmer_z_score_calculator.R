# Script to calculate z-scores from jackhmmer data
# Dependencies -----------------------------------------------------------------
library(tidyverse)

# Options ----------------------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)

jackhmmer_file <- "./data/jackhmmer/compiled_jackhmmer_data_corrected_names.csv"
output_dir     <- "./data/jackhmmer/"
side_meta      <- "./raw_data/siderophore_metadata.csv"
VFDB_meta      <- "./raw_data/VFDB_metadata_2.csv"
bact_meta      <- "./raw_data/bactibase_metadata.csv"
threads        <- 1

# Functions --------------------------------------------------------------------
# Load Data --------------------------------------------------------------------
# Jackhmmer data
writeLines(paste0(Sys.time(), " Loading jackhmmer data"))

jackhmmer_data <- read.csv(jackhmmer_file) %>%
  mutate(factor_id_name = ifelse(test = is.na(factor_id_name),
                                 yes  = factor_id,
                                 no   = factor_id_name),
         factor_category_name = ifelse(test = is.na(factor_category_name),
                                       yes  = factor_category,
                                       no   = factor_category_name),
         gene_id = query_gene)

writeLines(paste0(Sys.time(), " Done"))

# Load in protein metadata
writeLines(paste0(Sys.time(), " Loading metadata"))

VFDB_meta_trim <-  read.csv(file = VFDB_meta) %>%
  select(gene_id, VFID, VFCID, VF_Name, VFcategory ) %>%
  transmute(gene_id        = gene_id,
            factor_id      = VFID,
            factor_id_name = VF_Name,
            factor_category = VFCID,
            factor_category_name = VFcategory,
            type           = "Virulence_Factor")

bacti_meta_trim <- read.csv(file = bact_meta) %>%
  select(id, Name, Class) %>%
  transmute(gene_id   = id,
            factor_id = id,
            factor_id_name = Name,
            factor_category = ifelse(is.na(Class), "Unknown", Class),
            factor_category_name = ifelse(is.na(Class), "Unknown", Class),
            type     = "Bacteriocin")

sidero_meta_trim <-  read.csv(file = side_meta) %>%
  select(Protein.Referenece, Gene.name) %>%
  transmute(gene_id       = Protein.Referenece,
            factor_id     = Protein.Referenece,
            factor_id_name  = Gene.name,
            factor_category = Gene.name,
            factor_category_name = Gene.name,
            type     = "Siderophore_Biosynthesis")

jackhmmer_metadata <- do.call(rbind, list(VFDB_meta_trim, bacti_meta_trim, sidero_meta_trim))  %>%
  mutate(factor_id_name = ifelse(test = is.na(factor_id_name),
                                 yes  = factor_id,
                                 no   = factor_id_name),
         factor_category_name = ifelse(test = is.na(factor_category_name),
                                       yes  = factor_category,
                                       no   = factor_category_name))

writeLines(paste0(Sys.time(), " Done"))

# Calculate z-scores -----------------------------------------------------------
levels <- c("gene_id", "factor_id_name", "factor_category_name")
for(level in levels){
  
  writeLines(paste0(Sys.time(), " Calcualting z-scores for ", level))
  
  # Make table of all bacteria-gene combinations
  temp_function <- function(x){
    
    output <- jackhmmer_metadata %>%
      select(UQ(level), type) %>%
      distinct() %>%
      mutate(bacteria = x)
    
    return(output)
  }
  
  tmp_full_table <- mclapply(X        = as.list(unique(jackhmmer_data$target_organism)),
                             FUN      = temp_function,
                             mc.cores = threads)
  
  full_table <- do.call(rbind, tmp_full_table)
  
  
  # Calculate mean count of each gene per genome
  mean_counts <- jackhmmer_data %>%
    group_by(pick(target_organism, UQ(level))) %>%
    summarise(n = n(),
              .groups = "drop") %>%
    right_join(full_table,
               by = setNames(c(level, "bacteria"),
                             c(level, "target_organism"))) %>%
    replace_na(list(n = 0)) %>%
    group_by(pick(UQ(level))) %>%
    summarise(mean = mean(n),
              sd   = sd(n),
              .groups = "drop") %>%
    mutate(sd = ifelse(is.na(sd), 0, sd))
  
  # Calculate z-scores 
  z_scores <- jackhmmer_data %>%
    group_by(pick(target_organism, UQ(level))) %>%
    summarise(n       = n(),
              .groups = "drop") %>%
    right_join(full_table,
               by = setNames(c(level, "bacteria"),
                             c(level, "target_organism"))) %>%
    replace_na(list(n = 0)) %>%
    left_join(mean_counts,
              by = setNames(level, level)) %>%
    mutate(z_score = ifelse(is.nan((n - mean)/sd),
                            0,
                            (n - mean)/sd))
  
  # Save out z-scores
  date <- gsub("-", "", Sys.Date())
  
  writeLines(paste0(Sys.time(), " Saving out..."))
  
  write.csv(x         = z_scores,
            file      = paste0(out_dir,"/", date, "_", level, "_jackhmmer_z_scores.csv"),
            row.names = FALSE)
  
  rm(mean_counts, z_scores)
  gc()
  
  writeLines(paste0(Sys.time(), " Done\n"))
  
}









