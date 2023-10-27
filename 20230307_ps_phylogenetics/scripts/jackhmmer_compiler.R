# Script to calcualte z-scores of genes of interest for my bacteria
# and export z-scores for plotting on heatmap
# Dependencies -----------------------------------------------------------------
library(tidyverse)
library(parallel)

# Options ----------------------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)

jackhmmer_dir <- args[1]
out_dir       <- args[2]
side_meta     <- args[3]
VFDB_meta     <- args[4]
bact_meta     <- args[5]
threads       <- as.numeric(args[6])

# Functions --------------------------------------------------------------------
read_jackhmmer <- function(file){
  
  print(paste0("Reading ", file))
  tmp_file   <- readLines(con = file)[-c(1:4)]
  
  raw_jackhmmer_data <- read.delim(file   = textConnection(tmp_file),
                                   quote  = "",
                                   sep    = "",
                                   header = FALSE)
  
  colnames(raw_jackhmmer_data) <- c("target_name", "accession_1", "query_name", "accession_2",
                                    "seq_E_value", "seq_score", "seq_bias",
                                    "seq_E_value",  "seq_score",  "seq_bias",
                                    "exp", "reg", "clu",  "ov", 
                                    "env", "dom", "rep", "inc", "description_of_target")
  
  output <- raw_jackhmmer_data %>%
    select(target_name, query_name)
  
  rm(raw_jackhmmer_data)
  gc()
  
  print(paste0("Read ", file))
  return(output)
  
}

# Data -------------------------------------------------------------------------
# Silenced this for now as I saved out the table earlier
if(FALSE){
print("Args are:")
print(args)
print("")

date     <- gsub("-", "", Sys.Date())
numcores <- detectCores()

print(paste0("Finding files in ", jackhmmer_dir))

jackhmmer_files <- list.files(path       = jackhmmer_dir,
                              pattern    = "*_tbl.txt",
                              full.names = TRUE)
							  
print(paste0("Found ", length(jackhmmer_files), " files"))
print("")
print(paste0("Reading in files with ", threads, " threads, out of ", numcores, " available"))

raw_jackhmmer_data <- mclapply(X        = jackhmmer_files,
							   FUN      = read_jackhmmer,
							   mc.cores = threads)

print("Done!")
print("")

jackhmmer_data <- do.call(rbind, raw_jackhmmer_data) %>%
  mutate(target_name = gsub("\\_cds_refseq.*", "", target_name),
         target_name = gsub("\\-.*", "", target_name),
         target_name = gsub("_", "", tolower(target_name)))

# Write out compiled table
print(paste0("Writing out compiled data to ", out_dir, "/", date, "_compiled_jackhmmer_data.csv"))
write.csv(x    = jackhmmer_data,
          file = paste0(out_dir, "/", date, "_compiled_jackhmmer_data.csv"))
print("Done!")

# Remove raw data to save space
rm("raw_jackhmmer_data")

VFDB_metadata        <- read.csv(file = VFDB_meta)
bactibase_metadata   <- read.csv(file = bact_meta)
siderophore_metadata <- read.csv(file = side_meta)
}

# Calculate z-scores -----------------------------------------------------------
print("Calculating mean counts...")
jackhmmer_data <- read.csv(file = "../data/jackhmmer/20230817_compiled_jackhmmer_data.csv")
  
mean_counts <- jackhmmer_data %>%
  group_by(target_name, query_name) %>%
  summarise(n_target = n()) %>%
  pivot_wider(names_from  = query_name,
              values_from = n_target,
              values_fill = 0) %>%
  pivot_longer(cols      = -target_name,
               names_to  = "query_name",
               values_to = "count") %>%
  group_by(query_name) %>%
  summarise(mean = mean(count),
            sd   = sd(count))
			
print("Saving out mean counts...")
write.csv(x    = mean_counts,
          file = paste0(out_dir, "/", date, "_jackhmmer_mean_counts.csv"))

print("Calculating z-scores...")
jackhmmer_z_scores <- jackhmmer_data %>%
  group_by(target_name, query_name) %>%
  summarise(n_target = n()) %>%
  pivot_wider(names_from = query_name,
              values_from = n_target,
              values_fill = 0) %>%
  pivot_longer(cols = -target_name,
               names_to = "query_name",
               values_to = "count") %>%
  left_join(mean_counts,
            by = "query_name") %>%
  mutate(z_score = count - mean / sd)

print("Done!")
print("")
# Split into metadata tables ---------------------------------------------------
bactibase_z_scores <- jackhmmer_z_scores %>%
  right_join(bactibase_metadata,
             by = c("query_name" = "id"))

VFDB_z_scores <- jackhmmer_z_scores %>%
  right_join(VFDB_metadata,
             by = c("query_name" = "VFID"))

siderophore_z_scores <- jackhmmer_z_scores %>%
  right_join(siderophore_metadata,
             by = c("query_name" = "Protein.Referenece"))

# Save outputs -----------------------------------------------------------------
print("Saving bactibase z-scores...")
write.csv(x    = bactibase_z_scores,
          file = paste0(out_dir, date, "_bactibase_z_scores.csv"))

print("Saving VFDB z-scores...")
write.csv(x    = VFDB_z_scores,
          file = paste0(out_dir, date, "_VFDB_z_scores.csv"))

print("Saving siderophore z-scores...")
write.csv(x    = siderophore_z_scores,
          file = paste0(out_dir, date, "_siderophore_z_scores.csv"))
		  
print("Done! Script complete")
