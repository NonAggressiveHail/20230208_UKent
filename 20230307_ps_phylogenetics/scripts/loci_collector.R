# Script to parse "./data/jackhmmer/compiled_jackhmmer_data_corrected_names.csv" ready for orthomcl
# Dependencies -----------------------------------------------------------------
library(tidyverse)
library(parallel)

# Options ----------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
jackhmmer_file <- args[1]
out_dir        <- args[2]
threads        <- as.numeric(args[3])

# Input error checking ---------------------------------------------------------
if(!file.exists(jackhmmer_file)){
  stop(paste0("jackhmmer_file (", jackhmmer_file, ") does not exist!"))
}

if(!dir.exists(out_dir)){
  stop(paste0("out_dir (", out_dir, ") does not exist!"))
}

# Load data --------------------------------------------------------------------
compiled_jackhmmer_data <- read.csv(file = jackhmmer_file)

# extract relevant data --------------------------------------------------------
## Filter to just siderophore data
sid_data <- compiled_jackhmmer_data %>%
  filter(type == "Siderophore_Biosynthesis")

## Get vector of organisms
organisms <- unique(sid_data$target_organism)

## Get vector of siderophore_genes
sid_genes <- unique(sid_data$query_gene)

# Extract loci             -----------------------------------------------------
counter <- 0
for(organism in organisms){
  
  # Make output dir
  organism_dir <- paste0(out_dir, organism)
  dir.create(path = organism_dir)
  
  # subset data
  organsim_siderophore_loci <- sid_data %>%
    filter(target_organism == UQ(organism)) %>%
    select(locus_id, query_gene)
  
  # save each subset out for seqkit to then be run on
  for(sid_gene in sid_genes){
    
    loci_of_sideropore <- organsim_siderophore_loci %>%
      filter(query_gene == UQ(sid_gene)) %>%
      mutate(fmt_id = sprintf(fmt = "%05d", locus_id),
             seq_id = paste0(organism, "_", fmt_id)) %>%
      select(seq_id) 
    
    write_delim(x    = loci_of_sideropore,
                file = paste0(organism_dir,"/", sid_gene, "_loci_ids.txt"),
                col_names = FALSE)

  }
  
  # Progress counter
  counter      <- counter + 1
  per_complete <- counter/length(organisms) * 100
  writeLines(paste0(round(per_complete, 2), "% complete"))
}
