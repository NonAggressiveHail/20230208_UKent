# Script to generate a file containing all of the siderophore genes for a each bacteria,
# in the same order accross all bacteria
# Dependencies -----------------------------------------------------------------
library(tidyverse)
library(parallel)

# Options ----------------------------------------------------------------------
args           <- commandArgs(trailingOnly = TRUE)
blast_file     <- args[1]
out_dir        <- sub("/$", "", args[2])
threads        <- as.numeric(args[3])

# options when running locally
#blast_file     <- "./data/blast/blastp_annotations.csv"
#out_dir        <- "./data/sid_gene_per_organism"
#threads        <- 1

# Input error checking ---------------------------------------------------------
if(!file.exists(blast_file)){
  stop(paste0("blast_file (", blast_file, ") does not exist!"))
}

if(!dir.exists(out_dir)){
  stop(paste0("out_dir (", out_dir, ") does not exist!"))
}

# Load data --------------------------------------------------------------------
blastp_annotations <- read.csv(file = blast_file)

# extract relevant data --------------------------------------------------------
## Get vector of organisms
organisms <- unique(blastp_annotations$bacteria)

writeLines(paste0(length(organisms), " bacteria found\n"))

## Get vector of siderophore_genes
sid_genes <- unique(blastp_annotations$locus_name_in_source)

writeLines(paste0(length(sid_genes), " siderophore genes present\n"))

# Extract loci             -----------------------------------------------------
if(FALSE){
counter <- 0
for(organism in organisms){
  
  # Make output dir
  organism_dir <- paste0(out_dir, organism)
  dir.create(path = organism_dir)
  
  # subset data
  organsim_siderophore_loci <- blastp_annotations %>%
    filter(bacteria == UQ(organism)) %>%
    select(locus_name_in_source, query)
  
  # save each subset out for seqkit to then be run on
  for(sid_gene in sid_genes){
    
    loci_of_sideropore <- organsim_siderophore_loci %>%
      filter(locus_name_in_source == UQ(sid_gene)) %>%
      select(query) 
    
    write_delim(x         = loci_of_sideropore,
                file      = paste0(organism_dir,"/", sid_gene, "_loci_ids.txt"),
                col_names = FALSE)

  }
  
  # Progress counter
  counter      <- counter + 1
  per_complete <- counter/length(organisms) * 100
  writeLines(paste0(round(per_complete, 2), "% complete"))
}
}

# Multithread
# extraction_function
x <- split(x = blastp_annotations,
           f = ~ bacteria)[[1]]
sid_loci_name_extract <- function(x){
  
  organism <- unique(x$bacteria)
  writeLines(paste0("Writing out IDs for ", organism))
  
  
  # Make output dir
  organism_dir <- paste0(out_dir, "/", organism)
  dir.create(path = organism_dir)
  
  # save each subset out for seqkit to then be run on
  for(sid_gene in sid_genes){
    
    loci_of_sideropore <- x %>%
      filter(locus_name_in_source == UQ(sid_gene)) %>%
      select(target) 
    
    write_delim(x         = loci_of_sideropore,
                file      = paste0(organism_dir,"/", sid_gene, "_loci_ids.txt"),
                col_names = FALSE)
    
  }
  
}

# Split into list
blastp_annotations_list <- split(x = blastp_annotations,
                                 f = ~ bacteria)

# Extract names 
mclapply(X   = blastp_annotations_list,
         FUN = sid_loci_name_extract,
         mc.cores = threads)
