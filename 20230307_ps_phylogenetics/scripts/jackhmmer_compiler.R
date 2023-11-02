# Script to calcualte z-scores of genes of interest for my bacteria
# and export z-scores for plotting on heatmap
library(tidyverse)
library(parallel)
library(moments)
library(rhmmer)

# Options ----------------------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)

jackhmmer_dir   <- args[1]
gff_dir         <- args[2]
out_dir         <- args[3]
side_meta       <- args[4]
VFDB_meta       <- args[5]
bact_meta       <- args[6]
full_names_file <- args[7]
threads         <- as.numeric(args[8])

writeLines(paste0(Sys.time(), ": Args are:"))
writeLines(args)
writeLines("")

# Functions --------------------------------------------------------------------
partial_join <- function(x, y, by_x, pattern_y){
  
  idx_x <- sapply(y[[pattern_y]], grep, x[[by_x]])
  idx_y <- sapply(seq_along(idx_x), function(i) rep(i, length(idx_x[[i]])))
  
  df <- dplyr::bind_cols(x[unlist(idx_x), , drop = F],
                         y[unlist(idx_y), , drop = F])
  return(df)
}

mcpartial_join <- function(x, y, by_x, pattern_y, cores){
  
  idx_x <- mclapply(y[[pattern_y]], grep, x[[by_x]], mc.cores = cores)
  idx_y <- mclapply(seq_along(idx_x), function(i) rep(i, length(idx_x[[i]])), mc.cores = cores)
  
  
  df <- dplyr::bind_cols(x[unlist(idx_x), , drop = F],
                         y[unlist(idx_y), , drop = F])
  return(df)
}

string_extract <- function(str, split, pattern){
  
  split_string <- unlist(str_split(str, split))
  
  string_output <- split_string[grep(pattern, split_string)]
  
  return(string_output)
  
}

read_gff <- function(x){
  
  # Load file 
  temp_file <- readLines(x)
  
  # Find line the fasta info starts at 
  fasta_line <- grep("^##FASTA$", temp_file)
  
  # Select text before the fasta line
  drop_fasta <- temp_file[1:fasta_line]
  
  # Drop remaining lines with ## in them
  drop_header <- drop_fasta[!grepl("^##", drop_fasta)]
  
  output <- read.delim(file   = textConnection(drop_header),
                       header = FALSE)
  
  colnames(output) <- c("seqname", "source", "type",
                        "start", "end", "score", "strand",
                        "phase", "attributes")
  
  output[, "file"] <- x
  
  return(output)
}

# Build metadata ---------------------------------------------------------------
writeLines(paste0(Sys.time(), ": Loading Metadata..."))

# Load in protein metadata
VFDB_meta_trim <-  read.csv(file = VFDB_meta) %>%
  select(gene_id, VFID, VFCID, VFID_name, VFCID_name ) %>%
  transmute(gene_id        = gene_id,
            factor_id      = VFID,
            factor_id_name = VFID_name,
            factor_category = VFCID_name,
            factor_category_name = VFCID_name,
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

all_meta_trim <- do.call(rbind, list(VFDB_meta_trim, bacti_meta_trim, sidero_meta_trim))

# Load gff files for loci co-ordinates -----------------------------------------
writeLines(paste0(Sys.time(), ": Loading gff files..."))

gff_files <- list.files(path       = gff_dir,
                        pattern    = ".gff",
                        full.names = TRUE,
                        recursive  = TRUE)


writeLines(paste0(Sys.time(), ": ", length(gff_files), " gff files found"))

raw_gff_data <- mclapply(X        = as.list(gff_files),
                         FUN      = read_gff,
                         mc.cores = threads)


compiled_gff_data <- do.call(rbind, raw_gff_data) %>%
  filter(type == "CDS") %>%
  mutate(organism = basename(dirname(file)),
         locus    = gsub(".*_(\\d+);.*", "\\1", attributes)) %>%
  select(organism, locus, start, end)


rm(raw_gff_data, gff_files)
gc()

# Load Jackhmmer Data ----------------------------------------------------------
date     <- gsub("-", "", Sys.Date())
numcores <- detectCores()

writeLines(paste0("Finding files in ", jackhmmer_dir))

target_organisms      <- read.csv(file = full_names_file)

jackhmmer_files <- list.files(path       = jackhmmer_dir,
                              pattern    = "*_tbl.txt",
                              full.names = TRUE)

writeLines(paste0(Sys.time(), ": Found ", length(jackhmmer_files), " files"))
writeLines("")
writeLines(paste0(Sys.time(), ": Reading in files with ", threads, " threads, out of ", numcores, " available"))

# If already compiled data exists, load it
if(file.exists(paste0(out_dir, "/",  "compiled_jackhmmer_data.csv"))){
  
  writeLines(paste0(Sys.time(), ": Previous compiled jackhmmer data found, using"))
  compiled_jackhmmer_data <- read.csv(paste0(out_dir, "/", "compiled_jackhmmer_data.csv"))
  
  # Otherwise compile it
} else {
  
  
  writeLines(paste0(Sys.time(), ": Previous compiled jackhmmer data not found, compiling"))
  
  raw_jackhmmer_data <- mclapply(X        = jackhmmer_files,
                                 FUN      = read_tblout,
                                 mc.cores = threads)
  
  compiled_jackhmmer_data <- do.call(rbind, raw_jackhmmer_data) %>%
    select(domain_name, query_name, sequence_evalue)
  
  writeLines(paste0(Sys.time(), ": Writing out compiled data to ", out_dir, "/", "compiled_jackhmmer_data.csv"))
  
  write.csv(x         = compiled_jackhmmer_data,
            file      = paste0(out_dir, "/", "compiled_jackhmmer_data.csv"),
            row.names = FALSE)
  
}

writeLines(paste0(Sys.time(), ": Done!"))
writeLines(paste0(Sys.time(), ": jackhmmer data has ", nrow(compiled_jackhmmer_data), " lines"))
writeLines("")

# Filter out bad hits
writeLines(paste0(Sys.time(), ": Removing hits with e-value <= 1e-03"))
sig_jackhmmer_hits <- compiled_jackhmmer_data %>%
  filter(sequence_evalue <= 1e-03)

writeLines(paste0(Sys.time(), ": ", nrow(compiled_jackhmmer_data) - nrow(sig_jackhmmer_hits), " were not significant"))

# Convert names and remove duplicated data
if(FALSE){
  jackhmmer_data <- sig_jackhmmer_hits %>%
    mutate(target_organism = sub("_\\d+~.*$", "", domain_name),                # extract the target PA strain
           target_gene     = sub("([A-Za-z0-9_]+)~", "", domain_name),
           locus_id        = gsub(".*_(\\d+)~.*", "\\1", domain_name)) %>%                # extract the locus in that strain
    mutate(query_gene = ifelse(test = grepl("^VFG", query_name),
                               yes  = sub("\\(.*", "", query_name),
                               no   = sub("\\~.*", "", query_name))) %>%
    left_join(all_meta_trim,
              by = c("query_gene" = "gene_id")) %>%
    group_by(target_organism, target_gene, type) %>%
    filter(sequence_evalue == min(sequence_evalue)) %>%                      # Only keep the sequences with the lowest e_value
    reframe(query_gene     = query_gene,
            locus_id       = locus_id,
            hits_for_locus = n()) %>%                                        # Calculate how many total queries hit that locus
    group_by(target_organism, target_gene, query_gene, type)  %>%                 
    reframe(hits_for_locus  = hits_for_locus,
            locus_id       = locus_id,
            match_quality   = n()/hits_for_locus * 100) %>%                  # Calculate the number of times each query hit
    group_by(target_organism, target_gene, type) %>%
    filter(match_quality == max(match_quality)) %>%
    group_by(target_organism, target_gene, type)  %>%    
    left_join(all_meta_trim,
              by = c("query_gene" = "gene_id",
                     "type"       = "type")) %>%
    reframe(target_organism = target_organism,
            locus_id        = locus_id,
            target_gene     = target_gene,
            hits_for_locus  = hits_for_locus,
            query_gene = ifelse(length(unique(query_gene)) == 1,
                                query_gene,
                                "Unclear"),
            factor_id = ifelse(length(unique(factor_id)) == 1,
                               factor_id,
                               "Unclear"),
            factor_id_name = ifelse(length(unique(factor_id_name)) == 1,
                                    factor_id_name,
                                    "Unclear"),
            factor_category = ifelse(length(unique(factor_category)) == 1,
                                     factor_category,
                                     "Unclear"),
            factor_category_name = ifelse(length(unique(factor_category_name)) == 1,
                                          factor_category_name,
                                          "Unclear"),
            type                 = type) %>%
    distinct() %>%
    left_join(compiled_gff_data,
              by = c("target_organism" = "organism",
                     "locus_id"        = "locus"))
}

# Try multithreading some of above
if(TRUE){
  temp_jackhmmer_data <- sig_jackhmmer_hits %>%
    mutate(target_organism = sub("_\\d+~.*$", "", domain_name),                # extract the target PA strain
           target_gene     = sub("([A-Za-z0-9_]+)~", "", domain_name),
           locus_id        = gsub(".*_(\\d+)~.*", "\\1", domain_name)) %>%                # extract the locus in that strain
    mutate(query_gene = ifelse(test = grepl("^VFG", query_name),
                               yes  = sub("\\(.*", "", query_name),
                               no   = sub("\\~.*", "", query_name))) %>%
    left_join(all_meta_trim,
              by = c("query_gene" = "gene_id"))
  
  # Function to use in lapply to speed up
  find_best_match <- function(x){
    
    output <- x  group_by(target_organism, target_gene, type) %>%
      filter(sequence_evalue == min(sequence_evalue)) %>%                      # Only keep the sequences with the lowest e_value
      reframe(query_gene     = query_gene,
              locus_id       = locus_id,
              hits_for_locus = n()) %>%                                        # Calculate how many total queries hit that locus
      group_by(target_organism, target_gene, query_gene, type)  %>%                 
      reframe(hits_for_locus  = hits_for_locus,
              locus_id       = locus_id,
              match_quality   = n()/hits_for_locus * 100) %>%                  # Calculate the number of times each query hit
      group_by(target_organism, target_gene, type) %>%
      filter(match_quality == max(match_quality)) %>%
      group_by(target_organism, target_gene, type)  %>%    
      left_join(all_meta_trim,
                by = c("query_gene" = "gene_id",
                       "type"       = "type")) %>%
      reframe(target_organism = target_organism,
              locus_id        = locus_id,
              target_gene     = target_gene,
              hits_for_locus  = hits_for_locus,
              query_gene = ifelse(length(unique(query_gene)) == 1,
                                  query_gene,
                                  "Unclear"),
              factor_id = ifelse(length(unique(factor_id)) == 1,
                                 factor_id,
                                 "Unclear"),
              factor_id_name = ifelse(length(unique(factor_id_name)) == 1,
                                      factor_id_name,
                                      "Unclear"),
              factor_category = ifelse(length(unique(factor_category)) == 1,
                                       factor_category,
                                       "Unclear"),
              factor_category_name = ifelse(length(unique(factor_category_name)) == 1,
                                            factor_category_name,
                                            "Unclear"),
              type                 = type) %>%
      distinct() %>%
      left_join(compiled_gff_data,
                by = c("target_organism" = "organism",
                       "locus_id"        = "locus"))
    
    return(output)
  }
  
  jackhmmer_data <- do.call(what = rbind, 
                            args = mclapply(X        = split(temp_jackhmmer_data, f = ~ target_organism),
                                            FUN      = find_best_match,
                                            mc.cores = threads))
  
  rm(temp_jackhmmer_data)
  gc()
}

n_query_orgs <- length(unique(jackhmmer_data$target_organism))
writeLines(paste0(Sys.time(), ": ", n_query_orgs, " query_organisms found"))

# Write out compiled table
writeLines(paste0(Sys.time(), ": Writing out compiled data with correct names to ", out_dir, "/", "compiled_jackhmmer_data_corrected_names.csv"))
write.csv(x    = jackhmmer_data,
          file = paste0(out_dir, "/", "compiled_jackhmmer_data_corrected_names.csv"),
          row.names = FALSE)
writeLines(paste0(Sys.time(), ": Done!"))