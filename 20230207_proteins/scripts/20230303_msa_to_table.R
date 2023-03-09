# Program which pulls out Uniprot headers from Uniprot fasta files and processes them into a table
# this assumes that the fasta headers conform to uniprot format https://www.uniprot.org/help/fasta-headers

# Dependencies ----
library(tidyverse)

# Functions ----
# wrapper around str_split for selection of string
str_select <- function(x, sep, sel){
  
  tmp <- unlist(str_split(string  = x,
                           pattern = sep))[sel]
  
  
}

# read in msa file as a tbl
msa_to_tbbl <- function(file){
  
  # read in all fasta headers
  raw_file <- readLines(con = file) 
  msa_tbl  <- raw_file[grep(pattern = "^>", 
                            x       = readLines(con = file))]
  msa_tbl  <- msa_tbl[-1]
  
  for(i in 1:length(msa_tbl)){
    
    strng <- unlist(str_split(string  = msa_tbl[i],
                              pattern = " "))
    
    # Process the first section of the string
    strt <- unlist(str_split(string  = strng,
                             pattern = "\\Q|\\E"))
    
    db  <- strt[1]
    acc <- strt[2] 
    nam <- strt[3]
    
    # Get the rest of the information in the header
    if(any(grepl("OS", strng))){
      
      os_strt <- grep("OS", strng)
      os_end  <- grep("OX", strng) - 1 
      
      os      <- paste(strng[os_strt:os_end], collapse = "_")
      os      <- str_select(os, "=", 2)
      
    } else {
      
      os_loc <- 0
      os     <- NA
      
    } 
    
    if(any(grepl("OX", strng))){
      
      ox_loc <- grep("OX", strng)
      ox     <- str_select(strng[ox_loc], "=", 2)
      
    } else {
      
      ox_loc <- 0
      ox     <- NA
      
    } 
    
    if(any(grepl("GN", strng))){
      
      gn_loc <- grep("GN", strng)
      gn     <- str_select(strng[gn_loc], "=", 2)
      
      
    } else {
      
      gn_loc <- 0
      gn     <- NA
      
    } 
    
    if(any(grepl("PE", strng))){
      
      pe_loc <- grep("PE", strng)
      pe     <- str_select(strng[pe_loc], "=", 2)
      
    } else {
      
      pe_loc <- 0
      pe     <- NA
      
    } 
    
    if(any(grepl("SV", strng))){
      
      sv_loc <- grep("SV", strng)
      sv     <- str_select(strng[sv_loc], "=", 2)
      
    } else {
      
      sv_loc <- 0
      sv     <- NA
      
    }
    
    # Get the protein name
    prt <- paste(strng[-c(1, os_strt:os_end, ox_loc, gn_loc, pe_loc, sv_loc)], collapse = "_")
    
    # Compile all the bits together in a data.frame row
    if(exists("output")){
      
      tmp <- data.frame(database       = db,
                        prot_accession = acc,
                        name           = nam,
                        protein        = prt,
                        organism       = os,
                        org_accession  = ox,
                        gene           = gn,
                        prot_existance = pe,
                        seq_ver        = sv)
      
      output <- rbind(output, tmp)
      
      
    } else {
      
      output <- data.frame(database       = db,
                           prot_accession = acc,
                           name           = nam,
                           protein        = prt,
                           organism       = os,
                           org_accession  = ox,
                           gene           = gn,
                           prot_existance = pe,
                           seq_ver        = sv)
      
    }
    
    
  }
  
  
  return(output)
  
}

# Options ----
args = commandArgs(trailingOnly = TRUE)

# msa file
msa_file <- args[1]

# output location
out_loc  <- paste0(args[2])

# Process data ----
msa_data <- msa_to_tbbl(msa_file)

# Save out
write.table(x    = msa_data,
            file = out_loc,
            sep  = "\t")
