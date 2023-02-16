# script to produce an itol annotation file from a fasta format msa file
# this assumes that the fasta headers conform to uniprot format https://www.uniprot.org/help/fasta-headers
#TODO redo with DATASET_COLORSTRIP type, using dataset_color_strip_template.txt

# Dependencies ----
library(tidyverse)
library(RColorBrewer)

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

# Generate qualitative pallett of n colours from RCOlorBrewer
gen_pal <- function(n, pie){
  
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  big_pal       <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  pall_n        <- sample(big_pal, n)
  
  if(pie){
    
    print(pie(rep(1, n), col = sample(big_pal, n)))
  
  }
  
  return(pall_n)
}

# Data ----
msa_data <- msa_to_tbbl("uniprot_matches_msa.afa")

# generate colour pallette
set.seed(1)
pall_n <- gen_pal(n   = length(unique(msa_data[, "organism"])),
                  pie = TRUE)

# Make Labels file
write_file(x    = "LABELS\nSEPARATOR SPACE\nDATA\n",
           file = "WSC4_Labels.txt")

labs <-  msa_data %>%
  replace(is.na(.), "Unknown") %>% 
  transmute(node_id  = gsub(" |>", "", paste(database, "|", prot_accession, "|",name)),
            label    = gsub(" ", "", paste(name, "_", protein)))

write.table(x         = labs,
            file      = "WSC4_Labels.txt",
            append    = TRUE,
            sep       = " ",
            quote     = FALSE,
            row.names = FALSE,
            col.names = FALSE)



# Make organism colours file
annot <- msa_data %>%
  replace(is.na(.), "Unknown") %>% 
  transmute(node_id  = gsub(" |>", "", paste(database, "|", prot_accession, "|",name)),
            color    = "",
            organism = organism) %>%
  mutate(color       = factor(organism,
                              labels = pall_n))

write_file(x    = paste("DATASET_COLORSTRIP",
                        "SEPARATOR SPACE",
                        "DATASET_LABEL Organism",
                        "COLOR #ff0000",
                        "COLOR_BRANCHES 1",
                        "STRIP_WIDTH 25",
                        "SHOW_STRIP_LABELS 0",
                        "LABEL_ROTATION 90",
                        "LEGEND_TITLE Organism",
                        paste0("LEGEND_LABELS ", paste(unique(annot[, "organism"]), collapse = " ")),
                        paste0("LEGEND_COLORS ", paste(unique(annot[, "color"]), collapse = " ")),
                        paste0("LEGEND_SHAPES ", paste(rep(1, length(unique(annot[, "organism"]))), collapse = " ")),
                        "DATA\n",
                        sep = "\n"),
           file = "WSC4_org_colors.txt")



write.table(x         = annot,
            file      = "WSC4_org_colors.txt",
            append    = TRUE,
            sep       = " ",
            quote     = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# make MSA file
msa <- readLines(con = "uniprot_matches_msa.afa")

write_file(x    = paste("DATASET_ALIGNMENT",
                        "SEPARATOR SPACE",
                        "DATASET_LABEL jackhmmer MSA",
                        "START_POSITION 1",
                        "END_POSITION 500",
                        "DATA",
                        paste(print(msa, quote = FALSE), collapse = "\n"),
                        sep ="\n"),
           file = "WSC4_MSA.txt")






















