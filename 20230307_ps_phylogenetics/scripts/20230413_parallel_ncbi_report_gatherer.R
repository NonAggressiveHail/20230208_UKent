# script to pull information out of NCBI assembly_report files and reformat as
# a nice table for easy access
# Dependencies ----
library(parallel)

# Options ----
args = commandArgs(trailingOnly = TRUE)

path     <- args[1]
out_file <- args[2]
cores    <- args[3]


# Functions ----
str_select <- function(string, sep, pos){
  
  return(unlist(strsplit(string, sep))[[pos]])
  
}

report_str_extract <- function(pattern, file){
  
  if(any(grepl(pattern, file))){
    
    tmp_string   <- str_select(file[grep(pattern, file)], ":", 2)
    clean_string <- gsub(" ", "_", trimws(tmp_string), "both")
    
  } else {
    
    clean_string <- NA
    
  }
  
  
  return(clean_string)
  
}

report_gather <- function(report_file){
  
  filename      <- basename(report_file)
  assembly_name <- gsub("_assembly_report.txt", "", filename)
  
  
  # Remove "#" from start of lines in file
  temp_file   <- gsub("#", "", readLines(con = report_file), perl = TRUE)
  
  # Get organism name
  org <- report_str_extract("Organism", temp_file)
  
  # Get Infra specific name
  inf <- report_str_extract("Infraspecific", temp_file)
  
  # Get refseq accession
  ref_acc <- report_str_extract("RefSeq assembly accession", temp_file)
  
  # Get Genbank accession
  gen_acc <- report_str_extract("GenBank assembly accession", temp_file)
  
  # TaxID 
  taxid <- report_str_extract("Taxid", temp_file)
  
  # Get strain info if present
  if(grepl("strain", inf)){
    
    strain <- str_select(inf, "=", 2)
    
  } else {
    
    strain <- NA
    
  }
  
  report_df <- data.frame("assembly"           = assembly_name,
                          "organsim"           = org,
                          "infraspecific_name" = inf,
                          "strain"             = strain,
                          "refseq_accession"   = ref_acc,
                          "genbank_accession"  = gen_acc,
                          "taxid"              = taxid)
  
  return(report_df)
  
  
  
}

# Data ----
report_files <- as.list(paste0(path, 
                       list.files(path      = path,
                                  pattern   = "*assembly_report.txt",
                                  recursive = TRUE)))

# Main Program ----
report_list <- mclapply(report_files, report_gather, mc.cores = cores)

report_df   <- as.data.frame(do.call(rbind, report_list))


# save out full table
write.table(x         = report_df,
            file      = out_file,
            sep       = "\t",
            row.names = FALSE)









