# Summarise all bakta .tsv's and summaries together, make graphs
# Dependencies -----------------------------------------------------------------
library(tidyverse) 
library(parallel)

# Options ----------------------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)

bakta_dir <- sub("/$", "", args[1])
out_dir   <- sub("/$", "", args[2])
threads   <- as.numeric(args[3])

# troubleshooting options
# bakta_dir <- "./data/bakta/"
# threads    <- 1

# Functions --------------------------------------------------------------------
# read bakta output statistics
read_bakta_stats <- function(x){
  
  # Get bacteria name
  bacteria <- gsub("^.*/(.*).txt", "\\1", x)
  
  # Load file
  temp_file <- read.delim(file      = x,
                          skip      = 1, 
                          sep       = ":",
                          header    = FALSE,
                          col.names = c("Feature", "Count")) %>%
    mutate(Bacteria = bacteria) %>%
    relocate("Bacteria", "Feature", "Count")
  
  return(temp_file)
  
}

# Check that each output folder contains the right outputs ---------------------
writeLines(paste0(Sys.time(), ": Checking .tsv files exist..."))

output_folders <- list.dirs(path       = bakta_dir,
                               full.names = TRUE,
                               recursive  = FALSE)

writeLines(paste0(Sys.time(), ": ", length(output_folders), " bakta output folders found"))

tsv_files <- list.files(path       = bakta_dir,
                        pattern    = ".tsv",
                        full.names = TRUE,
                        recursive  = TRUE)

## we don't want hypotheticals, as that is additinal info on the hypotheticals
## all of them are in the normal tsv but this gives extra data on them
protein_tsv_files <- tsv_files[!grepl("hypotheticals", tsv_files)]

writeLines(paste0(Sys.time(), ": ", length(protein_tsv_files), " bakta output tsv files found"))

missing_files <- data.frame(tsv_files = protein_tsv_files) %>%
  mutate(folder_tsv_in = sub("/$", "", sub("[^/]+$", "", tsv_files))) %>%
  full_join(y  = data_frame(tmp           = output_folders,
                            bakta_folders = output_folders),
            by = c("folder_tsv_in" = "tmp")) %>%
  filter(is.na(tsv_files)) %>%
  pull(bakta_folders)

if(length(missing_files) != 0){
  
  writeLines(paste0(Sys.time(), ": Folders missing tsv files :\n ", paste(missing_files, collapse = "\n")))
  
}

# Load tsv files ---------------------------------------------------------------
writeLines(paste0(Sys.time(), ": Loading .tsv files..."))




writeLines(paste0(Sys.time(), ": ", length(protein_tsv_files), " tsv files found"))

raw_tsv_data <- mclapply(X        = as.list(protein_tsv_files),
                         FUN      = read.delim,
                         skip     = 5,
                         mc.cores = threads)

# TODO move some of the  regex calculations into read_tsv so they are 
# multi threaded
compiled_tsv_data <- do.call(rbind, raw_tsv_data) %>%
  rename(Sequence  = X.Sequence.Id,
         Locus_Tag = Locus.Tag) %>%
  mutate(Gene = ifelse(Type == "cds" & Gene == "",
                       "Unknown",
                       Gene),
         Organism = sub(pattern     = "_\\d{5}$",
                        replacement = "",
                        x           = Locus_Tag),
         Locus    = sub(pattern     = ".*?_(\\d{5})$",
                        replacement = "\\1",
                        x           = Locus_Tag))

write.csv(x    = compiled_tsv_data,
          file = paste0(out_dir, "/compiled_tsvs.csv"),
          row.names = FALSE)

rm(raw_tsv_data, tsv_files)
gc()

writeLines(paste0(Sys.time(), ": tsv files colated"))

# Load summary files -----------------------------------------------------------
writeLines(paste0(Sys.time(), ": Loading summary files..."))

summary_files <- list.files(path       = bakta_dir,
                            pattern    = ".txt",
                            full.names = TRUE,
                            recursive  = TRUE)

writeLines(paste0(Sys.time(), ": ", length(summary_files), " summary files found"))

raw_summary_data <- mclapply(X        = as.list(summary_files),
                             FUN      = read_bakta_stats,
                             mc.cores = threads)

compiled_summary_data <- do.call(rbind, raw_summary_data) %>%
  filter(Count != "") %>%
  filter(Feature != "DOI") %>%
  filter(Feature != "URL") %>%
  filter(Feature != "Software") %>%
  filter(Feature != "Database") %>%
  mutate(Count = as.numeric(Count))

write.csv(x    = compiled_summary_data,
          file = paste0(out_dir, "/compiled_bakta_summarys.csv"),
          row.names = FALSE)

rm(raw_summary_data, summary_files)
gc()

writeLines(paste0(Sys.time(), ": summary files colated"))

# Do some summary graphs -------------------------------------------------------
feature_freq_plot <- compiled_summary_data %>%
  group_by(Feature, Count) %>%
  summarise(count_freq = n())

summary_file_plots <- ggplot(data    = feature_freq_plot,
                             mapping = aes(x = Count,
                                           y = count_freq)) +
  geom_col() +
  scale_x_continuous(name   = "\nCount") +
  scale_y_continuous(name   = "Frequency") +
  facet_wrap(facets = vars(Feature),
             scale  = "free") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0,
                                    vjust = 0.5))

ggsave(filename = paste0(gsub("-", "", Sys.Date()), "bacteria_feature_summaries.png"),
       path     = paste0(out_dir, "/graphs"),
       plot     = summary_file_plots,
       device   = "png",
       dpi      = 200,
       width    = 30,
       height   = 15,
       units    = "cm")
  





