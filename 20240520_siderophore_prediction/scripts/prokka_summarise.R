# Summarise all prokka .gff's and bacteria summaries together, including graphs
# Dependencies -----------------------------------------------------------------
library(tidyverse) 
library(parallel)

# Options ----------------------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)

prokka_dir <- sub("/$", "", args[1])
out_dir    <- sub("/$", "", args[2])
threads    <- as.numeric(args[3])

# troubleshooting options
# prokka_dir <- "./data/prokka/"
# threads    <- 1

# Functions --------------------------------------------------------------------
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

# read prokka output statistics
read_prokka_stats <- function(x){
  
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

# Load gff files ---------------------------------------------------------------
writeLines(paste0(Sys.time(), ": Loading gff files..."))

gff_files <- list.files(path       = prokka_dir,
                        pattern    = ".gff",
                        full.names = TRUE,
                        recursive  = TRUE)

writeLines(paste0(Sys.time(), ": ", length(gff_files), " gff files found"))

raw_gff_data <- mclapply(X        = as.list(gff_files),
                         FUN      = read_gff,
                         mc.cores = threads)

# TODO move some of the  regex calculations into read_gff so they are 
# multi threaded
compiled_gff_data <- do.call(rbind, raw_gff_data) %>%
  filter(type == "CDS") %>%
  mutate(organism = basename(dirname(file)),
         gene     = ifelse(test = grepl("gene=", attributes),
                           yes  = gsub("^.*;gene=(.*?);.*$", "\\1", attributes),
                           no   = "Unknown"),
         description = ifelse(test = grepl("product", attributes),
                              yes  = gsub("^.*;product=(.*?)$", "\\1", attributes),
                              no   = "Unknown"),
         ID       = gsub("^ID=(.*?);.*$", "\\1", attributes),
         start2   = ifelse(strand == "+",
                           start,
                           end),
         end2     = ifelse(strand == "+",
                           end,
                           start)) %>%
  rowwise() %>%
  mutate(locus = sub(paste0(organism, "_"), "", ID)) %>%
  transmute(bacteria = organism,
            locus    = locus,
            gene     = gene,
            start    = start2,
            end      = end2,
            strand   = strand,
            description = description)

write.csv(x    = compiled_gff_data,
          file = paste0(out_dir, "/compiled_gffs.csv"),
          row.names = FALSE)

rm(raw_gff_data, gff_files)
gc()

writeLines(paste0(Sys.time(), ": gff files colated"))

# Load summary files -----------------------------------------------------------
writeLines(paste0(Sys.time(), ": Loading summary files..."))

summary_files <- list.files(path       = prokka_dir,
                            pattern    = ".txt",
                            full.names = TRUE,
                            recursive  = TRUE)

writeLines(paste0(Sys.time(), ": ", length(summary_files), " summary files found"))

raw_summary_data <- mclapply(X        = as.list(summary_files),
                             FUN      = read_prokka_stats,
                             mc.cores = threads)

compiled_summary_data <- do.call(rbind, raw_summary_data)
write.csv(x    = compiled_summary_data,
          file = paste0(out_dir, "/compiled_prokka_summarys.csv"),
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
  





