# compute summary statistics and makes some graphs of prokka predictions search
# Dependencies -----------------------------------------------------------------
library(tidyverse)
library(parallel)

# Setup ------------------------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)

out_dir         <- args[1]
prokka_dir      <- args[2]
side_meta       <- args[3]
full_names_file <- args[4]
threads         <- as.numeric(args[5])


#out_dir         <- "../data/prokka/"
#prokka_dir      <- "../data/prokka/"
#side_meta       <- "../raw_data/siderophore_metadata.csv"
#full_names_file <- "../raw_data/all_names.txt"
#threads         <- 1

# Functions --------------------------------------------------------------------
partial_join <- function(x, y, by_x, pattern_y){
  
  idx_x <- sapply(y[[pattern_y]], grep, x[[by_x]])
  idx_y <- sapply(seq_along(idx_x), function(i) rep(i, length(idx_x[[i]])))
  
  df <- dplyr::bind_cols(x[unlist(idx_x), , drop = F],
                         y[unlist(idx_y), , drop = F])
  return(df)
}

string_extract <- function(str, split, pattern){
  
  split_string <- unlist(str_split(str, split))
  
  string_output <- split_string[grep(pattern, split_string)]
  
  return(string_output)
  
}

tsv_sample_read <- function(path){
  
  read.delim(file = path,
             sep  = "\t",
             quote = "") %>%
    mutate(filename = basename(path))
  
}

mc_tsv_compile <- function(path, pattern, threads){
  
  files <- as.list(list.files(path       = path,
                              full.names = TRUE,
                              pattern    = pattern,
                              recursive  = TRUE))
  
  print(paste0(Sys.time(), ": Found ", length(files), " files"))
  
  output <- do.call(rbind, mclapply(X        = files,
                                    FUN      = tsv_sample_read,
                                    mc.cores = threads))
  
  
} 

tsv_compile2 <- function(path, pattern){
  
  files <- as.list(list.files(path       = path,
                              full.names = TRUE,
                              pattern    = pattern))
  
  output <- do.call(rbind, lapply(files, tsv_sample_read))
  
  
} 

# Load data --------------------------------------------------------------------
print(paste0(Sys.time(), ": Args are:"))
print(args)
print("")

date       <- gsub("-", "", Sys.Date())
numcores   <- detectCores()
full_names <- read.csv(file = full_names_file)


print(paste0(Sys.Date(), ": Loading files in ", prokka_dir))


prokka_data <- mc_tsv_compile(path    = prokka_dir,
                              pattern = "*.tsv",
                              threads = threads) %>%
  filter(ftype == "CDS") %>%
  filter(gene != "") %>%
  mutate(full_name = gsub(".tsv", "", filename)) %>%
  select(full_name, gene)

print(paste0(Sys.time(), ": Done!"))
print(paste0(Sys.time(), ": prokka data has ", nrow(prokka_data), " lines"))
print("")


# Write out raw table
print(paste0(Sys.time(), ": Writing out raw data with correct names to ", out_dir, "/", "raw_prokka_data_corrected_names.csv"))
write.csv(x    = prokka_data,
          file = paste0(out_dir, "/", "raw_prokka_data_corrected_names.csv"),
          row.names = FALSE)
print(paste0(Sys.time(), ": Done!"))

# Load siderophore metadata ----------------------------------------------------
siderophore_metadata <- read.csv(file = side_meta) %>%
  select(Gene.name)

# Make table of all protein:genome combinations
temp_function <- function(x){
  
  output <- siderophore_metadata %>%
    mutate(bacteria = x)
  
  return(output)
}

tmp_full_table <- mclapply(X        = as.list(unique(prokka_data$full_name)),
                           FUN      = temp_function,
                           mc.cores = threads)

full_table <- do.call(rbind, tmp_full_table)

# Get mean counts of siderophore genes per genome --------------------------------------
print(paste0(Sys.time(), ": Calculating counts of siderophore genes per genome"))

mc_group <- function(x){
  
  output <- x %>%
    group_by(full_name, gene) %>%
    summarise(n = n(),
              .groups = "drop_last")
  
  return(output)
}

per_genome_mean_counts <- do.call(what = rbind, 
                                  args = mclapply(X        = split(prokka_data, f = ~ full_name),
                                                  FUN      = mc_group,
                                                  mc.cores = threads)) %>%
  filter(gene %in% siderophore_metadata$Gene.name) %>%
  right_join(full_table,
             by = c("gene"      = "Gene.name",
                    "full_name" = "bacteria")) %>%
  replace_na(list(n = 0))

write.csv(x         = per_genome_mean_counts,
          file      = paste0(out_dir, "/", "siderophore_gene_counts_per_genome.csv"),
          row.names = FALSE)

print(paste0(Sys.time(), ": Done!"))

# Plot results
## Calculate order of genes
gene_order <- per_genome_mean_counts %>%
  group_by(gene) %>%
  summarise(mean_count = mean(n)) %>%
  arrange(desc(mean_count)) %>%
  pull(gene)

per_genome_mean_counts[, "gene"] <- factor(per_genome_mean_counts$gene, levels = gene_order)

# Plot counts of gene
gene_count_graph <- ggplot(data    = per_genome_mean_counts,
                           mapping = aes(x     = gene,
                                         y     = n)) +
  geom_boxplot(show.legend = TRUE) +
  scale_x_discrete(guide = guide_axis(n.dodge = 1),
                   name  = "") +
  scale_y_continuous(name = "Count") +
  labs(title = paste0("Counts of siderophore metabolism genes"),
       caption = "Predicted by Prokka") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5,
                                   size  = 5),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

ggsave(filename = paste0(date, "_", "siderophore_metabolism_count_graph.png"),
       plot     = gene_count_graph,
       path     = paste0(out_dir, "/graphs"),
       width    = 15,
       height   = 15,
       units    = "cm",
       limitsize = FALSE)

print(paste0(Sys.time(), ": Done!"))

# Counts of things in PAO1
print(paste0(Sys.time(), ": Making PAO1 graph!"))

PAO1s <- grep("PAO1", unique(per_genome_mean_counts$full_name), value = TRUE)
print(paste0(Sys.time(), ": These Genomes are present:")) 
print(unique(per_genome_mean_counts$full_name))

print(paste0(Sys.time(), ": These PAO1s are present:")) 
print(PAO1s)

if(any(grepl("Pseudomonas_aeruginosa_PAO1_107", unique(per_genome_mean_counts$full_name)))){
  
  print(paste0(Sys.time(), ": Pseudomonas_aeruginosa_PAO1_107 present:")) 
  
}
PAO1_graph <- ggplot(data = per_genome_mean_counts %>%
                       filter(full_name == "Pseudomonas_aeruginosa_PAO1_107"),
                     mapping = aes(x = gene,
                                   y = n)) +
  geom_col()+
  labs(title = "Genes in PAO1") +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Count\n") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5,
                                   size  = 5),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

ggsave(filename = paste0(date, "_", "PAO1_count_graph.png"),
       plot     = PAO1_graph,
       path     = paste0(out_dir, "/graphs"),
       width    = 15,
       height   = 15,
       units    = "cm",
       limitsize = FALSE)

rm(per_genome_mean_counts)
gc()

print(paste0(Sys.time(), ": Script Complete!"))


