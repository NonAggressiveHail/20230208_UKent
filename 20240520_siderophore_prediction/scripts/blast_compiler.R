# compute summary statistics and maks some graphs of blast search
# Dependencies -----------------------------------------------------------------
library(tidyverse)
library(parallel)
library(moments)

# Options  ---------------------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)

out_dir         <- args[1]
blast_output    <- args[2]
side_meta       <- args[3]
VFDB_meta       <- args[4]
bact_meta       <- args[5]
full_names_file <- args[6]
e_val_cutoff    <- as.numeric(args[7])
threads         <- as.numeric(args[8])

#out_dir        <- "./data/blast/"
#blast_output   <- "./data/blast/blastp_out_old.csv"
#side_meta      <- "./raw_data/siderophore_metadata.csv"
#VFDB_meta      <- "./raw_data/VFDB_metadata.csv"
#bact_meta      <- "./raw_data/bactibase_metadata.csv"
#full_names_file <- "./raw_data/all_names.txt"
#e_val_cutoff   <- 1e-03
#threads        <- 1

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


# Load data --------------------------------------------------------------------
print(paste0(Sys.time(), ": Args are:"))
print(args)
print("")

date       <- gsub("-", "", Sys.Date())
numcores   <- detectCores()
full_names <- read.csv(file = full_names_file)


print(paste0(Sys.Date(), ": Loading ", blast_output))

raw_blast_data <- read.csv(file   = blast_output,
                           header = FALSE)

colnames(raw_blast_data) <- c("query", "target", "per_ident", "length",
                              "mismatch", "gapopen", "qstart", "qend",
                              "sstart", "send", "evalue", "bitscore")



print(paste0(Sys.time(), ": Done!"))
print(paste0(Sys.time(), ": blast data has ", nrow(raw_blast_data), " lines"))
print("")

# Filter out bad hits
print(paste0(Sys.time(), ": Removing hits with e-value <= 1e-03"))
sig_blast_hits <- raw_blast_data %>%
  filter(evalue <= e_val_cutoff)

print(paste0(Sys.time(), ": ", nrow(raw_blast_data) - nrow(sig_blast_hits), " were not significant"))



# Need to do this partial join shenannigans to convert genome names to what I usually use
print(paste0(Sys.time(), ": Converting names..."))
blast_names <- raw_blast_data %>%
  select(target) %>%
  distinct()


names_conversion <- mclapply(X         = split(blast_names, 1:threads),
                             FUN       = partial_join,
                             y         = full_names,
                             by_x      = "target",
                             pattern_y = "full_name",
                             mc.cores  = threads)

# Convet names and convert query to query ID for later matching to metadata
blast_data <- sig_blast_hits %>%
  left_join(do.call(rbind, names_conversion),
            by = "target") %>%
  mutate(query_id = ifelse(test = grepl("^VFG", query),
                           yes  = string_extract(query, "~", "VF[0-9]"),
                           no   = ifelse(test = grepl("^BAC", query),
                                         yes  = string_extract(query, "~", "^BAC"),
                                         no   = ifelse(test = grepl("^NP", query),
                                                       yes  = string_extract(query, "~", "^NP"),
                                                       no   = "Unknown"))),
         query_id = gsub("\\(", "", query_id),
         query_id = gsub("\\)", "", query_id))

print(paste0(Sys.time(), ": Done!"))

# Write out raw table
print(paste0(Sys.time(), ": Writing out raw data with correct names to ", out_dir, "/", "raw_blast_data_corrected_names.csv"))
write.csv(x    = blast_data,
          file = paste0(out_dir, "/", "raw_blast_data_corrected_names.csv"),
          row.names = FALSE)
print(paste0(Sys.time(), ": Done!"))

# Remove raw data to save space
rm("raw_blast_data")
gc()

# Build metadata ---------------------------------------------------------------
print(paste0(Sys.time(), ": Loading Metadata..."))

# Load in protein metadata
VFDB_metadata        <- read.csv(file = VFDB_meta)
bactibase_metadata   <- read.csv(file = bact_meta)
siderophore_metadata <- read.csv(file = side_meta)

VFDB_meta_trim <- VFDB_metadata %>%
  select(VFID, VFcategory) %>%
  transmute(id       = VFID,
            category = VFcategory,
            type     = "Virulence_Factor")

bacti_meta_trim <- bactibase_metadata %>%
  select(id, Class) %>%
  transmute(id       = id,
            category = ifelse(is.na(Class), "Unknown", Class),
            type     = "Bacteriocin")

sidero_meta_trim <- siderophore_metadata %>%
  select(Protein.Referenece, Gene.name) %>%
  transmute(id       = Protein.Referenece,
            category = Gene.name,
            type     = "Siderophore_Biosynthesis")

all_meta_trim <- do.call(rbind, list(VFDB_meta_trim, bacti_meta_trim, sidero_meta_trim))

# Make table of all protein:genome combinations
temp_function <- function(x){
  
  output <- all_meta_trim %>%
    mutate(bacteria = x)
  
  return(output)
}

tmp_full_table <- mclapply(X        = as.list(unique(blast_data$full_name)),
                           FUN      = temp_function,
                           mc.cores = threads)

full_table <- do.call(rbind, tmp_full_table)

# Get mean counts of each gene per genome --------------------------------------
print(paste0(Sys.time(), ": Calculating counts of each gene category per genome"))

mc_group <- function(x){
  
  output <- x %>%
    group_by(full_name, query_id) %>%
    summarise(n = n(),
              .groups = "drop_last")
  
  return(output)
}

per_genome_mean_counts <- do.call(what = rbind, 
                                  args = mclapply(X        = split(blast_data, f = ~ full_name),
                                                  FUN      = mc_group,
                                                  mc.cores = threads)) %>%
  group_by(query_id, full_name) %>%
  summarise(count = sum(n)) %>%
  right_join(full_table,
             by = c("query_id" = "id",
                    "full_name"  = "bacteria"))  %>%
  mutate(label = ifelse(type == "Siderophore_Biosynthesis",
                        category,
                        query_id)) %>%
  replace_na(list(count = 0))

write.csv(x         = per_genome_mean_counts,
          file      = paste0(out_dir, "/", "gene_catagory_counts_per_genome.csv"),
          row.names = FALSE)

print(paste0(Sys.time(), ": Done!"))

# Plot results
## Calculate order of labels
label_order <- per_genome_mean_counts %>%
  group_by(label) %>%
  summarise(mean_count = mean(count)) %>%
  arrange(desc(mean_count)) %>%
  pull(label)

per_genome_mean_counts[, "label"] <- factor(per_genome_mean_counts$label, levels = label_order)

# Plot each type separately
for(type in unique(per_genome_mean_counts$type)){
  
  print(paste0(Sys.time(), ": Making ", type, " graph..."))
  
  n_genes <- length(unique(unlist(per_genome_mean_counts[per_genome_mean_counts$type == type, "label"])))
  
  plot_height <- 15
  plot_width  <- n_genes*0.2
  
  if(plot_width <= plot_height){
    
    plot_width <- plot_height
    
  }
  
  
  gene_count_graph <- ggplot(data    = per_genome_mean_counts %>%
                               filter(type == UQ(type)),
                             mapping = aes(x     = label,
                                           y     = count,
                                           fill  = category,
                                           color = category)) +
    geom_boxplot(show.legend = TRUE) +
    scale_x_discrete(guide = guide_axis(n.dodge = 1),
                     name  = type) +
    scale_y_continuous(name = "Count") +
    labs(title = paste0("Counts of ", n_genes, " ", type, "s")) +
    # guides(colour = guide_legend(nrow = 2),
    #         fill   = guide_legend(nrow = 2)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5,
                                     size  = 5),
          plot.title = element_text(hjust = 0.5),
          legend.position = "bottom")
  
  
  ggsave(filename = paste0(date, "_", type, "_", "count_graph.png"),
         plot     = gene_count_graph,
         path     = paste0(out_dir, "/graphs"),
         width    = plot_width,
         height   = plot_height,
         units    = "cm",
         limitsize = FALSE)
  
  print(paste0(Sys.time(), ": Done!"))
}

# Counts of things in PAO1
PAO1_graph <- ggplot(data = per_genome_mean_counts %>%
                       filter(full_name == "Pseudomonas_aeruginosa_PAO1_107"),
                     mapping = aes(x = label,
                                   y = count)) +
  geom_col()+
  labs(title = "Genes in PAO1") +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Count\n") +
  facet_wrap(facets = vars(type),
             scales = "free",
             nrow   = 3,
             ncol   = 1) +
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
       width    = 130.6,
       height   = 130.6/2,
       units    = "cm",
       limitsize = FALSE)

rm(per_genome_mean_counts)
gc()


# Calculate z-scores based on my chosen categories -----------------------------
# Get mean count of each gene in all genomes
print(paste0(Sys.time(), ": Calculating mean counts of each gene catagory..."))


mean_counts <- blast_data %>%
  select(-target) %>%
  group_by(query_id, full_name) %>%
  summarise(n = n()) %>%
  right_join(full_table,
             by = c("query_id" = "id",
                    "full_name"  = "bacteria")) %>%
  replace_na(list(n = 0)) %>%
  group_by(category, type) %>%
  summarise(mean = mean(n),
            sd   = sd(n))


print(paste0(Sys.time(), ": Saving out mean counts..."))
write.csv(x    = mean_counts,
          file = paste0(out_dir, "/", date, "_blast_mean_counts.csv"),
          row.names = FALSE)

print(paste0(Sys.time(), ": Calculating z-scores..."))

blast_z_scores <- blast_data %>%
  select(-target) %>%
  group_by(query_id, full_name) %>%
  summarise(n = n()) %>%
  right_join(full_table,
             by = c("query_id" = "id",
                    "full_name"  = "bacteria")) %>%
  replace_na(list(n = 0)) %>%
  group_by(category, type, full_name) %>%
  summarise(count = sum(n)) %>%
  left_join(mean_counts,
            by = c("category" = "category",
                   "type"     = "type")) %>%
  mutate(z_score = ifelse(is.nan((count - mean)/sd),
                          0,
                          (count - mean)/sd))

print(paste0(Sys.time(), ": Done!"))
print("")
print(paste0(Sys.time(), "Saving out z-scores..."))
write.csv(x    = blast_z_scores,
          file = paste0(out_dir, "/", date, "_blast_z_scores.csv"),
          row.names = FALSE)

print(paste0(Sys.time(), ": Done!"))

# Test to see if each category has a skewed distribution
print(paste0(Sys.time(), ": Calculating skewness"))

skew_data <- blast_z_scores %>%
  group_by(category, type) %>%
  summarise(skew     = skewness(count),
            kurtosis = kurtosis(count))


write.csv(file = paste0(out_dir, "/", date, "_skew_data.csv"),
          x    = skew_data,
          row.names = FALSE)

print(paste0(Sys.time(), ": Script complete!"))

