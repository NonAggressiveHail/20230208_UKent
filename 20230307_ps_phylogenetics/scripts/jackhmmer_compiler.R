# Script to calcualte z-scores of genes of interest for my bacteria
# and export z-scores for plotting on heatmap
library(tidyverse)
library(parallel)
library(moments)
library(rhmmer)

# Options ----------------------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)

jackhmmer_dir   <- args[1]
out_dir         <- args[2]
side_meta       <- args[3]
VFDB_meta       <- args[4]
bact_meta       <- args[5]
full_names_file <- args[6]
threads         <- as.numeric(args[7])
graphs          <- args[8]

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

# Build metadata ---------------------------------------------------------------
print(paste0(Sys.time(), ": Loading Metadata..."))

# Load in protein metadata
VFDB_meta_trim <-  read.csv(file = VFDB_meta) %>%
  select(gene_id, VFID, VFCID, VF_Name, VFcategory ) %>%
  transmute(gene_id        = gene_id,
            factor_id      = VFID,
            factor_id_name = VF_Name,
            factor_category = VFCID,
            factor_category_name = VFcategory,
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

# Load Jackhmmer Data ----------------------------------------------------------
print(paste0(Sys.time(), ": Args are:"))
print(args)
print("")

date     <- gsub("-", "", Sys.Date())
numcores <- detectCores()

print(paste0("Finding files in ", jackhmmer_dir))

target_organisms      <- read.csv(file = full_names_file)

jackhmmer_files <- list.files(path       = jackhmmer_dir,
                              pattern    = "*_tbl.txt",
                              full.names = TRUE)

print(paste0(Sys.time(), ": Found ", length(jackhmmer_files), " files"))
print("")
print(paste0(Sys.time(), ": Reading in files with ", threads, " threads, out of ", numcores, " available"))

# If already compiled data exists, load it
if(file.exists(paste0(out_dir, "/",  "compiled_jackhmmer_data.csv"))){
  
  print(paste0(Sys.time(), ": Previous compiled jackhmmer data found, using"))
  compiled_jackhmmer_data <- read.csv(paste0(out_dir, "/", "compiled_jackhmmer_data.csv"))
  
  # Otherwise compile it
} else {
  
  
  print(paste0(Sys.time(), ": Previous compiled jackhmmer data not found, compiling"))
  
  raw_jackhmmer_data <- mclapply(X        = jackhmmer_files,
                                 FUN      = read_tblout,
                                 mc.cores = threads)
  
  compiled_jackhmmer_data <- do.call(rbind, raw_jackhmmer_data) %>%
    select(domain_name, query_name, sequence_evalue)
  
  print(paste0(Sys.time(), ": Writing out compiled data to ", out_dir, "/", "compiled_jackhmmer_data.csv"))
  
  write.csv(x         = compiled_jackhmmer_data,
            file      = paste0(out_dir, "/", "compiled_jackhmmer_data.csv"),
            row.names = FALSE)
  
}

print(paste0(Sys.time(), ": Done!"))
print(paste0(Sys.time(), ": jackhmmer data has ", nrow(compiled_jackhmmer_data), " lines"))
print("")

# Filter out bad hits
print(paste0(Sys.time(), ": Removing hits with e-value <= 1e-03"))
sig_jackhmmer_hits <- compiled_jackhmmer_data %>%
  filter(sequence_evalue <= 1e-03)

print(paste0(Sys.time(), ": ", nrow(compiled_jackhmmer_data) - nrow(sig_jackhmmer_hits), " were not significant"))


# Convert names and remove duplicated data
jackhmmer_data <- sig_jackhmmer_hits %>%
  mutate(target_organism = sub("_[^_]+$", "", domain_name),                # exract the target PA strain
         target_gene     = sub(".*_", "", domain_name)) %>%                # extract the locus in that strain
  mutate(query_gene = ifelse(test = grepl("^VFG", query_name),
                             yes  = sub("\\(.*", "", query_name),
                             no   = sub("\\~.*", "", query_name))) %>%                         # Tidy the ID number up
  group_by(target_organism, target_gene) %>%
  filter(sequence_evalue == min(sequence_evalue)) %>%                      # Only keep the sequences with the lowest e_value
  reframe(query_gene     = query_gene,
          hits_for_locus = n()) %>%                                        # Calculate how many total queries hit that locus
  group_by(target_organism, target_gene, query_gene)  %>%                 
  reframe(hits_for_locus  = hits_for_locus,
          match_quality   = n()/hits_for_locus * 100) %>%                  # Calculate the number of times each query hit
  group_by(target_organism, target_gene) %>%
  filter(match_quality == max(match_quality)) %>%
  left_join(all_meta_trim,
            by = c("query_gene" = "gene_id")) %>%
  group_by(target_organism, target_gene) %>%
  reframe(target_organism = target_organism,
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
                                        "Unclear")) %>%
  distinct()

# Write out compiled table
print(paste0(Sys.time(), ": Writing out compiled data with correct names to ", out_dir, "/", "compiled_jackhmmer_data_corrected_names.csv"))
write.csv(x    = jackhmmer_data,
          file = paste0(out_dir, "/", "compiled_jackhmmer_data_corrected_names.csv"),
          row.names = FALSE)
print(paste0(Sys.time(), ": Done!"))

# Remove raw data to save space
rm("raw_jackhmmer_data")
gc()

stop()
# Get mean counts of each gene per genome --------------------------------------
print(paste0(Sys.time(), ": Calculating counts of each gene category per genome"))

mc_group <- function(x){
  
  output <- x %>%
    group_by(target_organism, query_id) %>%
    summarise(n = n(),
              .groups = "drop_last")
  
  return(output)
}

per_genome_mean_counts <- do.call(what = rbind, 
                                  args = mclapply(X        = split(jackhmmer_data, f = ~ target_organism),
                                                  FUN      = mc_group,
                                                  mc.cores = threads)) %>%
  group_by(query_id, target_organism) %>%
  summarise(count = sum(n)) %>%
  right_join(full_table,
             by = c("query_id"         = "id",
                    "target_organism"  = "bacteria"))  %>%
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
  # plot_width  <- n_genes*0.2
  
  if(plot_width <= plot_height){
    
    plot_width <- plot_height
    
  }
  
  
  gene_count_graph <- ggplot(data    = per_genome_mean_counts %>%
                               filter(count != 0) %>%
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
         width    = 50,
         height   = plot_height,
         units    = "cm",
         limitsize = FALSE)
  
  print(paste0(Sys.time(), ": Done!"))
}

# Counts of things in PAO1
PAO1_graph <- ggplot(data = per_genome_mean_counts %>%
                       filter(target_organism == "Pa_PAO1_107"),
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


mean_counts <- jackhmmer_data %>%
  select(-domain_name) %>%
  group_by(query_id, target_organism) %>%
  summarise(n = n()) %>%
  right_join(full_table,
             by = c("query_id" = "id",
                    "target_organism"  = "bacteria")) %>%
  replace_na(list(n = 0)) %>%
  group_by(category, type) %>%
  summarise(mean = mean(n),
            sd   = sd(n))


print(paste0(Sys.time(), ": Saving out mean counts..."))
write.csv(x    = mean_counts,
          file = paste0(out_dir, "/", date, "_jackhmmer_mean_counts.csv"),
          row.names = FALSE)

print(paste0(Sys.time(), ": Calculating z-scores..."))

jackhmmer_z_scores <- jackhmmer_data %>%
  select(-domain_name) %>%
  group_by(query_id, target_organism) %>%
  summarise(n = n()) %>%
  right_join(full_table,
             by = c("query_id" = "id",
                    "target_organism"  = "bacteria")) %>%
  replace_na(list(n = 0)) %>%
  group_by(category, type, target_organism) %>%
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
write.csv(x    = jackhmmer_z_scores,
          file = paste0(out_dir, "/", date, "_jackhmmer_z_scores.csv"),
          row.names = FALSE)

print(paste0(Sys.time(), ": Done!"))

# Test to see if each category has a skewed distribution
print(paste0(Sys.time(), ": Calculating skewness"))

skew_data <- jackhmmer_z_scores %>%
  group_by(category, type) %>%
  summarise(skew     = skewness(count),
            kurtosis = kurtosis(count))


write.csv(file = paste0(out_dir, "/", date, "_skew_data.csv"),
          x    = skew_data,
          row.names = FALSE)

print(paste0(Sys.time(), ": Script complete!"))

