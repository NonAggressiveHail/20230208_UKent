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
read_jackhmmer <- function(file){
  
  tmp_file   <- readLines(con = file)[-c(1:4)]
  
  raw_jackhmmer_data <- read.delim(file   = textConnection(tmp_file),
                                   quote  = "",
                                   sep    = "",
                                   header = FALSE)
  
  colnames(raw_jackhmmer_data) <- c("domain_name", "accession_1", "query_name", "accession_2",
                                    "sequence_evalue", "seq_score", "seq_bias",
                                    "dom_E_value", "dom_score", "dom_bias",
                                    "exp", "reg", "clu",  "ov", 
                                    "env", "dom", "rep", "inc", "description_of_target")
  output <- raw_jackhmmer_data %>%
    select(domain_name, query_name, sequence_evalue)
  
  return(output)
  
}

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


# Load Data --------------------------------------------------------------------
print(paste0(Sys.time(), ": Args are:"))
print(args)
print("")

date     <- gsub("-", "", Sys.Date())
numcores <- detectCores()

print(paste0("Finding files in ", jackhmmer_dir))

full_names      <- read.csv(file = full_names_file)

jackhmmer_files <- list.files(path       = jackhmmer_dir,
                              pattern    = "*.txt",
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
  
  str(raw_jackhmmer_data)
  
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


# Convert names
jackhmmer_data <- sig_jackhmmer_hits %>%
  mutate(target_organism = sub("_[^_]+$", "", domain_name),
         target_gene     = sub(".*_", "", domain_name)) %>%
  mutate(query_id = ifelse(test = grepl("^VFG", query_name),
                           yes  = sub("\\(.*", "", query_name),
                           no   = sub("\\~.*", "", query_name)))
if(FALSE){
  # Need to do this partial join shenannigans to convert genome names to what I usually use
  print(paste0(Sys.time(), ": Converting names..."))
  jackhmmer_names <- compiled_jackhmmer_data %>%
    select(domain_name) %>%
    distinct()
  
  
  names_conversion <- mclapply(X         = split(jackhmmer_names, 1:threads),
                               FUN       = partial_join,
                               y         = full_names,
                               by_x      = "domain_name",
                               pattern_y = "full_name",
                               mc.cores  = threads)
  
  # Convet names and convert query_name to query ID for later matching to metadata
  jackhmmer_data <- sig_jackhmmer_hits %>%
    left_join(do.call(rbind, names_conversion),
              by = "domain_name") %>%
    mutate(query_id = ifelse(test = grepl("^VFG", query_name),
                             yes  = string_extract(query_name, "~", "VF[0-9]"),
                             no   = ifelse(test = grepl("^BAC", query_name),
                                           yes  = string_extract(query_name, "~", "^BAC"),
                                           no   = ifelse(test = grepl("^NP", query_name),
                                                         yes  = string_extract(query_name, "~", "^NP"),
                                                         no   = "Unknown"))),
           query_id = gsub("\\(", "", query_id),
           query_id = gsub("\\)", "", query_id))
  
  print(paste0(Sys.time(), ": Done!"))
}

# Write out compiled table
print(paste0(Sys.time(), ": Writing out compiled data with correct names to ", out_dir, "/", "compiled_jackhmmer_data_corrected_names.csv"))
write.csv(x    = jackhmmer_data,
          file = paste0(out_dir, "/", "compiled_jackhmmer_data_corrected_names.csv"),
          row.names = FALSE)
print(paste0(Sys.time(), ": Done!"))

# Remove raw data to save space
rm("raw_jackhmmer_data")
gc()

# Build metadata ---------------------------------------------------------------
print(paste0(Sys.time(), ": Loading Metadata..."))

# Load in protein metadata
VFDB_meta_trim <-  read.csv(file = VFDB_meta) %>%
  select(gene_id, VF_Name) %>%
  transmute(id       = gene_id,
            category = VF_Name,
            type     = "Virulence_Factor")

bacti_meta_trim <- read.csv(file = bact_meta) %>%
  select(id, Class) %>%
  transmute(id       = id,
            category = ifelse(is.na(Class), "Unknown", Class),
            type     = "Bacteriocin")

sidero_meta_trim <-  read.csv(file = side_meta) %>%
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

tmp_full_table <- mclapply(X        = as.list(unique(jackhmmer_data$target_organism)),
                           FUN      = temp_function,
                           mc.cores = threads)

full_table <- do.call(rbind, tmp_full_table)

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
if(graphs){
# Plot results
## Calculate order of labels
label_order <- per_genome_mean_counts %>%
  group_by(label) %>%
  summarise(mean_count = mean(count)) %>%
  arrange(desc(mean_count)) %>%
  pull(label)

per_genome_mean_counts[, "label"] <- factor(per_genome_mean_counts$label, levels = label_order)

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



# Plot each type separately
for(type in unique(per_genome_mean_counts$type)){
  
  print(paste0(Sys.time(), ": Making ", type, " graph..."))
  
  n_genes <- length(unique(unlist(per_genome_mean_counts[per_genome_mean_counts$type == type, "label"])))
  
  plot_height <- 15
  plot_width  <- 50
  
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
         width    = plot_width,
         height   = plot_height,
         units    = "cm",
         limitsize = FALSE)
  
  print(paste0(Sys.time(), ": Done!"))
}
}

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
