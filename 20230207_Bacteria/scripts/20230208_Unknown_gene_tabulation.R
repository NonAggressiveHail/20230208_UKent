#script to analyse proportions of each gene in a genome
#dependencies
library(tidyverse)

#get file
file <- "./data/GCF_000006765.1_ASM676v1.tsv"

#load in
gene_data <- read.table(file = file,
                        header = TRUE,
                        sep = "\t",
                        quote = "")

#process for plotting
gene_freq_data <- gene_data %>%
  mutate(classification = ifelse(test = grepl(pattern = "hypothetical|putitative",
                                              x       = product),
                                 yes  = "Unannotated",
                                 no   = "Annotated")) %>%
  group_by(classification) %>%
  summarise(n = n()) %>%
  summarise(classification = classification,
            n = n,
            percentage = n/sum(n)*100) %>%
  pivot_longer(cols = -classification,
               names_to = "metric",
               values_to = "value")
  
head(gene_freq_data)

#plot
gene_freq_plot <- ggplot(data = gene_freq_data,
                         mapping = aes(x = "",
                                       y = value,
                                       fill = classification)) +
  geom_col() +
  facet_wrap(facets = vars(metric),
             scale  = "free") +
  theme_bw()

gene_freq_plot

#lets see what the unannotated proteins are 
unkn_freq_data <- gene_data %>%
  mutate(classification = ifelse(test = grepl(pattern = "hypothetical|putative",
                                              x       = product),
                                 yes  = "Unannotated",
                                 no   = "Annotated")) %>%
  filter(classification == "Unannotated") %>%
  group_by(product) %>%
  summarise(n = n())

#plot 
unknown_freq_plot <- ggplot(data = unkn_freq_data,
                            mapping = aes(y = product,
                                          x = n)) +
  geom_col() +
  theme_bw()

unknown_freq_plot








