#script to analyse proportions of each gene in a genome
#dependencies
library(tidyverse)

# options
args = commandArgs(trailingOnly = TRUE)
print(args)

#accesion 
acc <- args[1]
#get directory
dir <- args[2]

#load in
gene_data <- read.table(file   = paste0(dir, "/prokka/", acc, ".tsv"), 
                        header = TRUE,
                        sep    = "\t",
                        quote  = "")


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

#write out table
write.csv(x         = gene_freq_data,
          file      = paste0("../data/", acc, "_unknown_frequency.tsv"),
	  row.names = FALSE)





