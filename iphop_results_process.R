library(dplyr)
library(stringr)

# iphop_genus_df <- read.table('/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/iphop/iphop_v1.4_no_MAGs/Host_prediction_to_genus_m90.csv',
#                              sep = ',',
#                              header = TRUE)
# length(unique(iphop_genus_df$Virus))
# import iphop output
iphop_genome_df <- read.table('/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/iphop/iphop_v1.4_w_MAGs/Host_prediction_to_genome_m90.csv',
                              sep = ',',
                              header = TRUE)
length(unique(iphop_genome_df$Virus))

# add a new column 'iphop_host_genue' with string values of the predicted host_genus
iphop_genome_df <- iphop_genome_df %>%
  mutate(
    iphop_genus = str_extract(Host.taxonomy, "g__[^;]+"),
    iphop_genus = str_replace(iphop_genus, "^g__", "")
  )

iphop_genome_df$iphop_genus[iphop_genome_df$iphop_genus == "Acidithiobacillus_B"] <- "Acidithiobacillus"
iphop_genome_df$iphop_genus[iphop_genome_df$iphop_genus == "Acidithiobacillus_A"] <- "Acidithiobacillus"

# select the columns of 'Virus' and 'iphop_genus' and remove the duplicates
iphop_genome_df_distinct <- iphop_genome_df %>%
  select('Virus', 'iphop_genus') %>%
  distinct()

colnames(iphop_genome_df_distinct) <- c('contig_id', 'iphop_genus')
# wirte the output to a text file
write.table(iphop_genome_df_distinct, '/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/iphop/iphop_v1.4_w_MAGs/sulfurcave_iphop_v1.4_w_MAGs_genome_output_distinct.txt', 
            sep = '\t',
            quote = F,
            row.names = F,
            col.names = T)

nrow(iphop_genome_df_distinct)
length(unique(iphop_genome_df_distinct$contig_id))
length(unique(iphop_genome_df_distinct$iphop_genus))

virus_multiple_host <- iphop_genome_df_distinct %>%
  filter(duplicated(contig_id) | duplicated(contig_id, fromLast = TRUE))

length(unique(virus_multiple_host$contig_id))
