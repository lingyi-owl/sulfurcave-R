library(dplyr)
# library(stringr)
library(tidyr)
library(ggplot2)
library(ggalluvial)
library(purrr)
# library(forcats)

# 1) Import data
# import the virus character data (Note: the iphop data is old and not usable in the file)
virus_character_df <- read.csv2('/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/virus_character/sulfur_cave_virus_character_abundance.txt', sep = '\t')

# select the vcat virus taxonomy data
virus_vcat_df <- virus_character_df %>%
  select(contig_id, vcat_class)

# select the virus abundance data
virus_abundance_df <- virus_character_df %>%
  select(contig_id, starts_with("Sample_"))

# convert the abundance data frame from wide to long format
virus_abundance_df_long <- virus_abundance_df %>% 
  pivot_longer(
    cols = `Sample_ERR10036468`:`Sample_ERR10036470`, 
    names_to = "sample",
    values_to = "abundance"
  )

# convert the abundance column from character to numeric
virus_abundance_df_long$abundance <- as.numeric(virus_abundance_df_long$abundance)
class(virus_abundance_df_long$abundance)

# import the cleaned new iphop data frame without MAGs
iphop_genome_df <- read.csv2('/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/iphop/iphop_v1.4_no_MAGs/sulfurcave_iphop_v1.4_no_MAGs_genome_output_distinct.txt', sep = '\t')

# find the virus with multiple predicted hosts
multiple_hosts_virus <- iphop_genome_df$contig_id[duplicated(iphop_genome_df$contig_id)]

# 2) Weight predicted hosts evenly by the number of predictions
# distribute abundance evenly across all hosts per contig_id
hosts_weight <- iphop_genome_df %>%
  group_by(contig_id) %>%
  mutate(wt = 1 / n()) %>%
  ungroup()

# 3) Join abundance + taxonomy + host weights and allocate flow
# each contig–host edge gets a share of the contig's abundance in each sample.
flows <- virus_abundance_df_long %>%
  left_join(virus_vcat_df,  by = "contig_id") %>%
  left_join(hosts_weight,  by = "contig_id") %>%
  mutate(
    wt = ifelse(is.na(wt), 1, wt),       # assign full weight to contigs without hosts
    flow = abundance * wt
  )

# 4) Replace NA with 'Unclassified' and collapse the taxa with <1% abundance
flows_collapsed <- flows %>%
  # replace NA with 'Unclassified'
  mutate(
    vcat_class  = if_else(is.na(vcat_class),  "Unclassified", vcat_class),
    iphop_genus = if_else(is.na(iphop_genus), "Unclassified", iphop_genus)
  ) %>%
  # compute total abundance across all samples
  mutate(total_abundance = sum(abundance, na.rm = TRUE)) %>%
  group_by(vcat_class) %>%
  mutate(vcat_total = sum(abundance, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(iphop_genus) %>%
  mutate(host_total = sum(abundance, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    vcat_class = if_else(vcat_total / total_abundance < 0.01, "Below 1%", vcat_class),
    iphop_genus = if_else(host_total / total_abundance < 0.01, "Below 1%", iphop_genus)
  )

# examine the results
summary_df <- flows_collapsed %>%
  group_by(sample) %>%
  summarise(
    total_abundance = sum(abundance, na.rm = TRUE),
    total_flow = sum(flow, na.rm = TRUE),
    .groups = "drop"
  )

# 5) Build lodes and facet by sample (three panels side by side) with legend shown
lodes_all <- flows_collapsed %>%
  transmute(contig_id, Virus = vcat_class, Host = iphop_genus, weight = flow, sample) %>%
  ggalluvial::to_lodes_form(
    axes  = c("Virus", "Host"),
    key   = "axis",
    value = "stratum",
    id    = "contig_id"
  ) %>%
  dplyr::mutate(
    stratum = as.character(stratum),
    stratum = ifelse(is.na(stratum) | stratum == "", "Unclassified", stratum),
    stratum = factor(stratum)
  )

## 6) Renmae samples
name_map <- c(
  "Sample_ERR10036468" = "Cave biofilm 1",
  "Sample_ERR10036469" = "Cave biofilm 2",
  "Sample_ERR10036470" = "Lab CH4"
)

## 7) Color code
my_colors <- c(
  # viruses (Taxon)
  "Caudoviricetes" = "#e7298a",
  "Faserviricetes" = "#2ca02c", 
  "Megaviricetes" = "#d62728",
  "Below 1%" = "#1f77b4", 
  "Unclassified" = 'gray',
  
  # hosts
  "Acidithiobacillus" = "#fdc086",
  "Below 1%" = "#1f77b4",
  "Cuniculiplasma" = "#b3cde3",
  "Escherichia" = "#C49C94FF",
  "Ferroplasma" = "#ccebc5",
  "Igneacidithiobacillus" = "#decbe4",
  "JAKAFX01" = "#9EDAE5FF",
  "Lachnospira" = "#e5d8bd",
  # "Marinobacter" = "#fddaec",
  "Mobilitalea" = "#DBDB8DFF",
  "Mycobacterium" = "#fbb4ae",
  # "Streptomyces" = "#ffffcc",
  # "Sulfobacillus" = "#98DF8AFF",
  "Unclassified" = 'gray'
)

p <- ggplot(lodes_all,
            aes(x = axis,
                stratum = stratum,
                alluvium = contig_id,
                y = weight*100,
                fill = stratum)) +
  geom_alluvium(alpha = 0.6, knot.pos = 0.4) +
  geom_stratum(width = 0.3, color = "grey30") +
  scale_x_discrete(limits = c("Virus","Host"), expand = c(.12, .12)) +
  scale_fill_manual(values = my_colors, name = "Virus / Host") +
  facet_wrap(~ sample, nrow = 1, labeller = labeller(sample = name_map)) +
  labs(
    x = "Viral Taxonomy --> Microbial Host Prediction",
    y = "Relative abundance (%)",
    title = "Virus and predicted microbial host abundance",
    # subtitle = "Categories with average abundance < 1% across samples are grouped as 'Below 1%'",
    fill = "Taxa / Hosts"
  ) +
  theme_minimal(base_size = 20) +
  theme(legend.position = "none")

p <- p + theme(
  axis.text.x = element_blank(),   # 👈 removes “Virus” and “Host” text on x-axis
  axis.ticks.x = element_blank(),   # 👈 removes x-axis ticks as well
  text = element_text(size = 20, colour = "black"),   # all text black, uniform size
  plot.title = element_text(size = 20, colour = "black"),
  axis.title = element_text(size = 20, colour = "black"),
  axis.text = element_text(size = 20, colour = "black"),
  strip.text = element_text(size = 20, colour = "black")         # facet labels
)

p
