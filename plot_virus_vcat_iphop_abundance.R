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

test_flows <- virus_abundance_df_long %>%
  left_join(hosts_weight, by = "contig_id") %>%
  mutate(
    wt = ifelse(is.na(wt), 1, wt),       # assign full weight to contigs without hosts
    flow = abundance * wt
  )

summary_df <- flows %>%
  group_by(sample) %>%
  summarise(
    total_abundance = sum(abundance, na.rm = TRUE),
    total_flow = sum(flow, na.rm = TRUE),
    .groups = "drop"
  )
# 4) Build lodes per sample (Virus → Host), keeping contig_id as alluvium
make_lodes <- function(df_one_sample) {
  df_one_sample %>%
    transmute(contig_id, Virus = vcat_class, Host = iphop_genus, weight = flow) %>%
    ggalluvial::to_lodes_form(
      axes  = c("Virus", "Host"),
      key   = "axis",
      value = "stratum",
      id    = "contig_id"
    )
}

lodes_by_sample <- flows %>%
  group_split(sample, .keep = TRUE) %>%
  setNames(unique(flows$sample)) %>%
  map(make_lodes)

# 5) A plotting helper to keep styling consistent across samples
plot_alluvial <- function(lodes_df, sample_name) {
  ggplot(lodes_df,
         aes(x = axis,
             stratum = stratum,
             alluvium = contig_id,
             y = weight,
             fill = stratum)) +
    geom_alluvium(alpha = 0.6, knot.pos = 0.4) +
    geom_stratum(width = 0.3, color = "#decbe4") +
    scale_x_discrete(limits = c("Virus","Host"), expand = c(.12, .12)) +
    labs(
      title = paste("Taxonomy → Host alluvial"),
      subtitle = paste("Sample:", sample_name),
      x = NULL, y = "Relative abundance (allocated)"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
}

plots <- imap(lodes_by_sample, plot_alluvial)
plots



# starting from `flows` (with columns: contig_id, virus, host, sample, flow)

lodes_all <- flows %>%
  transmute(contig_id, Virus = vcat_class, Host = iphop_genus, weight = flow, sample) %>%
  ggalluvial::to_lodes_form(
    axes  = c("Virus", "Host"),
    key   = "axis",
    value = "stratum",
    id    = "contig_id"
  )

ggplot(lodes_all,
       aes(x = axis,
           stratum = stratum,
           alluvium = contig_id,
           y = weight,
           fill = stratum)) +
  geom_alluvium(alpha = 0.6, knot.pos = 0.4) +
  geom_stratum(width = 0.3, color = "#decbe4") +
  scale_x_discrete(limits = c("Virus","Host"), expand = c(.12, .12)) +
  facet_wrap(~ sample, nrow = 1) +
  theme_minimal(base_size = 12) 

#####################################################################
# 3) Compute mean share across samples for each vcat_class and host_genus
#    Work with within-sample percentages so it's robust even if per-sample totals ≠ 1
flows_pct <- flows %>%
  group_by(sample) %>%
  mutate(sample_total = sum(flow),
         flow_pct = ifelse(sample_total > 0, flow / sample_total, 0)) %>%
  ungroup()

taxon_avg <- flows_pct %>%
  group_by(vcat_class, sample) %>%
  summarise(pct = sum(flow_pct), .groups = "drop") %>%
  group_by(vcat_class) %>%
  summarise(avg_pct = mean(pct), .groups = "drop")

host_avg <- flows_pct %>%
  group_by(host_genus, sample) %>%
  summarise(pct = sum(flow_pct), .groups = "drop") %>%
  group_by(host_genus) %>%
  summarise(avg_pct = mean(pct), .groups = "drop")

rare_taxa  <- taxon_avg  %>% filter(avg_pct < 0.01) %>% pull(vcat_class)
rare_hosts <- host_avg   %>% filter(avg_pct < 0.01) %>% pull(host_genus)

## 4) Collapse rare categories to "Below 1%"
flows_collapsed <- flows %>%
  mutate(
    vcat_class2 = fct_other(vcat_class,
                            keep = setdiff(unique(vcat_class), rare_taxa),
                            other_level = "Below 1%"),
    host_genus2 = fct_other(host_genus,
                            keep = setdiff(unique(host_genus), rare_hosts),
                            other_level = "Below 1%")
  )

## 5) Build lodes and facet by sample (three panels side by side) with legend shown
lodes_all <- flows_collapsed %>%
  transmute(contig_id, Virus = vcat_class2, Host = host_genus2, weight = flow, sample) %>%
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
  "Unclassified" = '#decbe4',
  
  # hosts
  "Acidithiobacillus" = "#fdc086",
  "Ferroplasma" = "#ccebc5",
  "Mycobacterium" = "#fbb4ae",
  "Cuniculiplasma" = "#b3cde3",
  "Below 1%" = "#1f77b4",
  "Unclassified" = '#decbe4'
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
    x = NULL,
    y = "Relative abundance (%)",
    title = "Virus and predicted host abundance",
    # subtitle = "Categories with average abundance < 1% across samples are grouped as 'Below 1%'",
    fill = "Taxa / Hosts"
  ) +
  theme_minimal(base_size = 20)
  # theme(legend.position = "none")

p <- p + theme(
  text = element_text(size = 20, colour = "black"),   # all text black, uniform size
  plot.title = element_text(size = 20, colour = "black"),
  axis.title = element_text(size = 20, colour = "black"),
  axis.text = element_text(size = 20, colour = "black"),
  strip.text = element_text(size = 20, colour = "black")         # facet labels
)

p

########################################################################
# summary
iphop_genome_df %>%
  filter(!is.na(host_genus) & host_genus != "") %>%
  summarise(unique_virus_with_hosts = n_distinct(Virus))
