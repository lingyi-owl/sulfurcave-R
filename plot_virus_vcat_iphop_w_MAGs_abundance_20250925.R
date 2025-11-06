library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggalluvial)
library(purrr)
library(forcats)

# vcat_df <- read.csv2('/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/vcat/UPDATED_filtered_virus_contigs_cb1_60_cb2_300_C_10_length_10kb_seq_fasta.tsv', sep = '\t')
# vcat_df <- vcat_df %>%
#   select(SequenceID,)
virus_character_df <- read.csv2('/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/virus_character/sulfur_cave_virus_character_abundance.txt', sep = '\t')

virus_vcat_df <- virus_character_df %>%
  select(contig_id, vcat_class)

virus_abundance_df <- virus_character_df %>%
  select(contig_id, starts_with("Sample_"))

virus_abundance_df_long <- virus_abundance_df %>% 
  pivot_longer(
    cols = `Sample_ERR10036468`:`Sample_ERR10036470`, 
    names_to = "sample",
    values_to = "abundance"
  )
class(virus_abundance_df_long$abundance)
virus_abundance_df_long$abundance <- as.numeric(virus_abundance_df_long$abundance)

iphop_genome_df <- read.csv2('/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/iphop/iphop_v1.4_w_MAGs/Host_prediction_to_genome_m90.csv', sep = ',')
iphop_genome_df <- iphop_genome_df %>%
  mutate(
    Host.taxonomy = str_trim(Host.taxonomy),
    # Prefer genus (even if empty), otherwise fall back to family, else NA
    host_genus_lineage = case_when(
      str_detect(Host.taxonomy, "g__") ~
        str_replace(Host.taxonomy, "^(.*?g__[^;]*).*", "\\1"),
      str_detect(Host.taxonomy, "f__") ~
        str_replace(Host.taxonomy, "^(.*?f__[^;]*).*", "\\1"),
      TRUE ~ NA_character_
    ),
    host_genus = case_when(
      str_detect(Host.taxonomy, "g__") ~ 
        str_replace(Host.taxonomy, ".*g__([^;]*).*", "\\1"),
      TRUE ~ NA_character_
    )
  )

iphop_genome_df2 <- iphop_genome_df %>%
  select(Virus, host_genus)

iphop_genome_df_distinct <- iphop_genome_df2 %>%
  distinct()

colnames(iphop_genome_df_distinct) <- c('contig_id', 'host_genus')

multiple_hosts_virus <- iphop_genome_df_distinct$contig_id[duplicated(iphop_genome_df_distinct$contig_id)]
multiple_hosts_virus_df <- iphop_genome_df_distinct[iphop_genome_df_distinct$contig_id %in% multiple_hosts_virus,]

host_df <- as.data.frame(virus_vcat_df$contig_id)
colnames(host_df) <- 'contig_id'
host_df <- left_join(host_df, iphop_genome_df_distinct, by = 'contig_id')
host_df$score <- NA
# 1) Normalize host weights per contig
hosts_w <- host_df %>%
  group_by(contig_id) %>%
  mutate(
    wt = if (all(is.na(score))) 1/n() else score / sum(score, na.rm = TRUE)
  ) %>%
  ungroup()

# 2) Join abundance + taxonomy + host weights and allocate flow
#    Each contig–host edge gets a share of the contig's abundance in each sample.
flows <- virus_abundance_df_long %>%
  inner_join(virus_vcat_df,  by = "contig_id") %>%
  inner_join(hosts_w,  by = "contig_id") %>%
  mutate(flow = abundance * wt)

# 3) Build lodes per sample (Virus → Host), keeping contig_id as alluvium
make_lodes <- function(df_one_sample) {
  df_one_sample %>%
    transmute(contig_id, Virus = vcat_class, Host = host_genus, weight = flow) %>%
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

# 4) A plotting helper to keep styling consistent across samples
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
  transmute(contig_id, Virus = vcat_class, Host = host_genus, weight = flow, sample) %>%
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
