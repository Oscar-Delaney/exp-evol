library(tidyverse)
library(patchwork)

# load the necessary csv's produced by breseq:
AB3 <- read.csv("./breseq_output/AB3.csv")
AB13 <- read.csv("./breseq_output/AB13.csv")
AB3_extra <- read.csv("./breseq_output/AB3_extra.csv")
AB13_extra <- read.csv("./breseq_output/AB13_extra.csv") %>%
  select(-multiple_polymorphic_SNPs_in_same_codon)
df <- rbind(AB3, AB13, AB3_extra, AB13_extra)

# load and process the sample key:
names_df <- read.csv("./data/names.csv")
# Separate the 'Name' column
names_df <- names_df %>%
  separate(Name, into = c("shift", "strain", "condition", "replicate"), sep = "_") %>%
  group_by(strain, condition) %>%
  mutate(rep = 1:n(),
    strain = case_when(
      strain == "AB3" ~ "com+",
      strain == "AB13" ~ "com-"
    )
  ) %>%
  ungroup() %>%
  select(-c(shift, replicate))

# Define the ancestral mutations
mutations <- c("rpsL_A128C", "rpsL_A262G", "rpsL_G275C",
               "rpoB_T443G", "rpoB_T1555C", "rpoB_C1556T", "rpoB_A1568C",
               "rpoB_C1712T")

# Create sensible names for the genes and systematise mutation names
df_filtered <- df %>%
  mutate(
    gene_name = case_when(
      gene_name == "ACIAD_RS04070" ~ "rpsL_",
      gene_name == "ACIAD_RS01460" ~ "rpoB_",
      TRUE ~ gene_name
    ),
    # amino acid labelling
    # mutation = str_c(gene_name, aa_ref_seq, aa_position, aa_new_seq, sep = ""),
    # nucleotide labelling
    mutation = str_c(gene_name, "_", ref_seq, gene_position, new_seq, sep = ""),
    Sample_ID = str_sub(title, start = 1, end = 6),
    reads = new_read_count + ref_read_count,
    freq = new_read_count / reads
  ) %>%
  filter(mutation %in% mutations) %>%
  left_join(names_df, by = "Sample_ID")

# Define a function to calculate the confidence interval from binom.test
get_conf_int <- function(wins, losses) {
  result <- binom.test(x = c(wins, losses))
  return(result$conf.int)
}

# Use mapply to apply the function to each pair of values from "new_read_count" and "ref_read_count"
conf_int <- mapply(get_conf_int, df_filtered$new_read_count, df_filtered$ref_read_count)
df_filtered$conf_low <- conf_int[1, ]
df_filtered$conf_high <- conf_int[2, ]

final_df <- df_filtered %>%
  select(strain, condition, rep, mutation, freq, conf_low, conf_high) %>%
  arrange(mutation, strain)

df_filtered_mix <- df_filtered |>
  filter(condition == "MIX") |>
  complete(strain, rep, mutation) |>
  mutate(strain = factor(strain, levels = c("com+", "com-"))) |>
  mutate(mutation = factor(mutation, levels = mutations))
    
colors <- c(hsv(seq(0.6, 0.75, length.out = 3), 0.7, 1),
            hsv(seq(0.15, 0, length.out = 5), 0.7, 1))

p <- ggplot(df_filtered_mix, aes(x = rep, y = freq, fill = mutation)) +
  geom_col(width = 0.8,
           position = position_dodge(0.8)) +
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high), 
                position = position_dodge(0.8), 
                width = 0.4) +
  facet_wrap(vars(strain)) +
  scale_x_continuous(breaks = 1:8) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = setNames(colors, mutations)) +
  labs(x = "Replicate", y = "Allele frequency", fill = "Mutation") +
  theme_bw()
  
ggsave("./plots/breseq_results_MIX.pdf", p, width = 8, height = 5)
  
# sum the frequency of mutations for each replicate, ignoring NAs:
mix_sum <- df_filtered_mix %>%
  group_by(strain, rep) %>%
  summarise(freq = sum(freq, na.rm = TRUE),
            conf_low = sum(conf_low, na.rm = TRUE),
            conf_high = sum(conf_high, na.rm = TRUE)) %>%
  ungroup()
mix_sum
