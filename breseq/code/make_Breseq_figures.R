library(tidyverse)
library(patchwork)

#############################################################################
### Load and wrangle breseq output data                                   ###
#############################################################################

# List all files in the directory
files <- list.files(path = "./breseq_output/", pattern = "^AB.*\\.csv$", full.names = TRUE)

# Read, select specific columns, and combine all the matching files
df <- do.call(rbind, lapply(files, function(file) read.csv(file) %>%
  select(title, seq_id, ref_seq, new_seq, new_read_count, ref_read_count,
         ref_read_count,  position_start, position_end, new_read_count,
         frequency, type, snp_type, gene_name, gene_position, gene_product,
         gene_strand, mutation_category, codon_ref_seq, codon_position,
         codon_new_seq, aa_ref_seq, aa_position, aa_new_seq)
         ))


# load the sample key:
sample_key <- read.csv("./data/sequencing_sample_key.csv")

# Define the ancestral mutations
mutations <- c("rpsL_A128C", "rpsL_A262G", "rpsL_G275C",
               "rpoB_T443G", "rpoB_T1555C", "rpoB_C1556T", "rpoB_A1568C",
               "rpoB_C1712T")

# Create sensible names for the genes and systematise mutation names
df <- df %>%
  mutate(
    Mutation = str_c(gene_name, "_", ref_seq, gene_position, new_seq, sep = ""),
    Reads = new_read_count + ref_read_count,
    Frequency = new_read_count / Reads
  ) %>%
  rename(Sample_ID = title,
         Reference_seq_id = seq_id) %>%
  left_join(sample_key, by = "Sample_ID")

#############################################################################
### Produce simplified table                                              ###
#############################################################################

df_simplified <- df |>
  select(Sample_ID, Reference_seq_id,
         Sample_type, Time_point, Competence, Treatment, Shift, Replicate,
         ref_read_count,  position_start, position_end, new_read_count, frequency,
         type, snp_type,
         gene_name, gene_position, gene_product, gene_strand, mutation_category,
         codon_ref_seq, codon_position, codon_new_seq,
         aa_ref_seq, aa_position, aa_new_seq)

write_csv(df_simplified, "./breseq_output/breseq_output_simplified.csv")

# retain only ancestral mutations:
df_filtered <- df |>
  filter(Mutation %in% mutations)

# Define a function to calculate the confidence interval from binom.test
get_conf_int <- function(wins, losses) {
  result <- binom.test(x = c(wins, losses))
  return(result$conf.int)
}

# Use mapply to apply the function to each pair of values from "new_read_count" and "ref_read_count"
conf_int <- mapply(get_conf_int, df_filtered$new_read_count, df_filtered$ref_read_count)
df_filtered <- df_filtered |>
  mutate(Frequency_conf_low = conf_int[1, ],
         Frequency_conf_high = conf_int[2, ])

#############################################################################
### Produce plot                                                          ###
#############################################################################

df_filtered_mix <- df_filtered |>
  filter(Treatment == "MIX", Sample_type == "Pop") |>
  complete(Competence, Replicate, Mutation) |>
  mutate(Competence = factor(Competence, levels = c("com+", "com-"))) |>
  mutate(Mutation = factor(Mutation, levels = mutations))
    
colors <- c(hsv(seq(0.6, 0.75, length.out = 3), 0.7, 1),
            hsv(seq(0.15, 0, length.out = 5), 0.7, 1))

p <- ggplot(df_filtered_mix, aes(x = Replicate, y = Frequency, fill = Mutation)) +
  geom_col(width = 0.8,
           position = position_dodge(0.8)) +
  geom_errorbar(aes(ymin = Frequency_conf_low, ymax = Frequency_conf_high), 
                position = position_dodge(0.8), 
                width = 0.4) +
  facet_wrap(vars(Competence)) +
  scale_x_continuous(breaks = 1:8) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = setNames(colors, mutations)) +
  labs(x = "Replicate", y = "Allele frequency", fill = "Mutation") +
  theme_bw()
  
ggsave("./plots/breseq_results_MIX.pdf", p, width = 8, height = 5)

#############################################################################
### Statistical analyses                                                  ###
#############################################################################

# sum the frequency of mutations for each replicate, ignoring NAs:
mix_sum <- df_filtered_mix %>%
  group_by(Competence, Replicate) %>%
  summarise(freq = sum(Frequency, na.rm = TRUE),
            conf_low = sum(Frequency_conf_low, na.rm = TRUE),
            conf_high = sum(Frequency_conf_high, na.rm = TRUE)) %>%
  ungroup()
mix_sum
