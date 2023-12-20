library(tidyverse)

# save the current working directory
parent_wd <- getwd()
setwd(paste0(parent_wd, "/Breseq_output"))

# load the necessary csv's
AB3 <- read.csv("AB3.csv")
AB13 <- read.csv("AB13.csv")
AB3_extra <- read.csv("AB3_extra.csv")
AB13_extra <- read.csv("AB13_extra.csv") %>%
  select(-multiple_polymorphic_SNPs_in_same_codon)
df <- rbind(AB3, AB13, AB3_extra, AB13_extra)
names_df <- read.csv("names.csv")

# Separate the 'Name' column
names_df <- names_df %>%
  separate(Name, into = c("shift", "strain", "condition", "replicate"), sep = "_") %>%
  group_by(strain, condition) %>%
  mutate(rep = 1:n(),
    strain = case_when(
      strain == "AB3" ~ "rec+",
      strain == "AB13" ~ "rec-"
    )
  ) %>%
  ungroup() %>%
  select(-c(shift, replicate))

# Define the ancestral mutations
mutations <- c("rpoB_T443G", "rpoB_T1555C", "rpoB_C1556T", "rpoB_A1568C",
              "rpoB_C1712T", "rpsL_A128C", "rpsL_A262G", "rpsL_G275C")

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
  # complete(strain, condition, rep = 1:6, mutation) %>%
  # replace_na(list(freq = 0)) %>%
  arrange(mutation, strain)

rec_plus <- final_df %>%
  filter(strain == "rec+") %>%
  select(-strain)

rec_minus <- final_df %>%
  filter(strain == "rec-") %>%
  select(-strain)

rec_plus_mix <- rec_plus %>%
  filter(condition == "MIX") %>%
  select(-condition)

rec_minus_mix <- rec_minus %>%
  filter(condition == "MIX") %>%
  select(-condition)

breseq_plot <- function(data, title, show_condition = TRUE) {
  # Define the custom color palette
  colors <- c(hsv(seq(0.15, 0, length.out = 5), 0.7, 1), 
                      hsv(seq(0.6, 0.75, length.out = 3), 0.7, 1))

  p <- ggplot(data, aes(x = mutation, y = freq, fill = mutation)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = conf_low, ymax = conf_high), width = 0.2) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_fill_manual(values = setNames(colors, mutations)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = "Frequency", fill = "Mutation", title = "Replicate") +
    theme(strip.text = element_text(size = 34),
          axis.text = element_text(size = 18),
          axis.title = element_text(size = 24),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18),
          plot.title = element_text(size = 24, hjust = 0.5))

  if(show_condition){
    p <- p + facet_grid(condition ~ rep)
  } else {
    p <- p + facet_grid(~ rep)
  }
  
  return(p)
}

setwd(paste0(parent_wd, "/Breseq_figures"))

pdf("rec+_all.pdf", width = 10, height = 16)
breseq_plot(rec_plus, "Strain: Rec+")
dev.off()

pdf("rec-_all.pdf", width = 10, height = 16)
breseq_plot(rec_minus, "Strain: Rec-")
dev.off()

pdf("rec+_mix.pdf", width = 10, height = 8)
breseq_plot(rec_plus_mix, "Strain: Rec+, Condition: MIX", show_condition = F)
dev.off()

pdf("rec-_mix.pdf", width = 10, height = 8)
breseq_plot(rec_minus_mix, "Strain: Rec-, Condition: MIX", show_condition = F)
dev.off()

setwd(parent_wd)

# sum the frequency of mutations for each replicate in rec_plus_mix, ignoring NAs
rec_plus_mix_sum <- rec_plus_mix %>%
  group_by(rep) %>%
  summarise(freq = sum(freq, na.rm = TRUE),
    conf_low = sum(conf_low, na.rm = TRUE)) %>%
  ungroup()
rec_plus_mix_sum
