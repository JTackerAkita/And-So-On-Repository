############################################################
# ICNALE Target Lexical Bundle Search
# Version 4: Adds functional categories
############################################################

packages <- c(
  "tidyverse",
  "stringr",
  "FactoMineR",
  "factoextra"
)

installed <- packages %in% rownames(installed.packages())

if(any(!installed)){
  install.packages(packages[!installed])
}

library(tidyverse)
library(stringr)
library(FactoMineR)
library(factoextra)

############################################################
# 1. Read TXT files from current working directory
############################################################

files <- list.files(pattern = "\\.txt$")

texts <- tibble(
  file = files,
  text = map_chr(files, ~ paste(readLines(.x, warn = FALSE), collapse = " "))
)

############################################################
# 2. Identify groups
############################################################

texts <- texts %>%
  mutate(
    group = case_when(
      str_detect(file, "^JPN") ~ "JPN",
      str_detect(file, "^CHN") ~ "CHN",
      str_detect(file, "^KOR") ~ "KOR",
      str_detect(file, "^ENS") ~ "ENS",
      TRUE ~ "OTHER"
    )
  )

############################################################
# 3. Clean text: remove POS tags and normalize
############################################################

texts <- texts %>%
  mutate(
    text = str_to_lower(text),
    text = str_replace_all(text, "_[a-z$]+", ""),
    text = str_replace_all(text, "[^a-z\\s']", " "),
    text = str_squish(text)
  )

############################################################
# 4. Target lexical bundles with categories
############################################################

target_tbl <- tribble(
  ~bundle, ~category,
  
  # stance / opinion
  "i think that", "stance",
  "i think it", "stance",
  "i believe that", "stance",
  "in my opinion", "stance",
  "i agree with", "stance",
  "i disagree with", "stance",
  "i do not", "stance",
  "we should not", "stance",
  "we have to", "stance",
  "it is important", "stance",
  "it is necessary", "stance",
  "it is better", "stance",
  "it is difficult", "stance",
  "there is no", "stance",
  "there are many", "stance",
  
  # organization / argument structure
  "first of all", "discourse_organizer",
  "on the other hand", "discourse_organizer",
  "as a result", "discourse_organizer",
  "to sum up", "discourse_organizer",
  "what is more", "discourse_organizer",
  
  # academic / explanatory phrasing
  "in order to", "referential",
  "be able to", "referential",
  "one of the", "referential",
  "the most important", "referential",
  "a lot of", "referential",
  "many kinds of", "referential",
  "more and more", "referential",
  "because of the", "referential",
  "due to the", "referential",
  
  # general extenders / vague language
  "and so on", "vague_extender",
  "and so forth", "vague_extender",
  "things like that", "vague_extender",
  "or something like", "vague_extender"
)

############################################################
# 5. Count total words by file and group
############################################################

word_counts <- texts %>%
  mutate(word_count = str_count(text, "\\S+")) %>%
  select(file, group, word_count)

group_word_counts <- word_counts %>%
  group_by(group) %>%
  summarise(total_words = sum(word_count), .groups = "drop")

############################################################
# 6. Count target bundle occurrences in each file
############################################################

count_bundle <- function(text, bundle){
  pattern <- paste0("\\b", str_replace_all(bundle, " ", "\\\\s+"), "\\b")
  str_count(text, regex(pattern))
}

bundle_counts_by_file <- texts %>%
  select(file, group, text) %>%
  crossing(target_tbl) %>%
  mutate(
    raw_frequency = map2_int(text, bundle, count_bundle)
  ) %>%
  select(file, group, category, bundle, raw_frequency)

############################################################
# 7. Summarise by group and bundle
############################################################

bundle_counts_by_group <- bundle_counts_by_file %>%
  group_by(group, category, bundle) %>%
  summarise(raw_frequency = sum(raw_frequency), .groups = "drop") %>%
  left_join(group_word_counts, by = "group") %>%
  mutate(
    normalized_per_100k = raw_frequency / total_words * 100000
  ) %>%
  arrange(group, category, desc(raw_frequency))

############################################################
# 8. Summarise by category
############################################################

category_counts_by_group <- bundle_counts_by_file %>%
  group_by(group, category) %>%
  summarise(raw_frequency = sum(raw_frequency), .groups = "drop") %>%
  left_join(group_word_counts, by = "group") %>%
  mutate(
    normalized_per_100k = raw_frequency / total_words * 100000
  ) %>%
  arrange(group, desc(raw_frequency))

############################################################
# 9. Save results
############################################################

if(!dir.exists("lexical_bundle_results")){
  dir.create("lexical_bundle_results")
}

write.csv(
  bundle_counts_by_file,
  "lexical_bundle_results/target_bundle_counts_by_file.csv",
  row.names = FALSE
)

write.csv(
  bundle_counts_by_group,
  "lexical_bundle_results/target_bundle_counts_by_group.csv",
  row.names = FALSE
)

write.csv(
  category_counts_by_group,
  "lexical_bundle_results/category_counts_by_group.csv",
  row.names = FALSE
)

write.csv(
  group_word_counts,
  "lexical_bundle_results/group_word_counts.csv",
  row.names = FALSE
)
############################################################
# 10. Inferential statistics: Category x Group
############################################################

# Create contingency table
category_table <- category_counts_by_group %>%
  select(group, category, raw_frequency) %>%
  pivot_wider(
    names_from = category,
    values_from = raw_frequency,
    values_fill = 0
  ) %>%
  column_to_rownames("group") %>%
  as.matrix()

# Chi-square test
chi_result <- chisq.test(category_table)

# Save chi-square result
chi_summary <- tibble(
  test = "Chi-square test of lexical bundle category by group",
  chi_square = as.numeric(chi_result$statistic),
  df = as.numeric(chi_result$parameter),
  p_value = chi_result$p.value
)

write.csv(
  chi_summary,
  "lexical_bundle_results/chi_square_category_by_group.csv",
  row.names = FALSE
)

# Standardized residuals
standardized_residuals <- as.data.frame(chi_result$stdres) %>%
  rownames_to_column("group") %>%
  pivot_longer(
    cols = -group,
    names_to = "category",
    values_to = "standardized_residual"
  ) %>%
  arrange(desc(abs(standardized_residual)))

write.csv(
  standardized_residuals,
  "lexical_bundle_results/standardized_residuals_category_by_group.csv",
  row.names = FALSE
)

# Cramer's V
n_total <- sum(category_table)
min_dim <- min(nrow(category_table) - 1, ncol(category_table) - 1)

cramers_v <- sqrt(as.numeric(chi_result$statistic) / (n_total * min_dim))

cramers_v_summary <- tibble(
  effect_size = "Cramer's V",
  value = cramers_v
)

write.csv(
  cramers_v_summary,
  "lexical_bundle_results/cramers_v_category_by_group.csv",
  row.names = FALSE
)

# Print results
cat("\nChi-square test: Category x Group\n")
print(chi_result)

cat("\nCramer's V:\n")
print(cramers_v)

cat("\nLargest standardized residuals:\n")
print(head(standardized_residuals, 10))
############################################################
# 10B. Pairwise chi-square tests with Holm correction
############################################################

groups <- unique(category_counts_by_group$group)

pairwise_results <- list()

counter <- 1

for(i in 1:(length(groups)-1)){
  
  for(j in (i+1):length(groups)){
    
    g1 <- groups[i]
    g2 <- groups[j]
    
    pair_data <- category_counts_by_group %>%
      filter(group %in% c(g1, g2)) %>%
      select(group, category, raw_frequency) %>%
      pivot_wider(
        names_from = category,
        values_from = raw_frequency,
        values_fill = 0
      ) %>%
      column_to_rownames("group") %>%
      as.matrix()
    
    pair_test <- chisq.test(pair_data)
    
    pairwise_results[[counter]] <- tibble(
      comparison = paste(g1, "vs", g2),
      chi_square = as.numeric(pair_test$statistic),
      df = as.numeric(pair_test$parameter),
      p_value = pair_test$p.value
    )
    
    counter <- counter + 1
  }
}

pairwise_results_df <- bind_rows(pairwise_results)

############################################################
# Holm correction
############################################################

pairwise_results_df <- pairwise_results_df %>%
  mutate(
    p_holm = p.adjust(p_value, method = "holm")
  ) %>%
  arrange(p_holm)

############################################################
# Save results
############################################################

write.csv(
  pairwise_results_df,
  "lexical_bundle_results/pairwise_chi_square_holm.csv",
  row.names = FALSE
)

############################################################
# Print results
############################################################

cat("\nPairwise chi-square tests with Holm correction:\n")
print(pairwise_results_df)

############################################################
# 11. Print summaries
############################################################

cat("\nTop bundles by group:\n")

top_target_bundles <- bundle_counts_by_group %>%
  filter(raw_frequency > 0) %>%
  group_by(group) %>%
  slice_max(raw_frequency, n = 20, with_ties = FALSE) %>%
  ungroup()

print(top_target_bundles)

cat("\nCategory totals by group:\n")
print(category_counts_by_group)

cat("\nDONE!\n")
cat("Files analyzed:", length(files), "\n")
cat("Results saved in: lexical_bundle_results\n")

############################################################
# 12. Create figures
############################################################

figures_folder <- file.path("lexical_bundle_results", "figures")

if(!dir.exists(figures_folder)){
  dir.create(figures_folder)
}

############################################################
# Figure 1: Category totals by group
############################################################

fig_category <- ggplot(
  category_counts_by_group,
  aes(x = group, y = normalized_per_100k, fill = category)
) +
  geom_col(position = "dodge", color = "black") +
  labs(
    title = "Lexical bundle categories by group",
    x = "Group",
    y = "Normalized frequency per 100,000 words",
    fill = "Category"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(figures_folder, "Figure_1_category_totals_by_group.png"),
  plot = fig_category,
  width = 9,
  height = 6,
  dpi = 300
)
############################################################
# Figure 1B: Heatmap of bundle categories by group
############################################################

fig_heatmap <- ggplot(
  category_counts_by_group,
  aes(x = category, y = group, fill = normalized_per_100k)
) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = round(normalized_per_100k, 1)),
    size = 4
  ) +
  labs(
    title = "Heatmap of lexical bundle categories by group",
    x = "Bundle category",
    y = "Group",
    fill = "Per 100,000 words"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "right"
  )

ggsave(
  filename = file.path(figures_folder, "Figure_1B_category_heatmap_by_group.png"),
  plot = fig_heatmap,
  width = 9,
  height = 5,
  dpi = 300
)
############################################################
# Figure 2: Top 15 bundles overall
############################################################

top_15_overall <- bundle_counts_by_group %>%
  group_by(bundle, category) %>%
  summarise(
    total_frequency = sum(raw_frequency),
    mean_norm_per_100k = mean(normalized_per_100k),
    .groups = "drop"
  ) %>%
  arrange(desc(total_frequency)) %>%
  slice_head(n = 15)

fig_top15 <- ggplot(
  top_15_overall,
  aes(x = reorder(bundle, mean_norm_per_100k), y = mean_norm_per_100k)
) +
  geom_col(color = "black", fill = "gray70") +
  coord_flip() +
  labs(
    title = "Top 15 lexical bundles overall",
    x = "Lexical bundle",
    y = "Mean normalized frequency per 100,000 words"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(figures_folder, "Figure_2_top_15_bundles_overall.png"),
  plot = fig_top15,
  width = 9,
  height = 6,
  dpi = 300
)

############################################################
# Figure 3: 'and so on' by group
############################################################

and_so_on_data <- bundle_counts_by_group %>%
  filter(bundle == "and so on")

fig_and_so_on <- ggplot(
  and_so_on_data,
  aes(x = group, y = normalized_per_100k)
) +
  geom_col(color = "black", fill = "gray70") +
  labs(
    title = "Use of 'and so on' by group",
    x = "Group",
    y = "Normalized frequency per 100,000 words"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(figures_folder, "Figure_3_and_so_on_by_group.png"),
  plot = fig_and_so_on,
  width = 7,
  height = 5,
  dpi = 300
)
############################################################
# Figure 4: Top 10 bundles by group
############################################################

top_10_by_group <- bundle_counts_by_group %>%
  filter(raw_frequency > 0) %>%
  group_by(group) %>%
  slice_max(normalized_per_100k, n = 10, with_ties = FALSE) %>%
  ungroup()

fig_top10_group <- ggplot(
  top_10_by_group,
  aes(x = reorder(bundle, normalized_per_100k), y = normalized_per_100k)
) +
  geom_col(color = "black", fill = "gray70") +
  coord_flip() +
  facet_wrap(~ group, scales = "free_y") +
  labs(
    title = "Top 10 lexical bundles by group",
    x = "Lexical bundle",
    y = "Normalized frequency per 100,000 words"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

ggsave(
  filename = file.path(figures_folder, "Figure_4_top_10_bundles_by_group.png"),
  plot = fig_top10_group,
  width = 11,
  height = 8,
  dpi = 300
)

cat("\nFigures saved in: lexical_bundle_results/figures\n")

############################################################
# 13. Correspondence Analysis
############################################################

# Correspondence Analysis
ca_result <- FactoMineR::CA(category_table, graph = FALSE)

############################################################
# Save coordinates
############################################################

write.csv(
  as.data.frame(ca_result$row$coord),
  "lexical_bundle_results/correspondence_analysis_group_coordinates.csv",
  row.names = TRUE
)

write.csv(
  as.data.frame(ca_result$col$coord),
  "lexical_bundle_results/correspondence_analysis_category_coordinates.csv",
  row.names = TRUE
)

############################################################
# Eigenvalues
############################################################

write.csv(
  as.data.frame(ca_result$eig),
  "lexical_bundle_results/correspondence_analysis_eigenvalues.csv",
  row.names = TRUE
)

############################################################
# Figure 5
############################################################

fig_ca <- factoextra::fviz_ca_biplot(
  ca_result,
  repel = TRUE,
  title = "Correspondence Analysis of Writer Groups and Lexical Bundle Categories"
)

ggsave(
  filename = file.path(
    figures_folder,
    "Figure_5_correspondence_analysis.png"
  ),
  plot = fig_ca,
  width = 8,
  height = 6,
  dpi = 300
)

cat("\nCorrespondence Analysis completed.\n")


############################################################
# 13. Correspondence Analysis
############################################################

# Correspondence Analysis provides a visual summary of the
# association between writer groups and lexical bundle categories.
# It uses the same group x category contingency table used for
# the chi-square analysis.

ca_result <- CA(category_table, graph = FALSE)

############################################################
# Save Correspondence Analysis outputs
############################################################

# Eigenvalues / explained inertia
ca_eigenvalues <- as.data.frame(ca_result$eig) %>%
  rownames_to_column("dimension") %>%
  rename(
    eigenvalue = eigenvalue,
    percent_inertia = `percentage of variance`,
    cumulative_percent_inertia = `cumulative percentage of variance`
  )

write.csv(
  ca_eigenvalues,
  "lexical_bundle_results/correspondence_analysis_eigenvalues.csv",
  row.names = FALSE
)

# Row coordinates: writer groups
ca_group_coordinates <- as.data.frame(ca_result$row$coord) %>%
  rownames_to_column("group")

write.csv(
  ca_group_coordinates,
  "lexical_bundle_results/correspondence_analysis_group_coordinates.csv",
  row.names = FALSE
)

# Column coordinates: lexical bundle categories
ca_category_coordinates <- as.data.frame(ca_result$col$coord) %>%
  rownames_to_column("category")

write.csv(
  ca_category_coordinates,
  "lexical_bundle_results/correspondence_analysis_category_coordinates.csv",
  row.names = FALSE
)

############################################################
# Figure 5: Correspondence Analysis biplot
############################################################

fig_ca <- fviz_ca_biplot(
  ca_result,
  repel = TRUE,
  title = "Correspondence analysis of writer groups and lexical bundle categories"
) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(figures_folder, "Figure_5_correspondence_analysis.png"),
  plot = fig_ca,
  width = 9,
  height = 6,
  dpi = 300
)

cat("\nCorrespondence Analysis complete.\n")
cat("CA files saved in: lexical_bundle_results\n")
cat("Figure 5 saved in: lexical_bundle_results/figures\n")

