############################################################
# ICNALE Target Lexical Bundle Search
# Version 3: Academic / argumentative / discourse bundles
############################################################

packages <- c("tidyverse", "stringr")

installed <- packages %in% rownames(installed.packages())

if(any(!installed)){
  install.packages(packages[!installed])
}

library(tidyverse)
library(stringr)

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
# 4. Target lexical bundles
############################################################

target_bundles <- c(
  # stance / opinion
  "i think that",
  "i think it",
  "i believe that",
  "in my opinion",
  "i agree with",
  "i disagree with",
  "i do not",
  "we should not",
  "we have to",
  "it is important",
  "it is necessary",
  "it is better",
  "it is difficult",
  "there is no",
  "there are many",
  
  # organization / argument structure
  "first of all",
  "on the other hand",
  "as a result",
  "for example",
  "for instance",
  "in conclusion",
  "to sum up",
  "in addition",
  "what is more",
  "moreover",
  "therefore",
  "however",
  
  # academic / explanatory phrasing
  "in order to",
  "be able to",
  "one of the",
  "the most important",
  "a lot of",
  "many kinds of",
  "more and more",
  "because of the",
  "due to the",
  "according to",
  
  # general extenders / vague language
  "and so on",
  "and so forth",
  "things like that",
  "or something like",
  "such as",
  "etc"
)

target_tbl <- tibble(bundle = target_bundles)

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
  select(file, group, bundle, raw_frequency)

############################################################
# 7. Summarise by group
############################################################

bundle_counts_by_group <- bundle_counts_by_file %>%
  group_by(group, bundle) %>%
  summarise(raw_frequency = sum(raw_frequency), .groups = "drop") %>%
  left_join(group_word_counts, by = "group") %>%
  mutate(
    normalized_per_100k = raw_frequency / total_words * 100000
  ) %>%
  arrange(group, desc(raw_frequency))

############################################################
# 8. Save results
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
  group_word_counts,
  "lexical_bundle_results/group_word_counts.csv",
  row.names = FALSE
)

############################################################
# 9. Print useful summary
############################################################

top_target_bundles <- bundle_counts_by_group %>%
  filter(raw_frequency > 0) %>%
  group_by(group) %>%
  slice_max(raw_frequency, n = 20, with_ties = FALSE) %>%
  ungroup()

print(top_target_bundles)

cat("\nDONE!\n")
cat("Files analyzed:", length(files), "\n")
cat("Results saved in: lexical_bundle_results\n")