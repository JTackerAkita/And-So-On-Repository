# =============================
# "and so on" analysis — v4.3 (XXL Times New Roman)
# - Font sizes per request:
#   base 28 pt, title 30 pt, axis/legend 26 pt, caption 22 pt
# - Keeps all small-sample–safe tests and logic from v4.2
# =============================

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(readr)
  library(scales)
  library(DescTools)
  library(broom)
  library(broom.mixed)
  library(lme4)
  library(MASS)
  library(binom)
  library(emmeans)
  library(showtext)
  library(sysfonts)
})

# ---------- CONFIG ----------
input_dir  <- "."
output_dir <- "and_so_on_results"
font_family <- "Times New Roman"

# >>> Font sizes (per your request) <<<
base_pt    <- 28
title_pt   <- 30
axis_pt    <- 26
legend_pt  <- 26
caption_pt <- 22

# Figure sizes (same as v4.2; already generous for these font sizes)
size_p1 <- c(12, 9)
size_p2 <- c(12, 9)
size_p3 <- c(13, 10)
size_p4 <- c(10, 8)
size_p5 <- c(11, 8)

rate_scale <- 100000
MC_B       <- 50000

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---------- FONTS ----------
try({
  if (!(font_family %in% sysfonts::font_families())) {
    if (.Platform$OS.type == "windows") {
      sysfonts::font_add(family = font_family, regular = "C:/Windows/Fonts/times.ttf")
    } else {
      sysfonts::font_add(family = font_family, regular = "/System/Library/Fonts/Supplemental/Times New Roman.ttf")
    }
  }
}, silent = TRUE)
showtext::showtext_auto(enable = TRUE)
available_families <- tryCatch(sysfonts::font_families(), error = function(e) character())
if (!(font_family %in% available_families)) {
  message("Times New Roman not found; falling back to 'serif' for plotting.")
  font_family <- "serif"
}

# ---------- OUTPUT NAMES ----------
output_csv         <- file.path(output_dir, "and_so_on_stats_final_vs_nonfinal.csv")
plot_pmi_final     <- file.path(output_dir, "and_so_on_histogram_pmi_FINAL.png")
plot_trigram_fin   <- file.path(output_dir, "and_so_on_histogram_trigram_FINAL.png")
plot_final_non     <- file.path(output_dir, "and_so_on_final_vs_nonfinal_grouped.png")
plot_stacked100    <- file.path(output_dir, "and_so_on_100pct_stacked_by_culture_FINAL.png")
plot_cult_finalnon <- file.path(output_dir, "and_so_on_final_vs_nonfinal_by_culture.png")

# Statistical outputs (unchanged)
csv_group_summary             <- file.path(output_dir, "group_summary_andso.csv")
csv_tests_freq_chisq          <- file.path(output_dir, "tests_frequency_main.csv")
csv_tests_finalshare_chisq    <- file.path(output_dir, "tests_finalshare_main.csv")
csv_tests_freq_chisq_cells    <- file.path(output_dir, "tests_frequency_chisq_stdresiduals.csv")
csv_tests_finalshare_cells    <- file.path(output_dir, "tests_finalshare_chisq_stdresiduals.csv")
csv_pairwise_freq             <- file.path(output_dir, "pairwise_frequency_proptests_or_fisher.csv")
csv_pairwise_finalshare       <- file.path(output_dir, "pairwise_finalshare_proptests_or_fisher.csv")
csv_pairwise_rate_poisson     <- file.path(output_dir, "pairwise_frequency_poisson_tests.csv")
csv_glmer_freq                <- file.path(output_dir, "glmer_frequency_summary.csv")
csv_glmer_finalshare          <- file.path(output_dir, "glmer_finalshare_summary.csv")
csv_glm_poisson_nb_freq       <- file.path(output_dir, "glm_rate_models_RQ1.csv")
csv_glm_binomial_finalshare   <- file.path(output_dir, "glm_binomial_RQ2.csv")
csv_emm_glm_rate_pairs        <- file.path(output_dir, "emm_glm_rate_pairs_RQ1.csv")
csv_emm_glm_final_pairs       <- file.path(output_dir, "emm_glm_final_pairs_RQ2.csv")
csv_emm_glmer_freq_pairs      <- file.path(output_dir, "emm_glmer_freq_pairs.csv")
csv_emm_glmer_final_pairs     <- file.path(output_dir, "emm_glmer_final_pairs.csv")

# ---------- HELPERS (same as v4.2) ----------
normalize_punct <- function(tok) {
  t <- tolower(tok)
  t <- str_replace_all(t, fixed("。"), ".")
  t <- str_replace_all(t, fixed("．"), ".")
  t <- str_replace_all(t, fixed("！"), "!")
  t <- str_replace_all(t, fixed("？"), "?")
  t <- str_replace_all(t, "[…‥⋯]", "…")
  t
}

read_annotated_tokens <- function(path) {
  txt <- read_file(path)
  raw_tokens <- str_split(txt, "\\s+")[[1]]
  raw_tokens <- raw_tokens[nzchar(raw_tokens)]
  surface <- ifelse(raw_tokens %in% c(".", "。", "．", "！", "？", "…", "‥", "⋯"),
                    raw_tokens,
                    str_replace(raw_tokens, "^(.*?)(?:_.+)?$", "\\1"))
  normalize_punct(surface)
}

count_and_so_on_total <- function(tokens) {
  n <- length(tokens)
  if (n < 3) return(0L)
  sum(tokens[1:(n-2)] == "and" & tokens[2:(n-1)] == "so" & tokens[3:n] == "on")
}

END_TOKENS <- c(".", "!", "?", "…")
CLOSERS    <- c("\"", "”", "’", ")", "]", "】", "〉", "》", "」", "』", "＞")

count_and_so_on_sentence_final <- function(tokens, max_closers = 3L) {
  n <- length(tokens)
  if (n < 4) return(0L)
  idx <- which(tokens[seq_len(n-2)] == "and" & tokens[seq(2, n-1)] == "so" & tokens[seq(3, n)] == "on")
  finals <- 0L
  for (i in idx) {
    j <- i + 3L
    k <- 0L
    while (j <= n && k < max_closers && tokens[j] %in% CLOSERS) {
      j <- j + 1L; k <- k + 1L
    }
    if (j <= n && tokens[j] %in% END_TOKENS) finals <- finals + 1L
  }
  finals
}

infer_culture <- function(fname) {
  f <- basename(fname)
  case_when(
    str_starts(f, regex("^JPN", ignore_case = TRUE)) ~ "JPN",
    str_starts(f, regex("^CHN", ignore_case = TRUE)) ~ "CHN",
    str_starts(f, regex("^KOR", ignore_case = TRUE)) ~ "KOR",
    str_starts(f, regex("^ENS", ignore_case = TRUE)) ~ "ENS",
    TRUE ~ "OTHER"
  )
}

# ---------- PER-FILE STATS ----------
stats_for_file <- function(path) {
  message("Processing: ", basename(path))
  tokens <- read_annotated_tokens(path)
  N <- length(tokens)
  if (N == 0) warning("No tokens found in file: ", basename(path))
  
  c_and <- sum(tokens == "and")
  c_so  <- sum(tokens == "so")
  c_on  <- sum(tokens == "on")
  
  c_total  <- count_and_so_on_total(tokens)
  c_final  <- count_and_so_on_sentence_final(tokens)
  c_nonfin <- c_total - c_final
  
  denom_trigram <- max(N - 2, 0)
  
  p_and <- if (N > 0) c_and / N else NA_real_
  p_so  <- if (N > 0) c_so  / N else NA_real_
  p_on  <- if (N > 0) c_on  / N else NA_real_
  
  p_total  <- if (denom_trigram > 0) c_total  / denom_trigram else NA_real_
  p_final  <- if (denom_trigram > 0) c_final  / denom_trigram else NA_real_
  p_nonfin <- if (denom_trigram > 0) c_nonfin / denom_trigram else NA_real_
  
  denom_unigrams <- p_and * p_so * p_on
  pmi_final <- if (!is.na(p_final) && !is.na(denom_unigrams) && denom_unigrams > 0) p_final / denom_unigrams else NA_real_
  pmi_total <- if (!is.na(p_total) && !is.na(denom_unigrams) && denom_unigrams > 0) p_total / denom_unigrams else NA_real_
  
  fname <- basename(path)
  tibble(
    file = fname,
    culture = infer_culture(fname),
    N_tokens = N,
    trigram_positions = denom_trigram,
    # counts
    count_and = c_and,
    count_so  = c_so,
    count_on  = c_on,
    count_and_so_on_total    = c_total,
    count_and_so_on_final    = c_final,
    count_and_so_on_nonfinal = c_nonfin,
    # probabilities
    p_and = p_and, p_so = p_so, p_on = p_on,
    p_and_so_on_total    = p_total,
    p_and_so_on_final    = p_final,
    p_and_so_on_nonfinal = p_nonfin,
    # PMI-style ratios
    pmi_final = pmi_final,
    pmi_total = pmi_total,
    # rates per 100k tokens
    rate_total_per_100k  = ifelse(N_tokens > 0, c_total  / N_tokens * rate_scale, NA_real_),
    rate_final_per_100k  = ifelse(N_tokens > 0, c_final  / N_tokens * rate_scale, NA_real_),
    rate_nonfin_per_100k = ifelse(N_tokens > 0, c_nonfin / N_tokens * rate_scale, NA_real_)
  )
}

# ---------- RUN ----------
paths <- list.files(
  path = input_dir,
  pattern = "\\.(txt|csv)$",
  full.names = TRUE,
  ignore.case = TRUE
)
if (length(paths) == 0) stop("No .txt or .csv files found in input_dir: ", normalizePath(input_dir))

results <- purrr::map_dfr(paths, stats_for_file)
write_csv(results, output_csv)

# ---------- GROUP SUMMARY ----------
group_sum <- results %>%
  group_by(culture) %>%
  summarise(
    files                = n(),
    tokens               = sum(N_tokens, na.rm = TRUE),
    trigram_positions    = sum(trigram_positions, na.rm = TRUE),
    andso_total          = sum(count_and_so_on_total,    na.rm = TRUE),
    andso_final          = sum(count_and_so_on_final,    na.rm = TRUE),
    andso_nonfinal       = sum(count_and_so_on_nonfinal, na.rm = TRUE),
    rate_total_per_100k  = ifelse(tokens > 0, andso_total   / tokens * rate_scale, NA_real_),
    rate_final_per_100k  = ifelse(tokens > 0, andso_final   / tokens * rate_scale, NA_real_),
    rate_nonfin_per_100k = ifelse(tokens > 0, andso_nonfinal/ tokens * rate_scale, NA_real_),
    prop_final_share     = ifelse(andso_total > 0, andso_final / andso_total, NA_real_)
  ) %>%
  ungroup()

# Wilson 95% CI for final share
wilson_ci <- group_sum %>%
  mutate(ci = purrr::pmap(list(andso_final, andso_total), function(x, n){
    if (is.na(x) || is.na(n) || n == 0) return(c(NA_real_, NA_real_))
    w <- binom::binom.wilson(x, n)
    c(w$lower, w$upper)
  })) %>%
  mutate(ci_lower = purrr::map_dbl(ci, 1),
         ci_upper = purrr::map_dbl(ci, 2)) %>%
  dplyr::select(-ci)
group_sum <- group_sum %>% left_join(wilson_ci %>% dplyr::select(culture, ci_lower, ci_upper), by="culture")
write_csv(group_sum, csv_group_summary)

# ---------- SAFE TEST HELPERS ----------
safe_chisq <- function(tab, B = MC_B) {
  chisq_plain <- suppressWarnings(chisq.test(tab))
  exp_counts <- chisq_plain$expected
  need_mc <- any(exp_counts < 5) || mean(exp_counts < 5) > 0.2
  if (need_mc) {
    chisq_mc <- chisq.test(tab, simulate.p.value = TRUE, B = B)
    list(res = chisq_mc, method = paste0("Chi-square (Monte Carlo, B=", B, ")"), expected = exp_counts)
  } else {
    list(res = chisq_plain, method = "Chi-square (asymptotic)", expected = exp_counts)
  }
}

safe_pairwise_prop_or_fisher <- function(a_yes, a_no, b_yes, b_no) {
  mat <- matrix(c(a_yes, a_no, b_yes, b_no), nrow = 2, byrow = TRUE)
  if (any(mat < 5)) {
    ft <- fisher.test(mat)
    list(p = ft$p.value, stat = NA_real_, method = "Fisher exact")
  } else {
    pt <- prop.test(x = c(a_yes, b_yes), n = c(a_yes + a_no, b_yes + b_no), correct = FALSE)
    list(p = pt$p.value, stat = unname(pt$statistic), method = "prop.test (chi-square)")
  }
}

# ---------- RQ1: FREQUENCY ----------
freq_table <- group_sum %>% transmute(culture, yes = andso_total, no = pmax(trigram_positions - andso_total, 0))
sc <- safe_chisq(as.matrix(freq_table %>% dplyr::select(yes, no)))
chisq_freq <- sc$res
cramer_freq <- DescTools::CramerV(as.matrix(freq_table %>% dplyr::select(yes, no)))
freq_out <- tibble(test = "Frequency across groups",
                   method = sc$method,
                   statistic = unname(chisq_freq$statistic),
                   df = unname(chisq_freq$parameter),
                   p_value = unname(chisq_freq$p.value),
                   cramer_v = unname(cramer_freq))
write_csv(freq_out, csv_tests_freq_chisq)

stdres_freq <- as_tibble(chisq_freq$stdres, .name_repair = "minimal")
colnames(stdres_freq) <- c("yes_stdres","no_stdres")
stdres_freq <- bind_cols(freq_table %>% dplyr::select(culture), stdres_freq)
write_csv(stdres_freq, csv_tests_freq_chisq_cells)

# Pairwise prop/fisher with Holm; add effect sizes
pairwise_freq <- {
  combn(freq_table$culture, 2, simplify = FALSE, FUN = function(pair){
    g1 <- pair[1]; g2 <- pair[2]
    r1 <- freq_table %>% filter(culture == g1)
    r2 <- freq_table %>% filter(culture == g2)
    a_yes <- r1$yes; a_no <- r1$no; b_yes <- r2$yes; b_no <- r2$no
    sp <- safe_pairwise_prop_or_fisher(a_yes,a_no,b_yes,b_no)
    p1 <- a_yes/(a_yes + a_no); p2 <- b_yes/(b_yes + b_no)
    rr <- p1/p2; rd <- p1 - p2
    tibble(group1=g1, group2=g2, method=sp$method,
           statistic=sp$stat, p_value_raw=sp$p,
           p1=p1, p2=p2, rate_ratio=rr, risk_diff=rd)
  }) %>% bind_rows() %>%
  mutate(p_value_holm = p.adjust(p_value_raw, method = "holm")) %>%
  arrange(p_value_holm)
}
write_csv(pairwise_freq, csv_pairwise_freq)

# Pairwise Poisson tests for rates (per tokens)
rate_df <- results %>% group_by(file, culture) %>%
  summarise(tokens = sum(N_tokens, na.rm = TRUE),
            andso_total = sum(count_and_so_on_total, na.rm = TRUE), .groups = "drop") %>%
  filter(tokens > 0)
pair_poisson <- {
  groups <- unique(rate_df$culture)
  combn(groups, 2, simplify = FALSE, FUN = function(pair){
    a <- rate_df %>% filter(culture == pair[1])
    b <- rate_df %>% filter(culture == pair[2])
    x <- c(sum(a$andso_total), sum(b$andso_total))
    T <- c(sum(a$tokens), sum(b$tokens))
    pt <- poisson.test(x, T)
    rr <- (x[1]/T[1])/(x[2]/T[2])
    tibble(group1=pair[1], group2=pair[2], method="poisson.test",
           rate1 = x[1]/T[1], rate2 = x[2]/T[2],
           rate_ratio = rr,
           conf.low = exp(log(rr) - 1.96*sqrt(1/max(x[1],1) + 1/max(x[2],1))),
           conf.high = exp(log(rr) + 1.96*sqrt(1/max(x[1],1) + 1/max(x[2],1))),
           p_value_raw = pt$p.value)
  }) %>% bind_rows() %>%
  mutate(p_value_holm = p.adjust(p_value_raw, method = "holm")) %>%
  arrange(p_value_holm)
}
write_csv(pair_poisson, csv_pairwise_rate_poisson)

# Poisson/NB GLM with offset
pois_fit <- glm(andso_total ~ culture + offset(log(tokens)), data = rate_df, family = poisson)
dispersion <- sum(residuals(pois_fit, type="pearson")^2) / pois_fit$df.residual
if (is.finite(dispersion) && dispersion > 1.5) {
  nb_fit <- MASS::glm.nb(andso_total ~ culture + offset(log(tokens)), data = rate_df)
  rate_model <- nb_fit; model_type <- "Negative Binomial"
} else {
  rate_model <- pois_fit; model_type <- "Poisson"
}
emm_rate <- emmeans(rate_model, ~ culture, type = "response")
emm_pairs_rate <- pairs(emm_rate, adjust = "tukey")
rate_summary <- broom::tidy(rate_model, conf.int = TRUE, exponentiate = TRUE) %>%
  mutate(model = model_type, note = "Exponentiated coefficients are rate ratios (per tokens via offset)")
write_csv(rate_summary, csv_glm_poisson_nb_freq)
write_csv(as_tibble(emm_pairs_rate), csv_emm_glm_rate_pairs)

# ---------- RQ2: FINAL SHARE ----------
final_table <- group_sum %>% transmute(culture, final = andso_final, nonfinal = andso_nonfinal)
sc2 <- safe_chisq(as.matrix(final_table %>% dplyr::select(final, nonfinal)))
chisq_final <- sc2$res
cramer_final <- DescTools::CramerV(as.matrix(final_table %>% dplyr::select(final, nonfinal)))
final_out <- tibble(test = "Final share across groups",
                    method = sc2$method,
                    statistic = unname(chisq_final$statistic),
                    df = unname(chisq_final$parameter),
                    p_value = unname(chisq_final$p.value),
                    cramer_v = unname(cramer_final))
write_csv(final_out, csv_tests_finalshare_chisq)

stdres_final <- as_tibble(chisq_final$stdres, .name_repair = "minimal")
colnames(stdres_final) <- c("final_stdres","nonfinal_stdres")
stdres_final <- bind_cols(final_table %>% dplyr::select(culture), stdres_final)
write_csv(stdres_final, csv_tests_finalshare_cells)

# Pairwise for final share: Fisher when needed
pairwise_final <- {
  combn(final_table$culture, 2, simplify = FALSE, FUN = function(pair){
    g1 <- pair[1]; g2 <- pair[2]
    r1 <- final_table %>% filter(culture == g1)
    r2 <- final_table %>% filter(culture == g2)
    a_yes <- r1$final; a_no <- r1$nonfinal; b_yes <- r2$final; b_no <- r2$nonfinal
    sp <- safe_pairwise_prop_or_fisher(a_yes,a_no,b_yes,b_no)
    p1 <- a_yes/(a_yes + a_no); p2 <- b_yes/(b_yes + b_no)
    or <- (p1/(1-p1)) / (p2/(1-p2))
    rd <- p1 - p2
    tibble(group1=g1, group2=g2, method=sp$method,
           statistic=sp$stat, p_value_raw=sp$p,
           p1=p1, p2=p2, odds_ratio=or, risk_diff=rd)
  }) %>% bind_rows() %>%
  mutate(p_value_holm = p.adjust(p_value_raw, method = "holm")) %>%
  arrange(p_value_holm)
}
write_csv(pairwise_final, csv_pairwise_finalshare)

# Binomial GLM + emmeans
final_df <- results %>%
  transmute(file, culture, final = count_and_so_on_final, nonfinal = count_and_so_on_nonfinal) %>%
  filter((final + nonfinal) > 0)
glm_final <- glm(cbind(final, nonfinal) ~ culture, data = final_df, family = binomial)
emm_final <- emmeans(glm_final, ~ culture, type = "response")
emm_pairs_final <- pairs(emm_final, adjust = "tukey")
glm_final_tidy <- broom::tidy(glm_final, conf.int = TRUE, exponentiate = TRUE) %>%
  mutate(note = "Exponentiated coefficients are odds ratios (final share)")
write_csv(glm_final_tidy, csv_glm_binomial_finalshare)
write_csv(as_tibble(emm_pairs_final), csv_emm_glm_final_pairs)

# ---------- MIXED-EFFECTS LOGISTIC WITH EMMEANS ----------
glmer_data_freq <- results %>%
  transmute(file, culture, yes = count_and_so_on_total, no = pmax(trigram_positions - count_and_so_on_total, 0)) %>%
  filter((yes + no) > 0)
model_freq <- tryCatch(
  glmer(cbind(yes, no) ~ culture + (1 | file), data = glmer_data_freq, family = binomial),
  error = function(e) { message("GLMER freq failed: ", e$message); NULL }
)
if (!is.null(model_freq)) {
  tidy_freq <- broom.mixed::tidy(model_freq, effects = "fixed", conf.int = TRUE, conf.method = "Wald") %>%
    mutate(odds_ratio = exp(estimate), conf.low.OR = exp(conf.low), conf.high.OR = exp(conf.high))
  write_csv(tidy_freq, csv_glmer_freq)
  emm_mfreq <- emmeans(model_freq, ~ culture, type = "response")
  write_csv(as_tibble(pairs(emm_mfreq, adjust = "tukey")), csv_emm_glmer_freq_pairs)
}

glmer_data_final <- results %>%
  transmute(file, culture, final = count_and_so_on_final, nonfinal = count_and_so_on_nonfinal) %>%
  filter((final + nonfinal) > 0)
model_final <- tryCatch(
  glmer(cbind(final, nonfinal) ~ culture + (1 | file), data = glmer_data_final, family = binomial),
  error = function(e) { message("GLMER final failed: ", e$message); NULL }
)
if (!is.null(model_final)) {
  tidy_final <- broom.mixed::tidy(model_final, effects = "fixed", conf.int = TRUE, conf.method = "Wald") %>%
    mutate(odds_ratio = exp(estimate), conf.low.OR = exp(conf.low), conf.high.OR = exp(conf.high))
  write_csv(tidy_final, csv_glmer_finalshare)
  emm_mfinal <- emmeans(model_final, ~ culture, type = "response")
  write_csv(as_tibble(pairs(emm_mfinal, adjust = "tukey")), csv_emm_glmer_final_pairs)
}

# ---------- PLOTTING THEME (XXL TNR) ----------
theme_set(theme_minimal(base_size = base_pt, base_family = font_family))
big_theme <- theme(
  plot.title   = element_text(size = title_pt, face = "bold"),
  axis.title   = element_text(size = base_pt),
  axis.text    = element_text(size = axis_pt),
  legend.title = element_text(size = legend_pt),
  legend.text  = element_text(size = legend_pt),
  plot.caption = element_text(size = caption_pt)
)

# (1) PMI-style ratio (FINAL) per file
plot_df_pmi <- results %>% mutate(file_label = forcats::fct_reorder(file, pmi_final, .fun = identity, .na_rm = TRUE))
p1 <- ggplot(plot_df_pmi, aes(x = file_label, y = pmi_final, fill = culture)) +
  geom_col() + coord_flip() +
  labs(title = 'PMI-style ratio (FINAL): P(and so on FINAL) / [P(and)·P(so)·P(on)]', x = "File", y = "PMI-style ratio (final)") +
  big_theme
ggsave(plot_pmi_final, p1, width = size_p1[1], height = size_p1[2], dpi = 300, units = "in")

# (2) Raw trigram probability (FINAL) per file
plot_df_fin <- results %>% mutate(file_label = forcats::fct_reorder(file, p_and_so_on_final, .fun = identity, .na_rm = TRUE))
p2 <- ggplot(plot_df_fin, aes(x = file_label, y = p_and_so_on_final, fill = culture)) +
  geom_col() + coord_flip() +
  labs(title = 'Raw probability (FINAL): P("and so on") at sentence end', x = "File", y = "P(and so on FINAL)") +
  big_theme
ggsave(plot_trigram_fin, p2, width = size_p2[1], height = size_p2[2], dpi = 300, units = "in")

# (3) Final vs Non-final grouped bars per file
plot_long <- results %>%
  dplyr::select(file, culture, p_and_so_on_final, p_and_so_on_nonfinal) %>%
  pivot_longer(cols = c(p_and_so_on_final, p_and_so_on_nonfinal), names_to = "position", values_to = "prob") %>%
  mutate(position = recode(position, p_and_so_on_final = "Sentence-final", p_and_so_on_nonfinal = "Non-final"),
         file_label = forcats::fct_reorder(file, prob, .fun = max, .na_rm = TRUE))
p3 <- ggplot(plot_long, aes(x = file_label, y = prob, fill = position)) +
  geom_col(position = position_dodge(width = 0.75)) + coord_flip() +
  labs(title = 'Final vs. Non-final "and so on" (probability per file)', x = "File", y = "Probability") +
  big_theme
ggsave(plot_final_non, p3, width = size_p3[1], height = size_p3[2], dpi = 300, units = "in")

# (4) 100% stacked: culture share of ALL sentence-final occurrences
culture_summary <- results %>% group_by(culture) %>%
  summarise(total_and_so_on_final = sum(count_and_so_on_final, na.rm = TRUE),
            total_trigram_positions = sum(trigram_positions, na.rm = TRUE), .groups = "drop")
total_all_final <- sum(culture_summary$total_and_so_on_final, na.rm = TRUE)
culture_share <- culture_summary %>% mutate(share = if (total_all_final > 0) total_and_so_on_final / total_all_final else NA_real_)
p4 <- ggplot(culture_share, aes(x = "All files", y = share, fill = culture)) +
  geom_col() + scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = 'Culture share of sentence-final "and so on" (100% stacked, aggregated)', x = NULL, y = "Share of total (100%)") +
  big_theme
ggsave(plot_stacked100, p4, width = size_p4[1], height = size_p4[2], dpi = 300, units = "in")

# (5) Final vs Non-final by culture (probabilities)
culture_agg <- results %>%
  group_by(culture) %>%
  summarise(final_count = sum(count_and_so_on_final, na.rm = TRUE),
            nonfinal_count = sum(count_and_so_on_nonfinal, na.rm = TRUE),
            trigram_positions = sum(trigram_positions, na.rm = TRUE), .groups = "drop") %>%
  mutate(p_final = ifelse(trigram_positions > 0, final_count / trigram_positions, NA_real_),
         p_nonfinal = ifelse(trigram_positions > 0, nonfinal_count / trigram_positions, NA_real_)) %>%
  filter(culture %in% c("JPN", "CHN", "KOR", "ENS"))
culture_long <- culture_agg %>% dplyr::select(culture, p_final, p_nonfinal) %>%
  pivot_longer(cols = starts_with("p_"), names_to = "position", values_to = "prob") %>%
  mutate(position = recode(position, p_final = "Sentence-final", p_nonfinal = "Non-final"))
p5 <- ggplot(culture_long, aes(x = culture, y = prob, fill = position)) +
  geom_col(position = position_dodge(width = 0.75)) +
  labs(title = 'Final vs. Non-final "and so on" by culture (probability)', x = "Culture", y = "Probability") +
  big_theme
ggsave(plot_cult_finalnon, p5, width = size_p5[1], height = size_p5[2], dpi = 300, units = "in")

# ---------- DONE ----------
print(results)
message("Saved results and plots with XXL Times New Roman to: ", normalizePath(output_dir))
