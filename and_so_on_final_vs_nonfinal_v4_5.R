# =============================
# "and so on" analysis — v4.5 (extras-only layer)
# Adds ONLY: KW + Dunn, Wilson & Newcombe–Wilson pair CIs,
# Permutation MANOVA, Sensitivity (max_closers), GLMER + EMMEANS.
# Reads v4.3/4.4 results when available; otherwise recomputes per-file stats.
# Robust: writes a 'note' CSV if a step can't run (pkg/data).
# =============================

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(readr)
  library(broom)
  library(broom.mixed)
  library(binom)
  library(lme4)
  library(emmeans)
})

# ---------- helpers ----------
has_pkg   <- function(x) requireNamespace(x, quietly = TRUE)
emit_note <- function(path, note) {
  if (!dir.exists(dirname(path))) dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(tibble(note = note), path)
  message("Wrote note -> ", basename(path), ": ", note)
}

# ---------- config ----------
input_dir  <- "."
output_dir <- "and_so_on_results"
rate_scale <- 100000
ALPHA      <- 0.05
N_PERM     <- 4999

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---------- output filenames (v4.5 suffix) ----------
csv_kw_rate_total             <- file.path(output_dir, "kw_rate_total_per_file_v45.csv")
csv_dunn_rate_total           <- file.path(output_dir, "dunn_posthoc_rate_total_v45.csv")
csv_kw_final_share            <- file.path(output_dir, "kw_final_share_per_file_v45.csv")
csv_dunn_final_share          <- file.path(output_dir, "dunn_posthoc_final_share_v45.csv")
csv_pairwise_newcombe         <- file.path(output_dir, "pairwise_finalshare_newcombe_wilson_v45.csv")
csv_perm_manova               <- file.path(output_dir, "permutation_manova_rate_share_v45.csv")
csv_sensitivity_finalshare    <- file.path(output_dir, "sensitivity_finalshare_maxclosers_v45.csv")
csv_glmer_freq                <- file.path(output_dir, "glmer_v45_frequency_summary.csv")
csv_emm_glmer_freq_pairs      <- file.path(output_dir, "emm_glmer_v45_freq_pairs.csv")
csv_glmer_finalshare          <- file.path(output_dir, "glmer_v45_finalshare_summary.csv")
csv_emm_glmer_final_pairs     <- file.path(output_dir, "emm_glmer_v45_final_pairs.csv")

# ---------- token helpers (for sensitivity / recompute path) ----------
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
  txt <- readr::read_file(path)
  raw_tokens <- str_split(txt, "\\s+")[[1]]
  raw_tokens <- raw_tokens[nzchar(raw_tokens)]
  surface <- ifelse(raw_tokens %in% c(".", "。", "．", "！", "？", "…", "‥", "⋯"),
                    raw_tokens,
                    str_replace(raw_tokens, "^(.*?)(?:_.+)?$", "\\1"))
  normalize_punct(surface)
}
count_and_so_on_total <- function(tokens) {
  n <- length(tokens); if (n < 3) return(0L)
  sum(tokens[1:(n-2)] == "and" & tokens[2:(n-1)] == "so" & tokens[3:n] == "on")
}
END_TOKENS <- c(".", "!", "?", "…")
CLOSERS    <- c("\"", "”", "’", ")", "]", "】", "〉", "》", "」", "『", "』", "＞")
count_and_so_on_sentence_final <- function(tokens, max_closers = 3L) {
  n <- length(tokens); if (n < 4) return(0L)
  idx <- which(tokens[seq_len(n-2)] == "and" & tokens[seq(2, n-1)] == "so" & tokens[seq(3, n)] == "on")
  finals <- 0L
  for (i in idx) {
    j <- i + 3L; k <- 0L
    while (j <= n && k < max_closers && tokens[j] %in% CLOSERS) { j <- j + 1L; k <- k + 1L }
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

# ---------- load base per-file results (prefer v4.3/4.4 CSV) ----------
base_csv <- file.path(output_dir, "and_so_on_stats_final_vs_nonfinal.csv")
if (file.exists(base_csv)) {
  results <- readr::read_csv(base_csv, show_col_types = FALSE)
} else {
  # Recompute from raw inputs, excluding prior outputs
  all_candidates <- list.files(path = input_dir, pattern = "\\.(txt|csv)$", full.names = TRUE, ignore.case = TRUE)
  exclude_patterns <- c(
    paste0("^", normalizePath(output_dir, winslash = "/")),
    "group_summary_andso\\.csv$", "and_so_on_stats_final_vs_nonfinal\\.csv$",
    "_pairs?|pairwise_", "tests_.*\\.csv$", "glm_.*\\.csv$", "glmmTMB_.*\\.csv$",
    "kw_.*\\.csv$", "dunn_.*\\.csv$", "permutation_.*\\.csv$", "sensitivity_.*\\.csv$",
    "dharma_.*\\.csv$", "_BH\\.csv$", "_v45\\.csv$"
  )
  paths <- all_candidates[!Reduce(`|`, lapply(exclude_patterns, function(pat) str_detect(all_candidates, regex(pat, ignore_case = TRUE))))]
  if (length(paths) == 0) stop("v4.5 could not find inputs, and v4.3/4.4 CSV not present.")
  
  stats_for_file <- function(path) {
    tokens <- read_annotated_tokens(path)
    N <- length(tokens)
    c_total  <- count_and_so_on_total(tokens)
    c_final  <- count_and_so_on_sentence_final(tokens)
    c_nonfin <- c_total - c_final
    denom_trigram <- max(N - 2, 0)
    tibble(
      file = basename(path),
      culture = infer_culture(basename(path)),
      N_tokens = N,
      trigram_positions = denom_trigram,
      count_and_so_on_total    = c_total,
      count_and_so_on_final    = c_final,
      count_and_so_on_nonfinal = c_nonfin,
      rate_total_per_100k  = ifelse(N > 0, c_total  / N * rate_scale, NA_real_),
      rate_final_per_100k  = ifelse(N > 0, c_final  / N * rate_scale, NA_real_),
      rate_nonfin_per_100k = ifelse(N > 0, c_nonfin / N * rate_scale, NA_real_)
    )
  }
  results <- purrr::map_dfr(paths, stats_for_file)
}

# ---------- per-file frame for tests ----------
per_file <- results %>%
  transmute(file, culture,
            rate_total_per_100k,
            rate_final_per_100k,
            rate_nonfin_per_100k,
            trigram_positions,
            final_hits = count_and_so_on_final,
            total_hits = count_and_so_on_final + count_and_so_on_nonfinal) %>%
  mutate(final_share = ifelse(total_hits > 0, final_hits/total_hits, NA_real_))

# ========== 1) Kruskal–Wallis + Dunn ==========
# KW: rate_total_per_100k
d_rate <- per_file %>% filter(!is.na(rate_total_per_100k))
if (nrow(d_rate) >= 2 && n_distinct(d_rate$culture) >= 2) {
  kw_rate <- kruskal.test(rate_total_per_100k ~ culture, data = d_rate)
  write_csv(tibble(outcome="rate_total_per_100k",
                   statistic = unname(kw_rate$statistic),
                   df = n_distinct(d_rate$culture) - 1,
                   p_value = kw_rate$p.value),
            csv_kw_rate_total)
} else emit_note(csv_kw_rate_total, "Not enough non-NA groups for KW on rate_total_per_100k.")

# Dunn post-hoc (Holm)
if (has_pkg("FSA")) {
  tryCatch({
    dout <- FSA::dunnTest(rate_total_per_100k ~ culture, data = per_file, method = "holm")$res %>% as_tibble()
    write_csv(dout, csv_dunn_rate_total)
  }, error = function(e) emit_note(csv_dunn_rate_total, paste("FSA::dunnTest failed:", e$message)))
} else if (has_pkg("dunn.test")) {
  d <- per_file %>% filter(!is.na(rate_total_per_100k))
  if (nrow(d) >= 2 && n_distinct(d$culture) >= 2) {
    dd <- dunn.test::dunn.test(x = d$rate_total_per_100k, g = d$culture, method = "holm", list = TRUE)
    write_csv(tibble(comparison = dd$comparisons, Z = dd$Z, p_adjusted = dd$P.adjusted), csv_dunn_rate_total)
  } else emit_note(csv_dunn_rate_total, "Not enough data for Dunn post-hoc on rate_total_per_100k.")
} else emit_note(csv_dunn_rate_total, "Install FSA or dunn.test for Dunn post-hocs.")

# KW: final_share
d_final <- per_file %>% filter(!is.na(final_share))
if (nrow(d_final) >= 2 && n_distinct(d_final$culture) >= 2) {
  kw_final <- kruskal.test(final_share ~ culture, data = d_final)
  write_csv(tibble(outcome="final_share",
                   statistic = unname(kw_final$statistic),
                   df = n_distinct(d_final$culture) - 1,
                   p_value = kw_final$p.value),
            csv_kw_final_share)
} else emit_note(csv_kw_final_share, "Not enough non-NA groups for KW on final_share.")

# Dunn post-hoc: final_share
if (has_pkg("FSA")) {
  tryCatch({
    dout <- FSA::dunnTest(final_share ~ culture, data = d_final, method = "holm")$res %>% as_tibble()
    write_csv(dout, csv_dunn_final_share)
  }, error = function(e) emit_note(csv_dunn_final_share, paste("FSA::dunnTest failed:", e$message)))
} else if (has_pkg("dunn.test")) {
  d_fin <- per_file %>% filter(!is.na(final_share))
  if (nrow(d_fin) >= 2 && n_distinct(d_fin$culture) >= 2) {
    dd <- dunn.test::dunn.test(x = d_fin$final_share, g = d_fin$culture, method = "holm", list = TRUE)
    write_csv(tibble(comparison = dd$comparisons, Z = dd$Z, p_adjusted = dd$P.adjusted), csv_dunn_final_share)
  } else emit_note(csv_dunn_final_share, "Not enough data for Dunn post-hoc on final_share.")
} else emit_note(csv_dunn_final_share, "Install FSA or dunn.test for Dunn post-hocs.")

# ========== 2) Newcombe–Wilson pairwise CIs (final share by culture) ==========
final_totals <- results %>%
  group_by(culture) %>%
  summarise(final = sum(count_and_so_on_final, na.rm=TRUE),
            nonfinal = sum(count_and_so_on_nonfinal, na.rm=TRUE),
            total = final + nonfinal, .groups="drop") %>%
  filter(total > 0)
pair_rows <- list()
if (nrow(final_totals) >= 2) {
  comb <- combn(final_totals$culture, 2, simplify = FALSE)
  for (pair in comb) {
    g1 <- pair[1]; g2 <- pair[2]
    a <- final_totals %>% filter(culture == g1)
    b <- final_totals %>% filter(culture == g2)
    x1 <- a$final; n1 <- a$total; p1 <- x1/n1
    x2 <- b$final; n2 <- b$total; p2 <- x2/n2
    
    if (has_pkg("PropCIs")) {
      ci <- tryCatch(PropCIs::diffscoreci(x1, n1, x2, n2, conf.level = 0.95), error=function(e) NULL)
      if (!is.null(ci)) {
        pair_rows[[paste(g1,g2,sep="|")]] <- tibble(
          group1=g1, group2=g2, diff = p1 - p2,
          ci_low = ci$conf.int[1], ci_high = ci$conf.int[2],
          method = "Newcombe score (PropCIs::diffscoreci)"
        )
        next
      }
    }
    # fallback approx via Wilson components
    w1 <- binom::binom.wilson(x1, n1)[1, c("lower","upper")]
    w2 <- binom::binom.wilson(x2, n2)[1, c("lower","upper")]
    dl <- (p1 - as.numeric(w1["lower"]))^2 + (as.numeric(w2["upper"]) - p2)^2
    du <- (as.numeric(w1["upper"]) - p1)^2 + (p2 - as.numeric(w2["lower"]))^2
    ci_low <- (p1 - p2) - sqrt(dl)
    ci_high <- (p1 - p2) + sqrt(du)
    pair_rows[[paste(g1,g2,sep="|")]] <- tibble(
      group1=g1, group2=g2, diff = p1 - p2,
      ci_low = ci_low, ci_high = ci_high,
      method = "Newcombe (approx via Wilson components)"
    )
  }
  bind_rows(pair_rows) %>% arrange(group1, group2) %>% write_csv(csv_pairwise_newcombe)
} else emit_note(csv_pairwise_newcombe, "Need >=2 cultures with total > 0 for Newcombe–Wilson CIs.")

# ========== 3) Permutation MANOVA on [rate_per_100k, final_share] ==========
if (has_pkg("vegan")) {
  suppressPackageStartupMessages(library(vegan))
  
  file_level <- results %>%
    group_by(file, culture) %>%
    summarise(tokens = sum(N_tokens, na.rm = TRUE),
              hits = sum(count_and_so_on_total, na.rm = TRUE),
              final_hits = sum(count_and_so_on_final, na.rm = TRUE),
              nonfinal_hits = sum(count_and_so_on_nonfinal, na.rm = TRUE),
              .groups = "drop") %>%
    filter(tokens > 0) %>%
    mutate(rate_per_100k = hits / tokens * rate_scale,
           total_hits = final_hits + nonfinal_hits,
           final_share = ifelse(total_hits > 0, final_hits / total_hits, NA_real_)) %>%
    filter(!is.na(rate_per_100k) & !is.na(final_share)) %>%
    mutate(culture = factor(culture))
  
  if (nrow(file_level) > 0 && n_distinct(file_level$culture) > 1) {
    # Explicit response matrix (fixes eval(YVAR...) bug)
    Y <- as.matrix(file_level[, c("rate_per_100k", "final_share")])
    colnames(Y) <- make.names(colnames(Y))
    
    ad <- tryCatch(
      adonis2(Y ~ culture,
              data = file_level,
              permutations = N_PERM,
              method = "euclidean"),
      error = function(e) e
    )
    
    if (inherits(ad, "error")) {
      emit_note(csv_perm_manova,
                paste("Permutation MANOVA failed:", ad$message))
    } else {
      ad_df <- as.data.frame(ad)
      ad_df$term <- rownames(ad_df)
      write_csv(ad_df, csv_perm_manova)
    }
  } else {
    emit_note(csv_perm_manova,
              "Not enough file-level data for MANOVA (need >0 rows & >1 culture).")
  }
} else {
  emit_note(csv_perm_manova,
            "vegan not installed; skipped permutation MANOVA.")
}


# ========== 4) Sensitivity: max_closers = 0,1,3 ==========
# Uses raw inputs if available; else approximates from existing results (not ideal).
# Try raw path first:
raw_ok <- FALSE
all_candidates <- tryCatch(list.files(path = input_dir, pattern = "\\.(txt)$", full.names = TRUE, ignore.case = TRUE), error=function(e) character())
if (length(all_candidates) > 0) {
  raw_ok <- TRUE
  sens_vals <- c(0L, 1L, 3L)
  sens_rows <- list()
  for (mc in sens_vals) {
    tmp <- purrr::map_dfr(all_candidates, function(p){
      toks <- read_annotated_tokens(p)
      tibble(file = basename(p),
             culture = infer_culture(basename(p)),
             final = count_and_so_on_sentence_final(toks, max_closers = mc),
             total = count_and_so_on_total(toks))
    })
    sens <- tmp %>% group_by(culture) %>% summarise(final=sum(final), total=sum(total), .groups="drop") %>%
      mutate(final_share = ifelse(total>0, final/total, NA_real_))
    sens <- sens %>% mutate(ci = purrr::pmap(list(final, total), function(x, n){
      if (is.na(x) || is.na(n) || n == 0) return(c(NA_real_, NA_real_))
      w <- binom::binom.wilson(x, n); c(w$lower, w$upper)
    }))
    sens$ci_lower <- purrr::map_dbl(sens$ci, 1)
    sens$ci_upper <- purrr::map_dbl(sens$ci, 2)
    sens$max_closers <- mc; sens$ci <- NULL
    sens_rows[[as.character(mc)]] <- sens
  }
  bind_rows(sens_rows) %>% write_csv(csv_sensitivity_finalshare)
}
if (!raw_ok) {
  emit_note(csv_sensitivity_finalshare, "Raw .txt files not found; sensitivity requires tokenization.")
}

# ========== 5) Mixed-effects logistic (GLMER) + emmeans ==========
# (A) Frequency yes/no per file on trigram positions
glmer_data_freq <- results %>%
  transmute(file, culture,
            yes = count_and_so_on_total,
            no  = pmax(trigram_positions - count_and_so_on_total, 0)) %>%
  filter((yes + no) > 0)

model_freq <- tryCatch(
  lme4::glmer(cbind(yes, no) ~ culture + (1 | file), data = glmer_data_freq, family = binomial),
  error = function(e) { message("GLMER freq failed: ", e$message); NULL }
)
if (!is.null(model_freq)) {
  tidy_freq <- broom.mixed::tidy(model_freq, effects = "fixed", conf.int = TRUE, conf.method = "Wald") %>%
    mutate(odds_ratio = exp(estimate), conf.low.OR = exp(conf.low), conf.high.OR = exp(conf.high))
  write_csv(tidy_freq, csv_glmer_freq)
  emm_mfreq <- emmeans(model_freq, ~ culture, type = "response")
  write_csv(as_tibble(pairs(emm_mfreq, adjust = "tukey")), csv_emm_glmer_freq_pairs)
} else emit_note(csv_glmer_freq, "GLMER frequency model failed.")

# (B) Final vs nonfinal
glmer_data_final <- results %>%
  transmute(file, culture,
            final = count_and_so_on_final,
            nonfinal = count_and_so_on_nonfinal) %>%
  filter((final + nonfinal) > 0)

model_final <- tryCatch(
  lme4::glmer(cbind(final, nonfinal) ~ culture + (1 | file), data = glmer_data_final, family = binomial),
  error = function(e) { message("GLMER final failed: ", e$message); NULL }
)
if (!is.null(model_final)) {
  tidy_final <- broom.mixed::tidy(model_final, effects = "fixed", conf.int = TRUE, conf.method = "Wald") %>%
    mutate(odds_ratio = exp(estimate), conf.low.OR = exp(conf.low), conf.high.OR = exp(conf.high))
  write_csv(tidy_final, csv_glmer_finalshare)
  emm_mfinal <- emmeans(model_final, ~ culture, type = "response")
  write_csv(as_tibble(pairs(emm_mfinal, adjust = "tukey")), csv_emm_glmer_final_pairs)
} else emit_note(csv_glmer_finalshare, "GLMER final-share model failed.")

message("v4.5 extras completed. Outputs written to: ", normalizePath(output_dir))
