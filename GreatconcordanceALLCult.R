# ============================================================
# TXT-only FULL-SENTENCE Concordance for "and so on"
# - Robust string-level pre-clean (handles "._.I_PRP", ". . .", "--_ :")
# - ONE row per sentence that contains the phrase
# - Lists ALL occurrences per sentence (n_hits + positions)
# - Creates fresh cleaned/ each run
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(glue)
  library(readr)
})

# ----------------- Config -----------------
write_txt <- FALSE        # set TRUE if you want a TXT view
kwic_sep  <- "  ||  "

# ----------------- String-level pre-clean -----------------
preclean_text <- function(x) {
  s <- paste(x, collapse = " ")
  s <- gsub("\\s+", " ", s)  # normalize spaces early
  
  # 1) Fused punctuation with underscores (e.g., "._.", ",_,", "?_?", "!_!", ";_;", ":_:")
  #    Replace the whole fused sequence with the terminal mark (2nd capture), padded with spaces.
  s <- stringr::str_replace_all(s, stringr::regex("([\\.,;:\\?!])_+([\\.,;:\\?!])"), " \\2 ")
  
  # 2) Hyphen/underscore chains ending in a real mark, e.g. "--_ :" or "--_:"
  s <- stringr::str_replace_all(s, stringr::regex("[-]+_+\\s*([\\.,;:\\?!:])"), " \\1 ")
  
  # 3) Normalize Japanese/alt enders and ellipses
  s <- gsub("。|．|…", ".", s)
  s <- gsub("！", "!", s)
  s <- gsub("？", "?", s)
  s <- gsub("\\.\\.+", ".", s)  # "...." -> "."
  s <- gsub("(\\s*\\.\\s*){2,}", " . ", s, perl = TRUE)  # ". . ." -> " . "
  
  # 4) Penn Treebank markers -> real chars
  s <- gsub("-LRB-", "(", s, fixed = TRUE)
  s <- gsub("-RRB-", ")", s, fixed = TRUE)
  s <- gsub("-LSB-", "[", s, fixed = TRUE)
  s <- gsub("-RSB-", "]", s, fixed = TRUE)
  s <- gsub("-LCB-", "{", s, fixed = TRUE)
  s <- gsub("-RCB-", "}", s, fixed = TRUE)
  s <- gsub("``|''", "\"", s)
  
  # 5) Break punctuation glued to words (on both sides), AFTER fused/alt normalization
  #    e.g., ".I_PRP" -> ". I_PRP", "word,Next" -> "word , Next"
  s <- gsub("(?<=\\S)([\\.,;:\\?!])", " \\1", s, perl = TRUE)
  s <- gsub("([\\.,;:\\?!])(?=\\S)", "\\1 ", s, perl = TRUE)
  
  # 6) Any remaining all-punct/underscore tokens -> keep the first of . , ; : ? ! (else drop)
  s <- stringr::str_replace_all(
    s,
    stringr::regex("\\b[[:punct:]_]+\\b"),
    function(m) {
      p  <- gsub("_", "", m)                          # remove underscores
      m2 <- stringr::str_match(p, "([\\.,;:\\?!])")   # first terminal mark
      ifelse(!is.na(m2[,2]), paste0(" ", m2[,2], " "), " ")
    }
  )
  
  # 7) Final space normalization
  s <- gsub("\\s+", " ", s)
  trimws(s)
}

# ----------------- Token helpers -----------------

# Drop ONLY the final POS tag when present: "word_POS" -> "word"
drop_final_pos <- function(tok) {
  if (str_detect(tok, "_")) sub("^(.*)_[^_]+$", "\\1", tok) else tok
}

# Preserve case for rendering
norm_preserve <- function(tok) str_squish(tok)

# Lower for matching & segmentation; map enders to "."
norm_lower <- function(tok) {
  t <- str_to_lower(str_squish(tok))
  if (t %in% c(".", "!", "?")) return(".")
  t
}

# Sentence ids from lower tokens ('.' ends a sentence)
make_sentence_ids <- function(tokens_lower) {
  ids <- integer(length(tokens_lower)); sid <- 1L
  for (i in seq_along(tokens_lower)) {
    ids[i] <- sid
    if (tokens_lower[i] == ".") sid <- sid + 1L
  }
  ids
}

collapse_sentence <- function(tokens_vec) {
  s <- paste(tokens_vec, collapse = " ")
  s <- gsub(" ([\\.,;:\\?!])", "\\1", s, perl = TRUE)
  str_squish(s)
}

# Highlight all occurrences of the 3-word span starting at each rel_start
highlight_all_occurrences <- function(sentence_tokens, rel_starts) {
  if (length(sentence_tokens) == 0) return("")
  if (length(rel_starts) == 0 || all(is.na(rel_starts))) return(collapse_sentence(sentence_tokens))
  out <- sentence_tokens
  rel_starts <- sort(unique(rel_starts[!is.na(rel_starts)]))
  for (rs in rel_starts) {
    re <- rs + 2L
    if (rs >= 1L && re <= length(out)) {
      out[rs] <- paste0("[", out[rs])
      out[re] <- paste0(out[re], "]")
    }
  }
  collapse_sentence(out)
}

safe_slice <- function(v, i, j) {
  i <- max(1L, i); j <- min(length(v), j)
  if (i > j) character(0) else v[i:j]
}

# ----------------- Source files (TXT only; skip *_clean.txt) -----------------

all_txt <- list.files(pattern = "\\.txt$", ignore.case = TRUE)
src_txt <- all_txt[!grepl("_clean\\.txt$", all_txt, ignore.case = TRUE)]
if (length(src_txt) == 0) stop("No source .txt files found (excluding *_clean.txt).")

# Start fresh each run
if (dir.exists("cleaned")) unlink("cleaned", recursive = TRUE, force = TRUE)
dir.create("cleaned", showWarnings = FALSE)

# ----------------- Pre-clean to cleaned/*.txt -----------------

for (f in src_txt) {
  raw_text  <- readLines(f, warn = FALSE, encoding = "UTF-8")
  s         <- preclean_text(raw_text)
  
  # split to tokens, then drop final POS tags
  toks      <- str_split(s, "\\s+")[[1]]
  toks      <- vapply(toks, drop_final_pos, character(1))
  # drop pure numbers and standalone culture tags (if present)
  toks      <- toks[nzchar(toks)]
  toks      <- toks[!(str_detect(toks, "^\\d+$") | toks %in% c("JPN","CHN","KOR"))]
  
  out_path  <- file.path("cleaned", paste0(tools::file_path_sans_ext(basename(f)), "_clean.txt"))
  
  if (length(toks) > 0) {
    chunks <- split(toks, ceiling(seq_along(toks)/500))
    lines  <- vapply(chunks, function(x) paste(x, collapse = " "), character(1))
    writeLines(lines, out_path, useBytes = TRUE)
  } else {
    writeLines(character(0), out_path, useBytes = TRUE)
  }
}

# ----------------- Concordance on cleaned/*.txt -----------------

cleaned_files <- list.files("cleaned", pattern = "_clean\\.txt$", full.names = TRUE)
if (length(cleaned_files) == 0) stop("No cleaned files found in 'cleaned/'.")

rows <- list()

for (cf in cleaned_files) {
  base <- sub("_clean\\.txt$", "", basename(cf))
  # Reapply preclean to be safe, then tokenize
  s     <- preclean_text(readLines(cf, warn = FALSE, encoding = "UTF-8"))
  toks0 <- str_split(s, "\\s+")[[1]]
  if (length(toks0) == 0) next
  
  toks_pres <- vapply(toks0, norm_preserve, character(1))
  toks_low  <- vapply(toks0, norm_lower,   character(1))
  
  n <- length(toks_low); if (n < 3) next
  
  sid <- make_sentence_ids(toks_low)
  spans <- tibble(idx = seq_along(toks_low), sid = sid) %>%
    group_by(sid) %>% summarise(start = min(idx), end = max(idx), .groups = "drop")
  
  # all hits (exact three-word lower-case match)
  hits <- which(
    toks_low[seq_len(n - 2L)] == "and" &
      toks_low[seq_len(n - 2L) + 1L] == "so" &
      toks_low[seq_len(n - 2L) + 2L] == "on"
  )
  if (length(hits) == 0) next
  
  # group hits by sentence -> ONE row per sentence
  hit_df  <- tibble(hit = hits, sid = sid[hits]) %>% arrange(sid, hit)
  by_sent <- split(hit_df$hit, hit_df$sid)
  
  for (this_sid in names(by_sent)) {
    this_sid_num <- as.integer(this_sid)
    hit_idxs     <- by_sent[[this_sid]]
    span         <- spans %>% filter(sid == this_sid_num)
    if (nrow(span) == 0) next
    s_start <- span$start[[1]]; s_end <- span$end[[1]]
    
    sent_tokens <- safe_slice(toks_pres, s_start, s_end)
    rel_starts  <- hit_idxs - s_start + 1L
    
    is_final_any <- any((hit_idxs + 3L) <= n & toks_low[hit_idxs + 3L] == ".")
    
    rows[[length(rows) + 1L]] <- tibble(
      file                      = base,
      sentence_id               = this_sid_num,
      sentence_start_i          = s_start,
      sentence_end_i            = s_end,
      n_hits                    = length(rel_starts),
      hit_positions_in_sentence = paste(rel_starts, collapse = ","),
      is_final_any              = is_final_any,
      sentence_text             = collapse_sentence(sent_tokens),
      sentence_highlighted      = highlight_all_occurrences(sent_tokens, rel_starts)
    )
  }
}

concordance <- dplyr::bind_rows(rows) %>%
  distinct(file, sentence_id, sentence_text, .keep_all = TRUE) %>%
  arrange(file, sentence_id)

# ----------------- Outputs -----------------

readr::write_csv(concordance, "and_so_on_sentences.csv")

if (write_txt && nrow(concordance) > 0) {
  lines <- glue_data(
    concordance,
    "{sentence_highlighted}{kwic_sep}{file}{kwic_sep}sid={sentence_id}{kwic_sep}hits={n_hits}{kwic_sep}final_any={ifelse(is_final_any,'Y','N')}"
  )
  writeLines(as.character(lines), "and_so_on_sentences.txt", useBytes = TRUE)
}

message("Done! Cleaned files in 'cleaned/'. Wrote: and_so_on_sentences.csv",
        if (write_txt && nrow(concordance) > 0) " and and_so_on_sentences.txt" else "",
        "  (rows: ", nrow(concordance), ").")
