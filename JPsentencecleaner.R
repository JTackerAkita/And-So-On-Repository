# ==============================================
# Remove the special "._." punctuation from JPN*.txt files
# - Reads all JPN*.txt in input_dir
# - Replaces every occurrence of "._." with a single space
# - Tidies extra spaces while preserving line breaks
# - Writes cleaned .txt files to a new folder
# ==============================================

suppressPackageStartupMessages({
  library(stringr)
  library(readr)
  library(glue)
})

remove_dot_underscore_punct <- function(input_dir = ".",
                                        output_dir = "JPN_no_punct",
                                        pattern = "^JPN.*\\.txt$") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  files <- list.files(input_dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) {
    stop(glue("No files matching '{pattern}' found in: {normalizePath(input_dir)}"))
  }
  
  for (fp in files) {
    txt <- read_file(fp)
    
    # 1) Replace every '._.' with a single space (prevents word concatenation)
    cleaned <- str_replace_all(txt, fixed("._."), " ")
    
    # 2) Tidy spaces without collapsing line breaks
    #    - collapse multiple spaces
    cleaned <- str_replace_all(cleaned, " {2,}", " ")
    #    - trim spaces at line ends and beginnings
    cleaned <- str_replace_all(cleaned, "[ \\t]+(?=\\n)", "")   # trailing spaces before newline
    cleaned <- str_replace_all(cleaned, "(?<=\\n)[ \\t]+", "")  # leading spaces after newline
    cleaned <- str_replace_all(cleaned, "^[ \\t]+|[ \\t]+$", "")  # start/end of whole text
    
    out_path <- file.path(output_dir, basename(fp))
    write_file(cleaned, out_path)
    message(glue("Wrote: {out_path}"))
  }
  
  invisible(normalizePath(output_dir))
}

# --- Example usage ---
# setwd("path/to/your/folder/with/JPN_txt")
# remove_dot_underscore_punct(input_dir = ".", output_dir = "JPN_clean_txt")
