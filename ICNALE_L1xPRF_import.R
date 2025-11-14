# Load required libraries
library(stringr)
library(stringi)

delete_real_newlines <- function(text) {
  gsub("\n", "", text, fixed = TRUE)   # removes actual newline characters
}


# Define the function
load_combined_texts <- function(folder_path) {
  # List all files in the folder
  files <- list.files(folder_path, full.names = TRUE)
  
  # Split filenames into parts
  parts <- str_split(basename(files), "_")
  
  # Extract 2nd, 5th, and 6th groups
  second_group <- sapply(parts, function(x) x[2])
  fifth_group  <- sapply(parts, function(x) x[5])
  sixth_group  <- sapply(parts, function(x) str_remove(x[6], "\\..*$"))
  
  # Combine 5th and 6th groups
  combined_group <- paste0(fifth_group, sixth_group)
  
  # Build combinations (keys)
  combinations <- paste(second_group, combined_group, sep = "_")
  
  # Get unique combinations
  unique_combinations <- unique(combinations)
  
  # Initialize result list
  result <- list()
  
  # Loop over each unique combination
  for (combo in unique_combinations) {
    
    # Find files that match this combination
    matching_files <- files[combinations == combo]
    
    # Read each file as a full text string
    file_texts <- sapply(matching_files, function(f) {
      readChar(f, file.info(f)$size)
    }, USE.NAMES = FALSE)
    
    # Concatenate all texts directly, no separator
    joined_text <- paste(file_texts, collapse = "")
    joined_text <-delete_real_newlines(joined_text)
    
    # Store in result list
    result[[combo]] <- joined_text
  }
  
  return(result)
}

#usage:
#texts <- load_combined_texts("ICNALE_ALL")

#Export desired files (this uses 'stringi'):
#stri_write_lines(texts$CHN_B12, con = "CHN_B12.txt")
