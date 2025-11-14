library(stringr)
library(dplyr)


# MODULE 3: 
# Words you want to count (case-insensitive)
target_words <- c("whisky", "horse", "gold", "cup")

#MODULE 1:
#name of text file you will import
FileName <- "JJU.txt"
#this line reads the text file into R
text <- readLines(FileName, warn = FALSE, encoding = "UTF-8")

#MODULE 2:
#this line turns the entire text into a single line
text <- paste(text, collapse = " ")
#this line makes everything lower case
text <- tolower(text)
#this line removes punctuation
text <- str_replace_all(text, "[[:punct:]]", " ")
#this line normalizes spacing
text <- str_replace_all(text, "\\s+", " ")
#this line splits the cleaned up text into a list of words
words <- unlist(str_split(text, "\\s+"))

#MODULE 4:
# Create a frequency table
word_freq <- data.frame(table(words))

# Keep only the target words
result <- word_freq %>%
  filter(words %in% target_words) %>%
  arrange(desc(Freq))

print(result)
