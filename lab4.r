rm(list = ls()) # Clear the environment

require(xtable)
setwd("~/Documents/erasmus_2019/CSN/csn-lab4/data")
file_names <- list.files(".")

language_statistics <- data.frame()
languages <- list()

get_row_if_invalid <- function(row) {
    n <- row["vertices"]
    k2 <- row["degree_2nd_moment"]
    d <- row["mean_length"]
    l1 <- 4-6/n
    r1 <- n-1
    l2 <- n/(8*(n-1)) * k2 + 1/2
    r2 <- r1
    if(!(l1 <= k2 && k2 <= r1) || !(l2 <= d && d <= r2)) {
      return(row)
    }
}

get_statistics <- function(file_path, data) {
  N <- nrow(data)
  mi_n <- mean(data$vertices)
  mi_x <- mean(data$mean_length)
  sigma_n <-sd(data$vertices)
  sigma_x <- sd(data$mean_length)
  language <- get_language_from_filename(file_path)
  data.frame("Language" = language, "N" = N, "mi_n" = mi_n, "sigma_n" = sigma_n, "mi_x" = mi_x, "sigma_x" = sigma_x)
} 

get_language_from_filename <- function(file_name) {
  unlist(strsplit(file_name, "_"))[1]
}

read_language <- function(file_path) {
  language_data <- read.table(file_path)
  colnames(language_data) <- c("vertices", "degree_2nd_moment", "mean_length")
  language_data <- language_data[order(language_data$vertices), ]
  
  # Remove invalid rows from a language
  rows_to_remove <- apply(language_data, 1, get_row_if_invalid)
  rows_to_remove <- rows_to_remove[!sapply(rows_to_remove, is.null)]
  if (length(rows_to_remove) != 0) {
    language_data <- language_data[!rows_to_remove %in% language_data,]
  }
  return(language_data)
}

# Read the languages generate the general language statistics
for (i in seq_along(file_names)) {
  file_name <- file_names[[i]]
  data <- read_language(file_name)
  languages[[i]] <- data
  names(languages)[i] <- get_language_from_filename(file_name)
  
  new_row <- get_statistics(file_name, data)
  language_statistics <- rbind(language_statistics, new_row)
}

print(language_statistics)
str(languages)