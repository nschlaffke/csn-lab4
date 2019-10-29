setwd("~/Documentos/erasmus_2019/CSN/csn-lab4/data")
file_names <- list.files(".")
language_statistics <- data.frame(Language=character(),
                                  N=integer(), 
                                  mi_n=numeric(), 
                                  sigma_n=numeric(),
                                  mi_x=numeric(),
                                  sigma_x=numeric()) 

get_statistics <- function(file_path, data) {
  N <- nrow(data)
  mi_n <- mean(data$vertices)
  mi_x <- mean(data$mean_length)
  sigma_n <-sd(data$vertices)
  sigma_x <- sd(data$mean_length)
  data.frame("Language" = file_path, "N" = N, "mi_n" = mi_n, "sigma_n" = sigma_n, "mi_x" = mi_x, "sigma_x" = sigma_x)
} 

read_language <- function(file_path) {
  language_data <- read.table(file_path)
  colnames(language_data) <- c("vertices", "degree_2nd_moment", "mean_length")
  language_data[order(language_data$vertices), ]
}

for (file_name in file_names) {
  data <- read_language(file_name)
  new_row <- get_statistics(file_name, data)
  language_statistics <- rbind(language_statistics, new_row)
}
head(language_statistics)
