rm(list = ls()) # Clear the environment

require(xtable)

data_path <- "./data"
file_names <- list.files(data_path)

row_is_valid <- function(row) {
  n <- row["vertices"]
  k2 <- row["degree_2nd_moment"]
  d <- row["mean_length"]
  l1 <- 4-6/n
  r1 <- n-1
  l2 <- n/(8*(n-1)) * k2 + 1/2
  r2 <- r1
  return((l1 <= k2 && k2 <= r1) && (l2 <= d && d <= r2) && (((n - 1)**2 + (n/4 - 1) * n * k2) != 0))
}

get_statistics <- function(language, data) {
  N <- nrow(data)
  mi_n <- mean(data$vertices)
  mi_x <- mean(data$mean_length)
  sigma_n <-sd(data$vertices)
  sigma_x <- sd(data$mean_length)
  data.frame("Language" = language, "N" = N, "mi_n" = mi_n, "sigma_n" = sigma_n, "mi_d_z" = mi_x, "sigma_d_z" = sigma_x)
} 

get_language_from_filename <- function(file_name) {
  unlist(strsplit(file_name, "_"))[1]
}

E_rla_d <- function(n) {
  return((n + 1) / 3)
}

sigma_rla_d <- function(n, k2) {
  return(sqrt(V_rla_D(n, k2)/(n - 1)))
}

V_rla_D <- function (n, k2) {
  return((n + 1) / 45 * ((n - 1)**2 + (n/4 - 1) * n * k2))
}

read_language <- function(name, file_path) {
  print(paste(data_path, "/", file_path))
  language_data <- read.table(paste(data_path, "/", file_path, sep=""))
  colnames(language_data) <- c("vertices", "degree_2nd_moment", "mean_length")
  
  # Remove invalid rows from a language
  valid_rows <- apply(language_data, 1, row_is_valid)
  language_data <- language_data[valid_rows,]
  print(paste("Language:", name, "Invalid rows:", sum(!valid_rows), "Valid rows:", nrow(language_data)))
  
  language_data$mean_length_normalized <- (language_data$mean_length - E_rla_d(language_data$vertices)) / sigma_rla_d(language_data$vertices, language_data$degree_2nd_moment)
  language_data <- language_data[order(language_data$vertices), ]
  
  return(language_data)
}

language_statistics <- data.frame()
languages <- list()
languages_aggr <- list()

# Read the languages generate the general language statistics
for (i in seq_along(file_names)) {
  file_name <- file_names[[i]]
  language = get_language_from_filename(file_name)
  data <- read_language(language, file_name)
  languages[[language]] <- data
  languages_aggr[[language]] <- aggregate(data, list(data$vertices), mean)
  languages_aggr[[language]]$aggr_mean_length_normalized <- (languages_aggr[[language]]$mean_length - E_rla_d(languages_aggr[[language]]$vertices)) / sigma_rla_d(languages_aggr[[language]]$vertices, languages_aggr[[language]]$degree_2nd_moment)
  
  new_row <- get_statistics(language, data)
  language_statistics <- rbind(language_statistics, new_row)
}

print_basic <- function(language) {
  plot(languages[[language]]$vertices, languages[[language]]$mean_length, xlab = "vertices", ylab = "mean dependency length")
  plot(log(languages[[language]]$vertices), log(languages[[language]]$mean_length), xlab = "log(vertices)", ylab = "log(mean dependency length)")
}

print_aggregated <- function(language) {
  plot(languages_aggr[[language]]$vertices, languages_aggr[[language]]$mean_length, xlab = "vertices", ylab = "mean mean dependency length")
  plot(log(languages_aggr[[language]]$vertices), log(languages_aggr[[language]]$mean_length), xlab = "log(vertices)", ylab = "log(mean mean dependency length)")
}

print_d_basic_estimation <- function(language) {
  plot(log(languages[[language]]$vertices), log(languages[[language]]$mean_length), xlab = "vertices", ylab = "mean dependency length")
  lines(log(languages_aggr[[language]]$vertices),log(languages_aggr[[language]]$mean_length), col = "green")
  lines(log(languages_aggr[[language]]$vertices),log((languages_aggr[[language]]$vertices+1)/3), col = "red") 
}

print_d_normalized_basic_estimation <- function(language) {
  plot(log(languages[[language]]$vertices), log(languages[[language]]$mean_length_normalized), xlab = "vertices", ylab = "mean dependency length")
  lines(log(languages_aggr[[language]]$vertices),log(languages_aggr[[language]]$mean_length_normalized), col = "green")
  lines(log(languages_aggr[[language]]$vertices),log((languages_aggr[[language]]$vertices+1)/3), col = "red") 
}

print_k2_basic_estimation <- function(language) {
  plot(languages[[language]]$vertices, languages[[language]]$degree_2nd_moment, xlab = "vertices", ylab = "degree 2nd moment")
  lines(languages_aggr[[language]]$vertices,languages_aggr[[language]]$degree_2nd_moment, col = "green")
  lines(languages[[language]]$vertices, (1 - 1/languages[[language]]$vertices)*(5 - 6/languages[[language]]$vertices), col = "red")
  lines(languages[[language]]$vertices,4-6/languages[[language]]$vertices, col = "blue")
  lines(languages[[language]]$vertices,languages[[language]]$vertices-1, col = "blue")
}

init_vals <- list(
  Catalan=list(
    "2"=list(
      a=-1,
      b=0.5
    ),
    "2+"=list(
      a=-0.77,
      b=0.11,
      d=-2
    ),
    "3"=list(
      a=-0.5,
      c=0.1
    ),
    "4"=list(
      a=-1
    )
  ),
  English=list(
    "3"=list(
      a=-1,
      c=-0.1
    )
  ),
  Arabic=list(
    "2+"=list(
      a=2,
      b=-0.5,
      d=-2
    ),
    "3"=list(
      a=-1,
      c=-0.1
    )
  ),
  Basque=list(
    "2+"=list(
      a=2,
      b=-0.5,
      d=-2
    )
  ),
  Chinese=list(
    "2+"=list(
      a=2,
      b=-0.5,
      d=-2
    )
  ),
  Czech=list(
    "2+"=list(
      a=2,
      b=-0.5,
      d=-2
    ),
    "3"=list(
      a=-1,
      c=0.1
    )
  )
)

get_init_val <- function(language, method) {
  ret <- list()
  for (param in names(init_vals$Catalan[[method]])) {
    if(!(language %in% names(init_vals)) || !(method %in% names(init_vals[[language]])) || !(param %in% names(init_vals[[language]][[method]]))) {
      ret[[param]] <- init_vals$Catalan[[method]][[param]]
    } else {
      ret[[param]] <- init_vals[[language]][[method]][[param]]
    }
  }
  return(ret)
}

est_language <- function(language) {
  res <- list("0"=list(), "2"=list(), "2+"=list(), "3"=list(), "4"=list())
  n <- length(languages_aggr[[language]]$vertices)
  
  res$"0"$RSS <- sum((languages_aggr[[language]]$aggr_mean_length_normalized-(n+1)/3)^2)
  res$"0"$AIC <- n*log(2*pi) + n*log(res$"0"$RSS/n) + n + 2*(0 + 1)
  res$"0"$s <- sqrt(res$"0"$RSS/(n - 0))
  
  nonlinear_model <- nls(aggr_mean_length_normalized~a*vertices^b,data=languages_aggr[[language]], start = get_init_val(language, "2"), trace = TRUE)
  res$"2"$a <- coef(nonlinear_model)["a"]
  res$"2"$b <- coef(nonlinear_model)["b"]
  res$"2"$RSS <- deviance(nonlinear_model)
  res$"2"$AIC <- AIC(nonlinear_model)
  res$"2"$s <- sqrt(deviance(nonlinear_model)/df.residual(nonlinear_model))
  
  nonlinear_model <- nls(aggr_mean_length_normalized~a*vertices^b+d,data=languages_aggr[[language]], start = get_init_val(language, "2+"), trace = TRUE)
  res$"2+"$a <- coef(nonlinear_model)["a"]
  res$"2+"$b <- coef(nonlinear_model)["b"]
  res$"2+"$d <- coef(nonlinear_model)["d"]
  res$"2+"$RSS <- deviance(nonlinear_model)
  res$"2+"$AIC <- AIC(nonlinear_model)
  res$"2+"$s <- sqrt(deviance(nonlinear_model)/df.residual(nonlinear_model))
  
  res = tryCatch({
    nonlinear_model = nls(aggr_mean_length_normalized~a*exp(vertices*c),data=languages_aggr[[language]], start = get_init_val(language, "3"), trace = TRUE)
    res$"3"$a <- coef(nonlinear_model)["a"]
    res$"3"$c <- coef(nonlinear_model)["c"]
    res$"3"$RSS <- deviance(nonlinear_model)
    res$"3"$AIC <- AIC(nonlinear_model)
    res$"3"$s <- sqrt(deviance(nonlinear_model)/df.residual(nonlinear_model))
    res
  }, error=function(cond) {
    res$"3"$a <- NaN
    res$"3"$c <- NaN
    res$"3"$RSS <- NaN
    res$"3"$AIC <- NaN
    res$"3"$s <- NaN
    message("ERROR FITTING MODEL 3")
    message("Here's the original error message:")
    message(cond)
    return(res)
  })
  
  nonlinear_model = nls(aggr_mean_length_normalized~a*log(vertices),data=languages_aggr[[language]], start = get_init_val(language, "4"), trace = TRUE)
  res$"4"$a <- coef(nonlinear_model)["a"]
  res$"4"$RSS <- deviance(nonlinear_model)
  res$"4"$AIC <- AIC(nonlinear_model)
  res$"4"$s <- sqrt(deviance(nonlinear_model)/df.residual(nonlinear_model))
  
  print("min")
  min_aic = min(res$"0"$AIC, res$"2"$AIC, res$"2+"$AIC, res$"3"$AIC, res$"4"$AIC, na.rm=TRUE)
  print(min_aic)
  res$"0"$AIC_d = res$"0"$AIC - min_aic
  res$"2"$AIC_d = res$"2"$AIC - min_aic
  res$"2+"$AIC_d = res$"2+"$AIC - min_aic
  res$"3"$AIC_d = res$"3"$AIC - min_aic
  res$"4"$AIC_d = res$"4"$AIC - min_aic
  
  print(c(res$"0"$AIC_d, res$"2"$AIC_d, res$"2+"$AIC_d, res$"3"$AIC_d, res$"4"$AIC_d))
  
  return(res)
}

results <- list()
for (l in names(languages)) {
  print(l)
  results[[l]] <- est_language(l)
}

s = data.frame(matrix(ncol=6, nrow=0))
aics = data.frame(matrix(ncol=6, nrow=0))
aics_d  = data.frame(matrix(ncol=6, nrow=0))
names(aics) = c("language", "0_s", "2_s", "2+_s", "3_s", "4_s")
names(aics) = c("language", "0_aic", "2_aic", "2+_aic", "3_aic", "4_aic")
names(aics) = c("language", "0_aic_d", "2_aic_d", "2+_aic_d", "3_aic_d", "4_aic_d")
for (language in names(results)) {
  r <- results[[language]]
  s[nrow(s) + 1,] = c(language, r$`0`$s, r$`2`$s, r$`2+`$s, r$`3`$s, r$`4`$s)
  aics[nrow(aics) + 1,] = c(language, r$`0`$AIC, r$`2`$AIC, r$`2+`$AIC, r$`3`$AIC, r$`4`$AIC)
  aics_d[nrow(aics_d) + 1,] = c(language, r$`0`$AIC_d, r$`2`$AIC_d, r$`2+`$AIC_d, r$`3`$AIC_d, r$`4`$AIC_d)
}

print_est <- function(language) {
  plot(languages_aggr[[language]]$vertices, languages_aggr[[language]]$aggr_mean_length_normalized, xlab = "vertices", ylab = "mean dependency length")
  lines(languages_aggr[[language]]$vertices,results[[language]]$`2`$a*languages_aggr[[language]]$vertices^results[[language]]$`2`$b, col = "blue")
  lines(languages_aggr[[language]]$vertices,results[[language]]$`2+`$a*languages_aggr[[language]]$vertices^results[[language]]$`2+`$b+results[[language]]$`2+`$d, col = "red")
  lines(languages_aggr[[language]]$vertices,results[[language]]$`3`$a*languages_aggr[[language]]$vertices^results[[language]]$`3`$c, col = "green")
  lines(languages_aggr[[language]]$vertices,results[[language]]$`4`$a*log(languages_aggr[[language]]$vertices), col = "brown")
}

print_est_loglog <- function(language) {
  plot(log(languages_aggr[[language]]$vertices), -log(-languages_aggr[[language]]$aggr_mean_length_normalized), xlab = "vertices", ylab = "mean dependency length")
  lines(log(languages_aggr[[language]]$vertices), -log(-results[[language]]$`2`$a*languages_aggr[[language]]$vertices^results[[language]]$`2`$b), col = "blue")
  lines(log(languages_aggr[[language]]$vertices), -log(-results[[language]]$`2+`$a*languages_aggr[[language]]$vertices^results[[language]]$`2+`$b+results[[language]]$`2+`$d), col = "red")
  lines(log(languages_aggr[[language]]$vertices), -log(-results[[language]]$`3`$a*languages_aggr[[language]]$vertices^results[[language]]$`3`$c), col = "green")
  lines(log(languages_aggr[[language]]$vertices), -log(-results[[language]]$`4`$a*log(languages_aggr[[language]]$vertices)), col = "brown")
}

stopifnot(TRUE)

print_d_basic_estimation('Catalan')
print_d_basic_estimation('English')
print_aggregated('Catalan')

plot(languages_aggr$Catalan$vertices, languages_aggr$Catalan$mean_length, xlab = "vertices", ylab = "mean dependency length")
plot(languages_aggr$Catalan$vertices, E_rla_d(languages_aggr$Catalan$vertices), xlab = "vertices", ylab = "mean dependency length")