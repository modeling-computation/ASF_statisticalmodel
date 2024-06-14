### For rank test 

scan_pois_rank <- function(est, dat, total_iter, iter) {
  if (length(est) != length(dat)) {
    stop("Length of est and dat must be the same")
  }
  
  dat.rank <- rank(-dat, ties.method = "average", na.last = "keep")
  est.count.matrix <- matrix(NA, nrow = length(est), ncol = iter)
  est.rank.matrix <- matrix(NA, nrow = length(est), ncol = iter)
  p_list <- numeric(length(est))
  
  # Generate all Poisson random numbers at once
  all_rpois <- matrix(rpois(length(est) * iter * total_iter, lambda = rep(est, each = iter * total_iter)), nrow = length(est))
  
  for (i in 1:total_iter) {
    # Fill est.count.matrix and est.rank.matrix
    for (j in 1:iter) {
      idx <- ((i - 1) * iter + j)
      est.count.matrix[, j] <- all_rpois[, idx]
      est.rank.matrix[, j] <- rank(-est.count.matrix[, j], ties.method = "average", na.last = "keep")
    }
    # Perform the U-test
    for (k in 1:length(est)) {
      wil.test <- suppressWarnings(wilcox.test(est.rank.matrix[k,], dat.rank[k], alternative = "greater"))
      pvalue <- wil.test$p.value
      p <- ifelse(pvalue < 0.05, 0, 1)
      p_list[k] <- p_list[k] + p
    }
  }
  
  return(p_list / total_iter)
}



scan_negbin_rank <- function(est, dat, total_iter, iter) {
  if (length(est) != length(dat)) {
    stop("Length of est and dat must be the same")
  }
  
  dat.rank <- rank(-dat, ties.method = "average", na.last = "keep")
  est.count.matrix <- matrix(NA, nrow = length(est), ncol = iter)
  est.rank.matrix <- matrix(NA, nrow = length(est), ncol = iter)
  p_list <- numeric(length(est))
  
  # Generate all Negative Binomial random numbers at once
  all_rnbinom <- matrix(rnbinom(length(est) * iter * total_iter, mu = rep(est, each = iter * total_iter), size = 1), nrow = length(est))
  
  for (i in 1:total_iter) {
    # Fill est.count.matrix and est.rank.matrix
    for (j in 1:iter) {
      idx <- ((i - 1) * iter + j)
      est.count.matrix[, j] <- all_rnbinom[, idx]
      est.rank.matrix[, j] <- rank(-est.count.matrix[, j], ties.method = "average", na.last = "keep")
    }
    # Perform the U-test
    for (k in 1:length(est)) {
      wil.test <- suppressWarnings(wilcox.test(est.rank.matrix[k,], dat.rank[k], alternative = "greater"))
      pvalue <- wil.test$p.value
      p <- ifelse(pvalue < 0.05, 0, 1)
      p_list[k] <- p_list[k] + p
    }
  }
  
  return(p_list / total_iter)
}

library(VGAM)
scan_hurdle_rank <- function(est, dat, total_iter, iter) {
  if (length(est) != length(dat)) {
    stop("Length of est and dat must be the same")
  }
  
  dat.rank <- rank(-dat, ties.method = "average", na.last = "keep")
  est.count.matrix <- matrix(NA, nrow = length(est), ncol = iter)
  est.rank.matrix <- matrix(NA, nrow = length(est), ncol = iter)
  p_list <- numeric(length(est))
  
  # Generate all Zero-Inflated Poisson random numbers at once
  all_rzip <- matrix(rzipois(length(est) * iter * total_iter, lambda = rep(est, each = iter * total_iter), pstr0 = 0.9), nrow = length(est))
  
  for (i in 1:total_iter) {
    # Fill est.count.matrix and est.rank.matrix
    for (j in 1:iter) {
      idx <- ((i - 1) * iter + j)
      est.count.matrix[, j] <- all_rzip[, idx]
      est.rank.matrix[, j] <- rank(-est.count.matrix[, j], ties.method = "average", na.last = "keep")
    }
    # Perform the U-test
    for (k in 1:length(est)) {
      wil.test <- suppressWarnings(wilcox.test(est.rank.matrix[k,], dat.rank[k], alternative = "greater"))
      pvalue <- wil.test$p.value
      p <- ifelse(pvalue < 0.05, 0, 1)
      p_list[k] <- p_list[k] + p
    }
  }
  
  return(p_list / total_iter)
}

# data : 추정 기간에 대한 데이터, pred : 특정 기간에 대한 추정 폐사체 수
risk_score <- function(data, pred, p, infected) {
	new_data <- data %>%
		mutate(pred = pred, p_value = p) %>%
		filter(p_value < 0.05) %>%
		dplyr::select(SIG_CD, infected, pred, p_value) %>%
		filter(infected > pred)
	new_df <- data.frame(table(new_data$SIG_CD) %>% sort(decreasing = TRUE))
	colnames(new_df) <- c('SIG_CD', 'risk_score')
	new_df$risk_score <- new_df$risk_score 
	return(new_df)
}














