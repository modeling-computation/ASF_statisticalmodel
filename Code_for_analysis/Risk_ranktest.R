library(dplyr)
library(spdep)
library(scanstatistics)
library(ggplot2)
library(tidyr)
library(readxl)
library(sp)
library(magrittr)
library(terra)
library(sf)
library(MASS)
library(pscl)
setwd('/Users/gogyeongtae/Library/CloudStorage/Dropbox/2024/2024 ASF/ASF_figure_code/Analysis_code/')
source('./test_functions.r')

# Set path
apple_path <- "/Users/gogyeongtae/Library/CloudStorage/Dropbox/2022/2022 ASF/1차정리_202307/Data/시군구별 데이터/통합_shp/"
apple_path <- "D:/Dropbox/2022/2022 ASF/1차정리_202307/Data/시군구별 데이터/통합_shp/"
gpd_file_path = '../Data/'
surv_file <- st_read(paste0(gpd_file_path, 'korea_with_surveillance.shp'))
surv_file$SIG_CD <- as.numeric(surv_file$SIG_CD)

# Data
shp <- st_read(paste0(apple_path, 'shp통합4.shp'))
coord <- read_excel(paste0(apple_path, '시군구위도경도.xlsx'))
used_data <- read.csv('../Data/weekly_merge.csv', header = TRUE)
surv_data <- read.csv('../Data/add_surveilance.csv', header = TRUE)

used_data

# sort and add infected
used_data <- used_data %>% arrange(SIG_CD, time)
used_data['infected'] <- surv_data$Infected

# Convert some types
used_data$time <- as.numeric(used_data$time)
used_data$distance <- exp(-1 * used_data$distance)
used_data$sin <- sin(2*pi*used_data$time/52)
use_colnames <- c('time', 'time_num', 'sin', 'distance', 'forest', 'NUMPOINTS', 'infected', 'season', 'SIG_CD', 'x', 'y')
data <- used_data[, use_colnames]

# add distance^2
data$distance2 <- data$distance^2

# Poisson GLM for first period
mod_1 = glm(infected ~ offset(log(forest))  + distance + distance2 +  time_num + season,
          data = data %>% filter(time < 20220901),
          family = poisson(link = 'log'))

mod_2 = glm.nb(infected ~ offset(log(forest))  + distance + distance2 +  time_num + season,
              data = data %>% filter(time < 20220901))

mod_3 = hurdle(infected ~ offset(log(forest))  + distance + distance2 +  time_num + season,
              data = data %>% filter(time < 20220901),
              dist = 'poisson', zero.dist = 'binomial')

pred <- round(predict(mod_1, newdata = data %>% filter(time >= 20220901 & time <= 20221231), type = 'response'))
true <- data$infected[data$time >= 20220901 & data$time <= 20221231]

pred2 <- round(predict(mod_2, newdata = data %>% filter(time >= 20220901 & time <= 20221231), type = 'response'))
true2 <- data$infected[data$time >= 20220901 & data$time <= 20221231]

pred3 <- round(predict(mod_3, newdata = data %>% filter(time >= 20220901 & time <= 20221231), type = 'response'))
true3 <- data$infected[data$time >= 20220901 & data$time <= 20221231]

# Poisson rank test
start_time <- Sys.time()
pp <- scan_pois_rank(est = pred, dat = true, total_iter = 1, iter = 10000)
end_time <- Sys.time()
print(end_time - start_time)

# Poisson risk score
risk_score(data %>% filter(time >= 20220901 & time <= 20221231),
	pred, pp)

# Negbin rank test
start_time <- Sys.time()
pp <- scan_negbin_rank(est = pred2, dat = true2, total_iter = 1, iter = 10000)
end_time <- Sys.time()
print(end_time - start_time)

# Negbin risk score
risk_score(data %>% filter(time >= 20220901 & time <= 20221231),
	pred2, pp)

# Hurdle rank test
start_time <- Sys.time()
pp <- scan_hurdle_rank(est = pred3, dat = true3, total_iter = 1, iter = 10000)
end_time <- Sys.time()
print(end_time - start_time)

# ZIP risk score
risk_score(data %>% filter(time >= 20220901 & time <= 20221231),
	pred3, pp)

# Second period setting
counts2 <- data %>%
  filter(time >= 20230101 & time < 20230501) %>% 
  df_to_matrix(time_col = 'time', location_col = 'SIG_CD', value_col = 'infected')

# Poisson GLM for second period
mod_11 = glm(infected ~ offset(log(forest))  + distance + distance2 +  time_num + season,
          data = data %>% filter(time < 20230501),
          family = poisson(link = 'log'))

mod_22 = glm.nb(infected ~ offset(log(forest))  + distance + distance2 +  time_num +  season,
              data = data %>% filter(time < 20230501))

mod_33 = hurdle(infected ~ offset(log(forest))  + distance + distance2 +  time_num + season,
                data = data %>% filter(time < 20230501),
                dist = 'poisson', zero.dist = 'binomial')

pred11 <- round(predict(mod_11, newdata = data %>% filter(time >= 20230101 & time <= 20230431), type = 'response'))
true11 <- data$infected[data$time >= 20230101 & data$time <= 20230431]

pred22 <- round(predict(mod_22, newdata = data %>% filter(time >= 20230101 & time <= 20230431), type = 'response'))
true22 <- data$infected[data$time >= 20230101 & data$time <= 20230431]

pred33 <- round(predict(mod_33, newdata = data %>% filter(time >= 20230101 & time <= 20230431), type = 'response'))
true33 <- data$infected[data$time >= 20230101 & data$time <= 20230431]

# Check improvement
## comparing with only counts
df = data.frame(cbind(pred11, true11))
df = df %>% filter(true11 > 0)

rank_df = data.frame(pred_rank = rank(-df$pred11),
           true_rank = rank(-df$true11))
sum((rank_df$true_rank - rank_df$pred_rank)^2)

plot(rank_df$true_rank, rank_df$pred_rank, xlim = c(0, 100), ylim = c(0, 100))
abline(0, 1)

plot(df$true11, df$pred11, xlim = c(0, 12), ylim = c(0, 12))
abline(0, 1)

# Poisson rank test
start_time <- Sys.time()
pp <- scan_pois_rank(est = pred11, dat = true11, total_iter = 1, iter = 10000)
end_time <- Sys.time()

# Poisson risk score
risk_score(data %>% filter(time >= 20230101 & time <= 20230431),
	pred11, pp)

# Negbin rank test
start_time <- Sys.time()
pp <- scan_negbin_rank(est = pred22, dat = true22, total_iter = 1, iter = 10000)
end_time <- Sys.time()
print(end_time - start_time)

# Negbin risk score
risk_score(data %>% filter(time >= 20230101 & time <= 20230431),
	pred22, pp)


# Hurdle rank test
start_time <- Sys.time()
pp <- scan_hurdle_rank(est = pred33, dat = true33, total_iter = 1, iter = 10000)
end_time <- Sys.time()
print(end_time - start_time)

# ZIP risk score
risk_score(data %>% filter(time >= 20230101 & time <= 20230431),
	pred33, pp)