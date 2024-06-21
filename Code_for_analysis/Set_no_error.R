library(dplyr)
library(spdep)
library(ggplot2)
library(tidyr)
library(readxl)
library(sp)
library(magrittr)
library(terra)
library(sf)
library(MASS)
library(pscl)


source('./test_functions.r')
load('../Data/surv_intensity_week.Rdata')

# Set path
gpd_file_path = '../Data/'
surv_file <- st_read(paste0(gpd_file_path, 'korea_with_surveillance2.shp'))
surv_file$SIG_CD <- as.numeric(surv_file$SIG_CD)

# Data
shp <- st_read(paste0(gpd_file_path, 'shp통합4.shp'))
coord <- read_excel(paste0(gpd_file_path, '시군구위도경도.xlsx'))
used_data <- read.csv('../Data/weekly_merge.csv', header = TRUE)
surv_data <- read.csv('../Data/add_surveilance.csv', header = TRUE)

surv_file$SIG_CD <- as.numeric(surv_file$SIG_CD)

# sort and add infected
used_data <- used_data %>% arrange(SIG_CD, time)
used_data['infected'] <- surv_data$Infected

# Convert some types
used_data$time <- as.numeric(used_data$time)
used_data$distance <- exp(-1 * used_data$distance)
used_data$sin <- sin(2*pi*used_data$time/52)
use_colnames <- c('time', 'infected', 'time_num', 'sin', 'distance', 'forest', 'NUMPOINTS', 'season', 'SIG_CD', 'x', 'y')
data <- used_data[, use_colnames]

# add distance^2
data$distance2 <- data$distance^2

# Add observation error variables
data$intensity <- surveillance_intensity
data$obs_error <- (3-data$intensity)
data[data$intensity == 0, 'obs_error'] <- 0

# Add observation error with surveillance intensity
plot(data$NUMPOINTS, data$infected)
abline(0, 1)

# For error : 0.1
# Poisson GLM for first period
mod_1 = glm(NUMPOINTS ~ offset(log(forest))  + distance + distance2 +  time_num + season,
            data = data %>% filter(time < 20220901),
            family = poisson(link = 'log'))

mod_2 = glm.nb(NUMPOINTS ~ offset(log(forest))  + distance + distance2 +  time_num + season,
               data = data %>% filter(time < 20220901))

mod_3 = hurdle(NUMPOINTS ~ offset(log(forest))  + distance + distance2 +  time_num + season,
               data = data %>% filter(time < 20220901),
               dist = 'poisson', zero.dist = 'binomial')

pred <- round(predict(mod_1, newdata = data %>% filter(time >= 20220901 & time <= 20221231), type = 'response'))
true <- data$NUMPOINTS[data$time >= 20220901 & data$time <= 20221231]

# Most highest estimate region
# dat <- data.frame(true=true3, pred=pred3)
# ind <- dat[which(pred == max(pred3)),]
# new.dat <- data %>% filter(time >= 20220901 & time <= 20221231)
# new.dat$true <- dat$true
# new.dat$pred <- dat$pred
# new.dat %>% group_by(SIG_CD) %>% summarise(est_pred = sum(pred),
#                                            est_true = sum(true)) %>%
#   arrange(desc(est_pred))
# new.dat %>% group_by(SIG_CD) %>% summarise(est_pred = sum(pred),
#                                            est_true = sum(true)) %>%
#   arrange(desc(est_true))

pred2 <- round(predict(mod_2, newdata = data %>% filter(time >= 20220901 & time <= 20221231), type = 'response'))
true2 <- data$NUMPOINTS[data$time >= 20220901 & data$time <= 20221231]

pred3 <- round(predict(mod_3, newdata = data %>% filter(time >= 20220901 & time <= 20221231), type = 'response'))
true3 <- data$NUMPOINTS[data$time >= 20220901 & data$time <= 20221231]

# plot(true, pred, ylim = c(0, 12))
# points(true2, pred2, col = 'red', pch = 2)
# points(true3, pred3, col = 'blue', pch = 3)
# abline(0, 1)
# 
# plot(true2, pred2, ylim = c(0, 12))
# abline(0, 1)

# Poisson rank test
start_time <- Sys.time()
pp <- scan_pois_rank(est = pred, dat = true, total_iter = 1, iter = 10000)
end_time <- Sys.time()
print(end_time - start_time)

# Poisson risk score
risk_score(data %>% filter(time >= 20220901 & time <= 20221231),
           pred, pp, infected = 'NUMPOINTS')

# Negbin rank test
start_time <- Sys.time()
pp <- scan_negbin_rank(est = pred2, dat = true2, total_iter = 1, iter = 10000)
end_time <- Sys.time()
print(end_time - start_time)

# Negbin risk score
risk_score(data %>% filter(time >= 20220901 & time <= 20221231),
           pred2, pp, 'NUMPOINTS')

# Hurdle rank test
start_time <- Sys.time()
pp <- scan_hurdle_rank(est = pred3, dat = true3, total_iter = 1, iter = 10000)
end_time <- Sys.time()
print(end_time - start_time)

# ZIP risk score
risk_score(data %>% filter(time >= 20220901 & time <= 20221231),
           pred3, pp, 'NUMPOINTS')

# Second period setting
# Poisson GLM for second period
mod_11 = glm(NUMPOINTS ~ offset(log(forest)) + distance + distance2 +  time_num + season,
             data = data %>% filter(time < 20230101),
             family = poisson(link = 'log'))

mod_22 = glm.nb(NUMPOINTS ~ offset(log(forest))  + distance + distance2 +  time_num +  season,
                data = data %>% filter(time < 20230101))

mod_33 = hurdle(NUMPOINTS ~ offset(log(forest))  + distance + distance2 +  time_num + season,
                data = data %>% filter(time < 20230101),
                dist = 'poisson', zero.dist = 'binomial')

pred11 <- round(predict(mod_11, newdata = data %>% filter(time >= 20230101 & time <= 20230431), type = 'response'))
true11 <- data$NUMPOINTS[data$time >= 20230101 & data$time <= 20230431]

# Most highest estimate region
dat <- data.frame(true=true11, pred=pred11)
ind <- dat[which(pred == max(pred11)),]
new.dat <- data %>% filter(time >= 20230101 & time <= 20230431)
new.dat$true <- dat$true
new.dat$pred <- dat$pred
new.dat %>% group_by(SIG_CD) %>% summarise(est_pred = sum(pred),
                                           est_true = sum(true)) %>%
  arrange(desc(est_pred))
new.dat %>% group_by(SIG_CD) %>% summarise(est_pred = sum(pred),
                                           est_true = sum(true)) %>%
  arrange(desc(est_true))

pred22 <- round(predict(mod_22, newdata = data %>% filter(time >= 20230101 & time <= 20230431), type = 'response'))
true22 <- data$NUMPOINTS[data$time >= 20230101 & data$time <= 20230431]

pred33 <- round(predict(mod_33, newdata = data %>% filter(time >= 20230101 & time <= 20230431), type = 'response'))
true33 <- data$NUMPOINTS[data$time >= 20230101 & data$time <= 20230431]

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
           pred11, pp, 'NUMPOINTS')

# Negbin rank test
start_time <- Sys.time()
pp <- scan_negbin_rank(est = pred22, dat = true22, total_iter = 1, iter = 10000)
end_time <- Sys.time()
print(end_time - start_time)

# Negbin risk score
risk_score(data %>% filter(time >= 20230101 & time <= 20230431),
           pred22, pp, 'NUMPOINTS')


# Hurdle rank test
start_time <- Sys.time()
pp <- scan_hurdle_rank(est = pred33, dat = true33, total_iter = 1, iter = 10000)
end_time <- Sys.time()
print(end_time - start_time)

# ZIP risk score
risk_score(data %>% filter(time >= 20230101 & time <= 20230431),
           pred33, pp, 'NUMPOINTS')