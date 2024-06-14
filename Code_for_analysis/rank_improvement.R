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
surv_file <- st_read(paste0(gpd_file_path, 'korea_with_surveillance.shp'))
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

# Rank improvement
mod_1 = glm(infected ~ offset(log(forest))  + distance + distance2 +  time_num + season,
            data = data %>% filter(time < 20220901),
            family = poisson(link = 'log'))

pred <- round(predict(mod_1, newdata = data %>% filter(time >= 20220901 & time <= 20221231), type = 'response'))
true <- data$NUMPOINTS[data$time >= 20220901 & data$time <= 20221231]

rank.df <- data.frame(pred=pred, true=true, pred_rank=rank(-pred, ties.method='random', na.last='keep'),
           true_rank=rank(-true, ties.method='random', na.last='keep'))
rank.df <- rank.df %>% filter(true > 0) %>% filter(true > pred)

plot(rank.df$true, rank.df$pred, xlab = 'True', ylab = 'Predicted',
     cex = 1.2, pch = 20, ylim = c(0, 11), xlim = c(0,11))
abline(0, 1, col = 'red')

plot(rank.df$true_rank, rank.df$pred_rank, xlab = 'True Rank', ylab = 'Predicted Rank',
     cex = 1.2, pch = 20, ylim = c(0, 100), xlim = c(0, 100))
abline(0, 1, col = 'red')

mod_11 = glm(NUMPOINTS ~ offset(log(forest))  + distance + distance2 +  time_num + season,
             data = data %>% filter(time < 20230101),
             family = poisson(link = 'log'))
pred11 <- round(predict(mod_11, newdata = data %>% filter(time >= 20230101 & time <= 20230431), type = 'response'))
true11 <- data$NUMPOINTS[data$time >= 20230101 & data$time <= 20230431]

set.seed(123)
rank.df11 <- data.frame(pred=pred11, true=true11, pred_rank=rank(-pred11, ties.method='average', na.last='keep'),
           true_rank=rank(-true11, ties.method='average', na.last='keep'))
rank.df11 <- rank.df11 %>% filter(true11 > 0) %>% filter(true > pred)

mod_11 = glm(NUMPOINTS ~ offset(log(forest))  + distance + distance2 +  time_num + season,
             data = data %>% filter(time < 20210901),
             family = poisson(link = 'log'))
pred11 <- round(predict(mod_11, newdata = data %>% filter(time >= 20210901 & time <= 20211231), type = 'response'))
true11 <- data$NUMPOINTS[data$time >= 20210901 & data$time <= 20211231]

set.seed(123)
rank.df11 <- data.frame(pred=pred11, true=true11, pred_rank=rank(-pred11, ties.method='average', na.last='keep'),
                        true_rank=rank(-true11, ties.method='average', na.last='keep'))
rank.df11 <- rank.df11 %>% filter(true11 > 0) %>% filter(true > pred)
