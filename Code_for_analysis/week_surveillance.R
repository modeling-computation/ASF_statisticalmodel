library(dplyr)
library(spdep)
library(tidyr)
library(readxl)
library(sp)
library(magrittr)
library(terra)
library(sf)

setwd('D:/Dropbox/2024/2024 ASF/ASF_figure_code/Analysis_code/')

# 주별 감시강도 설정
gpd_file_path = '../Data/'
surv_file <- st_read(paste0(gpd_file_path, 'korea_with_surveillance.shp'))

total <- cbind(surv_file$SIG_CD, surv_file$X2019_2, surv_file$X2020_1, surv_file$X2020_2,
	surv_file$X2021_1, surv_file$X2021_2, surv_file$X2022_1, surv_file$X2022_2,
	surv_file$X2023_1)
total <- apply(total, 2, as.numeric) # 숫자로 변경


# 2019년 : 13개, 2023년 : 18개, 나머지 : 26개
surveillance_intensity <- c()
for (row in 1:nrow(total)) {
	# 2019 year
	for (j in 1:13) {
		surveillance_intensity <- c(surveillance_intensity, total[row, 2])
	}
	# 2020 ~ 2022 years
	year_lst <- seq(3,8,1)
	for (year in year_lst) {
		for (j in 1:26) {
			surveillance_intensity <- c(surveillance_intensity, total[row, year])
		}
	}
	# 2023 year
	for (j in 1:18) {
		surveillance_intensity <- c(surveillance_intensity, total[row, 9])
	}
}

# Save
file_path <- paste0(gpd_file_path, "surv_intensity_week.Rdata")
save(surveillance_intensity, file = file_path)

