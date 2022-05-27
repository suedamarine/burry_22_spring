library(rcompanion)
library(tidyverse)
library(lubridate)

## import file - in future create new data.frame from CSV, using vectors to import count & mass into individual columns

temp_cockle <- read.csv("data/temp_cockle.csv")

# define 250grid stations that overlap 100grid
north_stations <- c(187, 188, 189, 207, 208, 209, 227, 228, 229)
south_stations <- c(391, 392, 393, 395, 396, 400, 401, 402, 478)

# number of stations visited in each grid
sampled_250 <- sum(temp_cockle$Sampled == "Y" & temp_cockle$Grid =="250", na.rm = TRUE)
sampled_100n <- sum(temp_cockle$Sampled == "Y" & temp_cockle$Block =="ZN", na.rm = TRUE)
sampled_100s <- sum(temp_cockle$Sampled == "Y" & temp_cockle$Block =="ZS", na.rm = TRUE)
sampled_250_n <- sum(temp_cockle$Sampled == "Y" & temp_cockle$Stn %in% north_stations, na.rm = TRUE)
sampled_250_s <- sum(temp_cockle$Sampled == "Y" & temp_cockle$Stn %in% south_stations, na.rm = TRUE)

# number of samples not accessible and not saltmarsh
not_sampled <- sum(temp_cockle$Sampled == "N" & temp_cockle$Grid =="250" & temp_cockle$Substrata != "saltmarsh")

saltmarsh <- sum(temp_cockle$Substrata == "saltmarsh" & temp_cockle$Grid =="250", na.rm = TRUE)

# start date
start_date <- temp_cockle %>%
  mutate(Date = dmy(Date)) %>%
  filter(Date == min(Date)) %>%
  head()

# end date
end_date <- temp_cockle %>%
  mutate(Date = dmy(Date)) %>%
  filter(Date == max(Date)) %>%
  head()

# import and tidy count data
cockle_counts <- temp_cockle %>%
  select(Stn, Block, Grid, Sampled, Y0Count, Y1Count, Y2Count, Y3Count) %>% 
  mutate(total_count = Y0Count + Y1Count + Y2Count + Y3Count) %>% 
  gather(year_class, count, 'Y0Count', 'Y1Count', 'Y2Count', 'Y3Count', 'total_count')

# summary table count
summary_counts <- cockle_counts %>% 
  filter(Grid %in% c ("250", "100")) %>%
  mutate(count_sum = count*10*Grid^2) %>%
  group_by(year_class, Grid) %>%
  summarize(count_totals = sum(count_sum, na.rm = TRUE))

# write file to csv
write.csv(summary_counts, "tabs/summary_counts.csv")

# summary table count (selected blocks)
summary_counts_blocks <- cockle_counts %>% 
  filter(Block %in% c ("NB", "NC", "ND", "NE")) %>%
  mutate(count_sum = case_when(Grid == "250m" ~ count*10*250^2, Grid != "250m" ~ count*10*100^2)) %>%
  group_by(year_class, Block) %>%
  summarize(count_totals = sum(count_sum, na.rm = TRUE))

# import and tidy mass data
cockle_mass <- temp_cockle %>%
  select(Stn, Block, Grid, Sampled, Y0Weight, Y1Weight, Y2Weight, Y3Weight, fifteen) %>%
  mutate(total_weight = Y0Weight + Y1Weight + Y2Weight+ Y3Weight) %>%
  gather(year_class, mass, 'Y0Weight', 'Y1Weight', 'Y2Weight', 'Y3Weight', 'total_weight', 'fifteen')

# check that twenty is not greater than total
cockle_totals <- temp_cockle %>%
  select(Stn, Block, Y0Weight, Y1Weight, Y2Weight, Y3Weight, fifteen) %>%
  mutate(total = Y0Weight + Y1Weight + Y2Weight+ Y3Weight) %>%
  filter(fifteen > total)

# summary table mass
summary_mass <- cockle_mass %>%
  filter(Grid %in% c("250", "100")) %>%
  mutate(mass_sum = mass*10/1000*Grid^2) %>%
  group_by(year_class, Grid) %>%
  summarise(mass_totals = sum(mass_sum, na.rm = TRUE))

# write file to csv
write.csv(summary_mass, "tabs/summary_mass.csv")

# groupwise means for counts

# groupwise means for main survey grids, filtered to remove NAs - remove Y0 during spring

g_250 <- cockle_counts %>% 
  filter(year_class %in% c("Y1Count", "Y2Count", "Y3Count", "total_count")) %>%
  filter(Grid %in% c("250", "100")) %>%
  filter(!is.na(count))

# perform the groupwisemean selecting 10000 replicates and Bca
count_conf_g_250 <- groupwiseMean(count ~ year_class + Grid, data = g_250, conf = 0.95, digits = 3, R = 10000, boot = TRUE, traditional = FALSE, normal = FALSE, basic = FALSE, percentile = FALSE, bca = TRUE)

# calculate confidence intervals
count_intervals_250 <- count_conf_g_250 %>%  
  mutate(count_lower =  Bca.lower * 10 * Grid^2 * n) %>%
  mutate(count_upper = Bca.upper * 10 * Grid^2 * n) %>%
  mutate(total_mean = Mean * 10 * Grid^2*n)

# arrange for easier reading + filter required fields
count_intervals_250_filtered <- count_intervals_250 %>%
  select(-c(Mean, Boot.mean, Conf.level, Bca.lower, Bca.upper)) %>%
  print()

# write file to csv
write.csv(count_intervals_250_filtered, "tabs/count_intervals_250_filtered.csv")

# groupwise means for north counts, filtered to remove NAs
n_250_c <- cockle_counts %>%
  filter(Stn %in% north_stations) %>%
  filter(Sampled == "Y") %>%
  filter(year_class %in% c("Y1Count", "Y2Count", "Y3Count", "total_count")) %>%
  filter(!is.na(count))

# perform the groupwisemean selecting 10000 replicates and Bca
count_conf_n_250 <- groupwiseMean(count ~ year_class + Grid, data = n_250_c, conf = 0.95, digits = 3, R = 10000, boot = TRUE, traditional = FALSE, normal = FALSE, basic = FALSE, percentile = FALSE, bca = TRUE)

# calculate confidence intervals
count_intervals_n250 <- count_conf_n_250 %>%  
  mutate(count_lower =  Bca.lower * 10 * Grid^2 * n) %>%
  mutate(count_upper = Bca.upper * 10 * Grid^2 * n) %>%
  mutate(total_mean = Mean * 10 * Grid^2*n) %>%
  select(-c(Mean, Boot.mean, Conf.level, Bca.lower, Bca.upper)) %>%
  print()

# write file to csv
write.csv(count_intervals_n250, "tabs/count_intervals_n250.csv")

# groupwise means for south counts, filtered to remove NAs
s_250_c <- cockle_counts %>%
  filter(Stn %in% south_stations) %>%
  filter(Sampled == "Y") %>%
  filter(year_class %in% c("Y1Count", "Y2Count", "Y3Count", "total_count")) %>%
  filter(!is.na(count))

# perform the groupwisemean selecting 10000 replicates and Bca
count_conf_s_250 <- groupwiseMean(count ~ year_class + Grid, data = s_250_c, conf = 0.95, digits = 3, R = 10000, boot = TRUE, traditional = FALSE, normal = FALSE, basic = FALSE, percentile = FALSE, bca = TRUE)

# calculate confidence intervals
count_intervals_s250 <- count_conf_s_250 %>%  
  mutate(count_lower =  Bca.lower * 10 * Grid^2 * n) %>%
  mutate(count_upper = Bca.upper * 10 * Grid^2 * n) %>%
  mutate(total_mean = Mean * 10 * Grid^2*n) %>%
  select(-c(Mean, Boot.mean, Conf.level, Bca.lower, Bca.upper)) %>%
  print()

# write file to csv
write.csv(count_intervals_s250, "tabs/count_intervals_s250.csv")

# groupwise means for mass

# groupwise means for main survey grids, filtered to remove NAs
g_250_m <- cockle_mass %>% 
  filter(year_class %in% c("Y1Weight", "Y2Weight", "Y3Weight", "total_weight")) %>%
  filter(Grid %in% c("250", "100")) %>%
  filter(!is.na(mass))

# perform the groupwisemean selecting 10000 replicates and Bca
mass_conf_g_250 <- groupwiseMean(mass ~ year_class + Grid, data = g_250_m, conf = 0.95, digits = 3, R = 10000, boot = TRUE, traditional = FALSE, normal = FALSE, basic = FALSE, percentile = FALSE, bca = TRUE)

# calculate confidence intervals
mass_intervals_250 <- mass_conf_g_250 %>%  
  mutate(mass_lower =  Bca.lower * 0.01 * Grid^2 * n) %>%
  mutate(mass_upper = Bca.upper * 0.01 * Grid^2 * n) %>%
  mutate(total_mean = Mean * 0.01 * Grid^2 * n) %>%
  select(-c(Mean, Boot.mean, Conf.level, Bca.lower, Bca.upper)) %>%
  mutate(mass_lower = mass_lower/1000,
         mass_upper = mass_upper/1000,
         total_mean = total_mean/1000) %>%
  print()

# write file to csv
write.csv(mass_intervals_250, "tabs/mass_intervals_250.csv")

# groupwise means for north mass, filtered to remove NAs
n_250_m <- cockle_mass %>%
  filter(Stn %in% north_stations) %>%
  filter(Sampled == "Y") %>%
  filter(year_class %in% c("Y1Weight", "Y2Weight", "Y3Weight", "total_weight")) %>%
  filter(!is.na(mass))

# perform the groupwisemean selecting 10000 replicates and Bca
mass_conf_n_250 <- groupwiseMean(mass ~ year_class + Grid, data = n_250_m, conf = 0.95, digits = 3, R = 10000, boot = TRUE, traditional = FALSE, normal = FALSE, basic = FALSE, percentile = FALSE, bca = TRUE)

# calculate confidence intervals
mass_intervals_n250 <- mass_conf_n_250 %>%  
  mutate(mass_lower =  Bca.lower * 0.01 * Grid^2 * n) %>%
  mutate(mass_upper = Bca.upper * 0.01 * Grid^2 * n) %>%
  mutate(total_mean = Mean * 0.01 * Grid^2 * n) %>%
  select(-c(Mean, Boot.mean, Conf.level, Bca.lower, Bca.upper)) %>%
  mutate(mass_lower = mass_lower/1000,
         mass_upper = mass_upper/1000,
         total_mean = total_mean/1000) %>%
  print()

# write file to csv
write.csv(mass_intervals_n250, "tabs/mass_intervals_n250.csv")

# groupwise means for south mass, filtered to remove NAs
s_250_m <- cockle_mass %>%
  filter(Stn %in% south_stations) %>%
  filter(Sampled == "Y") %>%
  filter(year_class %in% c("Y1Weight", "Y2Weight", "Y3Weight", "total_weight")) %>%
  filter(!is.na(mass))

# perform the groupwisemean selecting 10000 replicates and Bca
mass_conf_s_250 <- groupwiseMean(mass ~ year_class + Grid, data = s_250_m, conf = 0.95, digits = 3, R = 10000, boot = TRUE, traditional = FALSE, normal = FALSE, basic = FALSE, percentile = FALSE, bca = TRUE)

# calculate confidence intervals
mass_intervals_s250 <- mass_conf_s_250 %>%  
  mutate(mass_lower =  Bca.lower * 0.01 * Grid^2 * n) %>%
  mutate(mass_upper = Bca.upper * 0.01 * Grid^2 * n) %>%
  mutate(total_mean = Mean * 0.01 * Grid^2 * n) %>%
  select(-c(Mean, Boot.mean, Conf.level, Bca.lower, Bca.upper)) %>%
  mutate(mass_lower = mass_lower/1000,
         mass_upper = mass_upper/1000,
         total_mean = total_mean/1000) %>%
  print()

# write file to csv
write.csv(mass_intervals_s250, "tabs/mass_intervals_s250.csv")

# mass percentage change
mass_22 <- c(9995000, 0, 5704375, 2933125, 1357500)
mass_21 <- c(6755000, 0, 3735625, 1328750, 1690625)
count_22 <- c(6241250000, 0, 5265625000, 738750000, 236875000)
count_21 <- c(5762500000, 0, 4.98E+09, 461875000, 320625000)

mass_changes <- (mass_22 - mass_21) / mass_21 * 100
count_changes <- (count_22 - count_21) / count_21 * 100

