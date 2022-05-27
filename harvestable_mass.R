library(tidyverse)

temp_sizes <- read.csv("data/temp_sizes.csv")


harvestable_mass <- temp_sizes %>%
  mutate(total = y0 + y1 + y2)

## manually enter summary_tab$mass_totals[15] if not loaded (total mass in kilos) - 15608750 in autumn 2019

total_mass_tonnes <- summary_mass[4,3]/1000

est_mass <-  exp((log(harvestable_mass$length[1:37])* 3.315)-8.939)
size.biomass <- harvestable_mass$total[1:37] * est_mass

## manually enter summary_mass_[3,3] if not loaded (total mass in kilos) - 15608750 in autumn 2019

total_mass_tonnes <- summary_mass[4,3]/1000

est_mass <-  exp((log(cu_length)* 3.315)-8.939)

size.biomass <- total_cut*est_mass

size.biomass.sum <- sum(size.biomass)

bio.relfreq <- size.biomass/size.biomass.sum

bio.cumfreq <- cumsum(bio.relfreq)

inv.bio <- 1-bio.cumfreq

inv.mass <- inv.bio*total_mass_tonnes$mass_totals[1]

cu_mass_df <- data.frame(harvestable_mass$length[1:37], cum.biomass, inv.mass)

cu_h_mass <- cu_mass_df$inv.mass[12]

# write file to csv
write.csv(cu_h_mass, "tabs/cu_h_mass.csv")

cu_hor_mass <- bio.cumfreq[12]

cu_harvest_mass <- ggplot(cu_mass_df, aes(harvestable_mass$length[1:37],inv.mass)) +
  geom_point(col = "#a6611a") +
  geom_line(col = "#a6611a") +
  geom_vline(xintercept = 15, col = "#dfc27d", linetype = "dashed", alpha = 0.9) +
  geom_hline(yintercept = cu_h_mass, col = "#dfc27d", linetype = "dashed", alpha = 0.9 ) +
  theme_minimal() +
  labs(x= "Length mm", y = "Biomass Tonnes")

cu_harvest_mass

# Open a pdf file
pdf("plots/harvest_mass.pdf") 

# 2. Create a plot
cu_harvest_mass

# Close the pdf file
dev.off() 

# print dataframe to csv

write.csv(cu_mass_df, "tabs/cu_mass_df.csv")
