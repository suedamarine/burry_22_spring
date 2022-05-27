library(tidyverse)

temp_sizes <- read.csv("data/temp_sizes.csv")

cumulative_mass <- temp_sizes %>%
  mutate(total = y0 + y1 + y2)


y0 <- temp_sizes$y0

y1 <- temp_sizes$y1

y2 <- temp_sizes$y2

total <- y0 + y1 + y2

cu_length <- cumulative_mass$length[1:37]

total_cut <- cumulative_mass$total[1:37]

## manually enter summary_tab$mass_totals[15] if not loaded (total mass in kilos) - 15608750 in autumn 2019

total_mass_tonnes <- summary_mass[4,3]/1000

est_mass <-  exp((log(cumulative_mass$length[1:37])* 3.315)-8.939)
size.biomass <- cumulative_mass$total[1:37] * est_mass

size.biomass.sum <- sum(size.biomass)

bio.relfreq <- size.biomass/size.biomass.sum

bio.cumfreq <- cumsum(bio.relfreq)

cum.biomass <- bio.cumfreq*total_mass_tonnes$mass_totals[1]

inv.bio <- 1-bio.cumfreq

inv.mass <- inv.bio*total_mass_tonnes$mass_totals[1]
cu_mass_df <- data.frame(cumulative_mass$length[1:37], cum.biomass, inv.mass)

cu_h_mass <- cu_mass_df$inv.mass[12]

cu_hor_mass <- bio.cumfreq[12]

# proportion of harvetable cockles
proportion_harvestable_mass <- 1 - cu_hor_mass

# write file to csv
write.csv(proportion_harvestable_mass, "tabs/proportion_harvestble_mass.csv")

cu_mass <- ggplot(cu_mass_df, aes(cumulative_mass$length[1:37],bio.cumfreq)) +
  geom_point(col = "#a6611a") +
  geom_line(col = "#a6611a") +
  geom_vline(xintercept = 15, col = "#dfc27d", linetype = "dashed", alpha = 0.9) +
  geom_hline(yintercept = cu_hor_mass, col = "#dfc27d", linetype = "dashed", alpha = 0.9 ) +
  theme_minimal() +
  labs(x= "Length mm", y = "Cumulative Relative Biomass")

cu_mass

# Open a pdf file
pdf("plots/cumulative_mass.pdf") 

# 2. Create a plot
cu_mass

# Close the pdf file
dev.off() 
