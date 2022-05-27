# load libraries
library(tidyverse)

# import data
temp_sizes <- read.csv("data/temp_sizes.csv")

# select sizes relavent to season
sizes <- temp_sizes %>%
  select("length",  "y1", "y2") %>%
  gather(year_class, quantity,  "y1", "y2") 

# produce new vectors for individual sizes
temp_length <- rep(sizes$length, times = sizes$quantity) 
temp_class <- rep(sizes$year_class, times = sizes$quantity)

# new data frame for plots
temp_sizes_df <- data.frame(temp_length, temp_class)

# calculate means for each year class 
mu <- temp_sizes_df %>% group_by(temp_class) %>% summarize(grp_med = median(temp_length))

head(mu)


# write file to csv
write.csv(mu, "tabs/mu.csv")

# plot 
p <- ggplot(temp_sizes_df, aes(x=temp_length, color=temp_class, fill=temp_class, position="dodge")) +
  geom_histogram(aes(y=..density..), binwidth = 1.0, alpha=0.4, position="identity") +
  geom_density(alpha=0.2, bw = 0.75) +
  geom_vline(data = mu, aes(xintercept=grp_med, colour = temp_class), linetype = "dashed") +
  scale_color_manual(values=c("#dfc27d","#80cdc1"), guide = "none")  +
  scale_fill_manual(name = "Class",labels = c("Year1", "Year2+"), values=c("#dfc27d","#80cdc1")) +
  theme(legend.position="top") +
  labs(title="Spring 2022 Cockle Size Distribution", x="Size mm", y = "Frequency") +
  theme_classic() +
  ylim(0,0.25) +
  xlim(4,45)

# Open a pdf file
pdf("plots/freq_dist.pdf") 

# 2. Create a plot
p

# Close the pdf file
dev.off()
