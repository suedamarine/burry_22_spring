library(tidyverse)

# load data
temp_sizes <- read.csv("data/temp_sizes.csv")

# mutate to get total count
cockle_sizes <- temp_sizes %>%
  mutate(total = y0 + y1 + y2)

# repeat the cockle sizes by their abundance
cockleT <- rep(c(cockle_sizes$length[1:37]), times = cockle_sizes$total[1:37])

# set some breaks (by one milimetre)
breaks = seq(3,40, by = 1)

# cut the cockle sizes into sizes defined by breaks
cockleT.cut <- cut(cockleT, breaks, right=FALSE)

# return number of cockles in each size bracket (3,4) (4,5)...
cockleT.frequency <- table(cockleT.cut)

# return total number of cockles counted
sum_t_freq <- sum(cockleT.frequency) 

# calculate relative frequency per size bracket
cockleT.relfreq <- cockleT.frequency / sum_t_freq

# cumulative sum of cockle abundance per size bracket
cumrelfreqT <- c(0,cumsum(cockleT.relfreq))

# find the horizontal value that meets 15mm cockle length
cu_h <- cumrelfreqT[13]

# harvestable proportion
harvestable_proportion <- 1 - cu_h

# write file to csv
write.csv(harvestable_proportion, "tabs/harvestable_proportion.csv")

# find the number of cockles 15mm and above
harvestable_number <- harvestable_proportion * summary_counts[3,3]

# write table to csv
write.csv(harvestable_number, "tabs/harvestable_number.csv")

# make df
rel_freq.df <- data.frame(breaks, cumrelfreqT)

#plot
p <- ggplot(rel_freq.df, aes(x= breaks, y = cumrelfreqT)) +
  geom_point(col =  "#a6611a") + geom_line(col =  "#a6611a") +
  geom_vline(xintercept = 15, col = "#dfc27d", linetype = "dashed", alpha = 0.9) +
  geom_hline(yintercept = cu_h, col = "#dfc27d", linetype = "dashed", alpha = 0.9 ) +
  theme_minimal() +
  labs(x= "Length mm", y = "Cumulative Frequency")


# Open a pdf file
pdf("plots/cumulative_size.pdf") 

# 2. Create a plot
p

# Close the pdf file
dev.off() 
