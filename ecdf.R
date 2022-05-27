library(tidyverse)

# get data
temp_sizes <- read.csv("data/temp_sizes.csv")

# mutate to get total count
cockle_sizes <- temp_sizes %>%
  mutate(total = y0 + y1 + y2)

cockleT <- rep(c(cockle_sizes$length[1:37]), times = cockle_sizes$total[1:37])

cockleT_df <- data.frame(cockleT)

cockle_ecdf <- data.frame(x=unique(cockleT_df$cockleT),  y=ecdf(cockleT_df$cockleT)(unique(cockleT_df$cockleT))*length(cockleT_df$cockleT))

cockle_ecdf$y <- scale(cockle_ecdf$y,center=min(cockle_ecdf$y),scale=diff(range(cockle_ecdf$y)))

hl <- cockle_ecdf$y[9]

print(1 - hl)

cu_p <- ggplot(cockleT_df, aes(cockleT)) +
  stat_ecdf(geom ="point", col = "#a6611a") +
  geom_vline(xintercept = 15, col = "#dfc27d", linetype = "dashed", alpha = 0.9) + geom_hline(yintercept = hl, col = "#dfc27d", linetype = "dashed", alpha = 0.9 ) + geom_line(data = cockle_ecdf, aes(x, y, col =  "#a6611a"),show.legend = FALSE) + theme_minimal() + labs(x= "Length mm", y = "Cumulative Frequency")

cu_p

# Open a pdf file
pdf("plots/ecdf.pdf") 

# 2. Create a plot
cu_p

# Close the pdf file
dev.off() 
