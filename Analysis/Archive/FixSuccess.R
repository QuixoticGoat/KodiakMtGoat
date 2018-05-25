################################################################################
# Creates a summary dataframe of mountain goat ATS Globalstar GPS collar data  #
# and plots fix success by collar id, includes a line of the average           #
# mountain goats on Kodiak Island, Alaska                                      #
#                                                                              #
# Author: McCrea Cobb <mccrea_cobb@fws.gov>                                    #
# Date created: 8/1/2017                                                       #
# Date edited: 4/16/2018                                                       #
################################################################################


## ---- FixSuccess

# Create summary data frame:
sum <- data.frame(table(df$CollarID))
colnames(sum) <- c("CollarID", "Fixes")
sum$Days <- max(df$Day) - min(df$Day)
sum$FixSucc <-  round((sum$Fixes/2)/(sum$Days) * 100, 2)
# Remove mortalities
sum <- subset(sum, CollarID != "IG32" & CollarID != "IG39" & CollarID != "IG44")
sum
summary(sum)

# Plot fix success by CollarID, show average line:
library(ggplot2)
p <- ggplot(sum, aes( x = CollarID, y = FixSucc))
p <- p +
  geom_point(size = 2) +
  geom_hline(yintercept = mean(sum$FixSucc), linetype = 2) +
  labs(x = "Collar ID", y = "Fix success rate (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
p
