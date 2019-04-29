
################################################################################
# Summary tables of GPS collar data for Kodiak mountain goats                  #
#                                                                              #
# Author: McCrea Cobb <mccrea_cobb@fws.gov>                                    #
# Last modified: 5/23/2018                                                     #
################################################################################


#-------------------------------------------------------------------------------
##---- GPS.Summary.R

GPS.Summary <- function(dat) {
  # Function to create a data frame of summary information for each GPS collars, including:
  # Collar ID (CollarID),
  # Number of fixes/collar (Fixes),
  # Date of first fix (Start),
  # Date of last fix (End),
  # Number of days between first and last fix (Days), and
  # Average fix interval, in days (Interval).

  require(dplyr)

  dat <- dat[order(dat$Date), ]
  dat <- dat[order(dat$CollarID), ]

  dfSummary = as.data.frame(dat %>%
                              arrange(CollarID, Date) %>%
                              group_by(CollarID, Collar) %>%
                              summarise(NumFixes = length(CollarID),
                                        StartDate = min(Date),
                                        EndDate = max(Date),
                                        Days = length(seq(min(Date), max(Date), by='day')),
                                        FixInt = round(mean(as.vector(Date - lag(Date)),
                                                            na.rm = T),
                                                       digits = 2)))
  dfSummary$StartDate = strptime(dfSummary$StartDate, format = "%Y-%m-%d")
  dfSummary$EndDate = strptime(dfSummary$EndDate, format = "%Y-%m-%d")
  dfSummary <- dfSummary[order(as.character(dfSummary$CollarID)), ]
  row.names(dfSummary) <- NULL

  return(dfSummary)
}


# Run it:
#dfSummary <- GPS.Summary(df)

#write.csv(dfSummary, "Output/Tables/SummaryTbl.csv")



#-------------------------------------------------------------------------------
## Produce tables of output from the top models

# Required packages:
library(mclogit)
library(broom) # For tidy()


################
###### Summer
################

# Load the best model
load(file = "./Data/Rdata/mod_summer_final.Rdata")

# Create data frame of the model output:
best.df <- tidy(best)
best.df[, c(2:5)] <- round(best.df[ , c(2:5)], 2)

# Write as a .csv:
write.csv(best.df, "./Manuscript/Tables/ModelOutput_Summer.csv", 
          row.names = FALSE)


###############
###### Winter
###############

# Load the best model:
load(file = "./Data/Rdata/mod_winter_final.Rdata")

# Create data frame of the model output:
best.df <- tidy(best)
best.df[, c(2:5)] <- round(best.df[ , c(2:5)], 2)

# Write as a .csv:
write.csv(best.df, "./Manuscript/Tables/ModelOutput_Winter.csv", 
          row.names = FALSE)

