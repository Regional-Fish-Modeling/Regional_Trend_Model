# Script to process data for use in the abundance model

#--------- Load Libraries------
library(readr)
library(dplyr)
library(tidyr)

#-------- Load Data -----------
visits <- read_csv("Data/BKT_COVS_SAMPLES.csv")

#---------- Get Survey Length ----------

unique(visits$Year)
# check if 2105 is supposed to be 2015 or 2005
visits %>%
  dplyr::filter(Year == 2105) # should be 2015 - convert

visits <- visits %>%
  dplyr::mutate(Year = ifelse(Year == 2105, 2015, Year))

# need in form (site x year) to match N
# Convert to standard 100 m so abundance inference is fish per 100 m
# Fill NA with 100 m (actually =1 since relative to 100 m)

survey_length <- visits %>%
  dplyr::mutate(length_100 = SiteLength_m / 100) %>%
  dplyr::select(SiteID, Year, length_100) %>%
  dplyr::group_by(SiteID) %>%
  tidyr::spread(key = Year, value = length_100, fill = 1)

write_csv(survey_length, "Data/BKT_SURVEY_LENGTH.csv")

#----------- End ----------
# other variables already processed to necessary form