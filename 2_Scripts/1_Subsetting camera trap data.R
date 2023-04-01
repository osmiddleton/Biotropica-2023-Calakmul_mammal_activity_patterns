#  
#  "Neo-tropical felid activity patterns in relation to potential prey and
#  intra-guild competitors in the Calakmul Biosphere Reserve, Mexico"                                                                                                #                                                                                                 #                                                                          #
#  Authors: Cristina Argudin-Violante, Owen Middleton, Kathy Slater, 
#           Esteban Dominguez-Bonilla & Patrick Doncaster                                                  
#  Published in: Biotropica   
#  Date (Accepted): April 2023

# Script 1: Tidying and reformatting raw camera trap data
#
# Written by: Owen Middleton
# Contact: o.middleton@sussex.ac.uk
# Last updated: 01/04/2023

# Notes on data and scripts:
# The camera trap data was processed using the software Timelapse V1.0. The 
# format of the camera trap data is the standard format of the software output.
# We have removed the coordinates of the camera trap locations because of 
# sensitivity for future publications (coming soon, contact the lead author CAV).

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# Housekeeping -----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
if(1){
  # Load packages ----
  library(tidyverse)
  library(lubridate)
  
  # Load data ----
  dat <- as_tibble(read.csv("./1_Data/Calakmul 2019 CTs.csv", header = TRUE))

  # Tidy data ----
  # Get bird and mammals into the same column 
  dat[dat$Mammal_Species == " ",] <- NA
  dat[dat$Bird_Species == " ",] <- NA
  
  dat$Species <- paste(dat$Mammal_Species, dat$Bird_Species, sep = "")
  dat <- dat %>% filter(!Species == "")
  
  dat$date_time <- paste(dat$Date, dat$Time)
  dat$date_time
  
  dat$date_time <- dmy_hms(dat$date_time)
  str(dat)
}
  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# Step 1. Remove duplicate records of the same camera trap record ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# Need to not remove rows if they're just the duplicate photos.
# Do not keep row if its the same species up to 5 seconds after the first
# record.
if(1){
new.dat <- data.frame()
for(i in 1:nrow(dat)) {
  if(i == 1) {
    # If it's the first row, just add this to the tidy dataframe
    new.dat <- rbind(new.dat, dat[i,])
    print(i)
  } else {
    # Get first record time
    previousrecord <- ymd_hms(dat[i-1,]$date_time)
    previousspecies <- dat[i-1,]$Species
    # Get the time-range for the same record (i.e. the photos that are duplicates)
    identifysamerecord <- interval(previousrecord, previousrecord + seconds(5))
    # Get the next record/row...
    newrecord <- ymd_hms(dat[i,]$date_time)
    newspecies <- dat[i,]$Species
    #...and see if it's the same species in the duplicates time frame:
      if(newrecord %within% identifysamerecord &
         newspecies == previousspecies) {
        # Same species in time frame
        print(i)
      } else {
        # Same or different species outside of the duplicates time window.
       new.dat <- rbind(new.dat, dat[i,])
       print(i)
    }
  }
}
}
# Note: I know there's probably a more efficient way to do this, I'm just
# reviewing this code right now slightly hungover in a cafe and do not want to 
# rewrite anything....

# Save this step, if it's the first time doing it...
if(0){
write.csv(new.dat, "./1_Data/Duplicate records removed.csv")
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# Step 2. Remove non-independent records ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# At each site (single or camera stations), go through species records and remove
# those that are not independent (i.e. within 30 minutes)
if(1){

# If the above is already done, load the dataset where duplicates are removed...
dat <- read.csv("./1_Data/Duplicate records removed.csv")
# Reformat data:
if(1){
# First, arrange whole dataframe by date and time
dat <- dat %>% arrange(date_time)
# Remove the first problematic record
dat <- dat[-1,]
# Make sure Sites is a factor:
dat$Site <- as.factor(dat$Site)
}

# Create an empty dataframe to save the non-pseudoreplicates to
tidy.data <- data.frame()

# For each site, go through each species and remove records
# within 30 minutes of one another.
for(c in 1:length(levels(dat$Site))) {

  # Go through each site separately...
  site <- levels(dat$Site)[c]
  site.data <- dat %>% filter(Site == site)
  
  # New dataframe to filter through each site
  new.dat <- data.frame()

  for(i in 1:length(levels(factor(site.data$Species)))){
      
      # Go through each species at each site separately
      spec <- levels(factor(site.data$Species))[i]
      spec.data <- site.data %>% filter(Species == spec)
      
      # Create new dataframe to save non-pseudos to...
      new.dat2 <- data.frame()
      
      for(j in 1:nrow(spec.data)) {
        
        if(j == 1) {
            # If it's the first record in the site of a species, it'll be 
            # independent
            new.dat2 <- rbind(new.dat2, spec.data[j,])
            print(j)
            } else {
              # Get the time stamp of the previous record.
              previousrecord <- ymd_hms(new.dat2[nrow(new.dat2),]$date_time)
              # Get the time frame for the next record not being independent
              # (i.e. captured within 30 minutes of the previous record).
              identifysamerecord <- interval(previousrecord,previousrecord + minutes(30))
              # Get the time stamp of the new record.
              newrecord <- ymd_hms(spec.data[j,]$date_time)
              if(newrecord %within% identifysamerecord) {
                # New record is within the time frame, don't keep
                print(j)
              } else {
                # New record is outside of the time frame, keep as independent
                new.dat2 <- rbind(new.dat2, spec.data[j,])
                print(j)
              }
            }
      }
      # Once all species from a site have been done, save it and go to the next site
      new.dat <- rbind(new.dat, new.dat2)
    }
  
  tidy.data <- rbind(tidy.data, new.dat)
  print(i)
    
}
}

# Save the tidied dataframe:
write.csv(tidy.data, "./1_Data/Final 2019 Calakmul dataset.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# End of Script ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----