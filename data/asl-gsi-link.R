#wrangling then merging GSI and ASL data into one object. 

library(tidyverse)
library(readxl)
library(here)
library(lubridate)

#read in and wrangle tables-------------------------------------------------------------- 
  #to wrangle is to make all names and types the same so they can be joined, and to drop
  #observations/columns that aren't needed. 

#gsi
gsi <- readxl::read_xlsx(here("data/gsi/ChinookYukon_Retro_Sept_22_2021.xlsx"), 
                         sheet = 'Individual Region IDs', 
                         skip = 3) |>
  filter(!(Comment %in% c("YukonRetro", "Number of missing Loci exceeded"))|is.na(Comment), 
         !is.na(Fish))

colnames(gsi)[13:16] <- c("Region 6", "Prob 6", "Region 7", "Prob 7")

gsi <- gsi |>
  cbind(as.data.frame(str_split_fixed(gsi$Fish, " ", 5))) |>
  select(-Fish, -V1) |>
  rename(Year = V2, 
         Gear = V3, 
         JulDate = V4, 
         Fish = V5) |>
  mutate(Year = str_sub(Year, 2,3)) |>
  mutate(Year = as.numeric(ifelse(Year < 85, paste0(20, Year), paste0(19, Year))), 
         Fish = as.numeric(Fish), 
         JulDate = as.numeric(JulDate)) |>
  relocate(Fish, Year, JulDate, Gear, .before=Comment)


#asl
asl <- read.csv(here("data/asl/ASL_Output_Chinook_FishWheels_1982-2012.csv")) |>
  filter(Species == "Chinook", Gear == "Fishwheel") |>
  rename(Fish = Fish.Number, 
         Year = Sample.Year) |>
  mutate(JulDate = yday(as.Date(sampleDate)), 
         Year = as.numeric(Year), 
         Fish = as.numeric(Fish))

#index
index <- read.table(here("data/gsi/ExtractionSheets1982-2012.txt"), 
                    fill = TRUE,
                    sep = "\t",
                    header=TRUE) |>
  mutate(Fish = as.numeric(Fish), 
         JulDate = as.numeric(JulDate)) |>
  select(Year, JulDate, Fish, Vial)   

#merge-----------------------------------------------------------------------------------
gsi_indexed <- left_join(gsi, index, by = c("Year", "JulDate", "Fish"))

if(nrow(distinct(gsi_indexed, Fish, JulDate, Year, Vial)) != nrow(gsi_indexed)){
  warning("something's buggered with the gsi-index")
}

nrow(distinct(asl, Fish, JulDate, Year, Vial))

asl_gsi <- full_join(asl, gsi_indexed, by = c("JulDate", "Year", "Fish", "Vial"))

nrow(gsi_indexed) + nrow(asl) #why is this only 1 longer than asl_gsi?
