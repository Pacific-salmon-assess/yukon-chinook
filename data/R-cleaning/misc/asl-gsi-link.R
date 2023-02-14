#wrangling then merging GSI and ASL data into one object. 

library(tidyverse)
library(readxl)
library(here)
library(lubridate)

#read in and wrangle tables-------------------------------------------------------------- 
  #to wrangle is to make all names and types the same so they can be joined, and to drop
  #observations/columns that aren't needed. 

#gsi
gsi <- readxl::read_xlsx(here("data/gsi/Copy of ChinookYukon_Retro_Sept_22_2021_2021-09-22.xlsx"), 
                         sheet = 'Individual Region IDs', 
                         skip = 3) |>
  filter(!(Comment %in% c("YukonRetro", "Number of missing Loci exceeded"))|is.na(Comment), 
         !is.na(Fish))

colnames(gsi)[13:16] <- c("Region 6", "Prob 6", "Region 7", "Prob 7")

gsi <- gsi |>
  cbind(as.data.frame(str_split_fixed(gsi$Fish, " ", 5))) |> #swap to seperate() for clarity?
  select(-Fish, -V1) |>
  rename(Year = V2, 
         Gear = V3, 
         JulDate = V4, 
         Fish = V5) |>
  mutate(Year = str_sub(Year, 2, 3)) |>
  mutate(Year = as.numeric(ifelse(Year < 85, 
                                  paste0(20, Year), 
                                  paste0(19, Year))), 
         Fish = as.numeric(Fish), 
         JulDate = as.numeric(JulDate)) |>
  relocate(Fish, Year, JulDate, Gear, 
           .before=Comment)

nrow(gsi) - nrow(distinct(gsi, Year, Fish, JulDate)) #how many dups? 0 = good 

#extraction sheet (i.e. to index vial)
extraction <- readxl::read_xlsx(here("data/gsi/Copy of ChinookYukon_Retro_Sept_22_2021_2021-09-22.xlsx"), 
                                sheet = 'Extraction Sheet') |>
  #guess to correct the funny vial names
  mutate(Vial = ifelse(grepl("SR", Vial), 
                       str_sub(Vial, 
                               start = str_locate(Vial, "-")[,1]+1, 
                               end = str_locate(Vial, "_")[,1]-1), 
                       Vial), 
         Year = as.numeric(gsub("Yukon River retrospective |Yukon River Retrospective FW", "", SampleName))) |>
  select(Year, JulDate, Fish, Vial) 

nrow(extraction) - nrow(distinct(extraction, Year, JulDate, Fish)) #how many dups? 0 = good

#asl
asl <- read.csv(here("data/asl/ASL_Output_Chinook_FishWheels_1982-2012.csv")) |>
  filter(Species == "Chinook", Gear == "Fishwheel") |>
  rename(Fish = Fish.Number, 
         Year = Sample.Year) |>
  mutate(JulDate = yday(as.Date(sampleDate)), 
         Year = as.numeric(Year), 
         Fish = as.numeric(Fish))

#merge and write--------------------------------------------------------------------------
gsi_extraction <- left_join(gsi, extraction, by = c("Year", "JulDate", "Fish")) 
nrow(filter(gsi_extraction, is.na(Vial))) #number of rows that didn't get a vial

asl_gsi <- full_join(asl, gsi_extraction, by = c("Year", "JulDate", "Vial"))

nrow(filter(asl_gsi, !is.na(Project) & is.na(`Region 1`)))  #rows with asl only
nrow(filter(asl_gsi, is.na(Project) & !is.na(`Region 1`))) #rows with gsi only
nrow(filter(asl_gsi, !is.na(Project) & !is.na(`Region 1`))) #rows with both

write.csv(asl_gsi, here("data/cleaned/asl-gsi.csv"))

#check the data by year-------------------------------------------------------------------
full_year_check <- filter(asl_gsi, !is.na(Project) & !is.na(`Region 1`)) |>
  group_by(Year) |>
  summarise(n())

asl_year_check <- filter(asl_gsi, !is.na(Project) & is.na(`Region 1`)) |>
  group_by(Year) |>
  summarise(n())

#misc checks------------------------------------------------------------------------------
asl_dup_fish <- asl |>
  group_by(Year, JulDate, Fish, Vial) |>
  summarise(n = n()) |>
  filter(n != 1) #looks like there was some error to make lots of fis/vials have -1 in them