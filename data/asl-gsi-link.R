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
index <- read_tsv(here("data/gsi/ExtractionSheets1982-2012.txt")) |>
  select(Year, JulDate, Fish, Vial) 

#merge and write--------------------------------------------------------------------------
gsi_indexed <- left_join(gsi, index, by = c("Year", "JulDate", "Fish"))

asl_gsi <- full_join(asl, gsi_indexed, by = c("Year", "JulDate", "Fish", "Vial"))

write.csv(asl_gsi, here("data/cleaned/asl-gsi.csv"))






#doing some checks------------------------------------------------------------------------
#what didn't have a match in asl from gsi?
gsi_no_match <- anti_join(gsi_indexed, asl, by = c("Year", "JulDate", "Fish", "Vial"))
#only 12 had matches?!
semi_join(asl, gsi_indexed, by = c("Year", "JulDate", "Fish", "Vial"))

#what's distinct in asl?
asl <- asl |>
  mutate(id = row_number()) #helper col for later

asl_join <- asl |> #what we would join on
  select(Year, JulDate, Fish, Vial) |>
  arrange(Year, JulDate, Fish, Vial)

asl_distinct <- asl |>
  distinct(Year, JulDate, Fish, Vial, .keep_all = TRUE) |>
  arrange(Year, JulDate, Fish, Vial)

nrow(asl_join)-nrow(distinct(asl_join)) #2420 duplicates 

asl_dups <- asl |> #the duplicates
  filter(!(id %in% asl_distinct$id))

asl_dup_fish <- asl |>
  group_by(Year, JulDate, Fish, Vial) |>
  summarise(n = n()) |>
  filter(n != 1) #looks like there was some error to make lots of fis/vials have -1 in them

gsi_distinct <- gsi_indexed |>
  select(Year, JulDate, Fish, Vial) |>
  arrange(Year, JulDate, Fish, Vial)
