#wrangling ASL data into one array

library(tidyverse)
library(readxl)
library(here)
library(lubridate)

#helper funs-----------------------------------------------------------------------------
#needed since using native pipe (i.e. referencing '.' doesn't work)
my_replace <- function(x){
  return(replace(
    x = x,
    list = is.na(x),
    values = 0))
}

#read in and wrangle tables-------------------------------------------------------------- 
asl <- read.csv(here("data/asl/yukon-chinook-ASL-ADFG-database-acessed-13Feb2023.csv")) |>
  select(Sample.Year, Sex, Fresh.Water.Age, Salt.Water.Age) |>
  filter(!is.na(Fresh.Water.Age) & !is.na(Salt.Water.Age)) |>
  mutate(Age = paste0(Fresh.Water.Age, ".", Salt.Water.Age)) |>
  rename(Year = Sample.Year) |>
  group_by(Year, Sex, Age) |>
  summarise(n = n()) |>
  arrange(as.numeric(Age)) |>
  pivot_wider(names_from = Age, 
              values_from = n) |>
  my_replace()|>
  arrange(Year, Sex) |>
  as.data.frame()

har_age <- read.csv(here("data/harvest/YkCk_Harvest_byDistrTypeStockAge.csv"), 
                    skip = 5) |>
  #WHAT ABOUT WHEN DISTRICT==0?
  mutate(River.Section = case_when(District <= 2 ~ "below_pilot",
                                   District>=3 & District<=6 ~ "above_pilot", 
                                  District == 7 ~ "Canada")) |> 
  arrange(Year, Fishery, River.Section)

colnames(har_age) <- gsub("Age", "", colnames(har_age))

har_age <- har_age |>
  select(Year, Fishery, River.Section, '1.1', '1.2', '1.3', '1.4', '1.5', '1.6', '2.2', 
         '2.3', '2.4', '2.5')

#write tables-----------------------------------------------------------------------------
write.csv(asl, here("data/cleaned-data/age-table.csv"))
write.csv(har_age, here("data/cleaned-data/harvested-age-table.csv"))
