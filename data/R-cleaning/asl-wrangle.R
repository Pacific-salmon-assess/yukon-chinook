#wrangling ASL data into 3 matrices

library(tidyverse)
library(here)

#helper funs-----------------------------------------------------------------------------
#needed since using native pipe (i.e. referencing '.' doesn't work)
my_replace <- function(x){
  return(replace(
    x = x,
    list = is.na(x),
    values = 0))
}

#read in and wrangle tables-------------------------------------------------------------- 
as <- read.csv(here("data/asl/yukon-chinook-ASL-ADFG-database-acessed-13Feb2023.csv")) |>
  select(Sample.Year, Gear, Sex, Fresh.Water.Age, Salt.Water.Age) |>
  filter(!is.na(Fresh.Water.Age) & !is.na(Salt.Water.Age)) |>
  mutate(Age = paste0(Fresh.Water.Age, ".", Salt.Water.Age)) |>
  rename(Year = Sample.Year) |>
  group_by(Year, Sex, Gear, Age) |>
  summarise(n = n()) |>
  arrange(as.numeric(Age)) |>
  pivot_wider(names_from = Age, 
              values_from = n) |>
  my_replace()|>
  arrange(Year, Sex) |>
  as.data.frame()

asl <- read.csv(here("data/asl/yukon-chinook-ASL-ADFG-database-acessed-13Feb2023.csv")) |>
  select(Sample.Year, Gear, Sex, Fresh.Water.Age, Salt.Water.Age, Length, Length.Measurement.Type) |>
  filter(!is.na(Fresh.Water.Age) & !is.na(Salt.Water.Age)) |>
  filter(Sex != "unknown") |>
  mutate(Age = paste0(Fresh.Water.Age, ".", Salt.Water.Age)) |>
  rename(Year = Sample.Year) |>
  group_by(Year, Sex, Gear, Age, Length.Measurement.Type) |>
  summarise(Mean.Length = mean(Length)) |>
  arrange(as.numeric(Age)) |>
  pivot_wider(names_from = Age, 
              values_from = Mean.Length) |>
  my_replace()|>
  arrange(Year, Sex) |>
  as.data.frame()

har_age <- read.csv(here("data/harvest/YkCk_Harvest_byDistrTypeStockAge.csv"), 
                    skip = 5) |>
  filter(Stock == "Upper") |>
  mutate(River.Section = case_when(District <= 2 ~ "below_pilot",
                                   District >=3 & District<6 ~ "above_pilot", 
                                  District == 7 ~ "Canada")) |> 
  mutate(Fishery.type = case_when(Fishery == "US_Subsistence" ~ "Subsistence",
                                  Fishery == "US_Comm" ~ "Commercial")) |> 
  arrange(Year, Fishery, River.Section)

colnames(har_age) <- gsub("Age", "", colnames(har_age))

har_age <- har_age |>
  select(Year, Fishery.type, River.Section, '1.1', '1.2', '1.3', '1.4', '1.5', '1.6', '2.2', 
         '2.3', '2.4', '2.5')

#write tables-----------------------------------------------------------------------------
write.csv(as, here("data/cleaned-data/border-age-table.csv"), row.names = FALSE)
write.csv(asl, here("data/cleaned-data/border-length-table.csv"), row.names = FALSE)
write.csv(har_age, here("data/cleaned-data/harvest-age-table.csv"), row.names = FALSE)
