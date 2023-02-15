#wrangling GSI data. Merge historical data with updated

library(tidyverse)
library(here)
library(readxl)

#read in and wrangle tables---------------------------------------------------------------
old_gsi <- read.csv(here("data/gsi/gsiSamplesAllProbs.csv")) |>
  rename(fish = sample_num)

regions <- old_gsi |>
  select(region, region_name)|>
  distinct()

new_gsi <- read_xlsx(here("data/gsi/PID20200058-20210137_EagleSonar_(20-21)_sc270-371_2022-12-20.xlsx"),
                     sheet = "custom_table_ids") |>
  select(-collection, -mixture_collection, -ID_Source, -PBT_brood_year) 

#checkNA <- filter(new_gsi, grepl("NA", indiv)) #weird row of data where Julian == NA

#horribly disgusting rbind way of reshaping data
gsi1 <- select(new_gsi, indiv, group.1, prob.1) |>
  rename(group = group.1, 
         prob = prob.1)
gsi2 <- select(new_gsi, indiv, group.2, prob.2) |>
  rename(group = group.2, 
         prob = prob.2)
gsi3 <- select(new_gsi, indiv, group.3, prob.3) |>
  rename(group = group.3, 
         prob = prob.3)
gsi4 <- select(new_gsi, indiv, group.4, prob.4) |>
  rename(group = group.4, 
         prob = prob.4)
gsi5 <- select(new_gsi, indiv, group.5, prob.5) |>
  rename(group = group.5, 
         prob = prob.5)

clean_new_gsi <- rbind(gsi1, gsi2, gsi3, gsi4, gsi5) |>
  filter(!is.na(prob)) |>
  separate(indiv, into = c(NA, "year", "julian", "fish"), sep="_") |>
  mutate_at(c("year", "julian", "fish"), as.integer) |>
  mutate(data_label = "EagleSonar", 
         gear = "Test Fishery") |>
  rename(region_name = group) |>
  left_join(regions, by = "region_name")

all_gsi <- bind_rows(old_gsi, clean_new_gsi)

#write it---------------------------------------------------------------------------------
write.csv(all_gsi, here("data/cleaned-data/border-gsi.csv"), row.names = FALSE)
