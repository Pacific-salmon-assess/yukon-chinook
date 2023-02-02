# yukon-chinook-data

| Folder | File | Description |
| ------- | -------- | ---------------------------------------------------- |
| harvest | `YkCk_Harvest_byDistrTypeStockAge.csv` | Harvest by age and fate (subsistence and commercial) by river section 1981-2021. River sections 1-5 are in Alaska and section 7 is Canadian mainstem|
| border-assessment | `YkCk_EagleSonar_Data.csv` | Eagle sonar estimate of border passage 2005-2021|
| border-assessment | `YkCk_BorderMR_Data.csv` | Border mark-recapture estimate of border passage 1982-2008|
| border-assessment | `borderCounts.csv` | Daily counts of Chinook border passage at either Bio Island (fish wheels; 1985-2005) or Eagle (sonar; 2005-2022)|
| border-assessment | `aslEagle.csv` | Raw age, sex, length, data from multi-panel gillnet test fishery at Eagle (2005-2019)|
| gsi | `gsiSamplesAllProbs.csv` | Genetic stock assignments to each Conservation Unit by gear (fish wheel, gillnet test fishery), year, and day (1985-2019)|
| gsi | `ChinookYukon_Retro_Sept_22_2021.xlsx` | retrospective Genetic stock assignments to each Conservation Unit by gear (fish wheel, gillnet test fishery), year, and day (1985-2012)|
| gsi | `ExtractionSheets1982-2012.txt` | table used to link the GSI and ASL data |
| asl | `ASL_Output_Chinook_FishWheels_1982-2012.csv` | Age, sex, and length (asl) data from fish wheel 1982-2012) |

Notes:
- ASL data pre 2005 to be added.
- ASL and GSI data need to be wrangled and merged. Requires discussion with DFO Molecular Genetics and sclerochronology labs and ADF&G ASL database folks. GSI data in this repo does not have unique ID that can be linked to origional ASL sample. 
- 2022 Eagle sonar and harvest by distric data will be available in early 2023.
