# Bee sex ratios and climate change
# RMBL bee and plant phenology data project
# Melanie Kazenel
# Created 7 October 2022
# ---------------------------------

# load relevant libraries
library(dplyr)
library(tidyverse)
library(effects)
library(lme4)
library(dplyr)

# set wd
setwd("~/Documents/*RMBLPhenology/Projects/BeeSexRatios/Analyses")

##### Bee data manipulation and formatting #####

# read in bee and climate data
beedata<-read.csv("bees_2022-09-06.csv")

# read in trait data
traits<-read.csv("traits_stemkovski_2019_10_15.csv")
# rename species column
names(traits)[1]<-"genus_species"

# join trait data to beedata
beedatafocal<-left_join(beedata,traits,by=c("genus_species","family"))

# standardize data format for sex column
unique(beedatafocal$sex)
beedatafocal$sex[beedatafocal$sex == "f"] <- "F"
beedatafocal$sex[beedatafocal$sex == "m"] <- "M"

# subset to just include M and F-designated solitary bees
beedatafocal<-subset(beedatafocal, sex=="M" | sex=="F")
unique(beedatafocal$sex)

# check which bees were collected in which year
year_sex_summary<-beedatafocal %>% group_by(genus_species,sex,year) %>% summarise(obs=n())
year_sex_wide<-pivot_wider(year_sex_summary,id_cols = c(1,3),names_from = sex,values_from = obs)
# remove blank rows and Bombus
year_sex_wide<-year_sex_wide[-c(1:10,176:230),]
# replace NAs with Os
year_sex_wide<-replace_na(year_sex_wide, replace=list(F=0,M=0))
# add M and F presence/absence columns
year_sex_wide$f_present <- ifelse(year_sex_wide$F>0, 1, 0)
year_sex_wide$m_present <- ifelse(year_sex_wide$M>0, 1, 0)
year_sex_wide$sex_present_sum<-year_sex_wide$f_present + year_sex_wide$m_present 
year_sex_focal<-filter(year_sex_wide, sex_present_sum>1)
# get counts of years with both males and females present for each species
year_count_species<-year_sex_focal%>%group_by(genus_species)%>%summarise(years_both_sexes=n())
# subset to just include species present in 5, 6, 7 etc. or more years
fiveplus<-filter(year_count_species, years_both_sexes>=5) # 23 species
fiveplus$genus_species
sixplus<-filter(year_count_species, years_both_sexes>=6) # 14 species
sixplus$genus_species
sevenplus<-filter(year_count_species, years_both_sexes>=7) # 10 species
sevenplus$genus_species

# present in 7+ years but not included in original focal species
# Halictus virgatellus
# Hylaeus annulatus

# present in 6+ years but not included in original focal species
# Ceratina nanula
# Osmia albolateralis

# present in 5+ years but not included in original focal species
# Megachile melanophaea
# Melissodes lutulentus
# Osmia albolateralis
# Osmia ednae
# Osmia simillima
# various Lasioglossum species (don't include in analyses because males can't concretely be identified to species)

# check focal species not revealed above
filter(year_count_species, genus_species %in% "Andrena algida") # only 1 year with box sexes present - EXCLUDE
filter(year_count_species, genus_species %in% "Dianthidium heterulkei") # only 2 years with box sexes present - EXCLUDE

# subset to just include focal species
target<-c(# Willow's original focal species
  "Agapostemon texanus", # 6+ years with both sexes
  #"Andrena algida", #EXCLUDE
  # "Ceratina neomexicana",	# 5+ years with both sexes EXCLUDE
  #"Dianthidium heterulkei",	# EXCLUDE
  "Dufourea harveyi",	# 7+ years with both sexes
  "Halictus rubicundus", # 7+ years with both sexes
  "Hoplitis fulgida",	# 7+ years with both sexes
  "Hoplitis robusta",	# 7+ years with both sexes
  "Panurginus cressoniellus",	# 7+ years with both sexes
  "Panurginus ineptus",	# 7+ years with both sexes
  "Pseudopanurgus bakeri", # 7+ years with both sexes
  "Pseudopanurgus didirupa", # 7+ years with both sexes
  
  # Additional focus species, present in 6 or more years of the dataset 
  #(2/3 or more of years)
  "Halictus virgatellus",
  "Hylaeus annulatus",
  "Ceratina nanula",
  "Osmia albolateralis"
)

beedatafocal<-filter(beedatafocal, genus_species %in% target)
unique(beedatafocal$genus_species)

# recode female/male as 1/0
beedatafocal$sex_binom <- ifelse(beedatafocal$sex == "F", 1, 0)

##### Add climate data #####
snow<-read.csv("billy_barr_bare_ground.csv")
