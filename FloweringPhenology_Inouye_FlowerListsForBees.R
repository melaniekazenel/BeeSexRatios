# Working with Inouye flowering phenology data for bee sex ratios project
# Flower lists for individual bee species
# Melanie Kazenel
# Created 2022-10-13
# -------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)

setwd("~/Documents/*RMBLPhenology/Projects/BeeSexRatios/Analyses")

# Initial data formatting ##### 

# load data
load(file='all-data-2022.RData')

# create a copy with a new name to work with
flor<-all.data.2022

# look at data and add family information
species_list<-flor %>% group_by(species) %>% summarise(count=n())
genus_list<-read.csv("Inouye_plant_species_list_2020.csv")
genus_list<-genus_list[,-c(3:4)]
species_list<-left_join(species_list,genus_list,by="species")
#write.csv(species_list,"Inouye_plant_species_list_2022.csv")
# family IDs added manually for each plant

# read in updated list
species_list<-read.csv("Inouye_plant_species_list_2022.csv")
species_list$count<-NULL
flor<-left_join(flor,species_list,by="species")

# subset to remove grasses and sedges, and any plants outside of plots
flor2<-filter(flor,family!="Poaceae" & family!="Cyperaceae")

# get plot-year summary to check that all plots were sampled in all years
plot_year<-flor2 %>% group_by(year,plot) %>% summarise(count=n())
plot_year2<-plot_year %>% group_by(year) %>% summarise(count=n()) %>% filter(year>=2008)

# check out which plot was missing in 2012
plot_year_12_13<-plot_year %>% filter(year==2012 | year==2013)
# MDW2 present in 2013 but not 2012
plot_year_12<-plot_year %>% filter(year==2012)

# subset to just include the plots sampled in each year from 2009-2021
flor_focal<-flor2 %>% filter(year>=2008 & year<= 2021) %>% filter(plot %in% plot_year_12$plot)
# check to make sure that all plots were sampled in all years
summary<-flor_focal %>% group_by(year,plot) %>% summarise(count=n())
summary$count<-1
summary2<-pivot_wider(summary,names_from = plot, values_from = count)
# looks good!


# Plant lists for focal bee species #####

# read in Jane's data
focalplants<-read.csv("RMBL_solitary_bee_data_2020_forMelanie.csv")

# rename columns
focalplants$pollinator <- sub("_", " ", focalplants$pollinator)
focalplants$plant <- sub("_", " ", focalplants$plant)

# subset to just include focal bees
focalplants<-focalplants %>% filter(pollinator=="Halictus rubicundus" | 
                                      pollinator=="Halictus virgatellus" |
                                      pollinator=="Hoplitis fulgida" |
                                      pollinator=="Panurginus ineptus")
unique(focalplants$pollinator)

# read in Becky's bee data
focalbees<-read.csv("BeeDataFocalSpecies_2023-03-30.csv")

# filter out non-plant visit records
unique(focalbees$details)
focalbees<-filter(focalbees, details!="in flight" & details!="In flight" & 
                    details!="Blue" & details!="Yellow" & details!="White" &
                    details!="unknown" & details!="Unknown" & 
                    details!="resting" & details!="Resting")
unique(focalbees$details)



# Halictus rubicundus #####

# check how datasets match up
jane<-focalplants %>% filter(pollinator=="Halictus rubicundus") %>%  group_by(pollinator,plant) %>% summarise(total.visits=sum(total.visits)) %>% mutate(dataset="jane")

becky<-focalbees %>% filter(genus_species == "Halictus rubicundus") %>% group_by(details) %>% summarise(total.visits=n()) %>% rename(plant=details) %>% mutate(pollinator="Halictus rubicundus") %>% select(pollinator,plant,total.visits) %>% mutate(dataset="becky")

# rename plants if needed
# Amerosedum lanceolatum to Sedum lanceolatum
becky$plant[becky$plant == "Amerosedum lanceolatum"] <- "Sedum lanceolatum"
# Chrysothamnus viscidiflorus to Chrysothamnus humilis
becky$plant[becky$plant == "Chrysothamnus viscidiflorus"] <- "Chrysothamnus humilis"
# Dugaldia hoopesii to Hymenoxys hoopesii
becky$plant[becky$plant == "Dugaldia hoopesii"] <- "Hymenoxys hoopesii"
# Oreocarya flavoculata to Cryptantha flavoculata
becky$plant[becky$plant == "Oreocarya flavoculata"] <- "Cryptantha flavoculata"
# Pediocactus simpsonii to Pediocactus nigrispinus
becky$plant[becky$plant == "Pediocactus simpsonii"] <- "Pediocactus nigrispinus"
# Pentaphylloides floribunda to Dasiphora fruticosa
becky$plant[becky$plant == "Pentaphylloides floribunda"] <- "Dasiphora fruticosa"

# Amerosedum lanceolatum to Sedum lanceolatum
jane$plant[jane$plant == "Amerosedum lanceolatum"] <- "Sedum lanceolatum"
# Chrysothamnus viscidiflorus to Chrysothamnus humilis
jane$plant[jane$plant == "Chrysothamnus viscidiflorus"] <- "Chrysothamnus humilis"
# Dugaldia hoopesii to Hymenoxys hoopesii
jane$plant[jane$plant == "Dugaldia hoopesii"] <- "Hymenoxys hoopesii"
# Oreocarya flavoculata to Cryptantha flavoculata
jane$plant[jane$plant == "Oreocarya flavoculata"] <- "Cryptantha flavoculata"
# Pediocactus simpsonii to Pediocactus nigrispinus
jane$plant[jane$plant == "Pediocactus simpsonii"] <- "Pediocactus nigrispinus"
# Pentaphylloides floribunda to Dasiphora fruticosa
jane$plant[jane$plant == "Pentaphylloides floribunda"] <- "Dasiphora fruticosa"

matches<-filter(jane, plant %in% becky$plant)
# jane: 23 plants
# becky: 30 plants
# shared: 13


# combine data frames
halrub<-focalbees %>% filter(genus_species == "Halictus rubicundus") %>% group_by(details) %>% summarise(total.visits=n()) %>% rename(plant=details) %>% mutate(pollinator="Halictus rubicundus") %>% select(pollinator,plant,total.visits)

halrub<-bind_rows(halrub,
                  filter(focalplants,pollinator=="Halictus rubicundus")[,c(1,3,10)]) %>% 
  group_by(pollinator,plant) %>% summarise(total.visits=sum(total.visits))

names(halrub)[2]<-"species"

# check that plant names are in species list
nomatch<-anti_join(halrub,species_list,by="species")

# rename plants in Irwin dataset
# Amerosedum lanceolatum to Sedum lanceolatum
halrub$species[halrub$species == "Amerosedum lanceolatum"] <- "Sedum lanceolatum"
# Chrysothamnus viscidiflorus to Chrysothamnus humilis
halrub$species[halrub$species == "Chrysothamnus viscidiflorus"] <- "Chrysothamnus humilis"
# Dugaldia hoopesii to Hymenoxys hoopesii
halrub$species[halrub$species == "Dugaldia hoopesii"] <- "Hymenoxys hoopesii"
# Oreocarya flavoculata to Cryptantha flavoculata
halrub$species[halrub$species == "Oreocarya flavoculata"] <- "Cryptantha flavoculata"
# Pediocactus simpsonii to Pediocactus nigrispinus
halrub$species[halrub$species == "Pediocactus simpsonii"] <- "Pediocactus nigrispinus"
# Pentaphylloides floribunda to Dasiphora fruticosa
halrub$species[halrub$species == "Pentaphylloides floribunda"] <- "Dasiphora fruticosa"

flor_focal<-inner_join(flor_focal,halrub[,c(1:2)],by="species")

# calculate floral sum and floral days for each year

# annual floral sum
# from Ogilvie et al. 2017:
# pooled flower counts across plant species and plots on every sampling date
# calculated the sum of flowers from first flower until 80% of flowers had accumulated

# calculate floral sum for each day
daily<-flor_focal %>% group_by(year,date,doy) %>% summarise(daily_floral_count=sum(floralcount))

# # check year with high floral counts (2015)
daily15<-filter(daily,year==2015)
ggplot(aes(x=doy,y=daily_floral_count),data=daily15) + geom_point()
# trend looks real!

# calculate floral sum for each year (across the whole sampling season)
flor_annual<-daily %>% group_by(year) %>% summarise(total_yearly_flowers=sum(daily_floral_count))

# for each year, calculate when 80% of flowers have accumulated
daily2<-transform(daily, accum_flowers = ave(daily_floral_count, year, FUN = cumsum))
daily2<-transform(daily2, total_flowers = ave(daily_floral_count, year, FUN = sum))
daily2$percent_flowers<-daily2$accum_flowers/daily2$total_flowers

daily3<-filter(daily2,percent_flowers>=0.8)
annual_80<-daily3 %>% group_by(year) %>% summarise(percent_flowers=min(percent_flowers))
annual_80<-left_join(annual_80,daily3,by=c("year","percent_flowers"))

flor_annual<-left_join(flor_annual,annual_80[,c(1,6)],by="year")

# annual floral days
# from Ogilvie et al. 2017:
# number of days above a low flower threshold (0.75 flowers/m2 or 3 flowers per 2 x 2 m plot) between the first flower date and the date on which 80% of the season’s flowers had accumulated

annual_80_2<-annual_80
names(annual_80_2)[2]<-"percent_flowers_80_day"
floral_days<-left_join(daily2,annual_80_2[,c(1:2)],by="year")
floral_days<-filter(floral_days,percent_flowers<=percent_flowers_80_day)
floral_days$flowers_per_plot<-floral_days$daily_floral_count/length(plot_year_12$plot)
floral_days_high<-filter(floral_days,flowers_per_plot>3)
total_days_high<-floral_days_high %>% group_by(year) %>% summarise(floral_days_80=n())


floral_days_high_all<-filter(daily,daily_floral_count>3)
total_days_high_all<-floral_days_high_all %>% group_by(year) %>% summarise(floral_days_all=n())

# join to flor_annual data frame
flor_annual<-left_join(flor_annual,total_days_high_all,by="year")
flor_annual<-left_join(flor_annual,total_days_high,by="year")

# rename columns
names(flor_annual)[2:5]<-c("total_flowers","total_flowers_80","floral_days","floral_days_80")

#write.csv(flor_annual,"floral_data_annual_summaries_forsexratios_halictus_rubicundus.csv",row.names = FALSE)

# Halictus virgatellus #####
halvir<-focalbees %>% filter(genus_species == "Halictus virgatellus") %>% group_by(details) %>% summarise(total.visits=n()) %>% rename(plant=details) %>% mutate(pollinator="Halictus virgatellus") %>% select(pollinator,plant,total.visits)

halvir<-bind_rows(halvir,
                  filter(focalplants,pollinator=="Halictus virgatellus")[,c(1,3,10)]) %>% 
  group_by(pollinator,plant) %>% summarise(total.visits=sum(total.visits))

names(halvir)[2]<-"species"

# check that plant names are in species list
nomatch<-anti_join(halvir,species_list,by="species")

# rename plants in Irwin dataset
# Amerosedum lanceolatum to Sedum lanceolatum
halvir$species[halvir$species == "Amerosedum lanceolatum"] <- "Sedum lanceolatum"
# Chrysothamnus viscidiflorus to Chrysothamnus humilis
halvir$species[halvir$species == "Chrysothamnus viscidiflorus"] <- "Chrysothamnus humilis"
# Dugaldia hoopesii to Hymenoxys hoopesii
halvir$species[halvir$species == "Dugaldia hoopesii"] <- "Hymenoxys hoopesii"
# Oreocarya flavoculata to Cryptantha flavoculata
halvir$species[halvir$species == "Oreocarya flavoculata"] <- "Cryptantha flavoculata"
# Pediocactus simpsonii to Pediocactus nigrispinus
halvir$species[halvir$species == "Pediocactus simpsonii"] <- "Pediocactus nigrispinus"
# Pentaphylloides floribunda to Dasiphora fruticosa
halvir$species[halvir$species == "Pentaphylloides floribunda"] <- "Dasiphora fruticosa"
# Adenolinum lewisii to Linum lewisii
halvir$species[halvir$species == "Adenolinum lewisii"] <- "Linum lewisii"
# Ligularia bigelovii to Senecio bigelovii 
halvir$species[halvir$species == "Ligularia bigelovii"] <- "Senecio bigelovii"

flor_focal<-inner_join(flor_focal,halvir[,c(1:2)],by="species")

# calculate floral sum and floral days for each year

# annual floral sum
# from Ogilvie et al. 2017:
# pooled flower counts across plant species and plots on every sampling date
# calculated the sum of flowers from first flower until 80% of flowers had accumulated

# calculate floral sum for each day
daily<-flor_focal %>% group_by(year,date,doy) %>% summarise(daily_floral_count=sum(floralcount))

# # check year with high floral counts (2015)
daily15<-filter(daily,year==2015)
ggplot(aes(x=doy,y=daily_floral_count),data=daily15) + geom_point()
# trend looks real!

# calculate floral sum for each year (across the whole sampling season)
flor_annual<-daily %>% group_by(year) %>% summarise(total_yearly_flowers=sum(daily_floral_count))

# for each year, calculate when 80% of flowers have accumulated
daily2<-transform(daily, accum_flowers = ave(daily_floral_count, year, FUN = cumsum))
daily2<-transform(daily2, total_flowers = ave(daily_floral_count, year, FUN = sum))
daily2$percent_flowers<-daily2$accum_flowers/daily2$total_flowers

daily3<-filter(daily2,percent_flowers>=0.8)
annual_80<-daily3 %>% group_by(year) %>% summarise(percent_flowers=min(percent_flowers))
annual_80<-left_join(annual_80,daily3,by=c("year","percent_flowers"))

flor_annual<-left_join(flor_annual,annual_80[,c(1,6)],by="year")

# annual floral days
# from Ogilvie et al. 2017:
# number of days above a low flower threshold (0.75 flowers/m2 or 3 flowers per 2 x 2 m plot) between the first flower date and the date on which 80% of the season’s flowers had accumulated

annual_80_2<-annual_80
names(annual_80_2)[2]<-"percent_flowers_80_day"
floral_days<-left_join(daily2,annual_80_2[,c(1:2)],by="year")
floral_days<-filter(floral_days,percent_flowers<=percent_flowers_80_day)
floral_days$flowers_per_plot<-floral_days$daily_floral_count/length(plot_year_12$plot)
floral_days_high<-filter(floral_days,flowers_per_plot>3)
total_days_high<-floral_days_high %>% group_by(year) %>% summarise(floral_days_80=n())


floral_days_high_all<-filter(daily,daily_floral_count>3)
total_days_high_all<-floral_days_high_all %>% group_by(year) %>% summarise(floral_days_all=n())

# join to flor_annual data frame
flor_annual<-left_join(flor_annual,total_days_high_all,by="year")
flor_annual<-left_join(flor_annual,total_days_high,by="year")

# rename columns
names(flor_annual)[2:5]<-c("total_flowers","total_flowers_80","floral_days","floral_days_80")

#write.csv(flor_annual,"floral_data_annual_summaries_forsexratios_halictus_virgatellus.csv",row.names = FALSE)

# Hoplitis fulgida #####
hopful<-focalbees %>% filter(genus_species == "Hoplitis fulgida") %>% group_by(details) %>% summarise(total.visits=n()) %>% rename(plant=details) %>% mutate(pollinator="Hoplitis fulgida") %>% select(pollinator,plant,total.visits)

hopful<-bind_rows(hopful,
                  filter(focalplants,pollinator=="Hoplitis fulgida")[,c(1,3,10)]) %>% 
  group_by(pollinator,plant) %>% summarise(total.visits=sum(total.visits))

names(hopful)[2]<-"species"

# check that plant names are in species list
nomatch<-anti_join(hopful,species_list,by="species")

# rename plants in Irwin dataset
# Amerosedum lanceolatum to Sedum lanceolatum
hopful$species[hopful$species == "Amerosedum lanceolatum"] <- "Sedum lanceolatum"
# Chrysothamnus viscidiflorus to Chrysothamnus humilis
hopful$species[hopful$species == "Chrysothamnus viscidiflorus"] <- "Chrysothamnus humilis"
# Dugaldia hoopesii to Hymenoxys hoopesii
hopful$species[hopful$species == "Dugaldia hoopesii"] <- "Hymenoxys hoopesii"
# Oreocarya flavoculata to Cryptantha flavoculata
hopful$species[hopful$species == "Oreocarya flavoculata"] <- "Cryptantha flavoculata"
# Pediocactus simpsonii to Pediocactus nigrispinus
hopful$species[hopful$species == "Pediocactus simpsonii"] <- "Pediocactus nigrispinus"
# Pentaphylloides floribunda to Dasiphora fruticosa
hopful$species[hopful$species == "Pentaphylloides floribunda"] <- "Dasiphora fruticosa"
# Adenolinum lewisii to Linum lewisii
hopful$species[hopful$species == "Adenolinum lewisii"] <- "Linum lewisii"
# Ligularia bigelovii to Senecio bigelovii 
hopful$species[hopful$species == "Ligularia bigelovii"] <- "Senecio bigelovii"
# Chamerion danielsii to Chamerion angustifolium
hopful$species[hopful$species == "Chamerion danielsii"] <- "Chamerion angustifolium"
# Erigeron glacialis to Erigeron peregrinus
hopful$species[hopful$species == "Erigeron glacialis"] <- "Erigeron peregrinus"


nomatch<-anti_join(hopful,species_list,by="species")


flor_focal<-inner_join(flor_focal,hopful[,c(1:2)],by="species")

# calculate floral sum and floral days for each year

# annual floral sum
# from Ogilvie et al. 2017:
# pooled flower counts across plant species and plots on every sampling date
# calculated the sum of flowers from first flower until 80% of flowers had accumulated

# calculate floral sum for each day
daily<-flor_focal %>% group_by(year,date,doy) %>% summarise(daily_floral_count=sum(floralcount))

# # check year with high floral counts (2015)
daily15<-filter(daily,year==2015)
ggplot(aes(x=doy,y=daily_floral_count),data=daily15) + geom_point()
# trend looks real!

# calculate floral sum for each year (across the whole sampling season)
flor_annual<-daily %>% group_by(year) %>% summarise(total_yearly_flowers=sum(daily_floral_count))

# for each year, calculate when 80% of flowers have accumulated
daily2<-transform(daily, accum_flowers = ave(daily_floral_count, year, FUN = cumsum))
daily2<-transform(daily2, total_flowers = ave(daily_floral_count, year, FUN = sum))
daily2$percent_flowers<-daily2$accum_flowers/daily2$total_flowers

daily3<-filter(daily2,percent_flowers>=0.8)
annual_80<-daily3 %>% group_by(year) %>% summarise(percent_flowers=min(percent_flowers))
annual_80<-left_join(annual_80,daily3,by=c("year","percent_flowers"))

flor_annual<-left_join(flor_annual,annual_80[,c(1,6)],by="year")

# annual floral days
# from Ogilvie et al. 2017:
# number of days above a low flower threshold (0.75 flowers/m2 or 3 flowers per 2 x 2 m plot) between the first flower date and the date on which 80% of the season’s flowers had accumulated

annual_80_2<-annual_80
names(annual_80_2)[2]<-"percent_flowers_80_day"
floral_days<-left_join(daily2,annual_80_2[,c(1:2)],by="year")
floral_days<-filter(floral_days,percent_flowers<=percent_flowers_80_day)
floral_days$flowers_per_plot<-floral_days$daily_floral_count/length(plot_year_12$plot)
floral_days_high<-filter(floral_days,flowers_per_plot>3)
total_days_high<-floral_days_high %>% group_by(year) %>% summarise(floral_days_80=n())


floral_days_high_all<-filter(daily,daily_floral_count>3)
total_days_high_all<-floral_days_high_all %>% group_by(year) %>% summarise(floral_days_all=n())

# join to flor_annual data frame
flor_annual<-left_join(flor_annual,total_days_high_all,by="year")
flor_annual<-left_join(flor_annual,total_days_high,by="year")

# rename columns
names(flor_annual)[2:5]<-c("total_flowers","total_flowers_80","floral_days","floral_days_80")

#write.csv(flor_annual,"floral_data_annual_summaries_forsexratios_hoplitis_fulgida.csv",row.names = FALSE)

# Panurginus ineptus #####

# check how datasets match up
jane<-focalplants %>% filter(pollinator=="Panurginus ineptus") %>%  group_by(pollinator,plant) %>% summarise(total.visits=sum(total.visits)) %>% mutate(dataset="jane")

becky<-focalbees %>% filter(genus_species == "Panurginus ineptus") %>% group_by(details) %>% summarise(total.visits=n()) %>% rename(plant=details) %>% mutate(pollinator="Panurginus ineptus") %>% select(pollinator,plant,total.visits) %>% mutate(dataset="becky")

# rename plants if needed
# Amerosedum lanceolatum to Sedum lanceolatum
becky$plant[becky$plant == "Amerosedum lanceolatum"] <- "Sedum lanceolatum"
# Chrysothamnus viscidiflorus to Chrysothamnus humilis
becky$plant[becky$plant == "Chrysothamnus viscidiflorus"] <- "Chrysothamnus humilis"
# Dugaldia hoopesii to Hymenoxys hoopesii
becky$plant[becky$plant == "Dugaldia hoopesii"] <- "Hymenoxys hoopesii"
# Oreocarya flavoculata to Cryptantha flavoculata
becky$plant[becky$plant == "Oreocarya flavoculata"] <- "Cryptantha flavoculata"
# Pediocactus simpsonii to Pediocactus nigrispinus
becky$plant[becky$plant == "Pediocactus simpsonii"] <- "Pediocactus nigrispinus"
# Pentaphylloides floribunda to Dasiphora fruticosa
becky$plant[becky$plant == "Pentaphylloides floribunda"] <- "Dasiphora fruticosa"
# Adenolinum lewisii to Linum lewisii
becky$plant[becky$plant == "Adenolinum lewisii"] <- "Linum lewisii"
# Ligularia bigelovii to Senecio bigelovii 
becky$plant[becky$plant == "Ligularia bigelovii"] <- "Senecio bigelovii"
# Chamerion danielsii to Chamerion angustifolium
becky$plant[becky$plant == "Chamerion danielsii"] <- "Chamerion angustifolium"
# Erigeron glacialis to Erigeron peregrinus
becky$plant[becky$plant == "Erigeron glacialis"] <- "Erigeron peregrinus"
# Galium septentrionale to Galium boreale
becky$plant[becky$plant == "Galium septentrionale"] <- "Galium boreale"
# Pneumonanthe parryi to Gentiana parryi
becky$plant[becky$plant == "Pneumonanthe parryi"] <- "Gentiana parryi"
# Potentilla gracilis to Potentilla pulcherrima
becky$plant[becky$plant == "Potentilla gracilis"] <- "Potentilla pulcherrima"

# Amerosedum lanceolatum to Sedum lanceolatum
jane$plant[jane$plant == "Amerosedum lanceolatum"] <- "Sedum lanceolatum"
# Chrysothamnus viscidiflorus to Chrysothamnus humilis
jane$plant[jane$plant == "Chrysothamnus viscidiflorus"] <- "Chrysothamnus humilis"
# Dugaldia hoopesii to Hymenoxys hoopesii
jane$plant[jane$plant == "Dugaldia hoopesii"] <- "Hymenoxys hoopesii"
# Oreocarya flavoculata to Cryptantha flavoculata
jane$plant[jane$plant == "Oreocarya flavoculata"] <- "Cryptantha flavoculata"
# Pediocactus simpsonii to Pediocactus nigrispinus
jane$plant[jane$plant == "Pediocactus simpsonii"] <- "Pediocactus nigrispinus"
# Pentaphylloides floribunda to Dasiphora fruticosa
jane$plant[jane$plant == "Pentaphylloides floribunda"] <- "Dasiphora fruticosa"
# Adenolinum lewisii to Linum lewisii
jane$plant[jane$plant == "Adenolinum lewisii"] <- "Linum lewisii"
# Ligularia bigelovii to Senecio bigelovii 
jane$plant[jane$plant == "Ligularia bigelovii"] <- "Senecio bigelovii"
# Chamerion danielsii to Chamerion angustifolium
jane$plant[jane$plant == "Chamerion danielsii"] <- "Chamerion angustifolium"
# Erigeron glacialis to Erigeron peregrinus
jane$plant[jane$plant == "Erigeron glacialis"] <- "Erigeron peregrinus"
# Galium septentrionale to Galium boreale
jane$plant[jane$plant == "Galium septentrionale"] <- "Galium boreale"
# Pneumonanthe parryi to Gentiana parryi
jane$plant[jane$plant == "Pneumonanthe parryi"] <- "Gentiana parryi"
# Potentilla gracilis to Potentilla pulcherrima
jane$plant[jane$plant == "Potentilla gracilis"] <- "Potentilla pulcherrima"

matches<-filter(jane, plant %in% becky$plant)
# jane: 19 plants
# becky: 11 plants
# shared: 6

# combine datasets
panin<-focalbees %>% filter(genus_species == "Panurginus ineptus") %>% group_by(details) %>% summarise(total.visits=n()) %>% rename(plant=details) %>% mutate(pollinator="Panurginus ineptus") %>% select(pollinator,plant,total.visits)

panin<-bind_rows(panin,
                  filter(focalplants,pollinator=="Panurginus ineptus")[,c(1,3,10)]) %>% 
  group_by(pollinator,plant) %>% summarise(total.visits=sum(total.visits))

names(panin)[2]<-"species"

# check that plant names are in species list
nomatch<-anti_join(panin,species_list,by="species")

# rename plants in Irwin dataset
# Amerosedum lanceolatum to Sedum lanceolatum
panin$species[panin$species == "Amerosedum lanceolatum"] <- "Sedum lanceolatum"
# Chrysothamnus viscidiflorus to Chrysothamnus humilis
panin$species[panin$species == "Chrysothamnus viscidiflorus"] <- "Chrysothamnus humilis"
# Dugaldia hoopesii to Hymenoxys hoopesii
panin$species[panin$species == "Dugaldia hoopesii"] <- "Hymenoxys hoopesii"
# Oreocarya flavoculata to Cryptantha flavoculata
panin$species[panin$species == "Oreocarya flavoculata"] <- "Cryptantha flavoculata"
# Pediocactus simpsonii to Pediocactus nigrispinus
panin$species[panin$species == "Pediocactus simpsonii"] <- "Pediocactus nigrispinus"
# Pentaphylloides floribunda to Dasiphora fruticosa
panin$species[panin$species == "Pentaphylloides floribunda"] <- "Dasiphora fruticosa"
# Adenolinum lewisii to Linum lewisii
panin$species[panin$species == "Adenolinum lewisii"] <- "Linum lewisii"
# Ligularia bigelovii to Senecio bigelovii 
panin$species[panin$species == "Ligularia bigelovii"] <- "Senecio bigelovii"
# Chamerion danielsii to Chamerion angustifolium
panin$species[panin$species == "Chamerion danielsii"] <- "Chamerion angustifolium"
# Erigeron glacialis to Erigeron peregrinus
panin$species[panin$species == "Erigeron glacialis"] <- "Erigeron peregrinus"
# Galium septentrionale to Galium boreale
panin$species[panin$species == "Galium septentrionale"] <- "Galium boreale"
# Pneumonanthe parryi to Gentiana parryi
panin$species[panin$species == "Pneumonanthe parryi"] <- "Gentiana parryi"
# Potentilla gracilis to Potentilla pulcherrima
panin$species[panin$species == "Potentilla gracilis"] <- "Potentilla pulcherrima"

nomatch<-anti_join(panin,species_list,by="species")

flor_focal<-inner_join(flor_focal,panin[,c(1:2)],by="species")

# calculate floral sum and floral days for each year

# annual floral sum
# from Ogilvie et al. 2017:
# pooled flower counts across plant species and plots on every sampling date
# calculated the sum of flowers from first flower until 80% of flowers had accumulated

# calculate floral sum for each day
daily<-flor_focal %>% group_by(year,date,doy) %>% summarise(daily_floral_count=sum(floralcount))

# # check year with high floral counts (2015)
daily15<-filter(daily,year==2015)
ggplot(aes(x=doy,y=daily_floral_count),data=daily15) + geom_point()
# trend looks real!

# calculate floral sum for each year (across the whole sampling season)
flor_annual<-daily %>% group_by(year) %>% summarise(total_yearly_flowers=sum(daily_floral_count))

# for each year, calculate when 80% of flowers have accumulated
daily2<-transform(daily, accum_flowers = ave(daily_floral_count, year, FUN = cumsum))
daily2<-transform(daily2, total_flowers = ave(daily_floral_count, year, FUN = sum))
daily2$percent_flowers<-daily2$accum_flowers/daily2$total_flowers

daily3<-filter(daily2,percent_flowers>=0.8)
annual_80<-daily3 %>% group_by(year) %>% summarise(percent_flowers=min(percent_flowers))
annual_80<-left_join(annual_80,daily3,by=c("year","percent_flowers"))

flor_annual<-left_join(flor_annual,annual_80[,c(1,6)],by="year")

# annual floral days
# from Ogilvie et al. 2017:
# number of days above a low flower threshold (0.75 flowers/m2 or 3 flowers per 2 x 2 m plot) between the first flower date and the date on which 80% of the season’s flowers had accumulated

annual_80_2<-annual_80
names(annual_80_2)[2]<-"percent_flowers_80_day"
floral_days<-left_join(daily2,annual_80_2[,c(1:2)],by="year")
floral_days<-filter(floral_days,percent_flowers<=percent_flowers_80_day)
floral_days$flowers_per_plot<-floral_days$daily_floral_count/length(plot_year_12$plot)
floral_days_high<-filter(floral_days,flowers_per_plot>3)
total_days_high<-floral_days_high %>% group_by(year) %>% summarise(floral_days_80=n())


floral_days_high_all<-filter(daily,daily_floral_count>3)
total_days_high_all<-floral_days_high_all %>% group_by(year) %>% summarise(floral_days_all=n())

# join to flor_annual data frame
flor_annual<-left_join(flor_annual,total_days_high_all,by="year")
flor_annual<-left_join(flor_annual,total_days_high,by="year")

# rename columns
names(flor_annual)[2:5]<-c("total_flowers","total_flowers_80","floral_days","floral_days_80")

#write.csv(flor_annual,"floral_data_annual_summaries_forsexratios_panurginus_ineptus.csv",row.names = FALSE)


# Panurginus cressoniellus #####
pancres<-focalbees %>% filter(genus_species == "Panurginus cressoniellus") %>% group_by(details) %>% summarise(total.visits=n()) %>% rename(plant=details) %>% mutate(pollinator="Panurginus cressoniellus") %>% select(pollinator,plant,total.visits)

pancres<-bind_rows(pancres,
                 filter(focalplants,pollinator=="Panurginus cressoniellus")[,c(1,3,10)]) %>% 
  group_by(pollinator,plant) %>% summarise(total.visits=sum(total.visits))

names(pancres)[2]<-"species"

# check that plant names are in species list
nomatch<-anti_join(pancres,species_list,by="species")

# rename plants in Irwin dataset
# Amerosedum lanceolatum to Sedum lanceolatum
pancres$species[pancres$species == "Amerosedum lanceolatum"] <- "Sedum lanceolatum"
# Chrysothamnus viscidiflorus to Chrysothamnus humilis
pancres$species[pancres$species == "Chrysothamnus viscidiflorus"] <- "Chrysothamnus humilis"
# Dugaldia hoopesii to Hymenoxys hoopesii
pancres$species[pancres$species == "Dugaldia hoopesii"] <- "Hymenoxys hoopesii"
# Oreocarya flavoculata to Cryptantha flavoculata
pancres$species[pancres$species == "Oreocarya flavoculata"] <- "Cryptantha flavoculata"
# Pediocactus simpsonii to Pediocactus nigrispinus
pancres$species[pancres$species == "Pediocactus simpsonii"] <- "Pediocactus nigrispinus"
# Pentaphylloides floribunda to Dasiphora fruticosa
pancres$species[pancres$species == "Pentaphylloides floribunda"] <- "Dasiphora fruticosa"
# Adenolinum lewisii to Linum lewisii
pancres$species[pancres$species == "Adenolinum lewisii"] <- "Linum lewisii"
# Ligularia bigelovii to Senecio bigelovii 
pancres$species[pancres$species == "Ligularia bigelovii"] <- "Senecio bigelovii"
# Chamerion danielsii to Chamerion angustifolium
pancres$species[pancres$species == "Chamerion danielsii"] <- "Chamerion angustifolium"
# Erigeron glacialis to Erigeron peregrinus
pancres$species[pancres$species == "Erigeron glacialis"] <- "Erigeron peregrinus"
# Galium septentrionale to Galium boreale
pancres$species[pancres$species == "Galium septentrionale"] <- "Galium boreale"
# Pneumonanthe parryi to Gentiana parryi
pancres$species[pancres$species == "Pneumonanthe parryi"] <- "Gentiana parryi"
# Potentilla gracilis to Potentilla pulcherrima
pancres$species[pancres$species == "Potentilla gracilis"] <- "Potentilla pulcherrima"
# Eriogonum umbellatum var. aureum to Eriogonum umbellatum
pancres$species[pancres$species == "Eriogonum umbellatum var. aureum"] <- "Eriogonum umbellatum"

nomatch<-anti_join(pancres,species_list,by="species")

flor_focal<-inner_join(flor_focal,pancres[,c(1:2)],by="species")

# calculate floral sum and floral days for each year

# annual floral sum
# from Ogilvie et al. 2017:
# pooled flower counts across plant species and plots on every sampling date
# calculated the sum of flowers from first flower until 80% of flowers had accumulated

# calculate floral sum for each day
daily<-flor_focal %>% group_by(year,date,doy) %>% summarise(daily_floral_count=sum(floralcount))

# # check year with high floral counts (2015)
daily15<-filter(daily,year==2015)
ggplot(aes(x=doy,y=daily_floral_count),data=daily15) + geom_point()
# trend looks real!

# calculate floral sum for each year (across the whole sampling season)
flor_annual<-daily %>% group_by(year) %>% summarise(total_yearly_flowers=sum(daily_floral_count))

# for each year, calculate when 80% of flowers have accumulated
daily2<-transform(daily, accum_flowers = ave(daily_floral_count, year, FUN = cumsum))
daily2<-transform(daily2, total_flowers = ave(daily_floral_count, year, FUN = sum))
daily2$percent_flowers<-daily2$accum_flowers/daily2$total_flowers

daily3<-filter(daily2,percent_flowers>=0.8)
annual_80<-daily3 %>% group_by(year) %>% summarise(percent_flowers=min(percent_flowers))
annual_80<-left_join(annual_80,daily3,by=c("year","percent_flowers"))

flor_annual<-left_join(flor_annual,annual_80[,c(1,6)],by="year")

# annual floral days
# from Ogilvie et al. 2017:
# number of days above a low flower threshold (0.75 flowers/m2 or 3 flowers per 2 x 2 m plot) between the first flower date and the date on which 80% of the season’s flowers had accumulated

annual_80_2<-annual_80
names(annual_80_2)[2]<-"percent_flowers_80_day"
floral_days<-left_join(daily2,annual_80_2[,c(1:2)],by="year")
floral_days<-filter(floral_days,percent_flowers<=percent_flowers_80_day)
floral_days$flowers_per_plot<-floral_days$daily_floral_count/length(plot_year_12$plot)
floral_days_high<-filter(floral_days,flowers_per_plot>3)
total_days_high<-floral_days_high %>% group_by(year) %>% summarise(floral_days_80=n())


floral_days_high_all<-filter(daily,daily_floral_count>3)
total_days_high_all<-floral_days_high_all %>% group_by(year) %>% summarise(floral_days_all=n())

# join to flor_annual data frame
flor_annual<-left_join(flor_annual,total_days_high_all,by="year")
flor_annual<-left_join(flor_annual,total_days_high,by="year")

# rename columns
names(flor_annual)[2:5]<-c("total_flowers","total_flowers_80","floral_days","floral_days_80")

#write.csv(flor_annual,"floral_data_annual_summaries_forsexratios_panurginus_cressoniellus.csv",row.names = FALSE)
