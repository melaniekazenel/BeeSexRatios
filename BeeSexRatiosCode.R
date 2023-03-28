# Bee sex ratios and climate change
# RMBL bee and plant phenology data project
# Melanie Kazenel
# Created 7 October 2022
# ---------------------------------

# load relevant libraries
library(tidyr)
library(dplyr)
library(effects)
library(lme4)
?library(corrplot)
library(car)
library(MuMIn)
library(piecewiseSEM)
library(lubridate)
library(visreg)
library(ggplot2)
library(DiagrammeR)
library(tidyverse)
library(semEff)


# set wd
setwd("~/Documents/*RMBLPhenology/Projects/BeeSexRatios/Analyses")

# DECISIONS ABOUT ANALYSES
# (3) SEM across species: just prior year's variables, as originally done
# (4) SEM for individual species: as above


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
year_sex_wide<-pivot_wider(year_sex_summary,names_from = sex,values_from = obs)
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

# fix site data names
beedatafocal$site[beedatafocal$site == "Almont curve"] <- "Almont Curve"
beedatafocal$site[beedatafocal$site == "Kettle ponds"] <- "Kettle Ponds"
beedatafocal$site[beedatafocal$site == "snodgrass"] <- "Snodgrass"
unique(beedatafocal$site)

##### Look at correlation between billy's precip data and interpolated precip data #####
climate_billy<-read.csv("weather_summaries_billybarr_forsexratios.csv")
precip_interp<-read.csv("may-sept_precip_pred_2023-03-14.csv")

# rename site column
names(precip_interp)[1]<-"site"
# make column for precip in mm
precip_interp<-precip_interp %>% mutate(accum_summer_precip_cm=Precip_mm/10)

names(climate_billy)[2]<-"billy_accum_summer_precip_cm"
names(precip_interp)[4]<-"interp_accum_summer_precip_cm"

precip_interp<-left_join(precip_interp,climate_billy,by="year")

ggplot(data=precip_interp, aes(x=billy_accum_summer_precip_cm, y=interp_accum_summer_precip_cm)) + geom_point() + theme_bw() + geom_smooth(method="lm")

ggplot(data=precip_interp, aes(x=billy_accum_summer_precip_cm, y=interp_accum_summer_precip_cm)) + geom_point() + theme_bw() + geom_smooth(method="lm") + facet_wrap(~site) + geom_abline(intercept = 0, slope = 1)


##### Add climate data #####

# # read in billy barr climate data
# climate_billy<-read.csv("weather_summaries_billybarr_forsexratios.csv")
# # remove variables included in Ian's dataset
# climate_billy<-climate_billy[,c(1,2)]

# read in Ian's temperature data
temp<-read.csv("bee_sites_temp_2006_2022.csv")
names(temp)[2:3]<-c("site","date")

# calculate GDD
temp$GDD_calc<-(temp$Tmax+temp$Tmin)/2

# create column to use in calculating accumulated degree days above 0, and column for below zero
temp$DD_abovezero<-ifelse(temp$GDD_calc>0, temp$GDD_calc, 0)
temp$DD_belowzero<-ifelse(temp$GDD_calc<0, temp$GDD_calc, 0)

# convert date to standard format
temp$date<-mdy(temp$date)

# create year, month, and day columns
temp$year<-year(temp$date)
temp$month<-month(temp$date)
temp$day<-day(temp$date)
temp$doy<-yday(temp$date)

# calculate May-September cumulative degree days above 0
climate<-temp %>% filter(month>=5 & month<=9) %>% group_by(site,year) %>% summarise(accum_summer_dd=sum(DD_abovezero))

# calculate October-April mean temperature
temp2<-temp %>% mutate(Tmean=(Tmax+Tmin)/2,
                      water_year=ifelse(month>=10,year+1,year))
wintertemp<-temp2 %>% filter(month>=10 | month<=4) %>% group_by(site,water_year) %>% summarise(winter_temp_mean=mean(Tmean))

snow1<-read.csv("bee_sites_snow_2006_2022.csv")
names(snow1)[2:4]<-c("site","water_year","doy_bare_ground")
wintertemp<-left_join(wintertemp, snow1[,c(2,3,4)], by=c("site","water_year"))

#ggplot(data=wintertemp, aes(x=winter_temp_mean, y=doy_bare_ground)) + geom_point() + theme_bw() + geom_smooth(method="lm") #+ facet_wrap(~site)

names(wintertemp)[2]<-"year"
climate<-left_join(climate,wintertemp[,c(1:3)],by=c("site","year"))

# add snow data to climate data frame
snow<-read.csv("bee_sites_snow_2006_2022.csv")
names(snow)[2:4]<-c("site","year","doy_bare_ground")

# calculate length of snow-free season
snow <- snow %>% mutate(snow_free_days = SnowOnsetDOY_calendar - doy_bare_ground)

# join data frames
climate<-left_join(climate,snow[,c(2:8)],by=c("site","year"))

# # add billy barr precip data to climate data frame
# climate<-left_join(climate,climate_billy,by="year")

# read in interpolated precip data
precip_interp<-read.csv("may-sept_precip_pred_2023-03-14.csv")
# rename site column
names(precip_interp)[1]<-"site"
# make column for precip in mm
precip_interp<-precip_interp %>% mutate(accum_summer_precip_cm=Precip_mm/10)

# add interpolated precip. data to climate data frame
climate<-left_join(climate,precip_interp[,c(1,3,4)],by=c("site","year"))

# add present year's climate data to bee dataset
beedatafocal<-left_join(beedatafocal,climate,by=c("site","year"))

# calculate fall freezing degree days without snowpack
temp<-left_join(temp, snow[,c(2,3,7)], by=c("site","year"))
temp <- temp %>% mutate(SnowOnsetDOY_calendar=round(SnowOnsetDOY_calendar))
temp <- temp %>% mutate(snow_day=SnowOnsetDOY_calendar-doy)
temp2<-temp %>% filter(snow_day>0 & month>=10) %>% group_by(site,year) %>% summarise(fall_freezing_dd_presnow=sum(DD_belowzero)*(-1))
# assign values from prior year to following year in dataset
temp2 <- temp2 %>% mutate(year=year+1) %>% rename(prev_yr_fall_freezing_dd_presnow=fall_freezing_dd_presnow)

# add prior year's climate data to bee dataset
climate_prioryear<-climate
climate_prioryear$winter_temp_mean<-NULL
names(climate_prioryear)[3:9]<-c("prev_yr_accum_summer_dd","prev_yr_doy_bare_ground","prev_yr_SnowOnsetDOY","prev_yr_SnowLengthDays","prev_yr_SnowOnsetDOY_calendar","prev_yr_snow_free_days","prev_yr_accum_summer_precip_cm")
climate_prioryear$year<-climate_prioryear$year+1

climate_prioryear<-left_join(climate_prioryear,temp2,by=c("site","year"))

beedatafocal<-left_join(beedatafocal,climate_prioryear,by=c("site","year"))

# z-score variables of interest
beedatafocal$accum_summer_precip_cm_z<-as.numeric(scale(beedatafocal$accum_summer_precip_cm, center = TRUE, scale = TRUE))
beedatafocal$doy_bare_ground_z<-as.numeric(scale(beedatafocal$doy_bare_ground, center = TRUE, scale = TRUE))
beedatafocal$snow_free_days_z<-as.numeric(scale(beedatafocal$snow_free_days, center = TRUE, scale = TRUE))
beedatafocal$prev_yr_accum_summer_precip_cm_z<-as.numeric(scale(beedatafocal$prev_yr_accum_summer_precip_cm, center = TRUE, scale = TRUE))
beedatafocal$prev_yr_doy_bare_ground_z<-as.numeric(scale(beedatafocal$prev_yr_doy_bare_ground, center = TRUE, scale = TRUE))
beedatafocal$prev_yr_snow_free_days_z<-as.numeric(scale(beedatafocal$prev_yr_snow_free_days, center = TRUE, scale = TRUE))
beedatafocal$prev_yr_fall_freezing_dd_presnow_z<-as.numeric(scale(beedatafocal$prev_yr_fall_freezing_dd_presnow, center = TRUE, scale = TRUE))
beedatafocal$winter_temp_mean_z<-as.numeric(scale(beedatafocal$winter_temp_mean, center = TRUE, scale = TRUE))

# add floral data
flowers<-read.csv("focalbeedata_envdata_RMBLsexratios_2022-11-07.csv")
beedatafocal2<-left_join(beedatafocal,flowers[,c(2,50:65)],by=c("unique_id_year"))

# create csv file
#write.csv(beedatafocal2,"focalbeedata_envdata_RMBLsexratios_2023-03-24.csv",row.names=FALSE)

##### Look at correlations between climate variables #####
beedatafocal<-read.csv("focalbeedata_envdata_RMBLsexratios_2023-03-24.csv")

ggplot(data=beedatafocal, aes(x=doy_bare_ground, y=SnowOnsetDOY)) + geom_point() + theme_bw() + geom_smooth(method="lm") #+ facet_wrap(~site)

ggplot(data=beedatafocal, aes(x=doy_bare_ground, y=SnowOnsetDOY_calendar)) + geom_point() + theme_bw() + geom_smooth(method="lm") #+ facet_wrap(~site)

ggplot(data=beedatafocal, aes(x=doy_bare_ground, y=SnowLengthDays)) + geom_point() + theme_bw() + geom_smooth(method="lm") #+ facet_wrap(~site)

ggplot(data=beedatafocal, aes(x=doy_bare_ground, y=SnowLengthDays)) + geom_point() + theme_bw() + geom_smooth(method="lm") #+ facet_wrap(~site)



##### Graph relationships between each climate variable and sex ratios #####
beedatafocal<-read.csv("focalbeedata_envdata_RMBLsexratios_2023-03-24.csv")

# subset to just include focal sites
# 16 sites currently sampled
sitelist<-c("Almont","Almont Curve","Beaver","CDOT","Copper","Davids","Elko","Gothic","Hill","Little","Lypps","Mexican Cut","Rustlers","Seans","Tuttle","Willey")
# 18 sites used in Michael's 2020 paper
#sitelist<-c("Almont","Almont Curve","Beaver","CDOT","Copper","Davids","Elko","Gothic","Hill","Little","Lypps","Mexican Cut","Rustlers","Seans","Tuttle","Willey","Kettle Ponds", "Snodgrass")
beedatafocal<-filter(beedatafocal, site %in% sitelist)

siteinfo<-read.csv("site_info.csv")
beedatafocal<-left_join(beedatafocal,siteinfo,by="site")

library(tidyverse)
beedatafocal <- beedatafocal %>%
  mutate(elev_group=factor(elev_group)) %>% 
  mutate(elev_group=fct_relevel(elev_group,c("Low","Mid","High")))

levels(beedatafocal$elev_group)

# variables with low pairwise correlations to focus on:
# accum_summer_precip_cm
# doy_bare_ground
# prev_yr_accum_summer_precip_cm
# prev_yr_doy_bare_ground

# examine graphs

# accum_summer_precip_cm
ggplot(beedatafocal, aes(y = sex_binom, x = accum_summer_precip_cm)) +
  geom_point() +
  theme_classic(base_size = 15) +
  geom_smooth(method = glm, method.args= list(family="binomial"))
ggplot(beedatafocal, aes(y = sex_binom, x = accum_summer_precip_cm)) +
  geom_point() +
  theme_classic(base_size = 15) +
  facet_wrap(~site) +
  geom_smooth(method = glm, method.args= list(family="binomial"))
ggplot(beedatafocal, aes(y = sex_binom, x = accum_summer_precip_cm)) +
  geom_point() +
  theme_classic(base_size = 15) +
  facet_wrap(~genus_species) +
  geom_smooth(method = glm, method.args= list(family="binomial"))
ggplot(beedatafocal, aes(y = sex_binom, x = accum_summer_precip_cm)) +
  geom_point() +
  theme_classic(base_size = 15) +
  facet_wrap(~elev_group) +
  geom_smooth(method = glm, method.args= list(family="binomial"))
# model
m1 <- glmer(sex_binom ~ accum_summer_precip_cm_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)
summary(m1)
library(effects)
plot(allEffects(m1))

# doy_bare_ground
ggplot(beedatafocal, aes(y = sex_binom, x = doy_bare_ground)) +
  geom_point() +
  theme_classic(base_size = 15) +
  geom_smooth(method = glm, method.args= list(family="binomial"))
ggplot(beedatafocal, aes(y = sex_binom, x = doy_bare_ground)) +
  geom_point() +
  theme_classic(base_size = 15) +
  facet_wrap(~site) +
  geom_smooth(method = glm, method.args= list(family="binomial"))
ggplot(beedatafocal, aes(y = sex_binom, x = doy_bare_ground)) +
  geom_point() +
  theme_classic(base_size = 15) +
  facet_wrap(~genus_species) +
  geom_smooth(method = glm, method.args= list(family="binomial"))
ggplot(beedatafocal, aes(y = sex_binom, x = doy_bare_ground)) +
  geom_point() +
  theme_classic(base_size = 15) +
  facet_wrap(~elev_group) +
  geom_smooth(method = glm, method.args= list(family="binomial"))
# model
m2 <- glmer(sex_binom ~ doy_bare_ground_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)
summary(m2)
plot(allEffects(m2))

# prev_yr_accum_summer_precip_cm
ggplot(beedatafocal, aes(y = sex_binom, x = prev_yr_accum_summer_precip_cm)) +
  geom_point() +
  theme_classic(base_size = 15) +
  geom_smooth(method = glm, method.args= list(family="binomial"))
ggplot(beedatafocal, aes(y = sex_binom, x = prev_yr_accum_summer_precip_cm)) +
  geom_point() +
  theme_classic(base_size = 15) +
  facet_wrap(~site) +
  geom_smooth(method = glm, method.args= list(family="binomial"))
ggplot(beedatafocal, aes(y = sex_binom, x = prev_yr_accum_summer_precip_cm)) +
  geom_point() +
  theme_classic(base_size = 15) +
  facet_wrap(~genus_species) +
  geom_smooth(method = glm, method.args= list(family="binomial"))
ggplot(beedatafocal, aes(y = sex_binom, x = prev_yr_accum_summer_precip_cm)) +
  geom_point() +
  theme_classic(base_size = 15) +
  facet_wrap(~elev_group) +
  geom_smooth(method = glm, method.args= list(family="binomial"))
# model
# prev_yr_accum_summer_precip_cm_z
m3 <- glmer(sex_binom ~ prev_yr_accum_summer_precip_cm_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)
summary(m3)
plot(allEffects(m3))

# prev_yr_doy_bare_ground
ggplot(beedatafocal, aes(y = sex_binom, x = prev_yr_doy_bare_ground)) +
  geom_point() +
  theme_classic(base_size = 15) +
  geom_smooth(method = glm, method.args= list(family="binomial"))
ggplot(beedatafocal, aes(y = sex_binom, x = prev_yr_doy_bare_ground)) +
  geom_point() +
  theme_classic(base_size = 15) +
  facet_wrap(~site) +
  geom_smooth(method = glm, method.args= list(family="binomial"))
ggplot(beedatafocal, aes(y = sex_binom, x = prev_yr_doy_bare_ground)) +
  geom_point() +
  theme_classic(base_size = 15) +
  facet_wrap(~genus_species) +
  geom_smooth(method = glm, method.args= list(family="binomial"))
ggplot(beedatafocal, aes(y = sex_binom, x = prev_yr_doy_bare_ground)) +
  geom_point() +
  theme_classic(base_size = 15) +
  facet_wrap(~elev_group) +
  geom_smooth(method = glm, method.args= list(family="binomial"))
# model
m4 <- glmer(sex_binom ~ prev_yr_doy_bare_ground + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)
summary(m4)
plot(allEffects(m4))

library(visreg)
dev.off()
p<-visreg(m4,"prev_yr_doy_bare_ground",type="conditional", scale="response", gg=TRUE) 
p<- p + xlab("Prior year's snowmelt date") + ylab("Sex \n(female=1, male=0") +
  theme_bw(base_size = 17)
p



##### GLMMs: Across species, do climate variables together predict female/male counts? #####
beedatafocal<-read.csv("focalbeedata_envdata_RMBLsexratios_2023-03-24.csv")

# subset to just include focal sites
# 16 sites currently sampled
sitelist<-c("Almont","Almont Curve","Beaver","CDOT","Copper","Davids","Elko","Gothic","Hill","Little","Lypps","Mexican Cut","Rustlers","Seans","Tuttle","Willey")
# 18 sites used in Michael's 2020 paper
#sitelist<-c("Almont","Almont Curve","Beaver","CDOT","Copper","Davids","Elko","Gothic","Hill","Little","Lypps","Mexican Cut","Rustlers","Seans","Tuttle","Willey","Kettle Ponds", "Snodgrass")
beedatafocal<-filter(beedatafocal, site %in% sitelist)

# variables with low pairwise correlations to focus on:
# accum_summer_precip_cm
# doy_bare_ground
# prev_yr_accum_summer_precip_cm
# prev_yr_doy_bare_ground

# construct models
# accum_summer_precip_cm_z
m1 <- glmer(sex_binom ~ accum_summer_precip_cm_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)
summary(m1)

# doy_bare_ground_z
m2 <- glmer(sex_binom ~ doy_bare_ground_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)

# prev_yr_accum_summer_precip_cm_z
m3 <- glmer(sex_binom ~ prev_yr_accum_summer_precip_cm_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)

# prev_yr_doy_bare_ground_z - BEST SINGLE VARIABLE MODEL
m4 <- glmer(sex_binom ~ prev_yr_doy_bare_ground_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)

# accum_summer_precip_cm_z + doy_bare_ground_z
m5 <- glmer(sex_binom ~ accum_summer_precip_cm_z + doy_bare_ground_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)

# accum_summer_precip_cm_z + prev_yr_accum_summer_precip_cm_z
m6 <- glmer(sex_binom ~ accum_summer_precip_cm_z + prev_yr_accum_summer_precip_cm_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)

# accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z
m7 <- glmer(sex_binom ~ accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)

# doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z
m8 <- glmer(sex_binom ~ doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)

# doy_bare_ground_z + prev_yr_doy_bare_ground_z
m9 <- glmer(sex_binom ~ doy_bare_ground_z + prev_yr_doy_bare_ground_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)

# prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z
m10 <- glmer(sex_binom ~ prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)

# accum_summer_precip_cm_z + doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z
m11 <- glmer(sex_binom ~ accum_summer_precip_cm_z + doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)

# accum_summer_precip_cm_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z
m12 <- glmer(sex_binom ~ accum_summer_precip_cm_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)

# accum_summer_precip_cm_z + doy_bare_ground_z + prev_yr_doy_bare_ground_z
m13 <- glmer(sex_binom ~ accum_summer_precip_cm_z + doy_bare_ground_z + prev_yr_doy_bare_ground_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)

# doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z
m14 <- glmer(sex_binom ~ doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)

# all climate variables - BEST MODEL OVERALL
m15 <- glmer(sex_binom ~ doy_bare_ground_z + accum_summer_precip_cm_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)

# # adding models that include winter temperature - WEREN'T AS PREDICTIVE
# # winter_temp_mean_z
# ma<-glmer(sex_binom ~ winter_temp_mean_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)
# 
# # winter_temp_mean_z + prev_yr_accum_summer_precip_cm_z
# mb<-glmer(sex_binom ~ winter_temp_mean_z + prev_yr_accum_summer_precip_cm_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)
# 
# # winter_temp_mean_z + prev_yr_doy_bare_ground_z
# mc<-glmer(sex_binom ~ winter_temp_mean_z + prev_yr_doy_bare_ground_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)
# 
# # winter_temp_mean_z + accum_summer_precip_cm_z
# md<-glmer(sex_binom ~ winter_temp_mean_z + accum_summer_precip_cm_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)
# 
# # winter_temp_mean_z + prev_yr_accum_summer_precip_cm_z + accum_summer_precip_cm_z
# me<-glmer(sex_binom ~ winter_temp_mean_z + prev_yr_accum_summer_precip_cm_z + accum_summer_precip_cm_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)
# 
# # winter_temp_mean_z + accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z
# mf<-glmer(sex_binom ~ winter_temp_mean_z + accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)
# 
# # winter_temp_mean_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z
# mg<-glmer(sex_binom ~ winter_temp_mean_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)
# 
# # winter_temp_mean_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z + accum_summer_precip_cm_z
# mh<-glmer(sex_binom ~ winter_temp_mean_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z + accum_summer_precip_cm_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)

aicc<-AICc(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15
           #,ma,mb,mc,md,me,mf,mg,mh
           )
aicc<-arrange(aicc,AICc)
aicc$deltaAICc<-aicc$AICc-min(aicc$AICc)
aicc

# get summary from best model
summary(m15)
rsquared(m15)
vif(m15)

df<-data.frame(summary(m15)$coefficients)
#write.csv(df,"output.csv")

# get summary from best single-variable model
summary(m4)
rsquared(m4)

# get vifs for other models - all are ok
vif(m5)
vif(m6)
vif(m7)
vif(m8)
vif(m9)
vif(m10)
vif(m11)
vif(m12)
vif(m13)
vif(m14)

##### GLMMs: For each species individually, how does climate relate to female/male counts? #####
beedatafocal<-read.csv("focalbeedata_envdata_RMBLsexratios_2023-03-24.csv")

# subset to just include focal sites
sitelist<-c("Almont","Almont Curve","Beaver","CDOT","Copper","Davids","Elko","Gothic","Hill","Little","Lypps","Mexican Cut","Rustlers","Seans","Tuttle","Willey")
beedatafocal<-filter(beedatafocal, site %in% sitelist) 

species_list<-unique(beedatafocal$genus_species)
output<-matrix(nrow=1,ncol=21,byrow=TRUE,dimnames=list(c("row1"),c("model","df","AICc","pvalue1","pvalue2","pvalue3","pvalue4","slope1","slope2","slope3","slope4","se1","se2","se3","se4","delta_AICc","vif1","vif2","vif3","vif4","genus_species")))

for (i in 1:length(species_list)){
  data<-filter(beedatafocal, genus_species==species_list[i])
  
  # construct models
  # accum_summer_precip_cm_z
  m1 <- glmer(sex_binom ~ accum_summer_precip_cm_z + (1|site) , data = data, family = binomial)
  
  # doy_bare_ground_z
  m2 <- glmer(sex_binom ~ doy_bare_ground_z + (1|site) , data = data, family = binomial)
  
  # prev_yr_accum_summer_precip_cm_z
  m3 <- glmer(sex_binom ~ prev_yr_accum_summer_precip_cm_z + (1|site) , data = data, family = binomial)
  
  # prev_yr_doy_bare_ground_z
  m4 <- glmer(sex_binom ~ prev_yr_doy_bare_ground_z + (1|site) , data = data, family = binomial)
  
  # accum_summer_precip_cm_z + doy_bare_ground_z
  m5 <- glmer(sex_binom ~ accum_summer_precip_cm_z + doy_bare_ground_z + (1|site) , data = data, family = binomial)
  
  # accum_summer_precip_cm_z + prev_yr_accum_summer_precip_cm_z
  m6 <- glmer(sex_binom ~ accum_summer_precip_cm_z + prev_yr_accum_summer_precip_cm_z + (1|site) , data = data, family = binomial)
  
  # accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z
  m7 <- glmer(sex_binom ~ accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z + (1|site) , data = data, family = binomial)
  
  # doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z
  m8 <- glmer(sex_binom ~ doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z + (1|site) , data = data, family = binomial)
  
  # doy_bare_ground_z + prev_yr_doy_bare_ground_z
  m9 <- glmer(sex_binom ~ doy_bare_ground_z + prev_yr_doy_bare_ground_z + (1|site) , data = data, family = binomial)
  
  # prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z
  m10 <- glmer(sex_binom ~ prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z + (1|site) , data = data, family = binomial)
  
  # accum_summer_precip_cm_z + doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z
  m11 <- glmer(sex_binom ~ accum_summer_precip_cm_z + doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z + (1|site) , data = data, family = binomial)
  
  # accum_summer_precip_cm_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z
  m12 <- glmer(sex_binom ~ accum_summer_precip_cm_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z + (1|site) , data = data, family = binomial)
  
  # accum_summer_precip_cm_z + doy_bare_ground_z + prev_yr_doy_bare_ground_z
  m13 <- glmer(sex_binom ~ accum_summer_precip_cm_z + doy_bare_ground_z + prev_yr_doy_bare_ground_z + (1|site) , data = data, family = binomial)
  
  # doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z
  m14 <- glmer(sex_binom ~ doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z + (1|site) , data = data, family = binomial)
  
  # all climate variables
  m15 <- glmer(sex_binom ~ doy_bare_ground_z + accum_summer_precip_cm_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z + (1|site) , data = data, family = binomial)
  
  aicc<-AICc(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15)
  aicc<- tibble::rownames_to_column(aicc, "model")
  aicc$pvalue1<-NA
  aicc$pvalue2<-NA
  aicc$pvalue3<-NA
  aicc$pvalue4<-NA
  aicc$slope1<-NA
  aicc$slope2<-NA
  aicc$slope3<-NA
  aicc$slope4<-NA
  aicc$se1<-NA
  aicc$se2<-NA
  aicc$se3<-NA
  aicc$se4<-NA
  aicc$vif1<-NA
  aicc$vif2<-NA
  aicc$vif3<-NA
  aicc$vif4<-NA
  
  df<-data.frame(summary(m1)[["coefficients"]])
  df<-df[-1,]
  aicc[1,4:sum(length(df$Estimate)+3)]<-df[,4]
  aicc[1,8:sum(length(df$Estimate)+7)]<-df[,1]
  aicc[1,12:sum(length(df$Estimate)+11)]<-df[,2]
  
  df<-data.frame(summary(m2)[["coefficients"]])
  df<-df[-1,]
  aicc[2,4:sum(length(df$Estimate)+3)]<-df[,4]
  aicc[2,8:sum(length(df$Estimate)+7)]<-df[,1]
  aicc[2,12:sum(length(df$Estimate)+11)]<-df[,2]
  
  df<-data.frame(summary(m3)[["coefficients"]])
  df<-df[-1,]
  aicc[3,4:sum(length(df$Estimate)+3)]<-df[,4]
  aicc[3,8:sum(length(df$Estimate)+7)]<-df[,1]
  aicc[3,12:sum(length(df$Estimate)+11)]<-df[,2]
  
  df<-data.frame(summary(m4)[["coefficients"]])
  df<-df[-1,]
  aicc[4,4:sum(length(df$Estimate)+3)]<-df[,4]
  aicc[4,8:sum(length(df$Estimate)+7)]<-df[,1]
  aicc[4,12:sum(length(df$Estimate)+11)]<-df[,2]
  
  df<-data.frame(summary(m5)[["coefficients"]])
  df<-df[-1,]
  aicc[5,4:sum(length(df$Estimate)+3)]<-df[,4]
  aicc[5,8:sum(length(df$Estimate)+7)]<-df[,1]
  aicc[5,12:sum(length(df$Estimate)+11)]<-df[,2]
  aicc[5,16:sum(length(df$Estimate)+15)]<-data.frame(vif(m5))[,1]

  df<-data.frame(summary(m6)[["coefficients"]])
  df<-df[-1,]
  aicc[6,4:sum(length(df$Estimate)+3)]<-df[,4]
  aicc[6,8:sum(length(df$Estimate)+7)]<-df[,1]
  aicc[6,12:sum(length(df$Estimate)+11)]<-df[,2]
  aicc[6,16:sum(length(df$Estimate)+15)]<-data.frame(vif(m6))[,1]
  
  df<-data.frame(summary(m7)[["coefficients"]])
  df<-df[-1,]
  aicc[7,4:sum(length(df$Estimate)+3)]<-df[,4]
  aicc[7,8:sum(length(df$Estimate)+7)]<-df[,1]
  aicc[7,12:sum(length(df$Estimate)+11)]<-df[,2]
  aicc[7,16:sum(length(df$Estimate)+15)]<-data.frame(vif(m7))[,1]
  
  df<-data.frame(summary(m8)[["coefficients"]])
  df<-df[-1,]
  aicc[8,4:sum(length(df$Estimate)+3)]<-df[,4]
  aicc[8,8:sum(length(df$Estimate)+7)]<-df[,1]
  aicc[8,12:sum(length(df$Estimate)+11)]<-df[,2]
  aicc[8,16:sum(length(df$Estimate)+15)]<-data.frame(vif(m8))[,1]
  
  df<-data.frame(summary(m9)[["coefficients"]])
  df<-df[-1,]
  aicc[9,4:sum(length(df$Estimate)+3)]<-df[,4]
  aicc[9,8:sum(length(df$Estimate)+7)]<-df[,1]
  aicc[9,12:sum(length(df$Estimate)+11)]<-df[,2]
  aicc[9,16:sum(length(df$Estimate)+15)]<-data.frame(vif(m9))[,1]
  
  df<-data.frame(summary(m10)[["coefficients"]])
  df<-df[-1,]
  aicc[10,4:sum(length(df$Estimate)+3)]<-df[,4]
  aicc[10,8:sum(length(df$Estimate)+7)]<-df[,1]
  aicc[10,12:sum(length(df$Estimate)+11)]<-df[,2]
  aicc[10,16:sum(length(df$Estimate)+15)]<-data.frame(vif(m10))[,1]
  
  df<-data.frame(summary(m11)[["coefficients"]])
  df<-df[-1,]
  aicc[11,4:sum(length(df$Estimate)+3)]<-df[,4]
  aicc[11,8:sum(length(df$Estimate)+7)]<-df[,1]
  aicc[11,12:sum(length(df$Estimate)+11)]<-df[,2]
  aicc[11,16:sum(length(df$Estimate)+15)]<-data.frame(vif(m11))[,1]
  
  df<-data.frame(summary(m12)[["coefficients"]])
  df<-df[-1,]
  aicc[12,4:sum(length(df$Estimate)+3)]<-df[,4]
  aicc[12,8:sum(length(df$Estimate)+7)]<-df[,1]
  aicc[12,12:sum(length(df$Estimate)+11)]<-df[,2]
  aicc[12,16:sum(length(df$Estimate)+15)]<-data.frame(vif(m12))[,1]
  
  df<-data.frame(summary(m13)[["coefficients"]])
  df<-df[-1,]
  aicc[13,4:sum(length(df$Estimate)+3)]<-df[,4]
  aicc[13,8:sum(length(df$Estimate)+7)]<-df[,1]
  aicc[13,12:sum(length(df$Estimate)+11)]<-df[,2]
  aicc[13,16:sum(length(df$Estimate)+15)]<-data.frame(vif(m13))[,1]
  
  df<-data.frame(summary(m14)[["coefficients"]])
  df<-df[-1,]
  aicc[14,4:sum(length(df$Estimate)+3)]<-df[,4]
  aicc[14,8:sum(length(df$Estimate)+7)]<-df[,1]
  aicc[14,12:sum(length(df$Estimate)+11)]<-df[,2]
  aicc[14,16:sum(length(df$Estimate)+15)]<-data.frame(vif(m14))[,1]
  
  df<-data.frame(summary(m15)[["coefficients"]])
  df<-df[-1,]
  aicc[15,4:sum(length(df$Estimate)+3)]<-df[,4]
  aicc[15,8:sum(length(df$Estimate)+7)]<-df[,1]
  aicc[15,12:sum(length(df$Estimate)+11)]<-df[,2]
  aicc[15,16:sum(length(df$Estimate)+15)]<-data.frame(vif(m15))[,1]
  
  aicc<-arrange(aicc,AICc)
  aicc$delta_AICc<-aicc$AICc-min(aicc$AICc)
  aicc$genus_species<-species_list[i]

  output<-rbind(output,aicc)
}

output<-output[-1,]

# add column listing predictors in each model to "output" df
aicc<-AICc(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15)
aicc<- tibble::rownames_to_column(aicc, "model")
aicc$predictors<-c("accum_summer_precip_cm_z", "doy_bare_ground_z", "prev_yr_accum_summer_precip_cm_z", "prev_yr_doy_bare_ground_z", "accum_summer_precip_cm_z + doy_bare_ground_z", "accum_summer_precip_cm_z + prev_yr_accum_summer_precip_cm_z", "accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z", "doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z", "doy_bare_ground_z + prev_yr_doy_bare_ground_z", "prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z", "accum_summer_precip_cm_z + doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z", "accum_summer_precip_cm_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z", "accum_summer_precip_cm_z + doy_bare_ground_z + prev_yr_doy_bare_ground_z", "doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z","all predictors")
aicc<-aicc[,-c(2:3)]
output<-left_join(output,aicc,by="model")

#write.csv(output,"species_glmms_climate_sex_ratios_2023-03-23.csv", row.names = FALSE)





##### For individual species: graphs of relationships with different climate variables #####
results<-read.csv("species_glmms_climate_sex_ratios_2023-03-24.csv")
bestmodels<-filter(results,delta_AICc==0)
slopes_long<-pivot_longer(bestmodels[,c(21,1,22,8:11)],cols=4:7,values_to="slope")
names(slopes_long)[4]<-"term"
slopes_long$term<-substr(slopes_long$term,6,6)

se_long<-pivot_longer(bestmodels[,c(21,1,22,12:15)],cols=4:7,values_to="se")
names(se_long)[4]<-"term"
se_long$term<-substr(se_long$term,3,3)

pvalues_long<-pivot_longer(bestmodels[,c(21,1,22,4:7)],cols=4:7,values_to="pvalue")
names(pvalues_long)[4]<-"term"
pvalues_long$term<-substr(pvalues_long$term,7,7)

slope_se_long<-left_join(slopes_long,se_long, by=c("genus_species","model","predictors","term"))
slope_se_long<-left_join(slope_se_long,pvalues_long, by=c("genus_species","model","predictors","term"))


slope_se_long<-drop_na(slope_se_long)
summary_terms<-slope_se_long %>% group_by(model, predictors, term) %>% summarise(count=n())
summary_terms$count<-NULL

summary_terms$term_id<-ifelse(summary_terms$model=="m1" & summary_terms$term=="1" |
                                summary_terms$model=="m12" & summary_terms$term=="1" |
                                summary_terms$model =="m13" & summary_terms$term=="1" |
                                summary_terms$model=="m5" & summary_terms$term=="1"|
                                summary_terms$model=="m6" & summary_terms$term=="1"|
                                summary_terms$model=="m15" & summary_terms$term=="2",
                              "accum_summer_precip_cm_z", " ")

summary_terms$term_id<-ifelse(summary_terms$model=="m3" & summary_terms$term=="1" |
                                summary_terms$model=="m6" & summary_terms$term=="2" |
                                summary_terms$model=="m8" & summary_terms$term=="2" |
                                summary_terms$model=="m12" & summary_terms$term=="2" |
                                summary_terms$model=="m10" & summary_terms$term=="1" |
                                summary_terms$model=="m14" & summary_terms$term=="2" |
                                summary_terms$model=="m15" & summary_terms$term=="3", 
                              "prev_yr_accum_summer_precip_cm_z", summary_terms$term_id)

summary_terms$term_id<-ifelse(summary_terms$model=="m13" & summary_terms$term=="2" |
                                summary_terms$model=="m2" & summary_terms$term=="1" |
                                summary_terms$model=="m8" & summary_terms$term=="1" |
                                summary_terms$model=="m9" & summary_terms$term=="1" |
                                summary_terms$model=="m14" & summary_terms$term=="1" |
                                summary_terms$model=="m15" & summary_terms$term=="1" |
                                summary_terms$model=="m5" & summary_terms$term=="2",
                              "doy_bare_ground_z", summary_terms$term_id)

summary_terms$term_id<-ifelse(summary_terms$model=="m4" & summary_terms$term=="1" |
                                summary_terms$model=="m9" & summary_terms$term=="2" |
                                summary_terms$model=="m12" & summary_terms$term=="3" |
                                summary_terms$model=="m13" & summary_terms$term=="3"| 
                                summary_terms$model=="m10" & summary_terms$term=="2" |
                                summary_terms$model=="m14" & summary_terms$term=="3" |
                                summary_terms$model=="m15" & summary_terms$term=="4", 
                              "prev_yr_doy_bare_ground_z", summary_terms$term_id)

slope_se_long<-left_join(slope_se_long, summary_terms, by=c("model","predictors","term"))


slope_se_long$term_id_labels<-slope_se_long$term_id
slope_se_long$term_id_renamed<-as.factor(slope_se_long$term_id)

levels(slope_se_long$term_id_renamed) <- c("Summer \nprecipitation","Snowmelt date","Prior summer \nprecipitation","Prior year \nsnowmelt date")

slope_se_long$pvalue_sig<-ifelse(slope_se_long$pvalue>=0.05,"Nonsignificant","Significant")

ggplot(slope_se_long, aes(y = slope, x = term_id_renamed, color=pvalue_sig)) +
  geom_point(size=2, position=position_dodge(width=0.5)) +
  xlab("Climate predictor") + ylab("Parameter estimate") +
  geom_errorbar(aes(ymin=slope-se, ymax=slope+se), width=.2,position=position_dodge(width=0.5)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  geom_hline(yintercept=0, linetype='dashed', col = 'darkgray') + facet_wrap(~genus_species,ncol=2) + theme(legend.position = "none") +
  scale_color_manual(values=c("darkgray","black"))


##### Adding floral data: Inouye data #####

# read in floral data
flor<-read.csv("floral_data_annual_summaries_forsexratios.csv")
# add current year's floral data to be dataset
beedatafocal<-left_join(beedatafocal,flor,by="year")
# add prior year's floral data to bee dataset
flor_prior<-flor
flor_prior$year<-flor_prior$year+1
names(flor_prior)[2:5]<-c("prev_yr_total_flowers", "prev_yr_total_flowers_80", "prev_yr_floral_days", "prev_yr_floral_days_80")
beedatafocal<-left_join(beedatafocal,flor_prior, by="year")

# z-score floral variables
beedatafocal$total_flowers_z<-as.numeric(scale(beedatafocal$total_flowers, center = TRUE, scale = TRUE))
beedatafocal$prev_yr_total_flowers_z<-as.numeric(scale(beedatafocal$prev_yr_total_flowers, center = TRUE, scale = TRUE))

beedatafocal$floral_days_z<-as.numeric(scale(beedatafocal$floral_days, center = TRUE, scale = TRUE))
beedatafocal$prev_yr_floral_days_z<-as.numeric(scale(beedatafocal$prev_yr_floral_days, center = TRUE, scale = TRUE))

beedatafocal$total_flowers_80_z<-as.numeric(scale(beedatafocal$total_flowers_80, center = TRUE, scale = TRUE))
beedatafocal$prev_yr_total_flowers_80_z<-as.numeric(scale(beedatafocal$prev_yr_total_flowers_80, center = TRUE, scale = TRUE))

beedatafocal$floral_days_80_z<-as.numeric(scale(beedatafocal$floral_days_80, center = TRUE, scale = TRUE))
beedatafocal$prev_yr_floral_days_80_z<-as.numeric(scale(beedatafocal$prev_yr_floral_days_80, center = TRUE, scale = TRUE))


##### Relationships between climate and floral variables: Inouye data #####

beedatafocal<-read.csv("focalbeedata_envdata_RMBLsexratios_2023-03-24.csv")

# subset to just include focal sites
sitelist<-c("Almont","Almont Curve","Beaver","CDOT","Copper","Davids","Elko","Gothic","Hill","Little","Lypps","Mexican Cut","Rustlers","Seans","Tuttle","Willey")
beedatafocal<-filter(beedatafocal, site %in% sitelist) 

# just look at mid-elevation sites
siteinfo<-read.csv("site_info.csv")
beedatafocal<-left_join(beedatafocal,siteinfo,by="site")
beedatafocal<-filter(beedatafocal, elev_group=="Mid") 

# optional: remove outlier floral resource year (remove 2016 for examining previous year's resources)
#beedatafocal<-subset(beedatafocal,year!=2016)

# exploratory graphs: relationships between climate and floral variables

ggplot(beedatafocal, aes(y = total_flowers, x = doy_bare_ground)) +
  geom_point() +
  theme_classic(base_size = 15) +
  geom_smooth(method = lm) + 
  xlab("Snowmelt date") + ylab("Floral sum") #+ facet_wrap(~site)
summary(lm(total_flowers~doy_bare_ground,data=beedatafocal))

ggplot(beedatafocal, aes(y = total_flowers, x = accum_summer_precip_cm)) +
  geom_point() +
  theme_classic(base_size = 15) +
  geom_smooth(method = lm) + 
  xlab("Summer precipitation") + ylab("Floral sum") #+ facet_wrap(~site)
summary(lm(total_flowers~accum_summer_precip_cm,data=beedatafocal))

ggplot(beedatafocal, aes(y = floral_days, x = doy_bare_ground)) +
  geom_point() +
  theme_classic(base_size = 15) +
  geom_smooth(method = lm) + 
  xlab("Snowmelt date") + ylab("Floral days") #+ facet_wrap(~site)
summary(lm(floral_days~doy_bare_ground,data=beedatafocal))

ggplot(beedatafocal, aes(y = floral_days, x = accum_summer_precip_cm)) +
  geom_point() +
  theme_classic(base_size = 15) +
  geom_smooth(method = lm) + 
  xlab("Summer precipitation") + ylab("Floral days") + facet_wrap(~site)
summary(lm(floral_days~accum_summer_precip_cm,data=beedatafocal))


##### SEMs: all species #####

beedatafocal<-read.csv("focalbeedata_envdata_RMBLsexratios_2023-03-24.csv")
# subset to just include focal sites
sitelist<-c("Almont","Almont Curve","Beaver","CDOT","Copper","Davids","Elko","Gothic","Hill","Little","Lypps","Mexican Cut","Rustlers","Seans","Tuttle","Willey")
beedatafocal<-filter(beedatafocal, site %in% sitelist) 

# # just look at mid-elevation sites
# siteinfo<-read.csv("site_info.csv")
# beedatafocal<-left_join(beedatafocal,siteinfo,by="site")
# beedatafocal <- beedatafocal %>%
#   mutate(elev_group=factor(elev_group)) %>% 
#   mutate(elev_group=fct_relevel(elev_group,c("Low","Mid","High")))
# levels(beedatafocal$elev_group)
# beedatafocal<-filter(beedatafocal, elev_group=="Mid") 

### (full) SEM: all prior year's climate and floral predictor variables #####
semFull<-psem(

  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_floral_days_z +
          prev_yr_total_flowers_z +
          (1|site) + (1|genus_species),
        data = beedatafocal, family = binomial),

  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~ prev_yr_total_flowers_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),

  lm(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)

)

summary(semFull, .progressBar = F)


# plotting the SEM! resources:
# https://rdrr.io/cran/DiagrammeR/man/create_graph.html
# also: https://rich-iannone.github.io/DiagrammeR/graphs.html

coef<-summary(semFull, .progressBar = F)$coef
coef
unique(coef$Response)
unique(coef$Predictor)

coef$to<-coef$Response
coef$from<-coef$Predictor

coef$to[coef$to=="sex_binom"] <- 1
coef$to[coef$to=="prev_yr_floral_days_z"] <- 2
coef$to[coef$to=="prev_yr_total_flowers_z"] <- 3
coef$to[coef$to=="prev_yr_doy_bare_ground_z"] <- 4
coef$to[coef$to=="prev_yr_accum_summer_precip_cm_z"] <- 5

coef$from[coef$from=="sex_binom"] <- 1
coef$from[coef$from=="prev_yr_floral_days_z"] <- 2
coef$from[coef$from=="prev_yr_total_flowers_z"] <- 3
coef$from[coef$from=="prev_yr_doy_bare_ground_z"] <- 4
coef$from[coef$from=="prev_yr_accum_summer_precip_cm_z"] <- 5

coef$style<-ifelse(coef$Estimate>0,"solid","dashed")
coef$color<-ifelse(coef$P.Value>=0.05,"gray","black")

ndf <-
  create_node_df(
    n = 5,
    fontsize=6,
    fixedsize=TRUE,
    color="black",
    fillcolor="lightgray",
    penwidth=0.5,
    shape="rectangle",
    width=0.7,
    height=0.3,
    label = c("Female/male \nratio","Floral days \n(prior year)", "Floral sum \n(prior year)","Snowmelt date \n(prior year)","Summer precip. \n(prior year)"))

graph <-
  create_graph(
    nodes_df = ndf)

render_graph(graph)

edf <-
  create_edge_df(
    from = coef$from,
    to = coef$to,
    rel = "leading_to",
    label = round(coef$Estimate,digits=2),
    penwidth = abs(coef$Estimate*2),
    fontsize=5.5,
    fontcolor="brown3",
    color=coef$color,
    style=coef$style)

graph <-
  create_graph(
    nodes_df = ndf,
    edges_df = edf)

render_graph(graph)

graph <-
  graph %>%
  set_node_position(
    node = 1, # sex ratio
    x = 3, y = 1.5) %>%
  set_node_position(
    node = 2, # floral days
    x = 2, y = 1) %>%
  set_node_position(
    node = 3, # floral sum
    x = 2, y = 2) %>%
  set_node_position(
    node = 4, # snowmelt date
    x = 1, y = 2) %>%
  set_node_position(
    node = 5, # summer precip
    x = 1, y = 1)

render_graph(graph)


### (1) SEM: snowmelt #####
sem1<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          (1|site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~ prev_yr_total_flowers_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),
  
  lm(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem1, .progressBar = F)
plot(sem1)


### (2) SEM: snowmelt, floral sum #####
sem2<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_total_flowers_z +
          (1|site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~ prev_yr_total_flowers_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),
  
  lm(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem2, .progressBar = F)
plot(sem2)

### (3) SEM: snowmelt, floral days #####
sem3<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_floral_days_z +
          (1|site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~ prev_yr_total_flowers_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),
  
  lm(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem3, .progressBar = F)
plot(sem3)

### (4) SEM: snowmelt, precip #####
sem4<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_accum_summer_precip_cm_z +
          (1|site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~ prev_yr_total_flowers_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),
  
  lm(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem4, .progressBar = F)
plot(sem4)

### (5) SEM: precip #####
sem5<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ 
          prev_yr_accum_summer_precip_cm_z +
          (1|site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~ prev_yr_total_flowers_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),
  
  lm(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem5, .progressBar = F)
plot(sem5)

### (6) SEM: precip, floral sum #####
sem6<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ 
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_total_flowers_z +
          (1|site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~ prev_yr_total_flowers_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),
  
  lm(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem6, .progressBar = F)
plot(sem6)

### (7) SEM: precip, floral days #####
sem7<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ 
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_floral_days_z +
          (1|site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~ prev_yr_total_flowers_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),
  
  lm(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem7, .progressBar = F)
plot(sem7)


### (8) SEM: floral sum, floral days #####
sem8<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~
          prev_yr_floral_days_z +
          prev_yr_total_flowers_z +
          (1|site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~ prev_yr_total_flowers_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),
  
  lm(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem8, .progressBar = F)
plot(sem8)


### (9) SEM: floral sum #####
sem9<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~
          prev_yr_total_flowers_z +
          (1|site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~ prev_yr_total_flowers_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),
  
  lm(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem9, .progressBar = F)
plot(sem9)

### (10) SEM: floral days #####
sem10<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~
          prev_yr_floral_days_z +
          (1|site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~ prev_yr_total_flowers_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),
  
  lm(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem10, .progressBar = F)
plot(sem10)






### (11) SEM: snowmelt, precip, floral sum #####
sem11<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_total_flowers_z +
          (1|site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~ prev_yr_total_flowers_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),
  
  lm(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem11, .progressBar = F)
plot(sem11)


### (12) SEM: snowmelt, precip, floral days #####
sem12<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_floral_days_z +
          (1|site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~ prev_yr_total_flowers_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),
  
  lm(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem12, .progressBar = F)
plot(sem12)

### (13) SEM: snowmelt, floral sum, floral days #####
sem13<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_floral_days_z +
          prev_yr_total_flowers_z +
          (1|site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~ prev_yr_total_flowers_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),
  
  lm(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem13, .progressBar = F)
plot(sem13)

### (14) SEM: precip, floral sum, floral days #####
sem14<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ 
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_floral_days_z +
          prev_yr_total_flowers_z +
          (1|site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~ prev_yr_total_flowers_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),
  
  lm(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem14, .progressBar = F)
plot(sem14)


### Get summaries for SEMs #####

df<-bind_rows(data.frame(summary(semFull, .progressBar = F)$AIC),
data.frame(summary(sem1, .progressBar = F)$AIC),
data.frame(summary(sem2, .progressBar = F)$AIC),
data.frame(summary(sem3, .progressBar = F)$AIC),
data.frame(summary(sem4, .progressBar = F)$AIC), 
data.frame(summary(sem5, .progressBar = F)$AIC),
data.frame(summary(sem6, .progressBar = F)$AIC),
data.frame(summary(sem7, .progressBar = F)$AIC),
data.frame(summary(sem8, .progressBar = F)$AIC),
data.frame(summary(sem9, .progressBar = F)$AIC),
data.frame(summary(sem10, .progressBar = F)$AIC),
data.frame(summary(sem11, .progressBar = F)$AIC),
data.frame(summary(sem12, .progressBar = F)$AIC),
data.frame(summary(sem13, .progressBar = F)$AIC),
data.frame(summary(sem14, .progressBar = F)$AIC))
df$model<-0:14
#write.csv(df,"AIC_value_SEM.csv", row.names=FALSE)

summary(semFull, .progressBar = F) # NA
summary(sem1, .progressBar = F)  # 0
summary(sem2, .progressBar = F)  # 0
summary(sem3, .progressBar = F) # 0
summary(sem4, .progressBar = F)  # 0
summary(sem5, .progressBar = F) # 0
summary(sem6, .progressBar = F) # 0
summary(sem7, .progressBar = F) # 0
summary(sem8, .progressBar = F) # 0
summary(sem9, .progressBar = F) # 0
summary(sem10, .progressBar = F) # 0 
summary(sem11, .progressBar = F) #0
summary(sem12, .progressBar = F) #0
summary(sem13, .progressBar = F) # 0.339 - BEST MODEL
summary(sem14, .progressBar = F) # 0











##### SEMs: individual species #####

beedatafocal<-read.csv("focalbeedata_envdata_RMBLsexratios_2023-03-24.csv")
# subset to just include focal sites
sitelist<-c("Almont","Almont Curve","Beaver","CDOT","Copper","Davids","Elko","Gothic","Hill","Little","Lypps","Mexican Cut","Rustlers","Seans","Tuttle","Willey")
beedatafocal<-filter(beedatafocal, site %in% sitelist) 

# filter down to individual species data
#beedatafocal<-filter(beedatafocal,genus_species=="Panurginus cressoniellus")
#beedatafocal<-filter(beedatafocal,genus_species=="Panurginus ineptus")
beedatafocal<-filter(beedatafocal,genus_species=="Pseudopanurgus bakeri")
#beedatafocal<-filter(beedatafocal,genus_species=="Pseudopanurgus didirupa")
#beedatafocal<-filter(beedatafocal,genus_species=="Ceratina nanula")
#beedatafocal<-filter(beedatafocal,genus_species=="Hylaeus annulatus")
#beedatafocal<-filter(beedatafocal,genus_species=="Agapostemon texanus")
#beedatafocal<-filter(beedatafocal,genus_species=="Dufourea harveyi")
#beedatafocal<-filter(beedatafocal,genus_species=="Halictus rubicundus")
#beedatafocal<-filter(beedatafocal,genus_species=="Halictus virgatellus")
#beedatafocal<-filter(beedatafocal,genus_species=="Hoplitis fulgida")
#beedatafocal<-filter(beedatafocal,genus_species=="Hoplitis robusta")
#beedatafocal<-filter(beedatafocal,genus_species=="Osmia albolateralis")

### (full) SEM: all prior year's climate and floral predictor variables #####
semFull<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_floral_days_z +
          prev_yr_total_flowers_z +
          (1|site) ,
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~ prev_yr_total_flowers_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),
  
  lm(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(semFull, .progressBar = F)


# plotting the SEM! resources:
# https://rdrr.io/cran/DiagrammeR/man/create_graph.html
# also: https://rich-iannone.github.io/DiagrammeR/graphs.html

coef<-summary(semFull, .progressBar = F)$coef
coef
unique(coef$Response)
unique(coef$Predictor)

coef$to<-coef$Response
coef$from<-coef$Predictor

coef$to[coef$to=="sex_binom"] <- 1
coef$to[coef$to=="prev_yr_floral_days_z"] <- 2
coef$to[coef$to=="prev_yr_total_flowers_z"] <- 3
coef$to[coef$to=="prev_yr_doy_bare_ground_z"] <- 4
coef$to[coef$to=="prev_yr_accum_summer_precip_cm_z"] <- 5

coef$from[coef$from=="sex_binom"] <- 1
coef$from[coef$from=="prev_yr_floral_days_z"] <- 2
coef$from[coef$from=="prev_yr_total_flowers_z"] <- 3
coef$from[coef$from=="prev_yr_doy_bare_ground_z"] <- 4
coef$from[coef$from=="prev_yr_accum_summer_precip_cm_z"] <- 5

coef$style<-ifelse(coef$Estimate>0,"solid","dashed")
coef$color<-ifelse(coef$P.Value>=0.05,"gray","black")

ndf <-
  create_node_df(
    n = 5,
    fontsize=6,
    fixedsize=TRUE,
    color="black",
    fillcolor="lightgray",
    penwidth=0.5,
    shape="rectangle",
    width=0.7,
    height=0.3,
    label = c("Female/male \nratio","Floral days \n(prior year)", "Floral sum \n(prior year)","Snowmelt date \n(prior year)","Summer precip. \n(prior year)"))

graph <-
  create_graph(
    nodes_df = ndf)

render_graph(graph)

edf <-
  create_edge_df(
    from = coef$from,
    to = coef$to,
    rel = "leading_to",
    label = round(coef$Estimate,digits=2),
    penwidth = abs(coef$Estimate*2),
    fontsize=5.5,
    fontcolor="brown3",
    color=coef$color,
    style=coef$style)

graph <-
  create_graph(
    nodes_df = ndf,
    edges_df = edf)

render_graph(graph)

graph <-
  graph %>%
  set_node_position(
    node = 1, # sex ratio
    x = 3, y = 1.5) %>%
  set_node_position(
    node = 2, # floral days
    x = 2, y = 1) %>%
  set_node_position(
    node = 3, # floral sum
    x = 2, y = 2) %>%
  set_node_position(
    node = 4, # snowmelt date
    x = 1, y = 2) %>%
  set_node_position(
    node = 5, # summer precip
    x = 1, y = 1)

render_graph(graph)


### (1) SEM: snowmelt #####
sem1<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          (1|site) ,
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~ prev_yr_total_flowers_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),
  
  lm(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem1, .progressBar = F)
plot(sem1)


### (2) SEM: snowmelt, floral sum #####
sem2<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_total_flowers_z +
          (1|site) ,
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~ prev_yr_total_flowers_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),
  
  lm(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem2, .progressBar = F)
plot(sem2)

### (3) SEM: snowmelt, floral days #####
sem3<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_floral_days_z +
          (1|site) ,
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~ prev_yr_total_flowers_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),
  
  lm(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem3, .progressBar = F)
plot(sem3)

### (4) SEM: snowmelt, precip #####
sem4<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_accum_summer_precip_cm_z +
          (1|site) ,
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~ prev_yr_total_flowers_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),
  
  lm(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem4, .progressBar = F)
plot(sem4)

### (5) SEM: precip #####
sem5<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ 
          prev_yr_accum_summer_precip_cm_z +
          (1|site) ,
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~ prev_yr_total_flowers_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),
  
  lm(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem5, .progressBar = F)
plot(sem5)

### (6) SEM: precip, floral sum #####
sem6<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ 
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_total_flowers_z +
          (1|site) ,
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~ prev_yr_total_flowers_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),
  
  lm(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem6, .progressBar = F)
plot(sem6)

### (7) SEM: precip, floral days #####
sem7<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ 
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_floral_days_z +
          (1|site) ,
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~ prev_yr_total_flowers_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),
  
  lm(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem7, .progressBar = F)
plot(sem7)


### (8) SEM: floral sum, floral days #####
sem8<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~
          prev_yr_floral_days_z +
          prev_yr_total_flowers_z +
          (1|site) ,
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~ prev_yr_total_flowers_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),
  
  lm(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem8, .progressBar = F)
plot(sem8)


### (9) SEM: floral sum #####
sem9<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~
          prev_yr_total_flowers_z +
          (1|site) ,
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~ prev_yr_total_flowers_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),
  
  lm(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem9, .progressBar = F)
plot(sem9)

### (10) SEM: floral days #####
sem10<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~
          prev_yr_floral_days_z +
          (1|site) ,
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~ prev_yr_total_flowers_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),
  
  lm(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem10, .progressBar = F)
plot(sem10)






### (11) SEM: snowmelt, precip, floral sum #####
sem11<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_total_flowers_z +
          (1|site) ,
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~ prev_yr_total_flowers_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),
  
  lm(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem11, .progressBar = F)
plot(sem11)


### (12) SEM: snowmelt, precip, floral days #####
sem12<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_floral_days_z +
          (1|site) ,
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~ prev_yr_total_flowers_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),
  
  lm(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem12, .progressBar = F)
plot(sem12)

### (13) SEM: snowmelt, floral sum, floral days #####
sem13<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_floral_days_z +
          prev_yr_total_flowers_z +
          (1|site) ,
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~ prev_yr_total_flowers_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),
  
  lm(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem13, .progressBar = F)
plot(sem13)

### (14) SEM: precip, floral sum, floral days #####
sem14<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ 
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_floral_days_z +
          prev_yr_total_flowers_z +
          (1|site) ,
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~ prev_yr_total_flowers_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),
  
  lm(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem14, .progressBar = F)
plot(sem14)




### Get summaries for SEMs #####

df1<-bind_rows(data.frame(summary(semFull, .progressBar = F)$AIC),
              data.frame(summary(sem1, .progressBar = F)$AIC),
              data.frame(summary(sem2, .progressBar = F)$AIC),
              data.frame(summary(sem3, .progressBar = F)$AIC),
              data.frame(summary(sem4, .progressBar = F)$AIC), 
              data.frame(summary(sem5, .progressBar = F)$AIC),
              data.frame(summary(sem6, .progressBar = F)$AIC),
              data.frame(summary(sem7, .progressBar = F)$AIC),
              data.frame(summary(sem8, .progressBar = F)$AIC),
              data.frame(summary(sem9, .progressBar = F)$AIC),
              data.frame(summary(sem10, .progressBar = F)$AIC),
              data.frame(summary(sem11, .progressBar = F)$AIC),
              data.frame(summary(sem12, .progressBar = F)$AIC),
              data.frame(summary(sem13, .progressBar = F)$AIC),
              data.frame(summary(sem14, .progressBar = F)$AIC))

df2<-bind_rows(bind_cols(LLchisq(semFull),fisherC(dSep(semFull))),
              bind_cols(LLchisq(sem1),fisherC(dSep(sem1))),
              bind_cols(LLchisq(sem2),fisherC(dSep(sem2))),
              bind_cols(LLchisq(sem3),fisherC(dSep(sem3))),
              bind_cols(LLchisq(sem4),fisherC(dSep(sem4))),
              bind_cols(LLchisq(sem5),fisherC(dSep(sem5))),
              bind_cols(LLchisq(sem6),fisherC(dSep(sem6))),
              bind_cols(LLchisq(sem7),fisherC(dSep(sem7))),
              bind_cols(LLchisq(sem8),fisherC(dSep(sem8))),
              bind_cols(LLchisq(sem9),fisherC(dSep(sem9))),
              bind_cols(LLchisq(sem10),fisherC(dSep(sem10))),
              bind_cols(LLchisq(sem11),fisherC(dSep(sem11))),
              bind_cols(LLchisq(sem12),fisherC(dSep(sem12))),
              bind_cols(LLchisq(sem13),fisherC(dSep(sem13))),
              bind_cols(LLchisq(sem14),fisherC(dSep(sem14)))
              )
names(df2)[c(2,3,5,6)]<-c("df_chisq","p_chisq","df_fisher","p_fisher")
df2$model<-0:14

df3<-bind_cols(df2[,c(7,1:6)],df1)
  
#write.csv(df3,"GoodnessOfFit_SEM_Panurginus_cressoniellus.csv", row.names=FALSE)
#write.csv(df3,"GoodnessOfFit_SEM_Panurginus_ineptus.csv", row.names=FALSE)
#write.csv(df3,"GoodnessOfFit_SEM_Pseudopanurgus_bakeri.csv", row.names=FALSE)
#write.csv(df3,"GoodnessOfFit_SEM_Pseudopanurgus_didirupa.csv", row.names=FALSE)
#write.csv(df3,"GoodnessOfFit_SEM_Ceratina_nanula.csv", row.names=FALSE)
#write.csv(df3,"GoodnessOfFit_SEM_Hylaeus_annulatus.csv", row.names=FALSE)
#write.csv(df3,"GoodnessOfFit_SEM_Agapostemon_texanus.csv", row.names=FALSE)
#write.csv(df3,"GoodnessOfFit_SEM_Dufourea_harveyi.csv", row.names=FALSE)
#write.csv(df3,"GoodnessOfFit_SEM_Halictus_rubicundus.csv", row.names=FALSE)
#write.csv(df3,"GoodnessOfFit_SEM_Halictus_virgatellus.csv", row.names=FALSE)
#write.csv(df3,"GoodnessOfFit_SEM_Hoplitis_fulgida.csv", row.names=FALSE)
#write.csv(df3,"GoodnessOfFit_SEM_Hoplitis_robusta.csv", row.names=FALSE)
#write.csv(df3,"GoodnessOfFit_SEM_Osmia_albolateralis.csv", row.names=FALSE)

coefs_df<-bind_rows(data.frame(coefs(semFull),model=0),
           data.frame(coefs(sem1),model=1),
           data.frame(coefs(sem2),model=2),
           data.frame(coefs(sem3),model=3),
           data.frame(coefs(sem4),model=4), 
           data.frame(coefs(sem5),model=5),
           data.frame(coefs(sem6),model=6),
           data.frame(coefs(sem7),model=7),
           data.frame(coefs(sem8),model=8),
           data.frame(coefs(sem9),model=9),
           data.frame(coefs(sem10),model=10),
           data.frame(coefs(sem11),model=11),
           data.frame(coefs(sem12),model=12),
           data.frame(coefs(sem13),model=13),
           data.frame(coefs(sem14),model=14))

#write.csv(coefs_df,"coefficients_SEM_Panurginus_cressoniellus.csv", row.names=FALSE)
#write.csv(coefs_df,"coefficients_SEM_Panurginus_ineptus.csv", row.names=FALSE)
#write.csv(coefs_df,"coefficients_SEM_Pseudopanurgus_bakeri.csv", row.names=FALSE)
#write.csv(coefs_df,"coefficients_SEM_Pseudopanurgus_didirupa.csv", row.names=FALSE)
#write.csv(coefs_df,"coefficients_SEM_Ceratina_nanula.csv", row.names=FALSE)
#write.csv(coefs_df,"coefficients_SEM_Hylaeus_annulatus.csv", row.names=FALSE)
#write.csv(coefs_df,"coefficients_SEM_Agapostemon_texanus.csv", row.names=FALSE)
#write.csv(coefs_df,"coefficients_SEM_Dufourea_harveyi.csv", row.names=FALSE)
#write.csv(coefs_df,"coefficients_SEM_Halictus_rubicundus.csv", row.names=FALSE)
#write.csv(coefs_df,"coefficients_SEM_Halictus_virgatellus.csv", row.names=FALSE)
#write.csv(coefs_df,"coefficients_SEM_Hoplitis_fulgida.csv", row.names=FALSE)
#write.csv(coefs_df,"coefficients_SEM_Hoplitis_robusta.csv", row.names=FALSE)
#write.csv(coefs_df,"coefficients_SEM_Osmia_albolateralis.csv", row.names=FALSE)


summary(semFull, .progressBar = F) 
summary(sem1, .progressBar = F)  
summary(sem2, .progressBar = F)   
summary(sem3, .progressBar = F) 
summary(sem4, .progressBar = F)  
summary(sem5, .progressBar = F) 
summary(sem6, .progressBar = F)  
summary(sem7, .progressBar = F) 
summary(sem8, .progressBar = F) 
summary(sem9, .progressBar = F) 
summary(sem10, .progressBar = F) 
summary(sem11, .progressBar = F)  
summary(sem12, .progressBar = F) 
summary(sem13, .progressBar = F) 
summary(sem14, .progressBar = F) 


# perform model averaging, if necessary

# Pseudopanurgus bakeri
# sem 8, sem 13, sem 14
coef_list<-list(coef(summary(sem8[[1]]))[, 1], coef(summary(sem13[[1]]))[, 1],coef(summary(sem14[[1]]))[, 1])
weightslist<-c(summary(sem8, .progressBar = F)$AIC$AICc,summary(sem13, .progressBar = F)$AIC$AICc,summary(sem14, .progressBar = F)$AIC$AICc)
avgEst(coef_list)
avgEst(coef_list,weights=weightslist)

coef_list2<-list(coef(summary(sem8[[2]]))[, 1], coef(summary(sem13[[2]]))[, 1],coef(summary(sem14[[2]]))[, 1])
weightslist<-c(summary(sem8, .progressBar = F)$AIC$AICc,summary(sem13, .progressBar = F)$AIC$AICc,summary(sem14, .progressBar = F)$AIC$AICc)
avgEst(coef_list2)
avgEst(coef_list2,weights=weightslist)

coef_list3<-list(coef(summary(sem8[[3]]))[, 1], coef(summary(sem13[[3]]))[, 1],coef(summary(sem14[[3]]))[, 1])
weightslist<-c(summary(sem8, .progressBar = F)$AIC$AICc,summary(sem13, .progressBar = F)$AIC$AICc,summary(sem14, .progressBar = F)$AIC$AICc)
avgEst(coef_list3)
avgEst(coef_list3,weights=weightslist)

coefs(sem8)
coefs(sem13)
coefs(sem14)

# plot the SEM
est<-data.frame(avgEst(coef_list,weights=weightslist))
est <- tibble::rownames_to_column(est, "Predictor")
est<-est[-1,]
est$Response<-"sex_binom"
names(est)[2]<-"Estimate_MdlAvg"

est2<-data.frame(avgEst(coef_list2,weights=weightslist))
est2 <- tibble::rownames_to_column(est2, "Predictor")
est2<-est2[-1,]
est2$Response<-"prev_yr_floral_days_z"
names(est2)[2]<-"Estimate_MdlAvg"

est3<-data.frame(avgEst(coef_list3,weights=weightslist))
est3 <- tibble::rownames_to_column(est3, "Predictor")
est3<-est3[-1,]
est3$Response<-"prev_yr_total_flowers_z"
names(est3)[2]<-"Estimate_MdlAvg"

est_all<-bind_rows(est,est2,est3)


coef<-summary(semFull, .progressBar = F)$coef
coef
unique(coef$Response)
unique(coef$Predictor)

coef<-left_join(coef,est_all,by=c("Predictor","Response"))

coef$to<-coef$Response
coef$from<-coef$Predictor

coef$to[coef$to=="sex_binom"] <- 1
coef$to[coef$to=="prev_yr_floral_days_z"] <- 2
coef$to[coef$to=="prev_yr_total_flowers_z"] <- 3
coef$to[coef$to=="prev_yr_doy_bare_ground_z"] <- 4
coef$to[coef$to=="prev_yr_accum_summer_precip_cm_z"] <- 5

coef$from[coef$from=="sex_binom"] <- 1
coef$from[coef$from=="prev_yr_floral_days_z"] <- 2
coef$from[coef$from=="prev_yr_total_flowers_z"] <- 3
coef$from[coef$from=="prev_yr_doy_bare_ground_z"] <- 4
coef$from[coef$from=="prev_yr_accum_summer_precip_cm_z"] <- 5

coef$style<-ifelse(coef$Estimate_MdlAvg>0,"solid","dashed")
coef$color<-ifelse(coef$P.Value>=0.05,"gray","black")

ndf <-
  create_node_df(
    n = 5,
    fontsize=6,
    fixedsize=TRUE,
    color="black",
    fillcolor="lightgray",
    penwidth=0.5,
    shape="rectangle",
    width=0.7,
    height=0.3,
    label = c("Female/male \nratio","Floral days \n(prior year)", "Floral sum \n(prior year)","Snowmelt date \n(prior year)","Summer precip. \n(prior year)"))

graph <-
  create_graph(
    nodes_df = ndf)

render_graph(graph)

edf <-
  create_edge_df(
    from = coef$from,
    to = coef$to,
    rel = "leading_to",
    #label = round(coef$Estimate_MdlAvg,digits=2),
    penwidth = abs(coef$Estimate_MdlAvg*2),
    fontsize=5.5,
    fontcolor="brown3",
    #color=coef$color,
    color="black",
    style=coef$style)

graph <-
  create_graph(
    nodes_df = ndf,
    edges_df = edf)

render_graph(graph)

graph <-
  graph %>%
  set_node_position(
    node = 1, # sex ratio
    x = 3, y = 1.5) %>%
  set_node_position(
    node = 2, # floral days
    x = 2, y = 1) %>%
  set_node_position(
    node = 3, # floral sum
    x = 2, y = 2) %>%
  set_node_position(
    node = 4, # snowmelt date
    x = 1, y = 2) %>%
  set_node_position(
    node = 5, # summer precip
    x = 1, y = 1)

render_graph(graph)




