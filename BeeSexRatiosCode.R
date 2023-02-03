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

##### Add climate data #####
# read in billy barr climate data
climate_billy<-read.csv("weather_summaries_billybarr_forsexratios.csv")
# remove variables included in Ian's dataset
climate_billy<-climate_billy[,c(1,2)]

# read in Ian's temperature data
temp<-read.csv("bee_sites_temp_2006_2022.csv")
names(temp)[2:3]<-c("site","date")

### Calculate May-September degree days above 0 C
# calculate GDD
temp$GDD_calc<-(temp$Tmax+temp$Tmin)/2

# create column to use in calculating accumulated degree days above 0
temp$DD_abovezero<-ifelse(temp$GDD_calc>0, temp$GDD_calc, 0)

# create year, month, and day columns
temp$year<-year(temp$date)
temp$month<-month(temp$date)
temp$day<-day(temp$date)

# calculate May-September cumulative degree days above 0
climate<-temp %>% filter(month>=5 & month<=9) %>% group_by(site,year) %>% summarise(accum_summer_dd=sum(DD_abovezero))

# add snow data to climate data frame
snow<-read.csv("bee_sites_snow_2006_2022.csv")
names(snow)[2:4]<-c("site","year","doy_bare_ground")
climate<-left_join(climate,snow[,c(2:7)],by=c("site","year"))

# add billy barr presip data to climate data frame
climate<-left_join(climate,climate_billy,by="year")

# add present year's climate data to bee dataset
beedatafocal<-left_join(beedatafocal,climate,by=c("site","year"))

# add prior year's climate data to bee dataset
climate_prioryear<-climate
names(climate_prioryear)[3:8]<-c("prev_yr_accum_summer_dd","prev_yr_doy_bare_ground","prev_yr_SnowOnsetDOY","prev_yr_SnowLengthDays","prev_yr_SnowOnsetDOY_calendar","prev_yr_accum_summer_precip_cm")
#names(climate_prioryear)[2:8]<-c("prev_yr_accum_summer_precip_cm","prev_yr_accum_summer_dd", "prev_yr_accum_winter_snowfall_cm","prev_yr_mean_snow_depth_cm","prev_yr_snow_cover_days","prev_yr_date_bare_ground", "prev_yr_doy_bare_ground")
climate_prioryear$year<-climate_prioryear$year+1
beedatafocal<-left_join(beedatafocal,climate_prioryear,by=c("site","year"))

# z-score variables of interest
beedatafocal$accum_summer_precip_cm_z<-as.numeric(scale(beedatafocal$accum_summer_precip_cm, center = TRUE, scale = TRUE))
beedatafocal$doy_bare_ground_z<-as.numeric(scale(beedatafocal$doy_bare_ground, center = TRUE, scale = TRUE))
beedatafocal$prev_yr_accum_summer_precip_cm_z<-as.numeric(scale(beedatafocal$prev_yr_accum_summer_precip_cm, center = TRUE, scale = TRUE))
beedatafocal$prev_yr_doy_bare_ground_z<-as.numeric(scale(beedatafocal$prev_yr_doy_bare_ground, center = TRUE, scale = TRUE))

# add floral data
flowers<-read.csv("focalbeedata_envdata_RMBLsexratios_2022-11-07.csv")
beedatafocal2<-left_join(beedatafocal,flowers[,c(2,50:65)],by=c("unique_id_year"))

# create csv file
write.csv(beedatafocal2,"focalbeedata_envdata_RMBLsexratios_2023-01-26.csv",row.names=FALSE)


##### GLMMs: Across species, do climate variables together predict female/male counts? #####
beedatafocal<-read.csv("focalbeedata_envdata_RMBLsexratios_2023-01-26.csv")

# subset to just include focal sites
sitelist<-c("Almont","Almont Curve","Beaver","CDOT","Copper","Davids","Elko","Gothic","Hill","Little","Lypps","Mexican Cut","Rustlers","Seans","Tuttle","Willey")
beedatafocal<-filter(beedatafocal, site %in% sitelist)

# # subset to exclude sites not sampled prior to 2015
# sitelist<-c("Almont","Beaver","Copper","Davids","Gothic","Hill","Little","Mexican Cut","Seans","Tuttle","Willey")
# beedatafocal<-filter(beedatafocal, site %in% sitelist) 

# variables with low pairwise correlations to focus on:
# accum_summer_precip_cm
# doy_bare_ground
# prev_yr_accum_summer_precip_cm
# prev_yr_doy_bare_ground

# try simple models
ma<-glm(sex_binom ~ prev_yr_doy_bare_ground*site, data=beedatafocal,family="binomial")
summary(ma)
anova(ma, test="Chisq")

mb<-glm(sex_binom ~ prev_yr_doy_bare_ground+site+genus_species, data=beedatafocal,family="binomial")
summary(mb)
anova(mb, test="Chisq")

mc<-glm(sex_binom ~ prev_yr_doy_bare_ground*genus_species, data=beedatafocal,family="binomial")
summary(mc)
anova(mc, test="Chisq")


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

ggplot(beedatafocal2, aes(y = sex_binom, x = prev_yr_doy_bare_ground)) +
  geom_point() +
  theme_classic(base_size = 15) +
  facet_wrap(~elev_group) +
  geom_smooth(method = glm, method.args= list(family="binomial"))
m4 <- glmer(sex_binom ~ prev_yr_doy_bare_ground_z*elev_group + (1|year) + (1|genus_species), data = beedatafocal2, family = binomial)
summary(m4)
Anova(m4,type=3)

ggplot(beedatafocal2, aes(y = sex_binom, x = year)) +
  geom_point() +
  theme_classic(base_size = 15) +
  facet_wrap(~elev_group) +
  #facet_wrap(~site) +
  geom_smooth(method = glm, method.args= list(family="binomial"))

# construct models
# accum_summer_precip_cm_z
m1 <- glmer(sex_binom ~ accum_summer_precip_cm_z + (1|year/site) + (1|genus_species), data = beedatafocal, family = binomial)

# doy_bare_ground_z
m2 <- glmer(sex_binom ~ doy_bare_ground_z + (1|year/site) + (1|genus_species), data = beedatafocal, family = binomial)

# prev_yr_accum_summer_precip_cm_z
m3 <- glmer(sex_binom ~ prev_yr_accum_summer_precip_cm_z + (1|year/site) + (1|genus_species), data = beedatafocal, family = binomial)

# prev_yr_doy_bare_ground_z
#m4 <- glmer(sex_binom ~ prev_yr_doy_bare_ground_z + (1|year/site) + (1|genus_species), data = beedatafocal, family = binomial)

#other models examined to determine which random effects to include
m4a <- glmer(sex_binom ~ prev_yr_doy_bare_ground_z + (1|genus_species), data = beedatafocal, family = binomial)

#m4b <- glmer(sex_binom ~ prev_yr_doy_bare_ground_z + (1|year/site), data = beedatafocal, family = binomial)

#m4c <- glmer(sex_binom ~ prev_yr_doy_bare_ground_z + (1|year) + (1|genus_species), data = beedatafocal, family = binomial)

# USE THIS MODEL!!
m4d <- glmer(sex_binom ~ prev_yr_doy_bare_ground_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)
summary(m4d)


#m4e <- glmer(sex_binom ~ prev_yr_doy_bare_ground_z + (1|year), data = beedatafocal, family = binomial)

m4f <- glmer(sex_binom ~ prev_yr_doy_bare_ground_z + (1|site), data = beedatafocal, family = binomial)

#m4g <- glmer(sex_binom ~ prev_yr_doy_bare_ground_z + (1|site/year) + (1|genus_species), data = beedatafocal, family = binomial)

#
#m4h <- glmer(sex_binom ~ prev_yr_doy_bare_ground_z + (1|year) + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)

modelaiccs<-AICc(m4,m4a,m4b,m4c,m4d,m4e,m4f,m4g,m4h)
modelaiccs<-arrange(modelaiccs,AICc)
modelaiccs








# accum_summer_precip_cm_z + doy_bare_ground_z
m5 <- glmer(sex_binom ~ accum_summer_precip_cm_z + doy_bare_ground_z + (1|year/site) + (1|genus_species), data = beedatafocal, family = binomial)

# accum_summer_precip_cm_z + prev_yr_accum_summer_precip_cm_z
m6 <- glmer(sex_binom ~ accum_summer_precip_cm_z + prev_yr_accum_summer_precip_cm_z + (1|year/site) + (1|genus_species), data = beedatafocal, family = binomial)

# accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z
m7 <- glmer(sex_binom ~ accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z + (1|year/site) + (1|genus_species), data = beedatafocal, family = binomial)

# doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z
m8 <- glmer(sex_binom ~ doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z + (1|year/site) + (1|genus_species), data = beedatafocal, family = binomial)

# doy_bare_ground_z + prev_yr_doy_bare_ground_z
m9 <- glmer(sex_binom ~ doy_bare_ground_z + prev_yr_doy_bare_ground_z + (1|year/site) + (1|genus_species), data = beedatafocal, family = binomial)

# prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z
m10 <- glmer(sex_binom ~ prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z + (1|year/site) + (1|genus_species), data = beedatafocal, family = binomial)

# accum_summer_precip_cm_z + doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z
m11 <- glmer(sex_binom ~ accum_summer_precip_cm_z + doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z + (1|year/site) + (1|genus_species), data = beedatafocal, family = binomial)

# accum_summer_precip_cm_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z
m12 <- glmer(sex_binom ~ accum_summer_precip_cm_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z + (1|year/site) + (1|genus_species), data = beedatafocal, family = binomial)

# accum_summer_precip_cm_z + doy_bare_ground_z + prev_yr_doy_bare_ground_z
m13 <- glmer(sex_binom ~ accum_summer_precip_cm_z + doy_bare_ground_z + prev_yr_doy_bare_ground_z + (1|year/site) + (1|genus_species), data = beedatafocal, family = binomial)

# doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z
m14 <- glmer(sex_binom ~ doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z + (1|year/site) + (1|genus_species), data = beedatafocal, family = binomial)

# all climate variables
m15 <- glmer(sex_binom ~ doy_bare_ground_z + accum_summer_precip_cm_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z + (1|year/site) + (1|genus_species), data = beedatafocal, family = binomial)

aicc<-AICc(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15)
aicc<-arrange(aicc,AICc)
aicc$deltaAICc<-aicc$AICc-min(aicc$AICc)
aicc

summary(m4)
rsquared(m4)
vif(m4)

#     df     AICc  deltaAICc
# m2   5 9324.820  0.0000000  doy_bare_ground_z: sig. negative
summary(m2)
rsquared(m2)

# m8   6 9325.534  0.7141427  doy_bare_ground_z: sig. negative, prev_yr_accum_summer_precip_cm_z: nonsig. 
summary(m8)
rsquared(m8)
vif(m8) # ok

# m5   6 9326.077  1.2571301  doy_bare_ground_z: sig. negative, accum_summer_precip_cm_z: nonsig.
summary(m5)
rsquared(m5)
vif(m5) # ok

# m11  7 9326.148  1.3284050  doy_bare_ground_z: sig. negative, accum_summer_precip_cm_z: nonsig, prev_yr_accum_summer_precip_cm_z: nonsig.
summary(m11)
rsquared(m11)
vif(m11) # ok

# m10  6 9326.456  1.6360193  prev_yr_doy_bare_ground_z: sig. negative, prev_yr_accum_summer_precip_cm_z: marginally nonsig. positive
summary(m10)
rsquared(m10)
vif(m10) # ok

# m9   6 9326.764  1.9442547  doy_bare_ground_z: sig. negative, prev_yr_doy_bare_ground_z: nonsig
summary(m9)
rsquared(m9)
vif(m9) # ok

# graph results from best model
ggplot(beedatafocal, aes(y = sex_binom, x = doy_bare_ground_z)) +
  geom_point() +
  theme_classic(base_size = 15) +
  geom_smooth(method = glm, method.args= list(family="binomial")) +
  xlab("Snowmelt date (day of year)") + ylab("Sex \n(female=1, male=0)")





##### GLMMs: For each species individually, how does climate relate to female/male counts? #####
beedatafocal<-read.csv("focalbeedata_envdata_RMBLsexratios_2023-01-26.csv")

# subset to just include focal sites
sitelist<-c("Almont","Almont Curve","Beaver","CDOT","Copper","Davids","Elko","Gothic","Hill","Little","Lypps","Mexican Cut","Rustlers","Seans","Tuttle","Willey")
beedatafocal<-filter(beedatafocal, site %in% sitelist) 

species_list<-unique(beedatafocal$genus_species)
output<-matrix(nrow=1,ncol=21,byrow=TRUE,dimnames=list(c("row1"),c("model","df","AICc","pvalue1","pvalue2","pvalue3","pvalue4","slope1","slope2","slope3","slope4","se1","se2","se3","se4","delta_AICc","vif1","vif2","vif3","vif4","genus_species")))

for (i in 1:length(species_list)){
  data<-filter(beedatafocal, genus_species==species_list[i])
  
  # construct models
  # accum_summer_precip_cm_z
  m1 <- glmer(sex_binom ~ accum_summer_precip_cm_z + (1|year/site) , data = data, family = binomial)
  
  # doy_bare_ground_z
  m2 <- glmer(sex_binom ~ doy_bare_ground_z + (1|year/site) , data = data, family = binomial)
  
  # prev_yr_accum_summer_precip_cm_z
  m3 <- glmer(sex_binom ~ prev_yr_accum_summer_precip_cm_z + (1|year/site) , data = data, family = binomial)
  
  # prev_yr_doy_bare_ground_z
  m4 <- glmer(sex_binom ~ prev_yr_doy_bare_ground_z + (1|year/site) , data = data, family = binomial)
  
  # accum_summer_precip_cm_z + doy_bare_ground_z
  m5 <- glmer(sex_binom ~ accum_summer_precip_cm_z + doy_bare_ground_z + (1|year/site) , data = data, family = binomial)
  
  # accum_summer_precip_cm_z + prev_yr_accum_summer_precip_cm_z
  m6 <- glmer(sex_binom ~ accum_summer_precip_cm_z + prev_yr_accum_summer_precip_cm_z + (1|year/site) , data = data, family = binomial)
  
  # accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z
  m7 <- glmer(sex_binom ~ accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z + (1|year/site) , data = data, family = binomial)
  
  # doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z
  m8 <- glmer(sex_binom ~ doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z + (1|year/site) , data = data, family = binomial)
  
  # doy_bare_ground_z + prev_yr_doy_bare_ground_z
  m9 <- glmer(sex_binom ~ doy_bare_ground_z + prev_yr_doy_bare_ground_z + (1|year/site) , data = data, family = binomial)
  
  # prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z
  m10 <- glmer(sex_binom ~ prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z + (1|year/site) , data = data, family = binomial)
  
  # accum_summer_precip_cm_z + doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z
  m11 <- glmer(sex_binom ~ accum_summer_precip_cm_z + doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z + (1|year/site) , data = data, family = binomial)
  
  # accum_summer_precip_cm_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z
  m12 <- glmer(sex_binom ~ accum_summer_precip_cm_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z + (1|year/site) , data = data, family = binomial)
  
  # accum_summer_precip_cm_z + doy_bare_ground_z + prev_yr_doy_bare_ground_z
  m13 <- glmer(sex_binom ~ accum_summer_precip_cm_z + doy_bare_ground_z + prev_yr_doy_bare_ground_z + (1|year/site) , data = data, family = binomial)
  
  # doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z
  m14 <- glmer(sex_binom ~ doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z + (1|year/site) , data = data, family = binomial)
  
  # all climate variables
  m15 <- glmer(sex_binom ~ doy_bare_ground_z + accum_summer_precip_cm_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z + (1|year/site) , data = data, family = binomial)
  
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

#write.csv(output,"species_glmms_climate_sex_ratios_smallsitelist_2023-02-02.csv", row.names = FALSE)





##### For individual species: graphs of relationships with different climate variables #####
results<-read.csv("species_glmms_climate_sex_ratios_2023-02-01.csv")
#results<-read.csv("species_glmms_climate_sex_ratios_smallsitelist_2023-02-02.csv")
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


# look at species x site (and elevation) bee counts
summary<-beedatafocal %>% group_by(site,genus_species)%>%summarise(count=n())

siteinfo<-read.csv("site_info.csv")
summary<-left_join(summary,siteinfo,by="site")

ggplot(summary,aes(x=site,y=count)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  facet_wrap(~genus_species, scales="free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 

ggplot(summary,aes(x=elevation_m,y=count,color=site)) + 
  geom_bar(stat="identity") +
  facet_wrap(~genus_species, scales="free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 


# count ~ year for each species
summary2<-beedatafocal %>% group_by(site,genus_species,year)%>%summarise(count=n())
ggplot(summary2,aes(x=year,y=count,color=site)) + 
  geom_point() +
  facet_wrap(~genus_species, scales="free_y")

almont<-filter(summary2,site=="Almont")
ggplot(summary2,aes(x=year,y=count)) + 
  geom_point() +
  facet_wrap(~genus_species, scales="free_y")

summary3<-beedatafocal %>% group_by(site,year)%>%summarise(count=n())


# sex ~ doy bare ground for each elevation (aka site)
siteinfo<-read.csv("site_info.csv")
beedatafocal2<-left_join(beedatafocal,siteinfo,by="site")
havi<-filter(beedatafocal2,genus_species=="Halictus virgatellus")

m4 <- glmer(sex_binom ~ prev_yr_doy_bare_ground_z*elev_group + (1|year/site), data = havi, family = binomial)
summary(m4)

m4 <- glmer(sex_binom ~ prev_yr_doy_bare_ground_z*elev_group + (1|year), data = havi, family = binomial)
summary(m4)
Anova(m4, type=3)


ggplot(havi, aes(y = sex_binom, x = prev_yr_doy_bare_ground)) +
  geom_point() +
  theme_classic(base_size = 15) +
  #facet_wrap(~elev_group) +
  geom_smooth(method = glm, method.args= list(family="binomial"))




sitelist_mid<-c("Beaver","Copper","Davids","Gothic","Hill","Little","Rustlers","Seans","Tuttle","Willey")
havi2<-filter(havi, site %in% sitelist_mid) 
m4 <- glmer(sex_binom ~ prev_yr_doy_bare_ground_z + (1|year/site), data = havi2, family = binomial)
summary(m4)
ggplot(havi2, aes(y = sex_binom, x = prev_yr_doy_bare_ground)) +
  geom_point() +
  theme_classic(base_size = 15) +
  #facet_wrap(~as.factor(elevation_m)) +
  geom_smooth(method = glm, method.args= list(family="binomial"))


##### Adding floral data #####

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

##### Relationships between climate and floral variables #####

beedatafocal<-read.csv("focalbeedata_envdata_RMBLsexratios_2023-01-26.csv")

# subset to just include focal sites
sitelist<-c("Almont","Almont Curve","Beaver","CDOT","Copper","Davids","Elko","Gothic","Hill","Little","Lypps","Mexican Cut","Rustlers","Seans","Tuttle","Willey")
beedatafocal<-filter(beedatafocal, site %in% sitelist) 

# optional: remove outlier floral resource year (remove 2016 for examining previous year's resources)
#beedatafocal<-subset(beedatafocal,year!=2016)

# exploratory graphs: relationships between climate and floral variables

# subset to just include mid-elevation sites
# sitelist_mid<-c("Beaver","Copper","Davids","Gothic","Hill","Little","Rustlers","Seans","Tuttle","Willey")
# beedatafocal<-filter(beedatafocal, site %in% sitelist_mid) 

ggplot(beedatafocal, aes(y = total_flowers, x = doy_bare_ground)) +
  geom_point() +
  theme_classic(base_size = 15) +
  geom_smooth(method = lm) + 
  xlab("Snowmelt date") + ylab("Floral sum") + facet_wrap(~site)
summary(lm(total_flowers~doy_bare_ground,data=beedatafocal))

ggplot(beedatafocal, aes(y = total_flowers, x = accum_summer_precip_cm)) +
  geom_point() +
  theme_classic(base_size = 15) +
  geom_smooth(method = lm) + 
  xlab("Summer precipitation") + ylab("Floral sum") + facet_wrap(~site)
summary(lm(total_flowers~accum_summer_precip_cm,data=beedatafocal))

ggplot(beedatafocal, aes(y = floral_days, x = doy_bare_ground)) +
  geom_point() +
  theme_classic(base_size = 15) +
  geom_smooth(method = lm) + 
  xlab("Snowmelt date") + ylab("Floral days") + facet_wrap(~site)
summary(lm(floral_days~doy_bare_ground,data=beedatafocal))

ggplot(beedatafocal, aes(y = floral_days, x = accum_summer_precip_cm)) +
  geom_point() +
  theme_classic(base_size = 15) +
  geom_smooth(method = lm) + 
  xlab("Summer precipitation") + ylab("Floral days") + facet_wrap(~site)
summary(lm(floral_days~accum_summer_precip_cm,data=beedatafocal))


##### GLMMs effect of floral metrics on sex ratio #####

beedatafocal<-read.csv("focalbeedata_envdata_RMBLsexratios_2022-11-07.csv")

# just floral metrics

m1<- glmer(sex_binom ~ prev_yr_floral_days_z +
          (1|year/site) + (1|genus_species),
        data = beedatafocal, family = binomial)
summary(m1)
plot(allEffects(m1))

ggplot(beedatafocal, aes(y = sex_binom, x = prev_yr_floral_days_z)) +
  geom_point() +
  theme_classic(base_size = 15) +
  geom_smooth(method = glm, method.args= list(family="binomial")) +
  xlab("Previous year's floral days") + ylab("Sex \n(female=1, male=0)")


m2<- glmer(sex_binom ~ prev_yr_total_flowers_z +
          (1|year/site) + (1|genus_species),
        data = beedatafocal, family = binomial)
summary(m2)
plot(allEffects(m2))

ggplot(beedatafocal, aes(y = sex_binom, x = prev_yr_total_flowers_z)) +
  geom_point() +
  theme_classic(base_size = 15) +
  geom_smooth(method = glm, method.args= list(family="binomial")) +
  xlab("Previous year's total flowers") + ylab("Sex \n(female=1, male=0)")


# Analyses looking at floral days and climate together in GLMMs

# sex ~ floral days and snowmelt date
m1<- glmer(sex_binom ~ prev_yr_floral_days_z*prev_yr_doy_bare_ground_z +
             (1|year/site) + (1|genus_species),
           data = beedatafocal, family = binomial)
summary(m1)
vif(m1) # ok!
plot(allEffects(m1))
rsquared(m1)

visreg(m1,"prev_yr_floral_days_z",type="conditional", scale="response",
       by="prev_yr_doy_bare_ground_z",
       points.par=list(cex=1.2,col="black"),
       cex.axis=1.4, line.par=list(col="black"),
       xlab=list("Prior year's floral days", cex=1.8),
       ylab=list("Sex", cex=1.8))

# compare this model to the one that just contained snowmelt date
msnow<- glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
             (1|year/site) + (1|genus_species),
           data = beedatafocal, family = binomial)

AICc(m1,msnow)
# including floral resources improves the model!

beedatafocal$facet_prev_yr_doy_bare_ground_z<-ifelse(beedatafocal$prev_yr_doy_bare_ground_z<=0,"Early snowmelt \nin previous year","Late snowmelt \nin previous year")
ggplot(beedatafocal, aes(y = sex_binom, x = prev_yr_floral_days)) +
  geom_point() +
  theme_classic(base_size = 15) +
  geom_smooth(method = glm, method.args= list(family="binomial")) +
  facet_wrap(~facet_prev_yr_doy_bare_ground_z) +
  xlab("Previous year's floral days") + ylab("Sex \n(female=1, male=0)")


# sex ~ floral sum and snowmelt date - VIF is too high
m2<- glmer(sex_binom ~ prev_yr_total_flowers_z*prev_yr_doy_bare_ground_z +
             (1|year/site) + (1|genus_species),
           data = beedatafocal, family = binomial)
summary(m2)
vif(m2) # TOO HIGH!



##### SEMs: all species #####

beedatafocal<-read.csv("focalbeedata_envdata_RMBLsexratios_2022-11-07.csv")

### piecewise SEM: all prior year's climate and floral predictor variables ###
sem1<-psem(

  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_floral_days_z +
          prev_yr_total_flowers_z +
          (1|year/site) + (1|genus_species),
        data = beedatafocal, family = binomial),

  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),

  lm(prev_yr_total_flowers_z ~ prev_yr_floral_days_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)

)

summary(sem1, .progressBar = F)


# plotting the SEM! resources:
# https://rdrr.io/cran/DiagrammeR/man/create_graph.html
# also: https://rich-iannone.github.io/DiagrammeR/graphs.html

coef<-summary(sem1, .progressBar = F)$coef
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




# removing floral sum

sem2<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_floral_days_z +
          (1|year/site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem2, .progressBar = F)

coef<-summary(sem2, .progressBar = F)$coef
coef
unique(coef$Response)
unique(coef$Predictor)

coef$to<-coef$Response
coef$from<-coef$Predictor

coef$to[coef$to=="sex_binom"] <- 1
coef$to[coef$to=="prev_yr_floral_days_z"] <- 2
coef$to[coef$to=="prev_yr_doy_bare_ground_z"] <- 3
coef$to[coef$to=="prev_yr_accum_summer_precip_cm_z"] <- 4

coef$from[coef$from=="sex_binom"] <- 1
coef$from[coef$from=="prev_yr_floral_days_z"] <- 2
coef$from[coef$from=="prev_yr_doy_bare_ground_z"] <- 3
coef$from[coef$from=="prev_yr_accum_summer_precip_cm_z"] <- 4

coef$style<-ifelse(coef$P.Value>=0.05,"dashed","solid")

ndf <-
  create_node_df(
    n = 4,
    fontsize=6,
    fixedsize=TRUE,
    color="black",
    fillcolor="lightgray",
    penwidth=0.5,
    shape="ellipse",
    width=0.6,
    label = c("Sex \nratio","Floral \ndays", "Snowmelt \ndate","Summer \nprecipitation"))

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
    color="black",
    style=coef$style)

graph <-
  create_graph(
    nodes_df = ndf,
    edges_df = edf)

render_graph(graph)




# removing precip

sem3<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_floral_days_z +
          (1|year/site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~
       prev_yr_doy_bare_ground_z,
     data = beedatafocal)
  
)

summary(sem3, .progressBar = F)

coef<-summary(sem3, .progressBar = F)$coef
coef
unique(coef$Response)
unique(coef$Predictor)

coef$to<-coef$Response
coef$from<-coef$Predictor

coef$to[coef$to=="sex_binom"] <- 1
coef$to[coef$to=="prev_yr_floral_days_z"] <- 2
coef$to[coef$to=="prev_yr_doy_bare_ground_z"] <- 3

coef$from[coef$from=="sex_binom"] <- 1
coef$from[coef$from=="prev_yr_floral_days_z"] <- 2
coef$from[coef$from=="prev_yr_doy_bare_ground_z"] <- 3

coef$style<-ifelse(coef$Estimate>0,"solid","dashed")
coef$color<-ifelse(coef$P.Value>=0.05,"gray","black")

ndf <-
  create_node_df(
    n = 3,
    fontsize=6,
    fixedsize=TRUE,
    color="black",
    fillcolor="lightgray",
    penwidth=0.5,
    shape="rectangle",
    width=0.7,
    height=0.3,
    label = c("Female/male \nratio","Floral days \n(prior year)", "Snowmelt date \n(prior year)"))

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
    x = 1, y = 2) %>%
  set_node_position(
    node = 3, # snowmelt date
    x = 1, y = 1)

render_graph(graph)





# ### piecewise SEM: version that includes present year's climate ###
# sem_all<-psem(
#   
#   # effects of prior year's climate on sex ratios
#   glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
#           prev_yr_accum_summer_precip_cm_z +
#           prev_yr_floral_days_z +
#           prev_yr_total_flowers_z +
#           doy_bare_ground_z +
#           accum_summer_precip_cm_z +
#           floral_days_z +
#           total_flowers_z +
#           (1|year/site) + (1|genus_species),
#         data = beedatafocal, family = binomial),
#   
#   # effects of prior year's climate on flowers
#   lm(prev_yr_floral_days_z ~
#        prev_yr_doy_bare_ground_z +
#        prev_yr_accum_summer_precip_cm_z,
#      data = beedatafocal),
#   
#   lm(prev_yr_total_flowers_z ~ prev_yr_floral_days_z +
#        prev_yr_doy_bare_ground_z +
#        prev_yr_accum_summer_precip_cm_z,
#      data = beedatafocal),
#   
#   # effects of current year's climate on flowers
#   lm(floral_days_z ~ prev_yr_floral_days_z +
#        prev_yr_doy_bare_ground_z +
#        prev_yr_accum_summer_precip_cm_z +
#        doy_bare_ground_z +
#        accum_summer_precip_cm_z,
#      data = beedatafocal),
#   
#   lm(total_flowers_z ~ floral_days_z + prev_yr_floral_days_z + prev_yr_total_flowers_z +
#        prev_yr_doy_bare_ground_z +
#        prev_yr_accum_summer_precip_cm_z +
#        doy_bare_ground_z +
#        accum_summer_precip_cm_z,
#      data = beedatafocal)
#   
# )
# 
# summary(sem_all, .progressBar = F)
# 
# 
# 
# 
# 
# sem_all2<-psem(
#   
#   # effects of prior year's climate on sex ratios
#   glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
#           prev_yr_accum_summer_precip_cm_z +
#           prev_yr_floral_days_z +
#           doy_bare_ground_z +
#           accum_summer_precip_cm_z +
#           total_flowers_z +
#           (1|year/site) + (1|genus_species),
#         data = beedatafocal, family = binomial),
#   
#   # effects of prior year's climate on flowers
#   lm(prev_yr_floral_days_z ~
#        prev_yr_doy_bare_ground_z +
#        prev_yr_accum_summer_precip_cm_z,
#      data = beedatafocal),
#   
#   # effects of current year's climate on flowers
#   lm(total_flowers_z ~ prev_yr_floral_days_z +
#        prev_yr_doy_bare_ground_z +
#        prev_yr_accum_summer_precip_cm_z +
#        doy_bare_ground_z +
#        accum_summer_precip_cm_z,
#      data = beedatafocal)
#   
# )
# 
# summary(sem_all2, .progressBar = F)





##### SEMs: individual species #####

results<-read.csv("species_glmms_climate_sex_ratios.csv")
bestmodels<-filter(results,delta_AICc==0)

### species set 1 ###

filter(bestmodels,genus_species=="Halictus rubicundus")$predictors
focal<-filter(beedatafocal,genus_species=="Halictus rubicundus")

filter(bestmodels,genus_species=="Halictus virgatellus")$predictors
focal<-filter(beedatafocal,genus_species=="Halictus virgatellus")

filter(bestmodels,genus_species=="Panurginus cressoniellus")$predictors
focal<-filter(beedatafocal,genus_species=="Panurginus cressoniellus")

filter(bestmodels,genus_species=="Dufourea harveyi")$predictors
focal<-filter(beedatafocal,genus_species=="Dufourea harveyi")

sem1<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_floral_days_z +
          prev_yr_total_flowers_z +
          (1|year/site),
        data = focal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = focal),
  
  lm(prev_yr_total_flowers_z ~ prev_yr_floral_days_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = focal)
  
)

summary(sem1, .progressBar = F)


coef<-summary(sem1, .progressBar = F)$coef
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
    width=0.5,
    height=0.3,
    label = c("Sex \nratio","Floral \ndays", "Floral \nsum","Snowmelt \ndate","Summer \nprecipitation"))

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


# removing floral sum

sem2<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_floral_days_z +
          (1|year/site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem2, .progressBar = F)

coef<-summary(sem2, .progressBar = F)$coef
coef
unique(coef$Response)
unique(coef$Predictor)

coef$to<-coef$Response
coef$from<-coef$Predictor

coef$to[coef$to=="sex_binom"] <- 1
coef$to[coef$to=="prev_yr_floral_days_z"] <- 2
coef$to[coef$to=="prev_yr_doy_bare_ground_z"] <- 3
coef$to[coef$to=="prev_yr_accum_summer_precip_cm_z"] <- 4

coef$from[coef$from=="sex_binom"] <- 1
coef$from[coef$from=="prev_yr_floral_days_z"] <- 2
coef$from[coef$from=="prev_yr_doy_bare_ground_z"] <- 3
coef$from[coef$from=="prev_yr_accum_summer_precip_cm_z"] <- 4

coef$style<-ifelse(coef$P.Value>=0.05,"dashed","solid")

ndf <-
  create_node_df(
    n = 4,
    fontsize=6,
    fixedsize=TRUE,
    color="black",
    fillcolor="lightgray",
    penwidth=0.5,
    shape="ellipse",
    width=0.6,
    label = c("Sex \nratio","Floral \ndays", "Snowmelt \ndate","Summer \nprecipitation"))

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
    color="black",
    style=coef$style)

graph <-
  create_graph(
    nodes_df = ndf,
    edges_df = edf)

render_graph(graph)


### species set 2 ###

filter(bestmodels,genus_species=="Dufourea harveyi")$predictors
focal<-filter(beedatafocal,genus_species=="Dufourea harveyi")

sem1<-psem(
  
  # effects of climate on sex ratios
  # glmer(sex_binom ~ prev_yr_accum_summer_precip_cm_z +
  #         prev_yr_floral_days_z +
  #         (1|year/site),
  #       data = focal, family = binomial),
  glmer(sex_binom ~ prev_yr_accum_summer_precip_cm_z +
          prev_yr_total_flowers_z +
          (1|year/site),
        data = focal, family = binomial),
  
  # effects of climate on flowers
  # lm(prev_yr_floral_days_z ~
  #      prev_yr_accum_summer_precip_cm_z,
  #    data = focal)
  lm(prev_yr_total_flowers_z ~
       prev_yr_accum_summer_precip_cm_z,
     data = focal)
  
)

summary(sem1, .progressBar = F)

coef<-summary(sem1, .progressBar = F)$coef
coef
unique(coef$Response)
unique(coef$Predictor)

coef$to<-coef$Response
coef$from<-coef$Predictor

coef$to[coef$to=="sex_binom"] <- 1
coef$to[coef$to=="prev_yr_total_flowers_z"] <- 2
coef$to[coef$to=="prev_yr_accum_summer_precip_cm_z"] <- 3

coef$from[coef$from=="sex_binom"] <- 1
coef$from[coef$from=="prev_yr_total_flowers_z"] <- 2
coef$from[coef$from=="prev_yr_accum_summer_precip_cm_z"] <- 3

coef$style<-ifelse(coef$Estimate>0,"solid","dashed")
coef$color<-ifelse(coef$P.Value>=0.05,"gray","black")

ndf <-
  create_node_df(
    n = 3,
    fontsize=6,
    fixedsize=TRUE,
    color="black",
    fillcolor="lightgray",
    penwidth=0.5,
    shape="rectangle",
    width=0.5,
    height=0.3,
    label = c("Sex \nratio","Floral \nsum","Summer \nprecipitation"))
# label = c("sex_binom","prev_yr_floral_days_z", "prev_yr_total_flowers_z","prev_yr_doy_bare_ground_z","prev_yr_accum_summer_precip_cm_z"))

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


### species set 3 ###

filter(bestmodels,genus_species=="Pseudopanurgus bakeri")$predictors
focal<-filter(beedatafocal,genus_species=="Pseudopanurgus bakeri")

filter(bestmodels,genus_species=="Pseudopanurgus didirupa")$predictors
focal<-filter(beedatafocal,genus_species=="Pseudopanurgus didirupa")

sem1<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ doy_bare_ground_z +
          floral_days_z +
          total_flowers_z +
          (1|year/site),
        data = focal, family = binomial),
  
  # effects of climate on flowers
  lm(floral_days_z ~ total_flowers_z +
       doy_bare_ground_z,
     data = focal),
  
  # effects of climate on flowers
  lm(total_flowers_z ~
       doy_bare_ground_z,
     data = focal)
  
)

summary(sem1, .progressBar = F)

coef<-summary(sem1, .progressBar = F)$coef
coef
unique(coef$Response)
unique(coef$Predictor)

coef$to<-coef$Response
coef$from<-coef$Predictor

coef$to[coef$to=="sex_binom"] <- 1
coef$to[coef$to=="floral_days_z"] <- 2
coef$to[coef$to=="total_flowers_z"] <- 3
coef$to[coef$to=="doy_bare_ground_z"] <- 4

coef$from[coef$from=="sex_binom"] <- 1
coef$from[coef$from=="floral_days_z"] <- 2
coef$from[coef$from=="total_flowers_z"] <- 3
coef$from[coef$from=="doy_bare_ground_z"] <- 4

coef$style<-ifelse(coef$Estimate>0,"solid","dashed")
coef$color<-ifelse(coef$P.Value>=0.05,"gray","black")

ndf <-
  create_node_df(
    n = 4,
    fontsize=6,
    fixedsize=TRUE,
    color="black",
    fillcolor="lightgray",
    penwidth=0.5,
    shape="rectangle",
    width=0.5,
    height=0.3,
    label = c("Sex \nratio","Floral \ndays","Floral \nsum","Snowmelt date"))
# label = c("sex_binom","prev_yr_floral_days_z", "prev_yr_total_flowers_z","prev_yr_doy_bare_ground_z","prev_yr_accum_summer_precip_cm_z"))

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

### species set 4 ###

filter(bestmodels,genus_species=="Agapostemon texanus")$predictors
focal<-filter(beedatafocal,genus_species=="Agapostemon texanus")

filter(bestmodels,genus_species=="Panurginus ineptus")$predictors
focal<-filter(beedatafocal,genus_species=="Panurginus ineptus")

sem1<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_floral_days_z +
          prev_yr_total_flowers_z +
          (1|year/site),
        data = focal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = focal),
  
  lm(prev_yr_total_flowers_z ~ prev_yr_floral_days_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = focal)
  
)

summary(sem1, .progressBar = F)


coef<-summary(sem1, .progressBar = F)$coef
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
    width=0.5,
    height=0.3,
    label = c("Sex \nratio","Floral \ndays", "Floral \nsum","Snowmelt \ndate","Summer \nprecipitation"))

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




sem1<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ doy_bare_ground_z +
          accum_summer_precip_cm_z +
          floral_days_z +
          total_flowers_z +
          (1|year/site),
        data = focal, family = binomial),
  
  # effects of climate on flowers
  lm(floral_days_z ~
       doy_bare_ground_z +
       accum_summer_precip_cm_z,
     data = focal),
  
  lm(total_flowers_z ~ floral_days_z +
       doy_bare_ground_z +
       accum_summer_precip_cm_z,
     data = focal)
  
)

summary(sem1, .progressBar = F)


coef<-summary(sem1, .progressBar = F)$coef
coef
unique(coef$Response)
unique(coef$Predictor)

coef$to<-coef$Response
coef$from<-coef$Predictor

coef$to[coef$to=="sex_binom"] <- 1
coef$to[coef$to=="floral_days_z"] <- 2
coef$to[coef$to=="total_flowers_z"] <- 3
coef$to[coef$to=="doy_bare_ground_z"] <- 4
coef$to[coef$to=="accum_summer_precip_cm_z"] <- 5

coef$from[coef$from=="sex_binom"] <- 1
coef$from[coef$from=="floral_days_z"] <- 2
coef$from[coef$from=="total_flowers_z"] <- 3
coef$from[coef$from=="doy_bare_ground_z"] <- 4
coef$from[coef$from=="accum_summer_precip_cm_z"] <- 5

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
    width=0.5,
    height=0.3,
    label = c("Sex \nratio","Floral \ndays", "Floral \nsum","Snowmelt \ndate","Summer \nprecipitation"))

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


### species set 5 ###

filter(bestmodels,genus_species=="Hoplitis fulgida")$predictors
focal<-filter(beedatafocal,genus_species=="Hoplitis fulgida")

filter(bestmodels,genus_species=="Hoplitis robusta")$predictors
focal<-filter(beedatafocal,genus_species=="Hoplitis robusta")

sem1<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_floral_days_z +
          prev_yr_total_flowers_z +
          (1|year/site),
        data = focal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = focal),
  
  lm(prev_yr_total_flowers_z ~ prev_yr_floral_days_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = focal)
  
)

summary(sem1, .progressBar = F)


coef<-summary(sem1, .progressBar = F)$coef
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
    width=0.5,
    height=0.3,
    label = c("Sex \nratio","Floral \ndays", "Floral \nsum","Snowmelt \ndate","Summer \nprecipitation"))

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




sem1<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ doy_bare_ground_z +
          accum_summer_precip_cm_z +
          floral_days_z +
          total_flowers_z +
          (1|year/site),
        data = focal, family = binomial),
  
  # effects of climate on flowers
  lm(floral_days_z ~
       doy_bare_ground_z +
       accum_summer_precip_cm_z,
     data = focal),
  
  lm(total_flowers_z ~ floral_days_z +
       doy_bare_ground_z +
       accum_summer_precip_cm_z,
     data = focal)
  
)

summary(sem1, .progressBar = F)


coef<-summary(sem1, .progressBar = F)$coef
coef
unique(coef$Response)
unique(coef$Predictor)

coef$to<-coef$Response
coef$from<-coef$Predictor

coef$to[coef$to=="sex_binom"] <- 1
coef$to[coef$to=="floral_days_z"] <- 2
coef$to[coef$to=="total_flowers_z"] <- 3
coef$to[coef$to=="doy_bare_ground_z"] <- 4
coef$to[coef$to=="accum_summer_precip_cm_z"] <- 5

coef$from[coef$from=="sex_binom"] <- 1
coef$from[coef$from=="floral_days_z"] <- 2
coef$from[coef$from=="total_flowers_z"] <- 3
coef$from[coef$from=="doy_bare_ground_z"] <- 4
coef$from[coef$from=="accum_summer_precip_cm_z"] <- 5

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
    width=0.5,
    height=0.3,
    label = c("Sex \nratio","Floral \ndays", "Floral \nsum","Snowmelt \ndate","Summer \nprecipitation"))

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



##### Models involving floral_days_80 and total_flowers_80 variables #####

m1<- glmer(sex_binom ~ prev_yr_floral_days_80_z +
             (1|year/site) + (1|genus_species),
           data = beedatafocal, family = binomial)
summary(m1)
plot(allEffects(m1))

m2<- glmer(sex_binom ~ prev_yr_total_flowers_80_z +
             (1|year/site) + (1|genus_species),
           data = beedatafocal, family = binomial)
summary(m2)
plot(allEffects(m2))

m1<- glmer(sex_binom ~ prev_yr_floral_days_80_z*prev_yr_doy_bare_ground_z +
             (1|year/site) + (1|genus_species),
           data = beedatafocal, family = binomial)
summary(m1)
vif(m1) # ok!
plot(allEffects(m1))
rsquared(m1)



### piecewise SEM: all prior year's climate and floral predictor variables ###
sem1<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_floral_days_80_z +
          prev_yr_total_flowers_80_z +
          (1|year/site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lm(prev_yr_floral_days_80_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal),
  
  lm(prev_yr_total_flowers_80_z ~ prev_yr_floral_days_80_z +
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     data = beedatafocal)
  
)

summary(sem1, .progressBar = F)


# plotting the SEM! resources:
# https://rdrr.io/cran/DiagrammeR/man/create_graph.html
# also: https://rich-iannone.github.io/DiagrammeR/graphs.html

coef<-summary(sem1, .progressBar = F)$coef
coef
unique(coef$Response)
unique(coef$Predictor)

coef$to<-coef$Response
coef$from<-coef$Predictor

coef$to[coef$to=="sex_binom"] <- 1
coef$to[coef$to=="prev_yr_floral_days_80_z"] <- 2
coef$to[coef$to=="prev_yr_total_flowers_80_z"] <- 3
coef$to[coef$to=="prev_yr_doy_bare_ground_z"] <- 4
coef$to[coef$to=="prev_yr_accum_summer_precip_cm_z"] <- 5

coef$from[coef$from=="sex_binom"] <- 1
coef$from[coef$from=="prev_yr_floral_days_80_z"] <- 2
coef$from[coef$from=="prev_yr_total_flowers_80_z"] <- 3
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




