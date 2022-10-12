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
library(corrplot)
library(car)
library(MuMIn)

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
climate<-read.csv("weather_summaries_billybarr_forsexratios.csv")
# add doy for bare ground
climate$doy_bare_ground<-yday(climate$date_bare_ground)
# add present year's climate data to bee dataset
beedatafocal<-left_join(beedatafocal,climate,by="year")
# add prior year's climate data to bee dataset
climate_prioryear<-climate
names(climate_prioryear)[2:8]<-c("prev_yr_accum_summer_precip_cm","prev_yr_accum_summer_dd", "prev_yr_accum_winter_snowfall_cm","prev_yr_mean_snow_depth_cm","prev_yr_snow_cover_days","prev_yr_date_bare_ground", "prev_yr_doy_bare_ground")
climate_prioryear$year<-climate_prioryear$year+1
beedatafocal<-left_join(beedatafocal,climate_prioryear,by="year")


##### Look at correlations between climate variables #####
correlations_climate<-cor(beedatafocal[,c(32:36,38:43,45)])
corrplot(correlations_climate)
# variables with low pairwise correlations to focus on:

##### GLMMs: Across species, how do female/male counts vary with individual climate? (lme4 package) #####

# previous year's snowmelt date -- positive relationship
m1 <- glmer(sex_binom ~ prev_yr_doy_bare_ground + (1|year/site) + (1|genus_species), data = beedatafocal, family = binomial)

summary(m1)
print(m1, corr = FALSE)
plot(allEffects(m1))

ggplot(beedatafocal, aes(y = sex_binom, x = prev_yr_doy_bare_ground)) +
  geom_point() +
  theme_classic(base_size = 15) +
  geom_smooth(method = glm, method.args= list(family="binomial"))

# previous year's summer precipitation -- nonsignificant
m2 <- glmer(sex_binom ~ prev_yr_accum_summer_precip_cm + (1|year/site) + (1|genus_species), data = beedatafocal, family = binomial)

summary(m2)
print(m2, corr = FALSE)
plot(allEffects(m2))

ggplot(beedatafocal, aes(y = sex_binom, x = prev_yr_accum_summer_precip_cm)) +
  geom_point() +
  theme_classic(base_size = 15) +
  geom_smooth(method = glm, method.args= list(family="binomial"))

# previous year's summer temperature -- nonsignificant
m3 <- glmer(sex_binom ~ prev_yr_accum_summer_dd + (1|year/site) + (1|genus_species), data = beedatafocal, family = binomial)

summary(m3)
print(m3, corr = FALSE)
plot(allEffects(m3))

ggplot(beedatafocal, aes(y = sex_binom, x = prev_yr_accum_summer_dd)) +
  geom_point() +
  theme_classic(base_size = 15) +
  geom_smooth(method = glm, method.args= list(family="binomial"))


# current year's winter snow cover (winter preceding current growing season) -- nonsignificant
m4 <- glmer(sex_binom ~ accum_winter_snowfall_cm + (1|year/site) + (1|genus_species), data = beedatafocal, family = binomial)

summary(m4)
print(m4, corr = FALSE)
plot(allEffects(m4))

ggplot(beedatafocal, aes(y = sex_binom, x = accum_winter_snowfall_cm)) +
  geom_point() +
  theme_classic(base_size = 15) +
  geom_smooth(method = glm, method.args= list(family="binomial"))


# current year's snowmelt date -- nonsignificant
m5 <- glmer(sex_binom ~ doy_bare_ground + (1|year/site) + (1|genus_species), data = beedatafocal, family = binomial)

summary(m5)
print(m5, corr = FALSE)
plot(allEffects(m5))

ggplot(beedatafocal, aes(y = sex_binom, x = doy_bare_ground)) +
  geom_point() +
  theme_classic(base_size = 15) +
  geom_smooth(method = glm, method.args= list(family="binomial"))


# current year's summer precipitation -- nonsignificant
m6 <- glmer(sex_binom ~ accum_summer_precip_cm + (1|year/site) + (1|genus_species), data = beedatafocal, family = binomial)

summary(m6)
print(m6, corr = FALSE)
plot(allEffects(m6))

ggplot(beedatafocal, aes(y = sex_binom, x = accum_summer_precip_cm)) +
  geom_point() +
  theme_classic(base_size = 15) +
  geom_smooth(method = glm, method.args= list(family="binomial"))

# current year's summer temperature -- nonsignificant
m7 <- glmer(sex_binom ~ accum_summer_dd + (1|year/site) + (1|genus_species), data = beedatafocal, family = binomial)

summary(m7)
print(m7, corr = FALSE)
plot(allEffects(m7))

ggplot(beedatafocal, aes(y = sex_binom, x = accum_summer_dd)) +
  geom_point() +
  theme_classic(base_size = 15) +
  geom_smooth(method = glm, method.args= list(family="binomial"))




##### GLMMs: Across species, do climate variables together predict female/male counts? #####

# variables with low pairwise correlations to focus on:
# accum_summer_precip_cm
# doy_bare_ground
# prev_yr_accum_summer_precip_cm
# prev_yr_doy_bare_ground

# z-score variables
beedatafocal$accum_summer_precip_cm_z<-scale(beedatafocal$accum_summer_precip_cm, center = TRUE, scale = TRUE)

beedatafocal$doy_bare_ground_z<-scale(beedatafocal$doy_bare_ground, center = TRUE, scale = TRUE)

beedatafocal$prev_yr_accum_summer_precip_cm_z<-scale(beedatafocal$prev_yr_accum_summer_precip_cm, center = TRUE, scale = TRUE)

beedatafocal$prev_yr_doy_bare_ground_z<-scale(beedatafocal$prev_yr_doy_bare_ground, center = TRUE, scale = TRUE)










# construct models
# accum_summer_precip_cm_z
m1 <- glmer(sex_binom ~ accum_summer_precip_cm_z + (1|year/site) + (1|genus_species), data = beedatafocal, family = binomial)

# doy_bare_ground_z
m2 <- glmer(sex_binom ~ doy_bare_ground_z + (1|year/site) + (1|genus_species), data = beedatafocal, family = binomial)

# prev_yr_accum_summer_precip_cm_z
m3 <- glmer(sex_binom ~ prev_yr_accum_summer_precip_cm_z + (1|year/site) + (1|genus_species), data = beedatafocal, family = binomial)

# prev_yr_doy_bare_ground_z
m4 <- glmer(sex_binom ~ prev_yr_doy_bare_ground_z + (1|year/site) + (1|genus_species), data = beedatafocal, family = binomial)

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

aicc<-AICc(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14)
aicc<-arrange(aicc,AICc)
aicc$deltaAICc<-aicc$AICc-AICc(m4)

# df     AICc  deltaAICc         Model
# m4   5 12089.19  0.0000000     prev_yr_doy_bare_ground_z
summary(m4)
# m7   6 12089.58  0.3859968     accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z 
summary(m7)
vif(m7) # ok
# m13  7 12090.74  1.5503916     accum_summer_precip_cm_z + doy_bare_ground_z + prev_yr_doy_bare_ground_z
summary(m13)
vif(m13) # ok
# m10  6 12091.15  1.9647074     prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z
summary(m10)
vif(m10) # ok
# m9   6 12091.16  1.9747274     doy_bare_ground_z + prev_yr_doy_bare_ground_z
summary(m9)
vif(m9) # ok
# m12  7 12091.50  2.3117760     accum_summer_precip_cm_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z
summary(m12)
vif(m12) # ok
# m14  7 12093.11  3.9231739     doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z +  prev_yr_doy_bare_ground_z
summary(m14)
vif(m14) # ok
# m3   5 12096.74  7.5477262
# m2   5 12096.76  7.5683016
# m8   6 12097.81  8.6201019
# m6   6 12098.29  9.1022211
# m1   5 12098.63  9.4375146
# m5   6 12098.76  9.5663272
# m11  7 12099.67 10.4777415




#### GLMMs: For each species individually, how does climate relate to female/male counts? #####
species_list<-unique(beedatafocal$genus_species)
output<-matrix(nrow=1,ncol=17,byrow=TRUE,dimnames=list(c("row1"),c("model","df","AICc","pvalue1","pvalue2","pvalue3","slope1","slope2","slope3","se1","se2","se3","delta_AICc","vif1","vif2","vif3","genus_species")))

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
  
  aicc<-AICc(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14)
  aicc<- tibble::rownames_to_column(aicc, "model")
  aicc$pvalue1<-NA
  aicc$pvalue2<-NA
  aicc$pvalue3<-NA
  aicc$slope1<-NA
  aicc$slope2<-NA
  aicc$slope3<-NA
  aicc$se1<-NA
  aicc$se2<-NA
  aicc$se3<-NA
  aicc$vif1<-NA
  aicc$vif2<-NA
  aicc$vif3<-NA
  
  df<-data.frame(summary(m1)[["coefficients"]])
  df<-df[-1,]
  aicc[1,4:sum(length(df$Estimate)+3)]<-df[,4]
  aicc[1,7:sum(length(df$Estimate)+6)]<-df[,1]
  aicc[1,10:sum(length(df$Estimate)+9)]<-df[,2]
  
  df<-data.frame(summary(m2)[["coefficients"]])
  df<-df[-1,]
  aicc[2,4:sum(length(df$Estimate)+3)]<-df[,4]
  aicc[2,7:sum(length(df$Estimate)+6)]<-df[,1]
  aicc[2,10:sum(length(df$Estimate)+9)]<-df[,2]
  
  df<-data.frame(summary(m3)[["coefficients"]])
  df<-df[-1,]
  aicc[3,4:sum(length(df$Estimate)+3)]<-df[,4]
  aicc[3,7:sum(length(df$Estimate)+6)]<-df[,1]
  aicc[3,10:sum(length(df$Estimate)+9)]<-df[,2]
  
  df<-data.frame(summary(m4)[["coefficients"]])
  df<-df[-1,]
  aicc[4,4:sum(length(df$Estimate)+3)]<-df[,4]
  aicc[4,7:sum(length(df$Estimate)+6)]<-df[,1]
  aicc[4,10:sum(length(df$Estimate)+9)]<-df[,2]
  
  df<-data.frame(summary(m5)[["coefficients"]])
  df<-df[-1,]
  aicc[5,4:sum(length(df$Estimate)+3)]<-df[,4]
  aicc[5,7:sum(length(df$Estimate)+6)]<-df[,1]
  aicc[5,10:sum(length(df$Estimate)+9)]<-df[,2]
  aicc[5,13:sum(length(df$Estimate)+12)]<-data.frame(vif(m5))[,1]

  df<-data.frame(summary(m6)[["coefficients"]])
  df<-df[-1,]
  aicc[6,4:sum(length(df$Estimate)+3)]<-df[,4]
  aicc[6,7:sum(length(df$Estimate)+6)]<-df[,1]
  aicc[6,10:sum(length(df$Estimate)+9)]<-df[,2]
  aicc[6,13:sum(length(df$Estimate)+12)]<-data.frame(vif(m6))[,1]
  
  df<-data.frame(summary(m7)[["coefficients"]])
  df<-df[-1,]
  aicc[7,4:sum(length(df$Estimate)+3)]<-df[,4]
  aicc[7,7:sum(length(df$Estimate)+6)]<-df[,1]
  aicc[7,10:sum(length(df$Estimate)+9)]<-df[,2]
  aicc[7,13:sum(length(df$Estimate)+12)]<-data.frame(vif(m7))[,1]
  
  df<-data.frame(summary(m8)[["coefficients"]])
  df<-df[-1,]
  aicc[8,4:sum(length(df$Estimate)+3)]<-df[,4]
  aicc[8,7:sum(length(df$Estimate)+6)]<-df[,1]
  aicc[8,10:sum(length(df$Estimate)+9)]<-df[,2]
  aicc[8,13:sum(length(df$Estimate)+12)]<-data.frame(vif(m8))[,1]
  
  df<-data.frame(summary(m9)[["coefficients"]])
  df<-df[-1,]
  aicc[9,4:sum(length(df$Estimate)+3)]<-df[,4]
  aicc[9,7:sum(length(df$Estimate)+6)]<-df[,1]
  aicc[9,10:sum(length(df$Estimate)+9)]<-df[,2]
  aicc[9,13:sum(length(df$Estimate)+12)]<-data.frame(vif(m9))[,1]
  
  df<-data.frame(summary(m10)[["coefficients"]])
  df<-df[-1,]
  aicc[10,4:sum(length(df$Estimate)+3)]<-df[,4]
  aicc[10,7:sum(length(df$Estimate)+6)]<-df[,1]
  aicc[10,10:sum(length(df$Estimate)+9)]<-df[,2]
  aicc[10,13:sum(length(df$Estimate)+12)]<-data.frame(vif(m10))[,1]
  
  df<-data.frame(summary(m11)[["coefficients"]])
  df<-df[-1,]
  aicc[11,4:sum(length(df$Estimate)+3)]<-df[,4]
  aicc[11,7:sum(length(df$Estimate)+6)]<-df[,1]
  aicc[11,10:sum(length(df$Estimate)+9)]<-df[,2]
  aicc[11,13:sum(length(df$Estimate)+12)]<-data.frame(vif(m11))[,1]
  
  df<-data.frame(summary(m12)[["coefficients"]])
  df<-df[-1,]
  aicc[12,4:sum(length(df$Estimate)+3)]<-df[,4]
  aicc[12,7:sum(length(df$Estimate)+6)]<-df[,1]
  aicc[12,10:sum(length(df$Estimate)+9)]<-df[,2]
  aicc[12,13:sum(length(df$Estimate)+12)]<-data.frame(vif(m12))[,1]
  
  df<-data.frame(summary(m13)[["coefficients"]])
  df<-df[-1,]
  aicc[13,4:sum(length(df$Estimate)+3)]<-df[,4]
  aicc[13,7:sum(length(df$Estimate)+6)]<-df[,1]
  aicc[13,10:sum(length(df$Estimate)+9)]<-df[,2]
  aicc[13,13:sum(length(df$Estimate)+12)]<-data.frame(vif(m13))[,1]
  
  df<-data.frame(summary(m14)[["coefficients"]])
  df<-df[-1,]
  aicc[14,4:sum(length(df$Estimate)+3)]<-df[,4]
  aicc[14,7:sum(length(df$Estimate)+6)]<-df[,1]
  aicc[14,10:sum(length(df$Estimate)+9)]<-df[,2]
  aicc[14,13:sum(length(df$Estimate)+12)]<-data.frame(vif(m14))[,1]
  
  
  aicc<-arrange(aicc,AICc)
  aicc$delta_AICc<-aicc$AICc-min(aicc$AICc)
  aicc$genus_species<-species_list[i]

  output<-rbind(output,aicc)
}

output<-output[-1,]

# add column listing predictors in each model to "output" df
aicc<-AICc(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14)
aicc<- tibble::rownames_to_column(aicc, "model")
aicc$predictors<-c("accum_summer_precip_cm_z", "doy_bare_ground_z", "prev_yr_accum_summer_precip_cm_z", "prev_yr_doy_bare_ground_z", "accum_summer_precip_cm_z + doy_bare_ground_z", "accum_summer_precip_cm_z + prev_yr_accum_summer_precip_cm_z", "accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z", "doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z", "doy_bare_ground_z + prev_yr_doy_bare_ground_z", "prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z", "accum_summer_precip_cm_z + doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z", "accum_summer_precip_cm_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z", "accum_summer_precip_cm_z + doy_bare_ground_z + prev_yr_doy_bare_ground_z", "doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z")
aicc<-aicc[,-c(2:3)]
output<-left_join(output,aicc,by="model")

write.csv(output,"species_glmms_climate_sex_ratios.csv", row.names = FALSE)




