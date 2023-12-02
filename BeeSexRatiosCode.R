# ---------------------------------
# RMBL bee sex ratios and climate change
# Analyses using RMBL Phenology Project bee survey and flowering phenology data
# Melanie R. Kazenel 
# Created 7 October 2022
# ---------------------------------


# load relevant packages
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
library(DiagrammeRsvg)
library(rsvg)
library(nlme)
library(effects)
library(patchwork)


##### Read in and format data #####

# set working directory
setwd("~/Documents/*RMBLPhenology/Projects/BeeSexRatios/Analyses")

# read in file of bee, floral, and climate data
beedatafocal<-read.csv("focalbeedata_envdata_RMBLsexratios_2023-03-29.csv")

# filter data to just include focal sites
sitelist<-c("Almont","Almont Curve","Beaver","CDOT","Copper","Davids","Elko","Gothic","Hill","Little","Lypps","Mexican Cut","Rustlers","Seans","Tuttle","Willey")
beedatafocal<-filter(beedatafocal, site %in% sitelist)


##### GLMMs: Across species, how do sex ratios vary with climate? #####

# Focal climate variables used below (all z-scored):
# accum_summer_precip_cm_z (cumulative summer precipitation)
# doy_bare_ground_z (snowmelt date)
# prev_yr_accum_summer_precip_cm_z (previous year's cumulative summer precipitation)
# prev_yr_doy_bare_ground_z (previous year's snowmelt date)

### Construct models ###

# accum_summer_precip_cm_z
m1 <- glmer(sex_binom ~ accum_summer_precip_cm_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)
summary(m1)

# doy_bare_ground_z
m2 <- glmer(sex_binom ~ doy_bare_ground_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)

# prev_yr_accum_summer_precip_cm_z
m3 <- glmer(sex_binom ~ prev_yr_accum_summer_precip_cm_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)

# prev_yr_doy_bare_ground_z
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

# all climate variables 
m15 <- glmer(sex_binom ~ doy_bare_ground_z + accum_summer_precip_cm_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)

# get AICc values for each model
aicc<-AICc(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15)
# assign a number to each model
aicc$model<-1:15
# sort the data by AICc value
aicc<-arrange(aicc,AICc)
# calculate delta AICc
aicc$deltaAICc<-aicc$AICc-min(aicc$AICc)
aicc

# get r-squared values for each model
r2<-bind_rows(data.frame(rsquared(m1)),
data.frame(rsquared(m2)),
data.frame(rsquared(m3)),
data.frame(rsquared(m4)),
data.frame(rsquared(m5)),
data.frame(rsquared(m6)),
data.frame(rsquared(m7)),
data.frame(rsquared(m8)),
data.frame(rsquared(m9)),
data.frame(rsquared(m10)),
data.frame(rsquared(m11)),
data.frame(rsquared(m12)),
data.frame(rsquared(m13)),
data.frame(rsquared(m14)),
data.frame(rsquared(m15)))
r2$model<-1:15

# get summary stats for best model
summary(m15)
rsquared(m15)
vif(m15)
allEffects(m15)
plot(allEffects(m15))

# create a data frame of summary stats for each model
summary_df<-bind_rows(data.frame(summary(m1)$coefficients,model=1),
                    data.frame(summary(m2)$coefficients,model=2),
                    data.frame(summary(m3)$coefficients,model=3),
                    data.frame(summary(m4)$coefficients,model=4), 
                    data.frame(summary(m5)$coefficients,model=5),
                    data.frame(summary(m6)$coefficients,model=6),
                    data.frame(summary(m7)$coefficients,model=7),
                    data.frame(summary(m8)$coefficients,model=8),
                    data.frame(summary(m9)$coefficients,model=9),
                    data.frame(summary(m10)$coefficients,model=10),
                    data.frame(summary(m11)$coefficients,model=11),
                    data.frame(summary(m12)$coefficients,model=12),
                    data.frame(summary(m13)$coefficients,model=13),
                    data.frame(summary(m14)$coefficients,model=14),
                    data.frame(summary(m15)$coefficients,model=15))

summary_df$Term <- row.names(summary_df)  
summary_df<-left_join(summary_df,aicc,by="model")
summary_df<-left_join(summary_df,r2,by="model")

#write.csv(summary_df,"glmm_output_allspecies_2023_04_04.csv")

# get vif values for all relevant models
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
vif(m15)


##### Figure: Across species, how do sex ratios vary with climate? #####

# load required packages
library(visreg)
library(ggtext)
dev.off()

# plot sex ~ snowmelt date
summary_snowmelt<-beedatafocal %>% group_by(doy_bare_ground_z, sex_binom) %>% summarise(count=n())

p1a<-visreg(m15,"doy_bare_ground_z",type="conditional", scale="response", rug=FALSE, partial=FALSE, gg=TRUE, line=list(col="black")) 
p1<- p1a + xlab("Snowmelt date") + ylab("Proportion of females") +
  theme_linedraw(base_size = 14) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  annotate("text", x = 1, y = .4, label = expression(atop(italic("\u03B2")==-0.32, italic(P)*" < 0.0001"))) +
  ggtitle("a") +
  geom_point(aes(x=doy_bare_ground_z, y=sex_binom,size=count), data=summary_snowmelt) +
  labs(size="Count")
p1

# plot sex ~ summer precipitation
summary_precip<-beedatafocal %>% group_by(accum_summer_precip_cm_z, sex_binom) %>% summarise(count=n())

p2a<-visreg(m15,"accum_summer_precip_cm_z",type="conditional", scale="response", rug=FALSE, partial=FALSE, gg=TRUE, line=list(col="black"))
p2<- p2a + xlab("Summer precipitation") + ylab("Proportion of females") +
  theme_linedraw(base_size = 14) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  annotate("text", x = 1, y = .4, label = expression(atop(italic("\u03B2")==0.16, italic(P)*" < 0.0001"))) +
  ggtitle("b") +
  geom_point(aes(x=accum_summer_precip_cm_z, y=sex_binom,size=count), data=summary_precip) +
  labs(size="Count")
p2

# plot sex ~ prior year's snowmelt date
summary_snowmeltprior<-beedatafocal %>% group_by(prev_yr_doy_bare_ground_z, sex_binom) %>% summarise(count=n())

p3a<-visreg(m15,"prev_yr_doy_bare_ground_z",type="conditional", scale="response", rug=FALSE, partial=FALSE, gg=TRUE, line=list(col="black")) 
p3<-p3a + xlab("Prior year's \nsnowmelt date") + ylab("Proportion of females") +
  theme_linedraw(base_size = 14) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  annotate("text", x = 1, y = .25, label = expression(atop(italic("\u03B2")==0.58, italic(P)*" < 0.0001"))) +
  ggtitle("c") +
  geom_point(aes(x=prev_yr_doy_bare_ground_z, y=sex_binom,size=count), data=summary_snowmeltprior) +
  labs(size="Count")
p3

# plot sex ~ prior year's summer precipitation
summary_precipprior<-beedatafocal %>% group_by(prev_yr_accum_summer_precip_cm_z, sex_binom) %>% summarise(count=n())

p4a<-visreg(m15,"prev_yr_accum_summer_precip_cm_z",type="conditional", scale="response", rug=FALSE, partial=FALSE, gg=TRUE, line=list(col="black")) 
p4<- p4a + xlab("Prior year's \nsummer precipitation") + ylab("Proportion of females") +
  theme_linedraw(base_size = 14) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  annotate("text", x = 1, y = .4, label = expression(atop(italic("\u03B2")==-0.09, italic(P)*" < 0.0001"))) +
  ggtitle("d") +
  geom_point(aes(x=prev_yr_accum_summer_precip_cm_z, y=sex_binom,size=count), data=summary_precipprior) +
  labs(size="Count")
p4

# aggregate panels into a single figure
plot<-p1 + p2 + p3 + p4 + plot_layout(ncol = 2, guides = "collect")

# remove axis text and tick marks where they aren't needed
plot[[2]] = plot[[2]] + theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank() )
plot[[4]] = plot[[4]] + theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.y = element_blank() )
plot

#ggsave("glmm_allspecies_parameters_2023-12-02.png", plot,width=7,height=7,units = c("in"),dpi = 600)


##### GLMMs: For individual species, how do sex ratios vary with climate? #####

# get bee species list
species_list<-unique(beedatafocal$genus_species)

# create output data frame to hold results from for loop below
output<-matrix(nrow=1,ncol=21,byrow=TRUE,dimnames=list(c("row1"),c("model","df","AICc","pvalue1","pvalue2","pvalue3","pvalue4","slope1","slope2","slope3","slope4","se1","se2","se3","se4","delta_AICc","vif1","vif2","vif3","vif4","genus_species")))

# for each species, compare models of sex as a function of climate variables, and put the results in the "output" data frame
for (i in 1:length(species_list)){
  # get data for a given species
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
  
  # get AICc values from models
  aicc<-AICc(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15)
  aicc<- tibble::rownames_to_column(aicc, "model")
  
  # create columns for other model output in aicc data frame, and fill with NAs as placeholders 
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
  
  # for each model, add relevant output statistics to the aicc data frame
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
  
  # sort the aicc data frame by AICc value
  aicc<-arrange(aicc,AICc)
  # calculate delta AICc
  aicc$delta_AICc<-aicc$AICc-min(aicc$AICc)
  # add genus and species information to the data frame
  aicc$genus_species<-species_list[i]

  # add results for the given species to the output data frame
  output<-rbind(output,aicc)
}

# remove unnecessary row
output<-output[-1,]

# add column listing predictors in each model to the "output" data frame
aicc<-AICc(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15)
aicc<- tibble::rownames_to_column(aicc, "model")
aicc$predictors<-c("accum_summer_precip_cm_z", "doy_bare_ground_z", "prev_yr_accum_summer_precip_cm_z", "prev_yr_doy_bare_ground_z", "accum_summer_precip_cm_z + doy_bare_ground_z", "accum_summer_precip_cm_z + prev_yr_accum_summer_precip_cm_z", "accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z", "doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z", "doy_bare_ground_z + prev_yr_doy_bare_ground_z", "prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z", "accum_summer_precip_cm_z + doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z", "accum_summer_precip_cm_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z", "accum_summer_precip_cm_z + doy_bare_ground_z + prev_yr_doy_bare_ground_z", "doy_bare_ground_z + prev_yr_accum_summer_precip_cm_z + prev_yr_doy_bare_ground_z","all predictors")
aicc<-aicc[,-c(2:3)]
output<-left_join(output,aicc,by="model")

#write.csv(output,"species_glmms_climate_sex_ratios_2023-03-29.csv", row.names = FALSE)


##### Figure: For individual species, how do sex ratios vary with climate? #####

# read in results from individual species GLMMs
results<-read.csv("species_glmms_climate_sex_ratios_2023-03-29.csv")

# filter results to only include the best model for each species
bestmodels<-filter(results,delta_AICc==0)

# create long-form data frame of slopes (parameter estimates) from the best model for each species, with one row for each parameter estimate in the model
slopes_long<-pivot_longer(bestmodels[,c(21,1,22,8:11)],cols=4:7,values_to="slope")
names(slopes_long)[4]<-"term"
slopes_long$term<-substr(slopes_long$term,6,6) # create term column listing the numeric order of the parameter in the model

# create similar long-form dataset of se values
se_long<-pivot_longer(bestmodels[,c(21,1,22,12:15)],cols=4:7,values_to="se")
names(se_long)[4]<-"term"
se_long$term<-substr(se_long$term,3,3)

# create similar long-form dataset of p-values
pvalues_long<-pivot_longer(bestmodels[,c(21,1,22,4:7)],cols=4:7,values_to="pvalue")
names(pvalues_long)[4]<-"term"
pvalues_long$term<-substr(pvalues_long$term,7,7)

# combine the data frames above
slope_se_long<-left_join(slopes_long,se_long, by=c("genus_species","model","predictors","term"))
slope_se_long<-left_join(slope_se_long,pvalues_long, by=c("genus_species","model","predictors","term"))

# drop rows containing missing values
slope_se_long<-drop_na(slope_se_long)

# get summary of terms in each model
summary_terms<-slope_se_long %>% group_by(model, predictors, term) %>% summarise(count=n())
summary_terms$count<-NULL

# create a data frame listing the name of each model parameter and its term number (numeric order) in the model
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

# add term_id information to the data frame of parameter estimates and their statistics
slope_se_long<-left_join(slope_se_long, summary_terms, by=c("model","predictors","term"))

# rename levels of term_id for graphing purposes
slope_se_long$term_id_renamed<-as.factor(slope_se_long$term_id)
levels(slope_se_long$term_id_renamed) <- c("Precip.","Snowmelt","Prior year's precip.","Prior year's snowmelt")
slope_se_long$term_id_renamed <- factor(slope_se_long$term_id_renamed, levels=c("Snowmelt","Precip.","Prior year's snowmelt","Prior year's precip."))

# create column coding each p-value as significant or nonsignificant
slope_se_long$pvalue_sig<-ifelse(slope_se_long$pvalue>=0.05,"Nonsignificant","Significant")

# filter out bee species for which climate did not predict sex ratios
slope_se_long<-filter(slope_se_long,genus_species!="Hoplitis fulgida" & genus_species!="Hoplitis robusta" & genus_species!="Hylaeus annulatus" & genus_species!="Osmia albolateralis")

# create figure
p<-ggplot(slope_se_long, aes(y = slope, x = term_id_renamed, fill=pvalue_sig)) +
  xlab("Climate predictor") + ylab("Parameter estimate") +
  geom_errorbar(aes(ymin=slope-se, ymax=slope+se), width=.2,position=position_dodge(width=0.5)) +
  geom_point(size=2, shape=21, position=position_dodge(width=0.5), color="black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  geom_hline(yintercept=0, linetype='dashed', col = 'darkgray') + facet_wrap(~genus_species,ncol=3) + theme(legend.position = "none") +
  scale_fill_manual(values=c("white","black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(strip.text = element_text(face = "italic"), strip.background = element_rect(fill = "white"))
p

#ggsave("glmm_indivspecies_parameters.png", p,width=5.5,height=6.5,units = c("in"),dpi = 600)


##### For individual species, is there phylogenetic signal in how sex ratios vary with climate? ######

### Format the phylogenetic tree ###

# read in and plot the tree
rmbl_genus_tree<-read.tree("rmbl_genus_tree_2023-04-11.tre")
plot(rmbl_genus_tree)

# get list of focal bee species
specieslist<-read.csv("focalbeedata_envdata_RMBLsexratios_2023-03-29.csv")
specieslist2<- specieslist %>% group_by(genus_species) %>% summarise(n=n()) %>% separate(genus_species, c('genus', 'species'), remove=FALSE)
specieslist2$genus_species<-paste(specieslist2$genus,specieslist2$species,sep="_")

# turn the phylogeny into a species level tree by adding fake species epithets, designated with "_zzz", to each genus
rmbl_genus_tree$tip.label = paste0(rmbl_genus_tree$tip.label,'_zzz')
# add species to the tree as polytomies
rmbl_species_tree<-congeneric.merge(specieslist2$genus_species, tree = rmbl_genus_tree, split = "_")
plot(rmbl_species_tree)

# drop tips not included in dataset
rmbl_species_tree_2<-keep.tip(phy=rmbl_species_tree,tip=specieslist2$genus_species)
plot(rmbl_species_tree_2)

# label nodes
rmbl_species_tree_2$node.label <- 1:12


### Format GLMM results data in preparation for testing for phylogenetic signal ###

# read in results from individual species GLMMs
results<-read.csv("species_glmms_climate_sex_ratios_2023-03-29.csv")

# filter results to only include the best model for each species
bestmodels<-filter(results,delta_AICc==0)

# create long-form data frame of slopes (parameter estimates) from the best model for each species, with one row for each parameter estimate in the model
slopes_long<-pivot_longer(bestmodels[,c(21,1,22,8:11)],cols=4:7,values_to="slope")
names(slopes_long)[4]<-"term"
slopes_long$term<-substr(slopes_long$term,6,6) # create term column listing the numeric order of the parameter in the model

# create similar long-form dataset of se values
se_long<-pivot_longer(bestmodels[,c(21,1,22,12:15)],cols=4:7,values_to="se")
names(se_long)[4]<-"term"
se_long$term<-substr(se_long$term,3,3)

# create similar long-form dataset of p-values
pvalues_long<-pivot_longer(bestmodels[,c(21,1,22,4:7)],cols=4:7,values_to="pvalue")
names(pvalues_long)[4]<-"term"
pvalues_long$term<-substr(pvalues_long$term,7,7)

# combine the data frames above
slope_se_long<-left_join(slopes_long,se_long, by=c("genus_species","model","predictors","term"))
slope_se_long<-left_join(slope_se_long,pvalues_long, by=c("genus_species","model","predictors","term"))

# drop rows containing missing values
slope_se_long<-drop_na(slope_se_long)

# get summary of terms in each model
summary_terms<-slope_se_long %>% group_by(model, predictors, term) %>% summarise(count=n())
summary_terms$count<-NULL

# create a data frame listing the name of each model parameter and its term number (numeric order) in the model
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

# add term_id information to the data frame of parameter estimates and their statistics
slope_se_long<-left_join(slope_se_long, summary_terms, by=c("model","predictors","term"))

# rename levels of term_id for graphing purposes
slope_se_long$term_id_renamed<-as.factor(slope_se_long$term_id)
levels(slope_se_long$term_id_renamed) <- c("Precip.","Snowmelt","Prior year's precip.","Prior year's snowmelt")
slope_se_long$term_id_renamed <- factor(slope_se_long$term_id_renamed, levels=c("Snowmelt","Precip.","Prior year's snowmelt","Prior year's precip."))

# create column coding each p-value as significant or nonsignificant
slope_se_long$pvalue_sig<-ifelse(slope_se_long$pvalue>=0.05,"Nonsignificant","Significant")

# create wide-form data frame with parameter estimate values for each bee species
slope_se_wide<-pivot_wider(slope_se_long[,c(1,8,5)], names_from=term_id, values_from=slope)

# replace NAs with 0s
slope_se_wide[is.na(slope_se_wide)]<-0

# replace space with underscore in species epithet
slope_se_wide$genus_species <- sub(" ", "_", slope_se_wide$genus_species)

# make slope_se_wide into a data frame
slope_se_wide<-data.frame(slope_se_wide)


### Test for phylogenetic signal in the slope associated with each climate variable ###

### prev_yr_accum_summer_precip_z ###
# create data frame of slopes associated with the focal climate variable
rownames(slope_se_wide) <- slope_se_wide$genus_species
slope_se_wide2<-slope_se_wide[2]

# make the phylogeny and slope data into a "phylo4d" object
p4d <- phylo4d(x=rmbl_species_tree_2, tip.data=slope_se_wide2, node.data=NULL)

# visualize slope data
barplot.phylo4d(p4d, tree.type = "phylo")

# test for phylogenetic signal
phyloSignal(p4d = p4d, method = "all")


### prev_yr_doy_bare_ground_z ###
# create data frame of slopes associated with the focal climate variable
rownames(slope_se_wide) <- slope_se_wide$genus_species
slope_se_wide2<-slope_se_wide[3]

# make the phylogeny and slope data into a "phylo4d" object
p4d <- phylo4d(x=rmbl_species_tree_2, tip.data=slope_se_wide2, node.data=NULL)

# visualize slope data
barplot.phylo4d(p4d, tree.type = "phylo")

# test for phylogenetic signal
phyloSignal(p4d = p4d, method = "all") # significant!


### doy_bare_ground_z ###
# create data frame of slopes associated with the focal climate variable
rownames(slope_se_wide) <- slope_se_wide$genus_species
slope_se_wide2<-slope_se_wide[4]

# make the phylogeny and slope data into a "phylo4d" object
p4d <- phylo4d(x=rmbl_species_tree_2, tip.data=slope_se_wide2, node.data=NULL)

# visualize slope data
barplot.phylo4d(p4d, tree.type = "phylo")

# test for phylogenetic signal
phyloSignal(p4d = p4d, method = "all") # marginal


### accum_summer_precip_z ###
# create data frame of slopes associated with the focal climate variable
rownames(slope_se_wide) <- slope_se_wide$genus_species
slope_se_wide2<-slope_se_wide[5]

# make the phylogeny and slope data into a "phylo4d" object
p4d <- phylo4d(x=rmbl_species_tree_2, tip.data=slope_se_wide2, node.data=NULL)

# visualize slope data
barplot.phylo4d(p4d, tree.type = "phylo")

# test for phylogenetic signal
phyloSignal(p4d = p4d, method = "all") 


##### Test for differences in sex ratio between trap and net data #####

# get unique values for collection method
unique(beedatafocal$method)
# recode certain values
beedatafocal$method[beedatafocal$method=="bowl"] <- "Bowl"
beedatafocal$method[beedatafocal$method=="Net, AM"] <- "Net"
beedatafocal$method[beedatafocal$method=="Net, PM"] <- "Net"

# GLMM of sex as a function of collection method
m1 <- glmer(sex_binom ~ method + (1|site) + (1|genus_species), data = beedatafocal, family = binomial)
summary(m1)
library(effects)
plot(allEffects(m1))


##### SEMs: across species, does climate most strongly influence sex ratios directly, or indirectly via floral resources? #####

# for analyses excluding 2016, run the following line of code
#beedatafocal<-filter(beedatafocal,year!=2016)

### Build series of SEMs representing all combinations of predictors ###

### (full) SEM: all prior year's climate and floral predictor variables
semFull<-psem(

  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_floral_days_z +
          prev_yr_total_flowers_z +
          (1|site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z, random=~1|site,
     data = beedatafocal, method="ML"),

  lme(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z, random=~1|site,
      data = beedatafocal, method="ML"), 
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z

)

summary(semFull, .progressBar = F)
plot(semFull)


### (1) SEM: snowmelt
sem1<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          (1|site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     random=~1|site, data = beedatafocal, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     random=~1|site, data = beedatafocal, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem1, .progressBar = F)
plot(sem1)


### (2) SEM: snowmelt, floral sum 
sem2<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_total_flowers_z +
          (1|site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     random=~1|site, data = beedatafocal, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     random=~1|site, data = beedatafocal, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem2, .progressBar = F)
plot(sem2)

### (3) SEM: snowmelt, floral days 
sem3<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_floral_days_z +
          (1|site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     random=~1|site, data = beedatafocal, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     random=~1|site, data = beedatafocal, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem3, .progressBar = F)
plot(sem3)

### (4) SEM: snowmelt, precip 
sem4<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_accum_summer_precip_cm_z +
          (1|site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     random=~1|site, data = beedatafocal, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     random=~1|site, data = beedatafocal, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem4, .progressBar = F)
plot(sem4)

### (5) SEM: precip 
sem5<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ 
          prev_yr_accum_summer_precip_cm_z +
          (1|site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     random=~1|site, data = beedatafocal, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     random=~1|site, data = beedatafocal, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem5, .progressBar = F)
plot(sem5)

### (6) SEM: precip, floral sum 
sem6<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ 
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_total_flowers_z +
          (1|site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     random=~1|site, data = beedatafocal, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     random=~1|site, data = beedatafocal, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem6, .progressBar = F)
plot(sem6)

### (7) SEM: precip, floral days 
sem7<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ 
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_floral_days_z +
          (1|site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     random=~1|site, data = beedatafocal, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     random=~1|site, data = beedatafocal, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem7, .progressBar = F)
plot(sem7)


### (8) SEM: floral sum, floral days 
sem8<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~
          prev_yr_floral_days_z +
          prev_yr_total_flowers_z +
          (1|site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~ 
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     random=~1|site, data = beedatafocal, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     random=~1|site, data = beedatafocal, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem8, .progressBar = F)
plot(sem8)


### (9) SEM: floral sum 
sem9<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~
          prev_yr_total_flowers_z +
          (1|site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~ 
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     random=~1|site, data = beedatafocal, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     random=~1|site, data = beedatafocal, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem9, .progressBar = F)
plot(sem9)

### (10) SEM: floral days 
sem10<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~
          prev_yr_floral_days_z +
          (1|site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~ 
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     random=~1|site, data = beedatafocal, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     random=~1|site, data = beedatafocal, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem10, .progressBar = F)
plot(sem10)


### (11) SEM: snowmelt, precip, floral sum 
sem11<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_total_flowers_z +
          (1|site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~ 
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     random=~1|site, data = beedatafocal, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     random=~1|site, data = beedatafocal, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem11, .progressBar = F)
plot(sem11)


### (12) SEM: snowmelt, precip, floral days 
sem12<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_floral_days_z +
          (1|site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     random=~1|site, data = beedatafocal, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     random=~1|site, data = beedatafocal, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem12, .progressBar = F)
plot(sem12)

### (13) SEM: snowmelt, floral sum, floral days 
sem13<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_floral_days_z +
          prev_yr_total_flowers_z +
          (1|site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     random=~1|site, data = beedatafocal, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     random=~1|site, data = beedatafocal, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem13, .progressBar = F)
plot(sem13)

### (14) SEM: precip, floral sum, floral days 
sem14<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ 
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_floral_days_z +
          prev_yr_total_flowers_z +
          (1|site) + (1|genus_species),
        data = beedatafocal, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     random=~1|site, data = beedatafocal, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
       prev_yr_doy_bare_ground_z +
       prev_yr_accum_summer_precip_cm_z,
     random=~1|site, data = beedatafocal, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem14, .progressBar = F)
plot(sem14)


##### SEMs: for across-species models, get summaries of results #####

# get summaries
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

# get chi-squared and Fisher's C test results
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

# rename columns
names(df2)[c(2,3,5,6)]<-c("df_chisq","p_chisq","df_fisher","p_fisher")
# add column of model ID number
df2$model<-0:14
# bind data frames
df3<-bind_cols(df2[,c(7,1:6)],df1)

#write.csv(df3,"GoodnessOfFit_SEM_AllSpecies.csv", row.names=FALSE)

# get model coefficients
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

#write.csv(coefs_df,"coefficients_SEM_AllSpecies.csv", row.names=FALSE)

# print summaries into console
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


##### SEMs: for individual species, does climate most strongly influence sex ratios directly, or indirectly via floral resources? #####
### Halictus rubicundus #####

# read in data
halrub<-read.csv("halictus_rubicundus_SEMdata_RMBLsexratios_2023-07-23.csv")

# for analyses excluding 2016, run the following line of code
#halrub<-filter(halrub,year!=2016)

### Build series of SEMs representing all combinations of predictors ###

### (full) SEM: all prior year's climate and floral predictor variables
semFull<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_floral_days_z +
          prev_yr_total_flowers_z +
          (1|site) ,
        data = halrub, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~ 
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halrub, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halrub, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
)

summary(semFull, .progressBar = F)
plot(semFull)

### (1) SEM: snowmelt
sem1<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          (1|site) ,
        data = halrub, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halrub, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halrub, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem1, .progressBar = F)
plot(sem1)


### (2) SEM: snowmelt, floral sum 
sem2<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_total_flowers_z +
          (1|site) ,
        data = halrub, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halrub, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halrub, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem2, .progressBar = F)
plot(sem2)

### (3) SEM: snowmelt, floral days 
sem3<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_floral_days_z +
          (1|site) ,
        data = halrub, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halrub, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halrub, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem3, .progressBar = F)
plot(sem3)

### (4) SEM: snowmelt, precip 
sem4<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_accum_summer_precip_cm_z +
          (1|site) ,
        data = halrub, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halrub, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halrub, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem4, .progressBar = F)
plot(sem4)

### (5) SEM: precip 
sem5<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ 
          prev_yr_accum_summer_precip_cm_z +
          (1|site) ,
        data = halrub, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halrub, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halrub, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem5, .progressBar = F)
plot(sem5)

### (6) SEM: precip, floral sum 
sem6<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ 
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_total_flowers_z +
          (1|site) ,
        data = halrub, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halrub, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halrub, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem6, .progressBar = F)
plot(sem6)

### (7) SEM: precip, floral days 
sem7<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ 
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_floral_days_z +
          (1|site) ,
        data = halrub, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halrub, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halrub, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem7, .progressBar = F)
plot(sem7)


### (8) SEM: floral sum, floral days 
sem8<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~
          prev_yr_floral_days_z +
          prev_yr_total_flowers_z +
          (1|site) ,
        data = halrub, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halrub, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halrub, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem8, .progressBar = F)
plot(sem8)


### (9) SEM: floral sum 
sem9<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~
          prev_yr_total_flowers_z +
          (1|site) ,
        data = halrub, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halrub, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halrub, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem9, .progressBar = F)
plot(sem9)

### (10) SEM: floral days 
sem10<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~
          prev_yr_floral_days_z +
          (1|site) ,
        data = halrub, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halrub, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halrub, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem10, .progressBar = F)
plot(sem10)

### (11) SEM: snowmelt, precip, floral sum 
sem11<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_total_flowers_z +
          (1|site) ,
        data = halrub, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halrub, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halrub, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem11, .progressBar = F)
plot(sem11)


### (12) SEM: snowmelt, precip, floral days 
sem12<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_floral_days_z +
          (1|site) ,
        data = halrub, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halrub, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halrub, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem12, .progressBar = F)
plot(sem12)

### (13) SEM: snowmelt, floral sum, floral days 
sem13<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_floral_days_z +
          prev_yr_total_flowers_z +
          (1|site) ,
        data = halrub, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halrub, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halrub, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem13, .progressBar = F)
plot(sem13)

### (14) SEM: precip, floral sum, floral days 
sem14<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ 
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_floral_days_z +
          prev_yr_total_flowers_z +
          (1|site) ,
        data = halrub, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halrub, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halrub, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem14, .progressBar = F)
plot(sem14)

# get summaries
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

# get chi-squared and Fisher's C test results
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

# rename columns
names(df2)[c(2,3,5,6)]<-c("df_chisq","p_chisq","df_fisher","p_fisher")
# add column of model ID number
df2$model<-0:14
# bind data frames
df3<-bind_cols(df2[,c(7,1:6)],df1)
# calculate delta AICc
df3<-df3 %>% mutate(deltaAICc=AICc-min(AICc))

#write.csv(df3,"GoodnessOfFit_SEM_Halictus_rubicundus.csv", row.names=FALSE)

# get model coefficients
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

#write.csv(coefs_df,"coefficients_SEM_Halictus_rubicundus_no2016.csv", row.names=FALSE)

# print summaries into console
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


### Halictus virgatellus #####

# read in data
halvir<-read.csv("halictus_virgatellus_SEMdata_RMBLsexratios_2023-07-23.csv")

# for analyses excluding 2016, run the following line of code
#halvir<-filter(halvir,year!=2016)

### Build series of SEMs representing all combinations of predictors ###

### (full) SEM: all prior year's climate and floral predictor variables
semFull<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_floral_days_z +
          prev_yr_total_flowers_z +
          (1|site) ,
        data = halvir, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halvir, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halvir, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(semFull, .progressBar = F)
plot(semFull)


### (1) SEM: snowmelt
sem1<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          (1|site) ,
        data = halvir, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halvir, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halvir, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem1, .progressBar = F)
plot(sem1)


### (2) SEM: snowmelt, floral sum 
sem2<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_total_flowers_z +
          (1|site) ,
        data = halvir, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halvir, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halvir, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem2, .progressBar = F)
plot(sem2)

### (3) SEM: snowmelt, floral days 
sem3<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_floral_days_z +
          (1|site) ,
        data = halvir, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halvir, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halvir, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem3, .progressBar = F)
plot(sem3)

### (4) SEM: snowmelt, precip 
sem4<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_accum_summer_precip_cm_z +
          (1|site) ,
        data = halvir, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halvir, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halvir, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem4, .progressBar = F)
plot(sem4)

### (5) SEM: precip 
sem5<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ 
          prev_yr_accum_summer_precip_cm_z +
          (1|site) ,
        data = halvir, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halvir, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halvir, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem5, .progressBar = F)
plot(sem5)

### (6) SEM: precip, floral sum 
sem6<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ 
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_total_flowers_z +
          (1|site) ,
        data = halvir, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halvir, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halvir, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem6, .progressBar = F)
plot(sem6)

### (7) SEM: precip, floral days 
sem7<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ 
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_floral_days_z +
          (1|site) ,
        data = halvir, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halvir, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halvir, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem7, .progressBar = F)
plot(sem7)


### (8) SEM: floral sum, floral days 
sem8<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~
          prev_yr_floral_days_z +
          prev_yr_total_flowers_z +
          (1|site) ,
        data = halvir, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halvir, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halvir, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem8, .progressBar = F)
plot(sem8)


### (9) SEM: floral sum 
sem9<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~
          prev_yr_total_flowers_z +
          (1|site) ,
        data = halvir, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halvir, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halvir, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem9, .progressBar = F)
plot(sem9)

### (10) SEM: floral days 
sem10<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~
          prev_yr_floral_days_z +
          (1|site) ,
        data = halvir, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halvir, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halvir, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem10, .progressBar = F)
plot(sem10)


### (11) SEM: snowmelt, precip, floral sum 
sem11<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_total_flowers_z +
          (1|site) ,
        data = halvir, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halvir, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halvir, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem11, .progressBar = F)
plot(sem11)


### (12) SEM: snowmelt, precip, floral days 
sem12<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_floral_days_z +
          (1|site) ,
        data = halvir, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halvir, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halvir, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem12, .progressBar = F)
plot(sem12)

### (13) SEM: snowmelt, floral sum, floral days 
sem13<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ prev_yr_doy_bare_ground_z +
          prev_yr_floral_days_z +
          prev_yr_total_flowers_z +
          (1|site) ,
        data = halvir, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halvir, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halvir, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem13, .progressBar = F)
plot(sem13)

### (14) SEM: precip, floral sum, floral days 
sem14<-psem(
  
  # effects of climate on sex ratios
  glmer(sex_binom ~ 
          prev_yr_accum_summer_precip_cm_z +
          prev_yr_floral_days_z +
          prev_yr_total_flowers_z +
          (1|site) ,
        data = halvir, family = binomial),
  
  # effects of climate on flowers
  lme(prev_yr_floral_days_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halvir, method="ML"),
  
  lme(prev_yr_total_flowers_z ~
        prev_yr_doy_bare_ground_z +
        prev_yr_accum_summer_precip_cm_z,
      random=~1|site, data = halvir, method="ML"),
  
  prev_yr_total_flowers_z %~~% prev_yr_floral_days_z
  
)

summary(sem14, .progressBar = F)
plot(sem14)


# get summaries
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

# get chi-squared and Fisher's C test results
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

# rename columns
names(df2)[c(2,3,5,6)]<-c("df_chisq","p_chisq","df_fisher","p_fisher")
# add column of model ID number
df2$model<-0:14
# bind data frames
df3<-bind_cols(df2[,c(7,1:6)],df1)
# calculate delta AICc
df3<-df3 %>% mutate(deltaAICc=AICc-min(AICc))

#write.csv(df3,"GoodnessOfFit_SEM_Halictus_virgatellus.csv", row.names=FALSE)

# get model coefficients
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

#write.csv(coefs_df,"coefficients_SEM_Halictus_virgatellus.csv", row.names=FALSE)

# print summaries into console
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


##### Figures: SEMs #####

### Hypothesis figure ###
grViz("
digraph boxes_and_circles {

rotate=90

  # a 'graph' statement
  graph [overlap = true, fontsize = 10]

  # several 'node' statements
  node [shape = square,
        fontname = Helvetica,
        label = ''
        ]
  1; 2; 3; 5; 4
  

  # several 'edge' statements
  edge [color = black, penwidth=1, style=solid, arrowsize=0.75]
  5->2 
  
  edge [color = black, penwidth=1,style=solid, arrowsize=0.75]
  4->2 
  
  edge [color = black, penwidth=1, style=solid, arrowsize=0.75]
  4->3 
  
  edge [color = black, penwidth=1, style=solid, arrowsize=0.75]
  5->3 
  
  edge [color = black, penwidth=1, style=solid, arrowsize=0.75]
  4->1 
  
  edge [color = black, penwidth=1, style=solid, arrowsize=0.75]
  3->1
  
  edge [color = black, penwidth=1, style=solid, arrowsize=0.75]
  2->1
  
  edge [color = black, penwidth=1, style=solid, arrowsize=0.75]  
  5->1
  
}
")


### Best model: across species, all years of data ###
grViz("
digraph boxes_and_circles {

rotate=90

  # a 'graph' statement
  graph [overlap = true, fontsize = 10]

  # several 'node' statements
  node [shape = square,
        fontname = Helvetica,
        label = ''
        ]
  1; 2; 3; 5; 4
  
    node [shape = point,
        label = ''
        ]
  y; z; a; b; c; d; e; f;

  # several 'edge' statements
  edge [color = black, penwidth=0.8856, style=solid, arrowsize=0.75]
  5->2 
  
  edge [color = black, penwidth=1.0040,style=dashed, arrowsize=0.75]
  4->2 
  
  edge [color = black, penwidth=0.1364, style=solid, arrowsize=0.75]
  4->3 
  
  edge [color = black, penwidth=1.6568, style=solid, arrowsize=0.75]
  5->3 
  
  edge [color = black, penwidth=1.4984, style=solid, arrowsize=0.75]
  4->1 
  
  edge [arrowhead = halfopen, color = black, penwidth=0.0846, style=dashed, arrowsize=0.75]
  3->1
  
  edge [arrowhead = halfopen, color = black, penwidth=0.0478, style=dashed, arrowsize=0.75]
  2->1
  
  edge [color = white, arrowhead = none, penwidth=0.0000000000000001]
  5->1
  
  edge [color = black, penwidth=2, style=solid, arrowhead=none]
  y->z
  
  edge [color = black, penwidth=1.5, style=solid, arrowhead=none]
  a->b
  
  edge [color = black, penwidth=1, style=solid, arrowhead=none]
  c->d
  
  edge [color = black, penwidth=0.5, style=solid, arrowhead=none]
  e->f
}
")


### Best model: across species, 2016 excluded ###
grViz("
digraph boxes_and_circles {

rotate=90

  # a 'graph' statement
  graph [overlap = true, fontsize = 10]

  # several 'node' statements
  node [shape = square,
        fontname = Helvetica,
        label = ''
        ]
  1; 2; 3; 5; 4
  
    node [shape = point,
        label = ''
        ]
  y; z; a; b; c; d; e; f;

  # several 'edge' statements
  edge [color = black, penwidth=0.5728, style=solid, arrowsize=0.75]
  5->2 
  
  edge [color = black, penwidth=0.8328,style=dashed, arrowsize=0.75]
  4->2 
  
  edge [color = black, penwidth=0.9856, style=solid, arrowsize=0.75]
  4->3 
  
  edge [color = black, penwidth=0.2592, style=solid, arrowsize=0.75]
  5->3 
  
  edge [color = black, penwidth=1.6116, style=solid, arrowsize=0.75]
  4->1 
  
  edge [arrowhead = halfopen, color = black, penwidth=0.234, style=dashed, arrowsize=0.75]
  3->1
  
  edge [arrowhead = halfopen, color = black, penwidth=0.084, style=dashed, arrowsize=0.75]
  2->1
  
  edge [color = white, arrowhead = none, penwidth=0.0000000000000001]
  5->1
  
  edge [color = black, penwidth=2, style=solid, arrowhead=none]
  y->z
  
  edge [color = black, penwidth=1.5, style=solid, arrowhead=none]
  a->b
  
  edge [color = black, penwidth=1, style=solid, arrowhead=none]
  c->d
  
  edge [color = black, penwidth=0.5, style=solid, arrowhead=none]
  e->f
}
")


### Best model: Halictus rubicundus, all years of data ###

grViz("
digraph boxes_and_circles {

rotate=90

  # a 'graph' statement
  graph [overlap = true, fontsize = 10]

  # several 'node' statements
  node [shape = square,
        fontname = Helvetica,
        label = ''
        ]
  1; 2; 3; 5; 4
  
    node [shape = point,
        label = ''
        ]
  y; z; a; b; c; d; e; f;

  # several 'edge' statements
  edge [color = black, penwidth=0.7912, style=solid, arrowsize=0.75]
  5->2 
  
  edge [color = black, penwidth=2.0868,style=dashed, arrowsize=0.75]
  4->2 
  
  edge [color = black, penwidth=0.6006, style=solid, arrowsize=0.75]
  4->3 
  
  edge [color = black, penwidth=1.7182, style=solid, arrowsize=0.75]
  5->3 
  
  edge [color = black, penwidth=2.4698, style=solid, arrowsize=0.75]
  4->1 
  
  edge [color = black, penwidth=1.0450, style=solid, arrowsize=0.75]
  3->1
  
  edge [color = black, penwidth=0.5386, style=solid, arrowsize=0.75]
  2->1
  
  edge [color = white, arrowhead = none, penwidth=0.0000000000000001]
  5->1
  
  edge [color = black, penwidth=2, style=solid, arrowhead=none]
  y->z
  
  edge [color = black, penwidth=1.5, style=solid, arrowhead=none]
  a->b
  
  edge [color = black, penwidth=1, style=solid, arrowhead=none]
  c->d
  
  edge [color = black, penwidth=0.5, style=solid, arrowhead=none]
  e->f
}
")


### Best model: Halicitus rubicundus, 2016 excluded ###

grViz("
digraph boxes_and_circles {

rotate=90

  # a 'graph' statement
  graph [overlap = true, fontsize = 10]

  # several 'node' statements
  node [shape = square,
        fontname = Helvetica,
        label = ''
        ]
  1; 2; 3; 5; 4
  
    node [shape = point,
        label = ''
        ]
  y; z; a; b; c; d; e; f;

  # several 'edge' statements
  edge [color = black, penwidth=0.8702, style=solid, arrowsize=0.75]
  5->2 
  
  edge [color = black, penwidth=2.2958,style=dashed, arrowsize=0.75]
  4->2 
  
  edge [color = black, penwidth=4.0090, style=solid, arrowsize=0.75]
  4->3 
  
  edge [color = black, penwidth=0.7758, style=solid, arrowsize=0.75]
  5->3 
  
  edge [color = black, penwidth=1.9396, style=solid, arrowsize=0.75]
  3->1
  
  edge [color = black, penwidth=0.7982, style=solid, arrowsize=0.75]
  2->1
  
  edge [color = black, penwidth=0.7962, style=dashed, arrowsize=0.75]
  5->1
  
  edge [color = white, arrowhead = none, penwidth=0.0000000000000001]
  4->1 
  
  edge [color = black, penwidth=2, style=solid, arrowhead=none]
  y->z
  
  edge [color = black, penwidth=1.5, style=solid, arrowhead=none]
  a->b
  
  edge [color = black, penwidth=1, style=solid, arrowhead=none]
  c->d
  
  edge [color = black, penwidth=0.5, style=solid, arrowhead=none]
  e->f
}
")


### Best model: Halictus virgatellus, all years of data ###
grViz("
digraph boxes_and_circles {

rotate=90

  # a 'graph' statement
  graph [overlap = true, fontsize = 10]

  # several 'node' statements
  node [shape = square,
        fontname = Helvetica,
        label = ''
        ]
  1; 2; 3; 5; 4
  
    node [shape = point,
        label = ''
        ]
  y; z; a; b; c; d; e; f;

  # several 'edge' statements
  edge [color = black, penwidth=1.2538, style=solid, arrowsize=0.75]
  5->2 
  
  edge [color = black, penwidth=0.7554,style=dashed, arrowsize=0.75]
  4->2 
  
  edge [color = black, penwidth=0.5544, style=solid, arrowsize=0.75]
  4->3 
  
  edge [color = black, penwidth=1.8558, style=solid, arrowsize=0.75]
  5->3 
  
  edge [color = black, penwidth=1.0946, style=solid, arrowsize=0.75]
  4->1 
  
  edge [color = black, penwidth=0.6290, style=solid, arrowsize=0.75]
  3->1
  
  edge [arrowhead = halfopen, color = black, penwidth=0.0842, style=dashed, arrowsize=0.75]
  2->1
  
  edge [color = white, arrowhead = none, penwidth=0.0000000000000001]
  5->1
  
  edge [color = black, penwidth=2, style=solid, arrowhead=none]
  y->z
  
  edge [color = black, penwidth=1.5, style=solid, arrowhead=none]
  a->b
  
  edge [color = black, penwidth=1, style=solid, arrowhead=none]
  c->d
  
  edge [color = black, penwidth=0.5, style=solid, arrowhead=none]
  e->f
}
")


### Best model: Halictus virgatellus, 2016 excluded ###
grViz("
digraph boxes_and_circles {

rotate=90

  # a 'graph' statement
  graph [overlap = true, fontsize = 10]

  # several 'node' statements
  node [shape = square,
        fontname = Helvetica,
        label = ''
        ]
  1; 2; 3; 5; 4
  
    node [shape = point,
        label = ''
        ]
  y; z; a; b; c; d; e; f;

  # several 'edge' statements
  edge [color = black, penwidth=1.1386, style=solid, arrowsize=0.75]
  5->2 
  
  edge [color = black, penwidth=0.7400,style=dashed, arrowsize=0.75]
  4->2 
  
  edge [color = black, penwidth=3.2654, style=solid, arrowsize=0.75]
  4->3 
  
  edge [color = black, penwidth=1.2582, style=solid, arrowsize=0.75]
  5->3 
  
  edge [color = black, penwidth=0.8902, style=solid, arrowsize=0.75]
  3->1
  
  edge [arrowhead = halfopen, color = black, penwidth=0.1786, style=solid, arrowsize=0.75]
  2->1
  
  edge [color = white, arrowhead = none, penwidth=0.0000000000000001]
  5->1
  
  edge [color = white, arrowhead = none, penwidth=0.0000000000000001]
  4->1 
  
  edge [color = black, penwidth=2, style=solid, arrowhead=none]
  y->z
  
  edge [color = black, penwidth=1.5, style=solid, arrowhead=none]
  a->b
  
  edge [color = black, penwidth=1, style=solid, arrowhead=none]
  c->d
  
  edge [color = black, penwidth=0.5, style=solid, arrowhead=none]
  e->f
}
")