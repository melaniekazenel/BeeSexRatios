# Bee sex ratios
# SEMs using floral lists for individual bee species
# 30 March 2023
# ---------------------------------------------------

# load relevant libraries #####
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


# set wd
setwd("~/Documents/*RMBLPhenology/Projects/BeeSexRatios/Analyses")




############## (9) Halictus rubicundus ################

beedatafocal<-read.csv("focalbeedata_envdata_RMBLsexratios_2023-03-29.csv")
# subset to just include focal sites
sitelist<-c("Almont","Almont Curve","Beaver","CDOT","Copper","Davids","Elko","Gothic","Hill","Little","Lypps","Mexican Cut","Rustlers","Seans","Tuttle","Willey")
beedatafocal<-filter(beedatafocal, site %in% sitelist) 

beedatafocal<-filter(beedatafocal,genus_species=="Halictus rubicundus")

# just look at Gothic sites
# unique(beedatafocal$block)
# test<-beedatafocal %>% group_by(site,block) %>% summarise(n=n())
# beedatafocal<-filter(beedatafocal,block=="A" | block=="B" | block=="D")
# unique(beedatafocal$site)

# # just look at mid-elevation sites
# siteinfo<-read.csv("site_info.csv")
# beedatafocal<-left_join(beedatafocal,siteinfo,by="site")
# beedatafocal<-filter(beedatafocal, elev_group=="Mid")

beedatafocal<-beedatafocal[,-c(56:71)]

##### Adding floral data: Inouye data #####

# read in floral data
flor<-read.csv("floral_data_annual_summaries_forsexratios_halictus_rubicundus.csv")
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

summary(semFull, .progressBar = F)
plot(semFull)

plot(beedatafocal$total_flowers_z,beedatafocal$floral_days_z)

coef<-summary(semFull, .progressBar = F)$coef
coef
unique(coef$Response)
unique(coef$Predictor)

coef<-subset(coef,Predictor!="~~prev_yr_floral_days_z")

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
    fontcolor="black",     fillcolor="white",     fontname="Baskerville",
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
df3<-df3 %>% mutate(deltaAICc=AICc-min(AICc))

#write.csv(df3,"GoodnessOfFit_SEM_Halictus_rubicundus.csv", row.names=FALSE)

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

#write.csv(coefs_df,"coefficients_SEM_Halictus_rubicundus.csv", row.names=FALSE)


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



### Plot the best SEM (lowest AICc) #####
coef<-summary(sem13, .progressBar = F)$coef
coef
unique(coef$Response)
unique(coef$Predictor)

coef<-subset(coef,Predictor!="~~prev_yr_floral_days_z")

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
    fontcolor="black",
    fillcolor="white",
    fontname="Baskerville",
    penwidth=0.5,
    shape="rectangle",
    width=0.7,
    height=0.3,
    label = c("Female/male \nratio","Floral days \n(prior year)", "Floral sum \n(prior year)","Snowmelt date \n(prior year)","Summer precip. \n(prior year)"))

graph <-
  create_graph(
    nodes_df = ndf)

render_graph(graph)

# edf <-
#   create_edge_df(
#     from = coef$from,
#     to = coef$to,
#     rel = "leading_to",
#     #label = round(coef$Estimate,digits=2),
#     penwidth = abs(coef$Estimate*2),
#     fontsize=5.5,
#     fontcolor="brown3",
#     #color=coef$color,
#     color="black",
#     style=coef$style)

edf <-
  create_edge_df(
    from = coef$from,
    to = coef$to,
    rel = "leading_to",
    #label = round(coef$Estimate,digits=2),
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

#graph %>% export_graph(file_name = "SEM_Haliictus_rubicundus_bestmodel.pdf")

### Perform model averaging #####

# sem full, sem 13
coef_list<-list(coef(summary(semFull[[1]]))[, 1], coef(summary(sem13[[1]]))[, 1])
weightslist<-c(summary(semFull, .progressBar = F)$AIC$AICc,summary(sem13, .progressBar = F)$AIC$AICc)
avgEst(coef_list)
avgEst(coef_list,weights=weightslist)

coef_list2<-list(coef(summary(semFull[[2]]))[, 1], coef(summary(sem13[[2]]))[, 1])
weightslist<-c(summary(semFull, .progressBar = F)$AIC$AICc,summary(sem13, .progressBar = F)$AIC$AICc)
avgEst(coef_list2)
avgEst(coef_list2,weights=weightslist)

coef_list3<-list(coef(summary(semFull[[3]]))[, 1], coef(summary(sem13[[3]]))[, 1])
weightslist<-c(summary(semFull, .progressBar = F)$AIC$AICc,summary(sem13, .progressBar = F)$AIC$AICc)
avgEst(coef_list3)
avgEst(coef_list3,weights=weightslist)


## plot the model-averaged SEM #####
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

coef<-subset(coef,Predictor!="~~prev_yr_floral_days_z")

coef<-left_join(coef,est_all,by=c("Predictor","Response"))

# remove variables not included in best model
coef<-coef %>% drop_na(Estimate_MdlAvg)

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
    fontcolor="black",     fillcolor="white",     fontname="Baskerville",
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

#graph %>% export_graph(file_name = "SEM_Halictus_rubicundus_modelavg.pdf")










############## (10) Halictus virgatellus ################

beedatafocal<-read.csv("focalbeedata_envdata_RMBLsexratios_2023-03-29.csv")
# subset to just include focal sites
sitelist<-c("Almont","Almont Curve","Beaver","CDOT","Copper","Davids","Elko","Gothic","Hill","Little","Lypps","Mexican Cut","Rustlers","Seans","Tuttle","Willey")
beedatafocal<-filter(beedatafocal, site %in% sitelist) 

beedatafocal<-filter(beedatafocal,genus_species=="Halictus virgatellus")

# Gothic sites
# beedatafocal<-filter(beedatafocal,block=="A" | block=="B" | block=="D")

beedatafocal<-beedatafocal[,-c(56:71)]

##### Adding floral data: Inouye data #####

# read in floral data
flor<-read.csv("floral_data_annual_summaries_forsexratios_halictus_virgatellus.csv")
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

summary(semFull, .progressBar = F)
plot(semFull)

plot(beedatafocal$total_flowers_z,beedatafocal$floral_days_z)

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
    fontcolor="black",     fillcolor="white",     fontname="Baskerville",
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
df3<-df3 %>% mutate(deltaAICc=AICc-min(AICc))

#write.csv(df3,"GoodnessOfFit_SEM_Halictus_virgatellus.csv", row.names=FALSE)

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


### Plot the best SEM (lowest AICc) #####
coef<-summary(sem2, .progressBar = F)$coef
coef
unique(coef$Response)
unique(coef$Predictor)

coef<-subset(coef,Predictor!="~~prev_yr_floral_days_z")


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
    fontcolor="black",
    fillcolor="white",
    fontname="Baskerville",
    penwidth=0.5,
    shape="rectangle",
    width=0.7,
    height=0.3,
    label = c("Female/male \nratio","Floral days \n(prior year)", "Floral sum \n(prior year)","Snowmelt date \n(prior year)","Summer precip. \n(prior year)"))

graph <-
  create_graph(
    nodes_df = ndf)

render_graph(graph)

# edf <-
#   create_edge_df(
#     from = coef$from,
#     to = coef$to,
#     rel = "leading_to",
#     #label = round(coef$Estimate,digits=2),
#     penwidth = abs(coef$Estimate*2),
#     fontsize=5.5,
#     fontcolor="brown3",
#     #color=coef$color,
#     color="black",
#     style=coef$style)

edf <-
  create_edge_df(
    from = coef$from,
    to = coef$to,
    rel = "leading_to",
    #label = round(coef$Estimate,digits=2),
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

#graph %>% export_graph(file_name = "SEM_Halictus_virgatellus_bestmodel.pdf")

### Perform model averaging #####

# sem full, sem 2, sem 11, sem 13
coef_list<-list(coef(summary(semFull[[1]]))[, 1], coef(summary(sem2[[1]]))[, 1],coef(summary(sem11[[1]]))[, 1],coef(summary(sem13[[1]]))[, 1])
weightslist<-c(summary(semFull, .progressBar = F)$AIC$AICc,summary(sem2, .progressBar = F)$AIC$AICc,summary(sem11, .progressBar = F)$AIC$AICc,summary(sem13, .progressBar = F)$AIC$AICc)
avgEst(coef_list)
avgEst(coef_list,weights=weightslist)

coef_list2<-list(coef(summary(semFull[[2]]))[, 1], coef(summary(sem2[[2]]))[, 1],coef(summary(sem11[[2]]))[, 1], coef(summary(sem13[[2]]))[, 1])
weightslist<-c(summary(semFull, .progressBar = F)$AIC$AICc,summary(sem2, .progressBar = F)$AIC$AICc,summary(sem11, .progressBar = F)$AIC$AICc,summary(sem13, .progressBar = F)$AIC$AICc)
avgEst(coef_list2)
avgEst(coef_list2,weights=weightslist)

coef_list3<-list(coef(summary(semFull[[3]]))[, 1], coef(summary(sem2[[3]]))[, 1],coef(summary(sem11[[3]]))[, 1],coef(summary(sem13[[3]]))[, 1])
weightslist<-c(summary(semFull, .progressBar = F)$AIC$AICc,summary(sem2, .progressBar = F)$AIC$AICc,summary(sem11, .progressBar = F)$AIC$AICc,summary(sem13, .progressBar = F)$AIC$AICc)
avgEst(coef_list3)
avgEst(coef_list3,weights=weightslist)


## plot the SEM #####
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

coef<-subset(coef,Predictor!="~~prev_yr_floral_days_z")

coef<-left_join(coef,est_all,by=c("Predictor","Response"))

# remove variables not included in best model
coef<-coef %>% drop_na(Estimate_MdlAvg)

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
    fontcolor="black",     fillcolor="white",     fontname="Baskerville",
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

#graph %>% export_graph(file_name = "SEM_Halictus_virgatellus_modelavg.pdf")







############## (11) Hoplitis fulgida ################

beedatafocal<-read.csv("focalbeedata_envdata_RMBLsexratios_2023-03-29.csv")
# subset to just include focal sites
sitelist<-c("Almont","Almont Curve","Beaver","CDOT","Copper","Davids","Elko","Gothic","Hill","Little","Lypps","Mexican Cut","Rustlers","Seans","Tuttle","Willey")
beedatafocal<-filter(beedatafocal, site %in% sitelist) 

beedatafocal<-filter(beedatafocal,genus_species=="Hoplitis fulgida")

beedatafocal<-beedatafocal[,-c(56:71)]

##### Adding floral data: Inouye data #####

# read in floral data
flor<-read.csv("floral_data_annual_summaries_forsexratios_hoplitis_fulgida.csv")
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

summary(semFull, .progressBar = F)
plot(semFull)

plot(beedatafocal$total_flowers_z,beedatafocal$floral_days_z)

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
    fontcolor="black",     fillcolor="white",     fontname="Baskerville",
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
df3<-df3 %>% mutate(deltaAICc=AICc-min(AICc))

#write.csv(df3,"GoodnessOfFit_SEM_hoplitis_fulgida.csv", row.names=FALSE)

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

#write.csv(coefs_df,"coefficients_SEM_hoplitis_fulgida.csv", row.names=FALSE)


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


### Plot the best SEM (lowest AICc) #####
coef<-summary(sem13, .progressBar = F)$coef
coef
unique(coef$Response)
unique(coef$Predictor)

coef<-subset(coef,Predictor!="~~prev_yr_floral_days_z")

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
    fontcolor="black",
    fillcolor="white",
    fontname="Baskerville",
    penwidth=0.5,
    shape="rectangle",
    width=0.7,
    height=0.3,
    label = c("Female/male \nratio","Floral days \n(prior year)", "Floral sum \n(prior year)","Snowmelt date \n(prior year)","Summer precip. \n(prior year)"))

graph <-
  create_graph(
    nodes_df = ndf)

render_graph(graph)

# edf <-
#   create_edge_df(
#     from = coef$from,
#     to = coef$to,
#     rel = "leading_to",
#     #label = round(coef$Estimate,digits=2),
#     penwidth = abs(coef$Estimate*2),
#     fontsize=5.5,
#     fontcolor="brown3",
#     #color=coef$color,
#     color="black",
#     style=coef$style)

edf <-
  create_edge_df(
    from = coef$from,
    to = coef$to,
    rel = "leading_to",
    #label = round(coef$Estimate,digits=2),
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

#graph %>% export_graph(file_name = "SEM_Hoplitis_fulgida_bestmodel.pdf")

### Perform model averaging #####

# sem full, sem 2, sem 11
coef_list<-list(coef(summary(semFull[[1]]))[, 1], coef(summary(sem2[[1]]))[, 1],coef(summary(sem11[[1]]))[, 1])
weightslist<-c(summary(semFull, .progressBar = F)$AIC$AICc,summary(sem2, .progressBar = F)$AIC$AICc,summary(sem11, .progressBar = F)$AIC$AICc)
avgEst(coef_list)
avgEst(coef_list,weights=weightslist)

coef_list2<-list(coef(summary(semFull[[2]]))[, 1], coef(summary(sem2[[2]]))[, 1],coef(summary(sem11[[2]]))[, 1])
weightslist<-c(summary(semFull, .progressBar = F)$AIC$AICc,summary(sem2, .progressBar = F)$AIC$AICc,summary(sem11, .progressBar = F)$AIC$AICc)
avgEst(coef_list2)
avgEst(coef_list2,weights=weightslist)

coef_list3<-list(coef(summary(semFull[[3]]))[, 1], coef(summary(sem2[[3]]))[, 1],coef(summary(sem11[[3]]))[, 1])
weightslist<-c(summary(semFull, .progressBar = F)$AIC$AICc,summary(sem2, .progressBar = F)$AIC$AICc,summary(sem11, .progressBar = F)$AIC$AICc)
avgEst(coef_list3)
avgEst(coef_list3,weights=weightslist)

coefs(sem8)
coefs(sem13)
coefs(sem14)

## plot the SEM #####
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

coef<-subset(coef,Predictor!="~~prev_yr_floral_days_z")

coef<-left_join(coef,est_all,by=c("Predictor","Response"))

# remove variables not included in best model
coef<-coef %>% drop_na(Estimate_MdlAvg)

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
    fontcolor="black",     fillcolor="white",     fontname="Baskerville",
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

#graph %>% export_graph(file_name = "SEM_Hoplitis_fulgida_modelavg.pdf")







############## (2) Panurginus ineptus ################

beedatafocal<-read.csv("focalbeedata_envdata_RMBLsexratios_2023-03-29.csv")
# subset to just include focal sites
sitelist<-c("Almont","Almont Curve","Beaver","CDOT","Copper","Davids","Elko","Gothic","Hill","Little","Lypps","Mexican Cut","Rustlers","Seans","Tuttle","Willey")
beedatafocal<-filter(beedatafocal, site %in% sitelist) 

beedatafocal<-filter(beedatafocal,genus_species=="Panurginus ineptus")

# Gothic sites
# beedatafocal<-filter(beedatafocal,block=="A" | block=="B" | block=="D")

beedatafocal<-beedatafocal[,-c(56:71)]

##### Adding floral data: Inouye data #####

# read in floral data
flor<-read.csv("floral_data_annual_summaries_forsexratios_panurginus_ineptus.csv")
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

summary(semFull, .progressBar = F)
plot(semFull)

plot(beedatafocal$total_flowers_z,beedatafocal$floral_days_z)

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
    fontcolor="black",     fillcolor="white",     fontname="Baskerville",
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
df3<-df3 %>% mutate(deltaAICc=AICc-min(AICc))

#write.csv(df3,"GoodnessOfFit_SEM_panurginus_ineptus.csv", row.names=FALSE)

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

#write.csv(coefs_df,"coefficients_SEM_panurginus_ineptus.csv", row.names=FALSE)


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


### Plot the best SEM (lowest AICc) #####
coef<-summary(sem11, .progressBar = F)$coef
coef
unique(coef$Response)
unique(coef$Predictor)

coef<-subset(coef,Predictor!="~~prev_yr_floral_days_z")

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
    fontcolor="black",
    fillcolor="white",
    fontname="Baskerville",
    penwidth=0.5,
    shape="rectangle",
    width=0.7,
    height=0.3,
    label = c("Female/male \nratio","Floral days \n(prior year)", "Floral sum \n(prior year)","Snowmelt date \n(prior year)","Summer precip. \n(prior year)"))

graph <-
  create_graph(
    nodes_df = ndf)

render_graph(graph)

# edf <-
#   create_edge_df(
#     from = coef$from,
#     to = coef$to,
#     rel = "leading_to",
#     #label = round(coef$Estimate,digits=2),
#     penwidth = abs(coef$Estimate*2),
#     fontsize=5.5,
#     fontcolor="brown3",
#     #color=coef$color,
#     color="black",
#     style=coef$style)

edf <-
  create_edge_df(
    from = coef$from,
    to = coef$to,
    rel = "leading_to",
    #label = round(coef$Estimate,digits=2),
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

#graph %>% export_graph(file_name = "SEM_Panurginus_ineptus_bestmodel.pdf")

### Perform model averaging #####

# sem full, sem 11
coef_list<-list(coef(summary(semFull[[1]]))[, 1], coef(summary(sem11[[1]]))[, 1])
weightslist<-c(summary(semFull, .progressBar = F)$AIC$AICc,summary(sem11, .progressBar = F)$AIC$AICc)
avgEst(coef_list)
avgEst(coef_list,weights=weightslist)

coef_list2<-list(coef(summary(semFull[[2]]))[, 1], coef(summary(sem11[[2]]))[, 1])
weightslist<-c(summary(semFull, .progressBar = F)$AIC$AICc,summary(sem11, .progressBar = F)$AIC$AICc)
avgEst(coef_list2)
avgEst(coef_list2,weights=weightslist)

coef_list3<-list(coef(summary(semFull[[3]]))[, 1], coef(summary(sem11[[3]]))[, 1])
weightslist<-c(summary(semFull, .progressBar = F)$AIC$AICc,summary(sem11, .progressBar = F)$AIC$AICc)
avgEst(coef_list3)
avgEst(coef_list3,weights=weightslist)


## plot the SEM #####
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

coef<-subset(coef,Predictor!="~~prev_yr_floral_days_z")

coef<-left_join(coef,est_all,by=c("Predictor","Response"))

# remove variables not included in best model
coef<-coef %>% drop_na(Estimate_MdlAvg)

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
    fontcolor="black",     fillcolor="white",     fontname="Baskerville",
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

#graph %>% export_graph(file_name = "SEM_Panurginus_ineptus_modelavg.pdf")







############## (1) Panurginus cressoniellus ################

beedatafocal<-read.csv("focalbeedata_envdata_RMBLsexratios_2023-03-29.csv")
# subset to just include focal sites
sitelist<-c("Almont","Almont Curve","Beaver","CDOT","Copper","Davids","Elko","Gothic","Hill","Little","Lypps","Mexican Cut","Rustlers","Seans","Tuttle","Willey")
beedatafocal<-filter(beedatafocal, site %in% sitelist) 

beedatafocal<-filter(beedatafocal,genus_species=="Panurginus cressoniellus")

beedatafocal<-beedatafocal[,-c(56:71)]

##### Adding floral data: Inouye data #####

# read in floral data
flor<-read.csv("floral_data_annual_summaries_forsexratios_panurginus_cressoniellus.csv")
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

summary(semFull, .progressBar = F)
plot(semFull)

plot(beedatafocal$total_flowers_z,beedatafocal$floral_days_z)

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
    fontcolor="black",     fillcolor="white",     fontname="Baskerville",
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
          (1|site) ,
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
df3<-df3 %>% mutate(deltaAICc=AICc-min(AICc))

#write.csv(df3,"GoodnessOfFit_SEM_panurginus_cressoniellus.csv", row.names=FALSE)

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

#write.csv(coefs_df,"coefficients_SEM_panurginus_cressoniellus.csv", row.names=FALSE)


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


### Plot the best SEM (lowest AICc) #####
coef<-summary(sem13, .progressBar = F)$coef
coef
unique(coef$Response)
unique(coef$Predictor)

coef<-subset(coef,Predictor!="~~prev_yr_floral_days_z")

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
    fontcolor="black",
    fillcolor="white",
    fontname="Baskerville",
    penwidth=0.5,
    shape="rectangle",
    width=0.7,
    height=0.3,
    label = c("Female/male \nratio","Floral days \n(prior year)", "Floral sum \n(prior year)","Snowmelt date \n(prior year)","Summer precip. \n(prior year)"))

graph <-
  create_graph(
    nodes_df = ndf)

render_graph(graph)

# edf <-
#   create_edge_df(
#     from = coef$from,
#     to = coef$to,
#     rel = "leading_to",
#     #label = round(coef$Estimate,digits=2),
#     penwidth = abs(coef$Estimate*2),
#     fontsize=5.5,
#     fontcolor="brown3",
#     #color=coef$color,
#     color="black",
#     style=coef$style)

edf <-
  create_edge_df(
    from = coef$from,
    to = coef$to,
    rel = "leading_to",
    #label = round(coef$Estimate,digits=2),
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

#graph %>% export_graph(file_name = "SEM_Panurginus_cressoniellus_bestmodel.pdf")

### Perform model averaging #####

# sem full, sem 2, sem 11
coef_list<-list(coef(summary(semFull[[1]]))[, 1], coef(summary(sem2[[1]]))[, 1],coef(summary(sem11[[1]]))[, 1])
weightslist<-c(summary(semFull, .progressBar = F)$AIC$AICc,summary(sem2, .progressBar = F)$AIC$AICc,summary(sem11, .progressBar = F)$AIC$AICc)
avgEst(coef_list)
avgEst(coef_list,weights=weightslist)

coef_list2<-list(coef(summary(semFull[[2]]))[, 1], coef(summary(sem2[[2]]))[, 1],coef(summary(sem11[[2]]))[, 1])
weightslist<-c(summary(semFull, .progressBar = F)$AIC$AICc,summary(sem2, .progressBar = F)$AIC$AICc,summary(sem11, .progressBar = F)$AIC$AICc)
avgEst(coef_list2)
avgEst(coef_list2,weights=weightslist)

coef_list3<-list(coef(summary(semFull[[3]]))[, 1], coef(summary(sem2[[3]]))[, 1],coef(summary(sem11[[3]]))[, 1])
weightslist<-c(summary(semFull, .progressBar = F)$AIC$AICc,summary(sem2, .progressBar = F)$AIC$AICc,summary(sem11, .progressBar = F)$AIC$AICc)
avgEst(coef_list3)
avgEst(coef_list3,weights=weightslist)

coefs(sem8)
coefs(sem13)
coefs(sem14)

## plot the SEM #####
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

coef<-subset(coef,Predictor!="~~prev_yr_floral_days_z")


coef<-left_join(coef,est_all,by=c("Predictor","Response"))

# remove variables not included in best model
coef<-coef %>% drop_na(Estimate_MdlAvg)

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
    fontcolor="black",     fillcolor="white",     fontname="Baskerville",
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

#graph %>% export_graph(file_name = "SEM_Panurginus_cressoniellus_modelavg.pdf")

