# Code for calculating daily weather summaries, correcting errors,
# and merging with earlier billy barr station data
# For billy bar weather station data, 2016-2021
# Melanie Kazenel
# Created 2022-10-07
# ----------------------------------------

rm(list=ls(all=TRUE)) #give R a blank slate
setwd("~/Documents/*RMBLPhenology/Projects/BeeSexRatios/Analyses")

library(car)
library(reshape2)
library(date)
library(zoo)
library(lubridate)
library(readxl)
library(measurements)

##### 2016-2021 data formatting and cleaning #####

# read in data
data <- read.csv("billybarr_2016-2021.csv")

# subset to just include variables of interest
data<-data[,c("date","time","max_air_temp_C","min_air_temp_C","avg_air_temp_C","precip_mm","snow_depth_mm")]

# convert date to standard format
data$date<-mdy(data$date)

# create year, month, and day columns
data$year<-year(data$date)
data$month<-month(data$date)
data$day<-day(data$date)

# correct issues identified in data summaries below initially

# min air temp outlier
# 37374    2016-09-17    8:30    min air temp = -39.99
outlier<-subset(data,date=="2016-09-17")
# replace -39.99 with average value of hour before and hour after
data$min_air_temp_C[data$date == "2016-09-17" & data$time == "8:30"] <- mean(c(6.4270,8.1900))

# some missing precip data for part of 2021-09-01 through 2021-09-02; filled in as -9999
outlier<-subset(data,date=="2021-09-01" | date=="2021-09-02")
# for 2021-09-01, fill in missing data with zeros, based on Crested Butte 6.2N SCENIC station records
data$precip_mm[data$date == "2021-09-01" & data$precip_mm<0]<-0
# for 2021-09-02, backfill daily average with value from prior day (below)

# some negative snow depth values in 2016 and 2021
# look into this issue
outlier<-subset(data,snow_depth_mm<0)
outlier<-subset(data,date=="2016-01-04")
# replace with average value of hour before and hour after
data$snow_depth_mm[data$date == "2016-01-04" & data$snow_depth_mm == -9999] <- mean(c(907.8,907.6))
outlier<-subset(data,date=="2016-01-27" | date=="2016-01-28")
data$snow_depth_mm[data$date == "2016-01-27" & data$snow_depth_mm == -9999] <- mean(c(997.7,996.2))
outlier<-subset(data,date=="2016-03-31" | date=="2016-04-01")
data$snow_depth_mm[data$date == "2016-03-31" & data$snow_depth_mm == -9999]<-1107
data$snow_depth_mm[data$date == "2016-04-01" & data$snow_depth_mm == -9999]<-1091
outlier<-subset(data,date=="2016-10-23")
data$snow_depth_mm[data$date == "2016-10-23" & data$snow_depth_mm == -9999]<-0

outlier<-subset(data,snow_depth_mm<0 & year==2016)
# sensor seemingly down from
# 2016-05-30 6:40 to
# 2016-09-27 16:20
# fill with NAs for now; see below

# check 2016 high outlier snow depth outlier
outlier<-subset(data,snow_depth_mm>3000 & year==2016)
# 9/17 through 9/27
# fill in with NAs below

# handle high precip outliers
# 2018-09-22
outlier<-subset(data,date=="2018-09-22")
ggplot(data=outlier,aes(x=time,y=precip_mm)) + geom_point()
# replace weird precip values with 0, based on Crested Butte 6.2N SCENIC station data
data$precip_mm[data$date == "2018-09-22" & data$precip_mm >50]<-0

# 2019-10-12
outlier<-subset(data,date=="2019-10-12")
ggplot(data=outlier,aes(x=time,y=precip_mm)) + geom_point()
# replace weird precip values with 0, based on Crested Butte 6.2N SCENIC station data
data$precip_mm[data$date == "2019-10-12" & data$precip_mm >19]<-0

# 2019-11-19
outlier<-subset(data,date=="2019-11-19")
ggplot(data=outlier,aes(x=time,y=precip_mm)) + geom_point()
# replace weird precip values with 0, based on Crested Butte 6.2N SCENIC station data
data$precip_mm[data$date == "2019-11-19" & data$precip_mm >10]<-0

# 2020-10-22
outlier<-subset(data,date=="2020-10-22")
ggplot(data=outlier,aes(x=time,y=precip_mm)) + geom_point()
# replace weird precip values with 0, based on Crested Butte 6.2N SCENIC station data
data$precip_mm[data$date == "2020-10-22" & data$precip_mm >50]<-0

# 2021-07-22
outlier<-subset(data,date=="2021-07-22")
ggplot(data=outlier,aes(x=time,y=precip_mm)) + geom_point()
# replace weird precip values with 0, based on Crested Butte 6.2N SCENIC station data
data$precip_mm[data$date == "2021-07-22" & data$precip_mm >20]<-0

# 2020-05-01
outlier<-subset(data,date=="2020-05-01")
data$precip_mm[data$date == "2020-05-01" & data$precip_mm >20]<-0

# 2021-06-20
outlier<-subset(data,date=="2021-06-20")
data$precip_mm[data$date == "2021-06-20" & data$precip_mm >20]<-0

# 2020-02-15 through 2020-02-27
outlier<-subset(data,year==2020 & month==2)
outlier<-subset(outlier,precip_mm>200)
data$precip_mm[data$year == 2020 & data$precip_mm >200]<-0


# format date

#Add column for Julian day: Given a month, day, and year, returns the number of days since January 1, 1960
data$jul_day<-mdy.date(data$month, data$day, data$year, nineteen= TRUE, fillday = FALSE, fillmonth = FALSE)

# calculate values for each day
billybarrmelt<-melt(data,id=c("jul_day","month","day","year","date","time"),variable.name="variable")
billybarrmelt$value <- as.numeric(billybarrmelt$value)

billybarrdaily<- data %>% group_by(date,month,day,year,jul_day) %>% 
  summarise(MaxAirTemp=max(max_air_temp_C),
            MinAirTemp=min(min_air_temp_C),
            AvgAirTemp=mean(avg_air_temp_C),
            Precip=sum(precip_mm),
            MaxSnowDepth=max(snow_depth_mm),
            AvgSnowDepth=mean(snow_depth_mm)
            )

# calculate GDD
billybarrdaily$GDD_calc<-(billybarrdaily$MaxAirTemp+billybarrdaily$MinAirTemp)/2

# add day of year column
billybarrdaily$doy<-yday(billybarrdaily$date)

# look at data
ggplot(data=billybarrdaily, aes(x=doy,y=MaxAirTemp))+
  geom_point() +
  facet_wrap(~year)

ggplot(data=billybarrdaily, aes(x=doy,y=MinAirTemp))+
  geom_point() +
  facet_wrap(~year)
# ISSUE BELOW FIXED
# # outlier of around -40 in 2016
# # look into this issue
# outlier<-subset(data,min_air_temp_C<=-30 & year==2016)
# # 37374    2016-09-17    8:30    min air temp = -39.99
# # look at all data from that date
# outlier<-subset(data,date=="2016-09-17")

ggplot(data=billybarrdaily, aes(x=doy,y=AvgAirTemp))+
  geom_point() +
  facet_wrap(~year)

ggplot(data=billybarrdaily, aes(x=doy,y=Precip))+
  geom_point() +
  facet_wrap(~year)
# some weird negative values in 2021
# look into this issue
outlier<-subset(data,precip_mm<0)
# some missing precip data for part of 2021-09-01 through 2021-09-02; filled in as -9999
# For 2021-09-02, backfill daily average with value from prior day
outlier<-subset(billybarrdaily,year==2021 & month==9)
billybarrdaily$Precip[billybarrdaily$date == "2021-09-02"]<-8.382

outlier<-subset(billybarrdaily,Precip>60)
# 2018-09-22
# 2019-10-12
# 2019-11-19
# 2020-10-22
# 2021-07-22
# 2020-05-01
# 2021-06-20
# 2020-02-15 through 2020-02-27


ggplot(data=billybarrdaily, aes(x=doy,y=MaxSnowDepth))+
  geom_point() +
  facet_wrap(~year)
# some negative snow depth values in 2016 and 2021
# look into this issue
outlier<-subset(data,snow_depth_mm<0)

# sensor seemingly down from
# 2016-05-30 6:40 to
# 2016-09-27 16:20
# fill with NAs for now; see below
billybarrdaily$MaxSnowDepth[billybarrdaily$year==2016 & billybarrdaily$MaxSnowDepth<0]<-NA
billybarrdaily$AvgSnowDepth[billybarrdaily$year==2016 & billybarrdaily$AvgSnowDepth<0]<-NA
billybarrdaily$MaxSnowDepth[billybarrdaily$year==2016 & billybarrdaily$MaxSnowDepth>3000]<-NA
# treat weird data from 2021 in same way
billybarrdaily$MaxSnowDepth[billybarrdaily$year==2021 & billybarrdaily$MaxSnowDepth<0]<-NA
billybarrdaily$MaxSnowDepth[billybarrdaily$year==2021 & billybarrdaily$MaxSnowDepth>3000]<-NA
billybarrdaily$AvgSnowDepth[billybarrdaily$year==2021 & billybarrdaily$AvgSnowDepth<0]<-NA

# check 2016 high outlier snow depth outlier
outlier<-subset(data,snow_depth_mm>3000 & year==2016)


ggplot(data=billybarrdaily, aes(x=doy,y=AvgSnowDepth))+
  geom_point() +
  facet_wrap(~year)

write.csv(billybarrdaily, "billybarr_2016-2021_corrected.csv")


##### Merging 2016-2021 data with earlier data from Stemkovski et al. 2020 paper #####
# modified code from Michael Stemkovski

# read in data
read.excel <- function(...) return(as.data.frame(read_excel(...)))
weather_raw <- read.excel("aawx_Michael.xlsx") #this came from billy barr's weather station
bare_ground <- read.csv("billy_barr_bare_ground.csv") #this came from billy barr's observations too

# basic cleaning
bare_ground$date_bare_ground<-paste(bare_ground$year,bare_ground$month_day_bare_ground,sep="-")
bare_ground$doy <- yday(bare_ground$date_bare_ground)
weather <- weather_raw[3:nrow(weather_raw),]
colnames(weather) <- c("year","month","day","temp_min","temp_max","snow_cm","snow_water_in","snow_depth_cm","rain_in")
weather$snow_water_in <- conv_unit(as.numeric(weather$snow_water_in),"inch","cm")
weather$rain_in <- conv_unit(as.numeric(weather$rain_in),"inch","cm")
colnames(weather) <- c("water_year","wy_month","day","temp_min","temp_max","snow_cm","snow_water_cm","snow_depth_cm","rain_cm")
weather$rain_cm <- sapply(weather$rain_cm,function(x) ifelse(is.na(x), return(0), return(x)))

# date reformatting
get.year <- function(water_year,wy_month){
  if(wy_month >= 9 & wy_month <= 12){
    year <- paste("20",substr(water_year,nchar(water_year)-3,nchar(water_year)-2),sep="")
  } else {
    year <- paste("20",substr(water_year,nchar(water_year)-1,nchar(water_year)),sep="")
  }
  return(as.numeric(year))
}

head(weather,20)
str(weather)

month_letters <- toupper(letters[1:12])
month_nums <- c(9:12,1:8)

month_letters
month_nums

weather$month <- sapply(weather$wy_month, function(x) return(month_nums[which(month_letters == x)]), USE.NAMES = F)
weather$year <- apply(weather[,c("water_year","month")], 1, function(x) get.year(x[1], as.numeric(x[2])))
weather$date <- ymd(apply(weather[,c("year","month","day")], 1, function(x) paste(x[1],x[2],x[3],sep="-")))

# check for presence of all days in 2009-2018
days_check<-weather %>% filter(year>=2009) %>% group_by(year) %>% summarise(days=n())
# December 2018 missing

# read in later data
billy_late<-read.csv("billybarr_2016-2021_corrected.csv")
# subset to just include data from December 2018 and 2019-2020
dec18<-filter(billy_late,year==2018 & month==12)
weather19_20<-filter(billy_late,year==2019 | year==2020)
weather_late<-bind_rows(dec18,weather19_20)

# add water year as column to billy_late data frame
# function to calculate water year
wtr_yr <- function(dates, start_month = 10) {
  # Convert possible character vector into date
  d1 = as.Date(dates)
  # Year offset
  offset = ifelse(as.integer(format(d1, "%m")) < start_month, 0, 1)
  # Water year
  adj.year = as.integer(format(d1, "%Y")) + offset
  # Return the water year
  return(adj.year)
}

# convert dates to standard format
weather_late$date<-mdy(weather_late$date)
# run function and add output to data frame
weather_late$water_year<-wtr_yr(weather_late$date)

# remove columns from late data not contained in early
weather_late_formerge<-weather_late[,-c(5,8:11,13)]
# rename columns
names(weather_late_formerge)[5:6]<-c("temp_max","temp_min")

# format older weather data frame for merging
weather$water_year<-substr(weather$water_year,1,4)
weather2<-weather[,-c(2)]
weather2$day<-as.integer(weather2$day)
weather2$water_year<-as.numeric(weather2$water_year)
weather2$temp_max<-as.numeric(weather2$temp_max)
weather2$temp_min<-as.numeric(weather2$temp_min)
# calculate GDD
weather2$GDD_calc<-(weather2$temp_max+weather2$temp_min)/2
# reorder columns
weather2<-weather2[,c(11,9,2,10,1,4,3,5:8,12)]

# format late dataset for merging
weather_late_formerge$snow_cm<-NA
weather_late_formerge$snow_water_cm<-NA
weather_late_formerge$snow_depth_cm<-NA
weather_late_formerge$rain_cm<-NA
weather_late_formerge<-weather_late_formerge[,c(1:4,8,5:6,9:12,7)]

# merge data frames
weather_all<-bind_rows(weather2,weather_late_formerge)

# check for presence of all days in target years
days_check<-weather_all %>% filter(year>=2007) %>% group_by(year) %>% summarise(days=n())
# all there!

# save data
write.csv(filter(weather_all,year>=2007),"weather_billybarr_2007-2020.csv",row.names = FALSE)





##### Create climate summaries needed for bee sex ratio analyses ##### 

# read in data
climate<-read.csv("weather_billybarr_2007-2020.csv")

# target variables for 2008-2018
# May-September cumulative precipitation
# May-September degree days above 0 C
# October-April winter snow cover
# Date of bare ground

clim_target<-filter(climate,year>=2008 & year<=2018)

# create column to use in calculating accumulated degree days above 0
clim_target$DD_abovezero<-ifelse(clim_target$GDD_calc>0, clim_target$GDD_calc, 0)

# calculate May-September cumulative rainfall and degree days above 0
clim_summary<-clim_target %>% filter(month>=5 & month<=9) %>% group_by(year) %>% summarise(accum_summer_precip_cm=sum(rain_cm),accum_summer_dd=sum(DD_abovezero))

# winter snow cover
snow<-clim_target %>% filter(month %in% c(10,11,12,1,2,3,4)) %>% group_by(water_year) %>% summarise(accum_winter_snowfall_cm=sum(snow_cm),mean_snow_depth_cm=mean(snow_depth_cm))
snow<-snow[-13,]
# change water year column to year for sake of merging with other data set
names(snow)[1]<-"year"
snow$year<-snow$year+1

# merge data frames
clim_summ_all<-left_join(clim_summary,snow,by="year")

# add bare ground data
bare_ground<-read.csv("billy_barr_bare_ground.csv")
bare_ground$date_bare_ground<-paste(bare_ground$year,bare_ground$month_day_bare_ground,sep="-")
clim_summ_all<-left_join(clim_summ_all,bare_ground[,c(2,5,8)])

write.csv(clim_summ_all,"weather_summaries_billybarr_forsexratios.csv",row.names = FALSE)
