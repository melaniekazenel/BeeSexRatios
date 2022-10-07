# Code for calculating daily weather summaries and correcting errors
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

#Upload all data
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



write.csv(billybarrdaily, "billybarrdaily_corrected.csv")
