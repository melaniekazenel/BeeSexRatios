# Working with Inouye flowering phenology data for bee sex ratios project
# Melanie Kazenel
# Created 2022-10-13
# -------------------------------------------------------

setwd("~/Documents/*RMBLPhenology/Projects/BeeSexRatios/Analyses")

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

# write.csv(flor_annual,"floral_data_annual_summaries_forsexratios.csv",row.names = FALSE)
