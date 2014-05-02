######################### load packages #####################################

library(stats, lib="D:/Rlibrary")
library(plyr, lib="D:/Rlibrary")
library(reshape, lib="D:/Rlibrary")
library(mvtnorm, lib="D:/Rlibrary")
library(mc2d, lib="D:/Rlibrary")
library(abind, lib="D:/Rlibrary")
library(ggplot2, lib="D:/Rlibrary")
library(ggthemes, lib="D:/Rlibrary")
library(eha, lib="D:/Rlibrary")
library(msm, lib="D:/Rlibrary")
library(tkrplot, lib="D:/Rlibrary")
library(rriskDistributions, lib="D:/Rlibrary")
library(ggmap, lib="D:/Rlibrary")
library(mapproj, lib="D:/Rlibrary")
library(rgeos, lib="D:/Rlibrary")
library(maptools, lib="D:/Rlibrary")
library(raster, lib="D:/Rlibrary")
library(rgdal, lib="D:/Rlibrary")
library(reshape2, lib="D:/Rlibrary")
library(plyr, lib="D:/Rlibrary")

########################import epi data####################################

# endemicity status
status<-read.csv(file="C:/Users/fitzpatrickc/Dropbox/yaws/yaws-verbose.csv")
status<-subset(status, select=c(YEAR..CODE.,COUNTRY..CODE.,COUNTRY..DISPLAY.,Display.Value))
# reported number of cases
cases<-read.csv(file="C:/Users/fitzpatrickc/Dropbox/yaws/yaws-verbose2.csv")
cases<-subset(cases, select=c(YEAR..CODE.,COUNTRY..CODE.,COUNTRY..DISPLAY.,Numeric))
#average number of cases per year
sum(cases$Numeric[cases$YEAR..CODE.==2012 | (cases$YEAR..CODE.==2008 & cases$COUNTRY..CODE.=="COD")], na.rm=T)

########################estimate population at risk##############################

# http://gecon.yale.edu/g-econ-project-yale-university-september-2009
gecon<-read.csv(file="C:/Users/fitzpatrickc/Dropbox/yaws/data/gecon.csv", header=TRUE, na.strings=c("#NA","#DIV/0"))
gecon$PPP2005_40<-as.numeric(as.character(gecon$PPP2005_40))*1000000000
gecon$MER2005_40<-as.numeric(as.character(gecon$MER2005_40))*1000000000
gecon$PPP2005_pp<-gecon$PPP2005_40/gecon$POPGPW_2005_40
gecon$MER2005_pp<-gecon$MER2005_40/gecon$POPGPW_2005_40
gecon$POPDEN_2005<-gecon$POPGPW_2005_40/as.numeric(as.character(gecon$AREA))
levels(gecon$COUNTRY)[levels(gecon$COUNTRY)=="Democratic Republic of Congo"] <- "Democratic Republic of the Congo"
levels(gecon$COUNTRY)[levels(gecon$COUNTRY)=="Bolivia"] <- "Bolivia (Plurinational State of)"
levels(gecon$COUNTRY)[levels(gecon$COUNTRY)=="Cook Is."] <- "Cook Islands"
levels(gecon$COUNTRY)[levels(gecon$COUNTRY)=="Guinea Bissau"] <- "Guinea-Bissau"
levels(gecon$COUNTRY)[levels(gecon$COUNTRY)=="Laos"] <- "Lao People's Democratic Republic"
levels(gecon$COUNTRY)[levels(gecon$COUNTRY)=="St. Lucia"] <- "Saint Lucia"
levels(gecon$COUNTRY)[levels(gecon$COUNTRY)=="St. Vincent and the Grenadines"] <- "Saint Vincent and the Grenadines"
levels(gecon$COUNTRY)[levels(gecon$COUNTRY)=="Tanzania"] <- "United Republic of Tanzania"
levels(gecon$COUNTRY)[levels(gecon$COUNTRY)=="Timor Leste"] <- "Timor-Leste"
levels(gecon$COUNTRY)[levels(gecon$COUNTRY)=="Venezuela"] <- "Venezuela (Bolivarian Republic of)"
levels(gecon$COUNTRY)[levels(gecon$COUNTRY)=="Vietnam"] <- "Viet Nam"
unique(gecon$COUNTRY)
# calculate gdp per capita and per square kilometer
sum(gecon$PPP2005_40[which(gecon$POPGPW_2005_40>0)], na.rm=T)/sum(gecon$POPGPW_2005_40, na.rm=T)
sum(gecon$PPP2005_40[which(gecon$POPGPW_2005_40>0)], na.rm=T)/sum(as.numeric(gecon$AREA[which(gecon$POPGPW_2005_40>0)]))
sum(gecon$MER2005_40[which(gecon$POPGPW_2005_40>0)], na.rm=T)/sum(gecon$POPGPW_2005_40, na.rm=T)
sum(gecon$MER2005_40[which(gecon$POPGPW_2005_40>0)], na.rm=T)/sum(as.numeric(gecon$AREA[which(gecon$POPGPW_2005_40>0)]))
sum(gecon$MER2005_40)
# subset based on environmental and economic variables
# We include the entire populations of Solomon Islands and Vanuatu
list_bigsub<-as.list(subset(status, select="COUNTRY..DISPLAY.", Display.Value=="Currently endemic" | Display.Value=="Previously endemic (current status unknown)"))
temp<-subset(gecon, lapply(COUNTRY, function (x) any(x %in% list_bigsub$COUNTRY..DISPLAY.))==TRUE)
unique(temp$COUNTRY)
gecon_bigsub<-subset(gecon, (lapply(COUNTRY, function (x) any(x %in% list_bigsub$COUNTRY..DISPLAY.))==TRUE & POPDEN_2005<100 & PREC_NEW>500 & TEMP_NEW>20 & (MATTVEG==1 | MATTVEG==2 | MATTVEG==7 | MATTVEG==9 | MATTVEG==15)) | COUNTRY=="Vanuatu" | COUNTRY=="Solomon Islands")
unique(gecon_bigsub$COUNTRY)
sum(gecon_bigsub$POPGPW_2005_40)
# endemic countries only
gecon_sub<-subset(gecon, COUNTRY=="Benin" | COUNTRY=="Cameroon" | COUNTRY=="Central African Republic" | COUNTRY=="Congo" | COUNTRY=="Cote d'Ivoire"| COUNTRY=="Democratic Republic of the Congo"| COUNTRY=="Ghana"| COUNTRY=="Indonesia"| COUNTRY=="Papua New Guinea"| COUNTRY=="Solomon Islands"| COUNTRY=="Togo"| COUNTRY=="Vanuatu")
sum(gecon_sub$PPP2005_40[which(gecon_sub$POPGPW_2005_40>0)], na.rm=T)/sum(gecon_sub$POPGPW_2005_40, na.rm=T)
sum(gecon_sub$PPP2005_40[which(gecon_sub$POPGPW_2005_40>0)], na.rm=T)/sum(as.numeric(gecon_sub$AREA[which(gecon_sub$POPGPW_2005_40>0)]))
sum(gecon_sub$MER2005_40[which(gecon_sub$POPGPW_2005_40>0)], na.rm=T)/sum(gecon_sub$POPGPW_2005_40, na.rm=T)
sum(gecon_sub$MER2005_40[which(gecon_sub$POPGPW_2005_40>0)], na.rm=T)/sum(as.numeric(gecon_sub$AREA[which(gecon_sub$POPGPW_2005_40>0)]))
sum(gecon_sub$MER2005_40)
gecon_subsub<-subset(gecon_sub, POPGPW_2005_40>0 & POPDEN_2005<100 & PREC_NEW>500 & TEMP_NEW>20 & (MATTVEG==1 | MATTVEG==2 | MATTVEG==7 | MATTVEG==9 | MATTVEG==15) | (COUNTRY=="Vanuatu" & POPGPW_2005_40>0) | (COUNTRY=="Solomon Islands" & POPGPW_2005_40>0))
sum(gecon_subsub$POPGPW_2005_40)
# we limit also to Hackett (1951) latitudes and longitudes for DRC
gecon_subsub<-subset(gecon_subsub, COUNTRY!="Democratic Republic of the Congo" | (COUNTRY=="Democratic Republic of the Congo" & ((LAT<=2 & LAT>=-1 & LONGITUDE<=26)|(LAT<=2 & LAT>=-9 & LONGITUDE>=26))))
sum(gecon_subsub$POPGPW_2005_40)
# PAR in countries of unknown endemicity 
par_unknown<-sum(gecon_bigsub$POPGPW_2005_40)-sum(gecon_subsub$POPGPW_2005_40)
par_unknown
# calculate gdp per capita and per square kilometer
sum(gecon_subsub$PPP2005_40)/sum(gecon_subsub$POPGPW_2005_40)
sum(gecon_subsub$PPP2005_40)/sum(as.numeric(gecon_subsub$AREA))
sum(gecon_subsub$MER2005_40)/sum(gecon_subsub$POPGPW_2005_40)
sum(gecon_subsub$MER2005_40)/sum(as.numeric(gecon_subsub$AREA))
sum(gecon_subsub$MER2005_40)
sum(gecon_subsub$MER2005_40)/sum(gecon_sub$MER2005_40)
# aggregate PAR by country
gecon_sub_byc<-aggregate(cbind(gecon_subsub$POPGPW_2005_40,gecon_subsub$PPP2005_40,gecon_subsub$AREA), by=list(gecon_subsub$COUNTRY), sum, na.rm=TRUE)
names(gecon_sub_byc)<-c("country","par","gdp","area")
gecon_sub_byc$par<-as.numeric(gecon_sub_byc$par)
# merge in iso3 codes
iso3<-read.csv(file="C:/Users/fitzpatrickc/Dropbox/yaws/data/iso3_list.csv", header=TRUE)
gecon_sub_byc<-merge(gecon_sub_byc,iso3,by="country")
gecon_sub_byc$year<-2005
# format PAR estimates for outsheet
yaw.par<-as.data.frame(cbind(as.character(gecon_sub_byc$country),as.character(gecon_sub_byc$iso3),gecon_sub_byc$year,round(gecon_sub_byc$par, digits=0),round(gecon_sub_byc$area, digits=0)))
names(yaw.par)<-c("country","iso3","year","par","area")
levels(yaw.par$par)
yaw.par$par<-as.numeric(as.character(yaw.par$par))
# adjust par to 2015 using rural population numbers
pop_rur<-read.csv(file="C:/Users/fitzpatrickc/Dropbox/yaws/data/pop_rur.csv", header=TRUE)
names(pop_rur)[names(pop_rur) == "unpop_name"] <- "country"
levels(pop_rur$country)[57] <- "Cote d'Ivoire"
pop_rur$pop_rur_gr<-as.numeric(pop_rur$X2015/pop_rur$X2005)
yaw.par<-merge(yaw.par, subset(pop_rur, select=c("country","pop_rur_gr")), by=c("country"))
yaw.par$par<-yaw.par$par*yaw.par$pop_rur_gr
yaw.par$year<-2015
sum(yaw.par$par)
# write to csv
write.csv(yaw.par, file="C:/Users/fitzpatrickc/Dropbox/yaws/yaw_par.csv", row.names = FALSE)

############## map population at risk using Google maps ##############################

parmap<-qmap(location="India",zoom=2,maptype = c("satellite"),color="color",source="google")
breaks<-as.vector(floor(quantile(gecon_subsub$POPGPW_2005_40,c(0.2,0.4,0.6,0.8,1.00))))
parmap + geom_point(aes(x=LONGITUDE, y=LAT, color=POPGPW_2005_40), data=gecon_subsub) + scale_colour_gradient(low="orange", high = "red", guide="legend", name="population at risk", breaks=breaks) + scale_y_continuous(limits = c(-45,45))

############## map population at risk using Natural Earth ##############################

# Ctrl+Shift+C to uncomment lines
# nat.earth<-stack(list(x='D:/HYP_50M_SR_W.tif'))
# ne_lakes <- readOGR('D:/ne_10m_lakes.shp','ne_10m_lakes')
# quick.subset <- function(x, longlat){
#   # longlat should be a vector of four values: c(xmin, xmax, ymin, ymax)
#   x@data$id <- rownames(x@data)
#   x.f = fortify(x, region="id")
#   x.join = join(x.f, x@data, by="id")
#     x.subset <- subset(x.join, x.join$long > longlat[1] & x.join$long < longlat[2] & x.join$lat > longlat[3] & x.join$lat < longlat[4])
#     x.subset
# }
# domain <- c(-15, 170, -25, 15)
# lakes.subset <- quick.subset(ne_lakes, domain)
# nat.crop<-crop(x=nat.earth,y=extent(domain))
# rast.table<-data.frame(xyFromCell(nat.crop, 1:ncell(nat.crop)), getValues(nat.crop/255))
# rast.table$rgb<-with(rast.table, rgb(x_1,x_2,x_3,1))                                                                                                                         
# breaks<-as.vector(floor(quantile(gecon_subsub$POPGPW_2005_40,c(0.2,0.4,0.6,0.8,1.00))))
# ggplot(data = rast.table, aes(x = x, y = y)) +
#   geom_tile(fill = rast.table$rgb) +
#   #geom_polygon(data=lakes.subset, aes(x = long, y = lat, group = group), fill = '#ADD8E6') +
#   scale_alpha_discrete(range=c(1,0)) +
#   # geom_path(data=lakes.subset, aes(x = long, y = lat, group = group), color = 'blue') +
#   scale_x_continuous(expand=c(0,0)) +
#   scale_y_continuous(expand=c(0,0)) +
#   xlab('') + ylab('') +
#   geom_point(aes(x=LONGITUDE, y=LAT, color=POPGPW_2005_40), data=gecon_subsub) +
#   scale_colour_gradient(low="orange", high = "red", guide="legend", name="population at risk", breaks=breaks)

#####################set up model parameters#########################

#Number of iterations for uncertainty and variability dimensions
ndunc<-1000
ndunc(ndunc)
# we don't use the variability dimension in this file
ndvar<-1
ndvar(ndvar)
# period of analysis
years<-2050-2015+1
# cycle length (in years)
cycle_len<-0.5
# burn in period
burn<-10/cycle_len
# number of cycles including burn-in
cycle_num<-years/cycle_len+burn
# countries
ctry_num<-nrow(yaw.par)
# create list of iso3 codes
iso3<-yaw.par$iso3
# make empty 4D array for later use
array4d<-array(NA,dim=c(ctry_num, cycle_num, ndunc, ndvar), dimnames=list(iso3=iso3,cycle=1:cycle_num,ndunc=1:ndunc,ndvar=1:ndvar))
array4dz<-array(0,dim=c(ctry_num, cycle_num, ndunc, ndvar), dimnames=list(iso3=iso3,cycle=1:cycle_num,ndunc=1:ndunc,ndvar=1:ndvar))
dim(array4d)
dimnames(array4d)

#############################import population data###########################

# population by age
pop_by_age<-read.csv(file="C:/Users/fitzpatrickc/Dropbox/yaws/data/pop_by_age.csv", header=TRUE)
names(pop_by_age)[names(pop_by_age) == "unpop_name"] <- "country"
levels(pop_by_age$country)[49] <- "Cote d'Ivoire"
# use 2015 populations
pop_by_age<-pop_by_age[pop_by_age$year==2015,]
pop_by_age$pop_0_4<-rowSums(pop_by_age[,c("age_0","age_1","age_2","age_3","age_4")])
pop_by_age$pop_5_9<-rowSums(pop_by_age[,c("age_5","age_6","age_7","age_8","age_9")])
pop_by_age$pop_10_14<-rowSums(pop_by_age[,c("age_10","age_11","age_12","age_13","age_14")])
pop_by_age$pop_15_99<-rowSums(pop_by_age[,which(colnames(pop_by_age)=="age_15"):which(colnames(pop_by_age)=="age_80.")])
pop_by_age$pop_0_pc<-pop_by_age$age_0
pop_by_age<-pop_by_age[,c("country", "pop_0_4","pop_5_9","pop_10_14","pop_15_99","pop_0_pc")]
pop_by_age$pop_total<-pop_by_age$pop_0_4+pop_by_age$pop_5_9+pop_by_age$pop_10_14+pop_by_age$pop_15_99
# eligible population (more than 6 months)
pop_by_age$pop_elig<-1-(pop_by_age$pop_0_pc/2)/pop_by_age$pop_total
# merge and drop
yaw.pop<-merge(yaw.par, subset(pop_by_age, select=c("country","pop_0_4","pop_5_9","pop_10_14","pop_15_99","pop_total","pop_elig")), by=c("country"))

#########################estimate population at risk and treatment targets ####################

# population at risk % (country specific, not cycle specific, and uncertain)
par_p<-array(NA, dim=c(ctry_num,ndunc), dimnames=list(iso3=iso3,ndunc=1:ndunc))
# minimum is 5%, maximum is 10% or the value obtained by G-Econ environmental conditions
yaw.pop$par_pc<-yaw.pop$par/(yaw.pop$pop_total*1000)
for (i in 1:ctry_num) {par_p[i,]<-runif(ndunc,min=0.05, max=max(0.10, yaw.pop$par_pc[i]))}
par_p["SLB",]<-runif(ndunc,min=1.00, max=1.00)
par_p["VUT",]<-runif(ndunc,min=1.00, max=1.00)
par_p["PNG",1:10]
# turn into full 4D array
temp<-array4d
for (i in 1:ctry_num) {
  for (j in 1:ndunc) {temp[i,,j,]<-par_p[i,j]}
}
par_p<-temp
# population at risk, by age group
par_0_4<-yaw.pop$pop_0_4*par_p*1000
par_5_9<-yaw.pop$pop_5_9*par_p*1000
par_10_14<-yaw.pop$pop_10_14*par_p*1000
par_15_99<-yaw.pop$pop_15_99*par_p*1000

##############################enter baseline number of cases###############################

# case detection rate at baseline (country but not year specific)
temp<-array(runif(ndunc, min=0.50, max=0.80),dim=c(ctry_num, ndunc))
cdr<-array4d
for (i in 1:ctry_num) {
  for (j in 1:ndunc){
    cdr[i,,j,]<-temp[i,j]
  }
}

# new cases
# use latest reported value
new<-array4dz
new["BEN",1:2,,]<-max(cases[which(cases$COUNTRY..CODE.=="BEN"),]$Numeric, na.rm=T)*cycle_len/cdr["BEN",1:2,,]
new["CMR",1:2,,]<-max(cases[which(cases$COUNTRY..CODE.=="CMR"),]$Numeric, na.rm=T)*cycle_len/cdr["CMR",1:2,,]
new["CAF",1:2,,]<-max(cases[which(cases$COUNTRY..CODE.=="CAF"),]$Numeric, na.rm=T)*cycle_len/cdr["CAF",1:2,,]
new["COG",1:2,,]<-max(cases[which(cases$COUNTRY..CODE.=="COG"),]$Numeric, na.rm=T)*cycle_len/cdr["COG",1:2,,]
new["CIV",1:2,,]<-max(cases[which(cases$COUNTRY..CODE.=="CIV"),]$Numeric, na.rm=T)*cycle_len/cdr["CIV",1:2,,]
new["COD",1:2,,]<-max(cases[which(cases$COUNTRY..CODE.=="COD"),]$Numeric, na.rm=T)*cycle_len/cdr["COD",1:2,,]
new["GHA",1:2,,]<-max(cases[which(cases$COUNTRY..CODE.=="GHA"),]$Numeric, na.rm=T)*cycle_len/cdr["GHA",1:2,,]
new["IDN",1:2,,]<-max(cases[which(cases$COUNTRY..CODE.=="IDN"),]$Numeric, na.rm=T)*cycle_len/cdr["IDN",1:2,,]
new["PNG",1:2,,]<-max(cases[which(cases$COUNTRY..CODE.=="PNG"),]$Numeric, na.rm=T)*cycle_len/cdr["PNG",1:2,,]
new["SLB",1:2,,]<-max(cases[which(cases$COUNTRY..CODE.=="SLB"),]$Numeric, na.rm=T)*cycle_len/cdr["SLB",1:2,,]
new["TGO",1:2,,]<-max(cases[which(cases$COUNTRY..CODE.=="TGO"),]$Numeric, na.rm=T)*cycle_len/cdr["TGO",1:2,,]
new["VUT",1:2,,]<-max(cases[which(cases$COUNTRY..CODE.=="VUT"),]$Numeric, na.rm=T)*cycle_len/cdr["VUT",1:2,,]
pri<-array4dz
sec<-array4dz
lat<-array4dz
ter<-array4dz
dea<-array4dz

################################ define other effect model parameters#############################

#discount rate on effects (country but not year specific)
temp<-array(1-exp(-runif(ndunc,min=0.00, max=0.03)*cycle_len),dim=c(ctry_num, ndunc))
disc_e<-array4d
for (i in 1:ctry_num) {
  for (j in 1:ndunc){
    disc_e[i,,j,]<-temp[i,j]
  }
}

# 45-year probability of death (between 15 and 60 years of age)
admort<-array(NA, dim=c(ctry_num), dimnames=list(iso3=iso3))
admort["BEN"]<-319/1000
admort["CMR"]<-413/1000
admort["CAF"]<-464/1000
admort["COG"]<-365/1000
admort["CIV"]<-495/1000
admort["COD"]<-387/1000
admort["GHA"]<-332/1000
admort["IDN"]<-190/1000
admort["PNG"]<-248/1000
admort["SLB"]<-145/1000
admort["TGO"]<-307/1000
admort["VUT"]<-180/1000
# adult mortality within population at risk
temp<-runif(ndunc,min=1, max=1.2)
admort_par<-array4d
for (i in 1:ndunc) {admort_par[,,i,]<-temp[i]*admort}
# yearly rate of death
rate_death<--log(1-admort_par)/(60-15)
# 1-cycle probability of death.
prob_death<-1-exp(-rate_death*cycle_len)

# exit probability (from population at risk / susceptible population)
temp<-1-exp(-runif(ndunc,min=0.02,max=0.07)*cycle_len)
exit_prob<-array4d
for (i in 1:ndunc) {exit_prob[,,i,]<-temp[i]}

# reproduction number
temp<-runif(ndunc, min=0.9, max=1.0)
repro_num<-array4d
for (i in 1:ndunc) {repro_num[,,i,]<-temp[i]}
# during the burn-in use a reproduction number of 1
repro_num[,1:burn,,]<-1

# generation time in cycles
temp<-runif(ndunc, min=0.08, max=5)/cycle_len 
gen_time<-array4d
for (i in 1:ndunc) {gen_time[,,i,]<-temp[i]}

# long-term proportion of new cases that enter primary stage in absence of treatment
new_pri_p<-0.999
# long-term proportion of primary cases that enter secondary stage in absence of treatment
pri_sec_p<-0.999
# long-term proportion of secondary cases that enter tertiary stage
temp<-runif(ndunc, min=.08, max=.12)
sec_ter_p<-array4d
for (i in 1:ndunc) {sec_ter_p[,,i,]<-temp[i]}

# duration in early stages (to establish baseline values to start burn-in)
temp<-runif(ndunc,min=0.25,max=0.5)/cycle_len
dur_pri_sec<-array4d
for (i in 1:ndunc) {dur_pri_sec[,,i,]<-temp[i]}
temp<-runif(ndunc,min=5,max=10)/cycle_len
dur_pri_ter<-array4d
for (i in 1:ndunc) {dur_pri_ter[,,i,]<-temp[i]}
dur_sec_ter<-dur_pri_ter-dur_pri_sec

# baseline proportion of non-tertiary cases that are secondary
sec_p<-dur_sec_ter/dur_pri_ter
# ratio of latent to non-latent among secondary stage cases
temp<-runif(ndunc, min=2, max=6)
lat_r<-array4d
for (i in 1:ndunc) {lat_r[,,i,]<-temp[i]}
# proportion of secondary stage cases that are latent
lat_p<-lat_r/(lat_r+1)

# transition probabilities
# new is a compartment to keep track of transmission
# new cases are produced for each primary cases, but move into the primary compartment over the generation time
new_death<-0
rate_new_pri<--log(1-new_pri_p)/gen_time
new_pri<-(1-exp(-rate_new_pri))*(1-new_death)
new_new<-1-new_death-new_pri
# pri is pimary cases  
pri_death<-prob_death
rate_pri_sec<--log(1-pri_sec_p)/dur_pri_sec
pri_sec<-(1-exp(-rate_pri_sec))*(1-pri_death)
pri_pri<-1-pri_death-pri_sec
# sec is secondary cases
sec_death<-prob_death
rate_sec_ter<--log(1-sec_ter_p)/dur_sec_ter
sec_ter<-(1-exp(-rate_sec_ter))*(1-sec_death)
sec_lat<-(1-sec_death-sec_ter)*lat_p
sec_sec<-1-sec_death-sec_ter-sec_lat
# lat is latent cases
lat_death<-prob_death
lat_sec<-(1-sec_death-sec_ter)*(1-lat_p)
lat_lat<-1-lat_death-lat_sec
# ter is tertiary cases
ter_death<-prob_death
ter_ter<-1-ter_death

#######################run Markov model of effects (no eradication)############################

# run model
for (i in 2:cycle_num) {
  new[,i,,]<-new[,i-1,,]*new_new[,i-1,,]+repro_num[,i-1,,]*(1-exit_prob[,i-1,,])^max(i-1-burn,0)*pri[,i-1,,]
  pri[,i,,]<-pri[,i-1,,]*pri_pri[,i-1,,]+new_pri[,i-1,,]*new[,i-1,,]
  sec[,i,,]<-sec[,i-1,,]*sec_sec[,i-1,,]+pri_sec[,i-1,,]*pri[,i-1,,]+lat_sec[,i-1,,]*lat[,i-1,,]
  lat[,i,,]<-lat[,i-1,,]*lat_lat[,i-1,,]+sec_lat[,i-1,,]*sec[,i-1,,]
  ter[,i,,]<-ter[,i-1,,]*ter_ter[,i-1,,]+sec_ter[,i-1,,]*sec[,i-1,,]
}

# we collapse for world total (by cycle, mainting uncertainty)
# new
new_collapse<-apply(new,c(2,3,4),sum)
new_collapse_best<-apply(new_collapse,1,function(x) quantile(x, probs=c(0.50),na.rm=TRUE))
plot(new_collapse_best)
# primary
pri_collapse<-apply(pri,c(2,3,4),sum)
pri_collapse_best<-apply(pri_collapse,1,function(x) quantile(x, probs=c(0.5),na.rm=TRUE))
plot(pri_collapse_best)
# secondary
sec_collapse<-apply(sec,c(2,3,4),sum)
sec_collapse_best<-apply(sec_collapse,1,function(x) quantile(x, probs=c(0.5),na.rm=TRUE))
plot(sec_collapse_best)
# latent
lat_collapse<-apply(lat,c(2,3,4),sum)
lat_collapse_best<-apply(lat_collapse,1,function(x) quantile(x, probs=c(0.5),na.rm=TRUE))
plot(lat_collapse_best)
# tertiary
ter_collapse<-apply(ter,c(2,3,4),sum)
ter_collapse_best<-apply(ter_collapse,1,function(x) quantile(x, probs=c(0.5),na.rm=TRUE))
plot(ter_collapse_best)
# early=primary and secondary
early<-pri+sec
early_collapse<-pri_collapse+sec_collapse
early_collapse_best<-apply(early_collapse,1,function(x) quantile(x, probs=c(0.5),na.rm=TRUE))
plot(early_collapse_best)

#######################run Markov model of effects (eradication)#########################

# programmatic parameters
# these we allow to vary by country and cycle, with no correlation (see "tool box" below if we wish to impose correlation)

# coverage
cover<-array4dz
# in first (TCT) round
for (i in 1:ctry_num) {cover[i,(burn+1):cycle_num,,]<-runif(ndunc,min=0.90, max=0.99)}
# therafter in (TTT) round(s)
cover[,(burn+1+1):cycle_num,,]<-1

# eligibility 
eligible<-array4dz
# in first (TCT) round
for (i in 1:ctry_num) {eligible[i,(burn+1):cycle_num,,]<-runif(ndunc,min=0.98, max=0.99)}
# therafter in (TTT) round(s)
eligible[,(burn+1+1):cycle_num,,]<-1

# cure
cure<-array4d
for (i in 1:ctry_num) {cure[i,,,]<-runif(ndunc,min=0.96, max=1.00)}

# use same starting numbers as in baseline 
new_erad<-new
pri_erad<-pri
sec_erad<-sec
lat_erad<-lat
ter_erad<-ter

# run the model starting from the TCT (total community treatment) round in the first period (after burn-in)
for (i in (burn+1):cycle_num) {
  new_erad[,i,,]<-new_erad[,i-1,,]*new_new[,i-1,,]*(1-cover[,i-1,,]*eligible[,i-1,,]*cure[,i-1,,])+repro_num[,i-1,,]*(1-exit_prob[,i-1,,])^max(i-1-burn,0)*pri_erad[,i-1,,]*(1-cover[,i-1,,]*eligible[,i-1,,]*cure[,i-1,,])
  pri_erad[,i,,]<-pri_erad[,i-1,,]*pri_pri[,i-1,,]*(1-cover[,i-1,,]*eligible[,i-1,,]*cure[,i-1,,])+new_pri[,i-1,,]*new_erad[,i-1,,]*(1-cover[,i-1,,]*eligible[,i-1,,]*cure[,i-1,,])
  sec_erad[,i,,]<-sec_erad[,i-1,,]*sec_sec[,i-1,,]*(1-cover[,i-1,,]*eligible[,i-1,,]*cure[,i-1,,])+pri_sec[,i-1,,]*pri_erad[,i-1,,]*(1-cover[,i-1,,]*eligible[,i-1,,]*cure[,i-1,,])+lat_sec[,i-1,,]*lat_erad[,i-1,,]*(1-cover[,i-1,,]*eligible[,i-1,,]*cure[,i-1,,])
  lat_erad[,i,,]<-lat_erad[,i-1,,]*lat_lat[,i-1,,]*(1-cover[,i-1,,]*eligible[,i-1,,]*cure[,i-1,,])+sec_lat[,i-1,,]*sec_erad[,i-1,,]*(1-cover[,i-1,,]*eligible[,i-1,,]*cure[,i-1,,])
  ter_erad[,i,,]<-ter_erad[,i-1,,]*ter_ter[,i-1,,]+sec_ter[,i-1,,]*sec_erad[,i-1,,]
  }

# collapse
# primary
pri_erad_collapse<-apply(pri_erad,c(2,3,4),sum)
pri_erad_collapse_best<-apply(pri_erad_collapse,1,function(x) quantile(x, probs=c(0.5),na.rm=TRUE))
plot(pri_erad_collapse_best)
# secondary
sec_erad_collapse<-apply(sec_erad,c(2,3,4),sum)
sec_erad_collapse_best<-apply(sec_erad_collapse,1,function(x) quantile(x, probs=c(0.5),na.rm=TRUE))
plot(sec_erad_collapse_best)
# latent
lat_erad_collapse<-apply(lat_erad,c(2,3,4),sum)
lat_erad_collapse_best<-apply(lat_erad_collapse,1,function(x) quantile(x, probs=c(0.5),na.rm=TRUE))
plot(lat_erad_collapse_best)
# tertiary
ter_erad_collapse<-apply(ter_erad,c(2,3,4),sum)
ter_erad_collapse_best<-apply(ter_erad_collapse,1,function(x) quantile(x, probs=c(0.5),na.rm=TRUE))
plot(ter_erad_collapse_best)
# early
early_erad<-pri_erad+sec_erad
early_erad_collapse<-pri_erad_collapse+sec_erad_collapse
early_erad_collapse_best<-apply(early_erad_collapse,1,function(x) quantile(x, probs=c(0.5),na.rm=TRUE))
plot(early_erad_collapse_best)

####################calculate incremental effects########################

# DALYs
# based on mean and 95th CI from IHME GBD
dw_early<-get.beta.par(p=c(0.025,0.5,0.975),q=c(0.016,0.029,0.048))
dw_late<-get.beta.par(p=c(0.025,0.5,0.975),q=c(0.271,0.398,0.543))
temp<-rbeta(ndunc,shape1=dw_early[1], shape2=dw_early[2])
early_daly<-array4dz
for (i in 1:ndunc) {early_daly[,,i,]<-temp[i]}
temp<-rbeta(ndunc,shape1=dw_late[1], shape2=dw_late[2])
late_daly<-array4dz
for (i in 1:ndunc) {late_daly[,,i,]<-temp[i]}

# difference in effects between eradication and no eradication
pri_diff<-pri_collapse-pri_erad_collapse
sec_diff<-sec_collapse-sec_erad_collapse
ter_diff<-ter_collapse-ter_erad_collapse
early_diff<-early_collapse-early_erad_collapse
# plot mean differences for visual check
pri_diff_best<-apply(pri_diff,1,mean)
sec_diff_best<-apply(sec_diff,1,mean)
ter_diff_best<-apply(ter_diff,1,mean)
early_diff_best<-apply(early_diff,1,mean)
plot(pri_diff_best[1:cycle_num]/1000, x=1:cycle_num, xlab="", ylab="thousands")
plot(sec_diff_best[1:cycle_num]/1000, x=1:cycle_num, xlab="", ylab="thousands")
plot(ter_diff_best[1:cycle_num]/1000, x=1:cycle_num, xlab="", ylab="thousands")
plot(early_diff_best[1:cycle_num]/1000, x=1:cycle_num, xlab="", ylab="thousands")
# highs and low
early_diff_lo<-apply(early_diff,1,function(x) quantile(x, probs=c(0.05),na.rm=TRUE))
early_diff_hi<-apply(early_diff,1,function(x) quantile(x, probs=c(0.95),na.rm=TRUE))
ter_diff_lo<-apply(ter_diff,1,function(x) quantile(x, probs=c(0.05),na.rm=TRUE))
ter_diff_hi<-apply(ter_diff,1,function(x) quantile(x, probs=c(0.95),na.rm=TRUE))

# for presentation of disease dynamics and impact, collapse into years not cylces
early_diff_best_yearly<-rep(NA, years)
early_diff_lo_yearly<-rep(NA, years)
early_diff_hi_yearly<-rep(NA, years)
ter_diff_best_yearly<-rep(NA, years)
ter_diff_lo_yearly<-rep(NA, years)
ter_diff_hi_yearly<-rep(NA, years)
num=1
for (i in seq(from=burn+1, to=cycle_num, by=2)) {
  early_diff_best_yearly[num]<-(early_diff_best[i]+early_diff_best[i+1])/2
  early_diff_lo_yearly[num]<-(early_diff_lo[i]+early_diff_lo[i+1])/2
  early_diff_hi_yearly[num]<-(early_diff_hi[i]+early_diff_hi[i+1])/2
  ter_diff_best_yearly[num]<-(ter_diff_best[i]+ter_diff_best[i+1])/2
  ter_diff_lo_yearly[num]<-(ter_diff_lo[i]+ter_diff_lo[i+1])/2
  ter_diff_hi_yearly[num]<-(ter_diff_hi[i]+ter_diff_hi[i+1])/2
  num=num+1
}

plot(early_diff_best_yearly[1:years]/1000, x=2015:2050, xlab="", ylab="thousands")
p <- ggplot(as.data.frame(cbind(year=c(2015:2050),early_diff_best_yearly,early_diff_lo_yearly,early_diff_hi_yearly)), aes(year,early_diff_best_yearly/1000, ymin=early_diff_lo_yearly/1000,ymax=early_diff_hi_yearly/1000)) 
p + geom_pointrange() + scale_y_continuous(limits=c(), "thousands")
p + geom_pointrange() + expand_limits(y = 0) + scale_y_continuous(limits=c(), "thousands", expand = c(0, 0)) + theme_tufte() + geom_rangeframe()

plot(ter_diff_best_yearly[1:years]/1000, x=2015:2050, xlab="", ylab="thousands")
p <- ggplot(as.data.frame(cbind(year=c(2015:2050),ter_diff_best_yearly,ter_diff_lo_yearly,ter_diff_hi_yearly)), aes(year,ter_diff_best_yearly/1000, ymin=ter_diff_lo_yearly/1000,ymax=ter_diff_hi_yearly/1000)) 
p + geom_pointrange() + scale_y_continuous(limits=c(), "thousands") + scale_x_continuous(limits=c(), "year")
p + geom_pointrange() +  expand_limits(y = 0) + scale_y_continuous(limits=c(), "thousands", expand = c(0, 0)) + theme_tufte() + geom_rangeframe()

# discounted effects (no eradication)
early_disc<-early
ter_disc<-ter
early_a_late_disc<-early+ter
for (i in 2:cycle_num) {
	early_disc[,i,,]<-early_disc[,i,,]/(1+disc_e[,i,,])^max(i-1-burn,0)
	ter_disc[,i,,]<-ter_disc[,i,,]/(1+disc_e[,i,,])^max(i-1-burn,0)
	early_a_late_disc[,i,,]<-early_a_late_disc[,i,,]/(1+disc_e[,i,,])^max(i-1-burn,0)
	}
early_disc_collapse<-apply(early_disc, c(1,3,4), sum)
ter_disc_collapse<-apply(ter_disc, c(1,3,4), sum)
early_a_late_disc_collapse<-apply(early_a_late_disc, c(1,3,4), sum)

# summary table of (discounted period) effects
# early
early_disc_collapse_t<-abind(early_disc_collapse,apply(early_disc_collapse, c(2,3), sum), along=1)
early_disc_collapse_summary<-as.data.frame(apply(early_disc_collapse_t, 1, mean))
colnames(early_disc_collapse_summary)<-"best"
early_disc_collapse_summary$lo<-apply(early_disc_collapse_t, 1, function(x) quantile(x, probs=c(0.05),na.rm=TRUE))
early_disc_collapse_summary$hi<-apply(early_disc_collapse_t, 1, function(x) quantile(x, probs=c(0.95),na.rm=TRUE))
early_disc_collapse_summary
early_disc_collapse_summary$country<-row.names(early_disc_collapse_summary)
early_disc_collapse_summary<-early_disc_collapse_summary[,c(4,1,2,3)]
early_disc_collapse_summary
# late
ter_disc_collapse_t<-abind(ter_disc_collapse,apply(ter_disc_collapse, c(2,3), sum), along=1)
ter_disc_collapse_summary<-as.data.frame(apply(ter_disc_collapse_t, 1, mean))
colnames(ter_disc_collapse_summary)<-"best"
ter_disc_collapse_summary$lo<-apply(ter_disc_collapse_t, 1, function(x) quantile(x, probs=c(0.05),na.rm=TRUE))
ter_disc_collapse_summary$hi<-apply(ter_disc_collapse_t, 1, function(x) quantile(x, probs=c(0.95),na.rm=TRUE))
ter_disc_collapse_summary
ter_disc_collapse_summary$country<-row.names(ter_disc_collapse_summary)
ter_disc_collapse_summary<-ter_disc_collapse_summary[,c(4,1,2,3)]
ter_disc_collapse_summary
# early and late
early_a_late_disc_collapse_t<-abind(early_a_late_disc_collapse,apply(early_a_late_disc_collapse, c(2,3), sum), along=1)
early_a_late_disc_collapse_summary<-as.data.frame(apply(early_a_late_disc_collapse_t, 1, mean))
colnames(early_a_late_disc_collapse_summary)<-"best"
early_a_late_disc_collapse_summary$lo<-apply(early_a_late_disc_collapse_t, 1, function(x) quantile(x, probs=c(0.05),na.rm=TRUE))
early_a_late_disc_collapse_summary$hi<-apply(early_a_late_disc_collapse_t, 1, function(x) quantile(x, probs=c(0.95),na.rm=TRUE))
early_a_late_disc_collapse_summary
early_a_late_disc_collapse_summary$country<-row.names(early_a_late_disc_collapse_summary)
early_a_late_disc_collapse_summary<-early_a_late_disc_collapse_summary[,c(4,1,2,3)]
early_a_late_disc_collapse_summary

# using DALY weights
early_disc_daly<-early_disc*early_daly
ter_disc_daly<-ter_disc*late_daly
early_a_late_disc_daly<-early_disc_daly+ter_disc_daly
early_disc_daly_collapse<-apply(early_disc_daly, c(1,3,4), sum)
ter_disc_daly_collapse<-apply(ter_disc_daly, c(1,3,4), sum)
early_a_late_disc_daly_collapse<-apply(early_a_late_disc_daly, c(1,3,4), sum)
early_a_late_disc_daly_collapse_t<-abind(early_a_late_disc_daly_collapse,apply(early_a_late_disc_daly_collapse, c(2,3), sum), along=1)
apply(early_a_late_disc_daly_collapse_t, 1, mean)
apply(early_a_late_disc_daly_collapse_t, 1, function(x) quantile(x, probs=c(0.05),na.rm=TRUE))
apply(early_a_late_disc_daly_collapse_t, 1, function(x) quantile(x, probs=c(0.95),na.rm=TRUE))

# discounted effects (eradication)
early_erad_disc<-early_erad
ter_erad_disc<-ter_erad
early_a_late_erad_disc<-early_erad+ter_erad
for (i in 2:cycle_num) {
	early_erad_disc[,i,,]<-early_erad_disc[,i,,]/(1+disc_e[,i,,])^max(i-1-burn,0)
	ter_erad_disc[,i,,]<-ter_erad_disc[,i,,]/(1+disc_e[,i,,])^max(i-1-burn,0)
	early_a_late_erad_disc[,i,,]<-early_a_late_erad_disc[,i,,]/(1+disc_e[,i,,])^max(i-1-burn,0)
	}
early_erad_disc_collapse<-apply(early_erad_disc, c(1,3,4), sum)
ter_erad_disc_collapse<-apply(early_erad_disc, c(1,3,4), sum)
early_a_late_erad_disc_collapse<-apply(early_a_late_erad_disc, c(1,3,4), sum)
early_a_late_erad_disc_collapse_t<-abind(early_a_late_erad_disc_collapse,apply(early_a_late_erad_disc_collapse, c(2,3), sum), along=1)
apply(early_a_late_erad_disc_collapse_t, 1, mean)
apply(early_a_late_erad_disc_collapse_t, 1, function(x) quantile(x, probs=c(0.05),na.rm=TRUE))
apply(early_a_late_erad_disc_collapse_t, 1, function(x) quantile(x, probs=c(0.95),na.rm=TRUE))
# using DALY weights
early_erad_disc_daly<-early_erad_disc*early_daly
ter_erad_disc_daly<-ter_erad_disc*late_daly
early_a_late_erad_disc_daly<-early_erad_disc_daly+ter_erad_disc_daly
early_erad_disc_daly_collapse<-apply(early_erad_disc_daly, c(1,3,4), sum)
ter_erad_disc_daly_collapse<-apply(ter_erad_disc_daly, c(1,3,4), sum)
early_a_late_erad_disc_daly_collapse<-apply(early_a_late_erad_disc_daly, c(1,3,4), sum)
early_a_late_erad_disc_daly_collapse_t<-abind(early_a_late_erad_disc_daly_collapse,apply(early_a_late_erad_disc_daly_collapse, c(2,3), sum), along=1)
apply(early_a_late_erad_disc_daly_collapse_t, 1, mean)
apply(early_a_late_erad_disc_daly_collapse_t, 1, function(x) quantile(x, probs=c(0.05),na.rm=TRUE))
apply(early_a_late_erad_disc_daly_collapse_t, 1, function(x) quantile(x, probs=c(0.95),na.rm=TRUE))

#########################################calculate costs#################################

# target percentage per year
tar_p<-array4dz
# for simplicity we assumed all population targeted in first period
tar_p[,burn+1,,]<-1
tar_0_4<-par_0_4*tar_p
tar_5_9<-par_5_9*tar_p
tar_10_14<-par_10_14*tar_p
tar_15_99<-par_15_99*tar_p

# contacts
# in (TTT) round(s) only
contacts<-array4dz
for (i in 1:ctry_num) {contacts[i,(burn+1+1):cycle_num,,]<-runif(ndunc,min=8, max=12)}

# actual number targeted for treatment will include "mop up" of remaining cases and their contacts in TTT only
mop<-array4dz
mop[,(burn+1+1):cycle_num,,]<-(pri_erad[,(burn+1+1):cycle_num,,]+sec_erad[,(burn+1+1):cycle_num,,])*(1+contacts[,(burn+1+1):cycle_num,,])
# we assume contacts average contact is aged 10-14
txd_0_4<-tar_0_4
txd_5_9<-tar_5_9
txd_10_14<-tar_10_14+mop
txd_15_99<-tar_15_99

# drug costs (including 10% buffer)
azi_mgs<-array4dz
azi_mgs[,,,]<-(txd_0_4[,,,]*500+txd_5_9[,,,]*1000+txd_10_14[,,,]*1500+txd_15_99[,,,]*2000)*eligible[,,,]*cover[,,,]
uc_azi<-0.17
azi_usd<-uc_azi*azi_mgs/500*1.1

# clinical misdiagnosis
misdx<-array4dz
for (i in 1:ctry_num) {misdx[i,,,]<-runif(ndunc,min=0.1, max=0.3)}
# diagnostic costs
# population to be tested serogically 
sero_num<-array4dz
sero_num[,,,]<-(pri_erad+sec_erad)*(1+misdx)
uc_sero<-2
sero_usd<-uc_sero*sero_num

# implementation
# based on benchmark values
uc_mda_bmk<-read.csv(file="C:/Users/fitzpatrickc/Dropbox/yaws/mda_uc_yaw.csv")
uc_mda_bmk<-uc_mda_bmk[with(uc_mda_bmk, order(country)),]
# sort by country, as in other databases
uc_mda<-array4dz
for (i in 1:ctry_num) {
  temp<-exp(rnorm(ndunc,mean=uc_mda_bmk$ln_uc_econ_iu0_ac1_xb[i], sd=uc_mda_bmk$ln_uc_econ_iu0_ac1_se[i]))
  for (j in 1:ndunc) {uc_mda[i,,j,]<-temp[j]}
}

# TCT
implem_num<-array4dz
implem_num[,burn+1,,]<-txd_0_4[,burn+1,,]+txd_5_9[,burn+1,,]+txd_10_14[,burn+1,,]+txd_15_99[,burn+1,,]
# TTT, if required, 30-50% of cost of TCT based on Lihir island
ttt_pc<-array4dz
for (i in 1:ctry_num) {ttt_pc[i,,,]<-runif(ndunc,min=0.3, max=0.5)}
implem_num[,burn+1+1,,]<-implem_num[,burn+1,,]*ttt_pc[,burn+1+1,,]
implem_usd<-uc_mda*implem_num

# surveillance in the subsequent period, for unknown duration
surv_years<-array4dz
for (i in 1:ctry_num) {surv_years[i,(burn+1+1+1),,]<-runif(ndunc,min=1, max=3)}
uc_surv<-array4dz
for (i in 1:ctry_num) {uc_surv[i,,,]<-runif(ndunc,min=2000/100000, max=30000/100000)}
surv_num<-array4dz
surv_num[,,,]<-(par_0_4[,,,]+par_5_9[,,,]+par_10_14[,,,]+par_15_99[,,,])*surv_years[,,,]
surv_usd<-uc_surv*surv_num

# quantities by line item
mean(apply(azi_mgs,c(3,4),sum))
mean(apply(azi_mgs,c(3,4),sum))/375000000000
quantile(apply(azi_mgs,c(3,4),sum),c(0.5,0.05,0.95))
mean(apply(sero_num,c(3,4),sum))
quantile(apply(sero_num,c(3,4),sum),c(0.5,0.05,0.95))

# cost by line item
mean(apply(azi_usd,c(3,4),sum))
quantile(apply(azi_usd,c(3,4),sum),c(0.5,0.05,0.95))
mean(apply(sero_usd,c(3,4),sum))
quantile(apply(sero_usd,c(3,4),sum),c(0.5,0.05,0.95))
mean(apply(implem_usd,c(3,4),sum))
quantile(apply(implem_usd,c(3,4),sum),c(0.5,0.05,0.95))
mean(apply(surv_usd,c(3,4),sum))
quantile(apply(surv_usd,c(3,4),sum),c(0.5,0.05,0.95))

# total cost (undiscounted) 
total_usd<-azi_usd+sero_usd+implem_usd+surv_usd
mean(apply(total_usd,c(3,4),sum))
quantile(apply(total_usd,c(3,4),sum),c(0.5,0.05,0.95))
mean(apply(total_usd,c(3,4),sum))/sum(gecon_subsub$MER2005_40)
mean(apply(total_usd,c(3,4),sum))/sum(gecon_sub$MER2005_40)
total_ntx_usd<-sero_usd+implem_usd+surv_usd
mean(apply(total_ntx_usd,c(3,4),sum))
quantile(apply(total_ntx_usd,c(3,4),sum),c(0.5,0.05,0.95))

# as % of GDP
mean(apply(total_usd,c(3,4),sum)/sum(gecon_subsub$MER2005_40))
mean(apply(total_usd,c(3,4),sum)/sum(gecon_sub$MER2005_40))
  
# financial cost excluding drugs and diagnostic tests
fin_pc<-array4dz
for (i in 1:ctry_num) {fin_pc[i,,,]<-runif(ndunc,min=0.28, max=0.86)}
total_fin_usd<-(sero_usd+implem_usd+surv_usd)*fin_pc
mean(apply(total_fin_usd,c(3,4),sum))
quantile(apply(total_fin_usd,c(3,4),sum),c(0.5,0.05,0.95))

# surveillance in the 76 countries of unknown endemicity
surv_usd_unknown<-runif(ndunc,min=2000/100000, max=30000/100000)*par_unknown
mean(surv_usd_unknown)
quantile(surv_usd_unknown,c(0.5,0.05,0.95))

####################### discount total costs for the period#################

#discount rate on costs (country but not year specific)
temp<-array(1-exp(-runif(ndunc,min=0.03, max=0.06)*cycle_len),dim=c(ctry_num, ndunc))
disc_c<-array4d
for (i in 1:ctry_num) {
  for (j in 1:ndunc){
    disc_c[i,,j,]<-temp[i,j]
  }
}

# discounted costs (eradication) 
total_usd_disc<-total_usd
total_ntx_usd_disc<-total_ntx_usd
for (i in (burn+1):cycle_num) {
  total_usd_disc[,i,,]<-total_usd_disc[,i,,]/(1+disc_c[,i,,])^max(i-1-burn,0)
  total_ntx_usd_disc[,i,,]<-total_ntx_usd_disc[,i,,]/(1+disc_c[,i,,])^max(i-1-burn,0)
}

# summary table of costs (discounted) 
total_usd_disc_collapse<-apply(total_usd_disc, c(1,3,4), sum)
total_usd_disc_t<-abind(total_usd_disc_collapse,apply(total_usd_disc_collapse, c(2,3), sum), along=1)
total_usd_disc_summary<-as.data.frame(apply(total_usd_disc_t, 1, mean))
colnames(total_usd_disc_summary)<-"best"
total_usd_disc_summary$lo<-apply(total_usd_disc_t, 1, function(x) quantile(x, probs=c(0.05),na.rm=TRUE))
total_usd_disc_summary$hi<-apply(total_usd_disc_t, 1, function(x) quantile(x, probs=c(0.95),na.rm=TRUE))
total_usd_disc_summary$country<-row.names(total_usd_disc_summary)
total_usd_disc_summary<-total_usd_disc_summary[,c(4,1,2,3)]
total_usd_disc_summary

# summary table of costs excluding drugs (discounted) 
total_ntx_usd_disc_collapse<-apply(total_ntx_usd_disc, c(1,3,4), sum)
total_ntx_usd_disc_t<-abind(total_ntx_usd_disc_collapse,apply(total_ntx_usd_disc_collapse, c(2,3), sum), along=1)
total_ntx_usd_disc_summary<-as.data.frame(apply(total_ntx_usd_disc_t, 1, mean))
colnames(total_ntx_usd_disc_summary)<-"best"
total_ntx_usd_disc_summary$lo<-apply(total_ntx_usd_disc_t, 1, function(x) quantile(x, probs=c(0.05),na.rm=TRUE))
total_ntx_usd_disc_summary$hi<-apply(total_ntx_usd_disc_t, 1, function(x) quantile(x, probs=c(0.95),na.rm=TRUE))
total_ntx_usd_disc_summary$country<-row.names(total_ntx_usd_disc_summary)
total_ntx_usd_disc_summary<-total_ntx_usd_disc_summary[,c(4,1,2,3)]
total_ntx_usd_disc_summary

########################estimate cost-effectiveness###########################

# incremental cost effectiveness
incr_cost<-apply(total_usd_disc_collapse, c(2,3), sum)
incr_effect<-apply(early_a_late_disc_collapse,c(2,3),sum)-apply(early_a_late_erad_disc_collapse, c(2,3), sum)
icer<-incr_cost/incr_effect
mean(icer)
quantile(icer,c(0.05,0.95))

# using disability weights 
incr_effect_daly<-apply(early_a_late_disc_daly_collapse,c(2,3),sum)-apply(early_a_late_erad_disc_daly_collapse,c(2,3),sum)
icer_daly<-incr_cost/incr_effect_daly
mean(icer_daly)
quantile(icer_daly,c(0.05,0.95))

######################## produce summary outputs ###########################

# cost-effectiveness acceptability curve (CEAC)
quantile(icer,c(0.50,0.90))
plot(ecdf(icer), main="", xlab="acceptability threshold (US$)", ylab="probability that cost per life-year < threshold", xlim=c(0,100), do.points = FALSE)
df <- as.data.frame(icer)
names(df)<-c("x")
ggplot(df, aes(x=x)) + stat_ecdf(size=1)  + theme_tufte() + geom_rangeframe(stat="ecdf", sides="bt") + coord_cartesian(xlim = c(0, 50)) +
  geom_segment(aes(x=0,y=0.5,xend=quantile(icer,c(0.50)),yend=0.50), colour="red") + 
  geom_segment(aes(x=0,y=0.9,xend=quantile(icer,c(0.90)),yend=0.90), colour="red") + 
  geom_segment(aes(x=quantile(icer,c(0.50)),y=0,xend=quantile(icer,c(0.50)),yend=0.50), colour="red") + 
  geom_segment(aes(x=quantile(icer,c(0.90)),y=0,xend=quantile(icer,c(0.90)),yend=0.90), colour="red") + 
  scale_x_continuous("acceptability threshold (US$)") + 
  scale_y_continuous(breaks=c(seq(0,1,by=0.1)), "probability that cost per life-year < threshold", expand = c(0, 0)) 

#### optional (not verified)

# cost-effectiveness plane
scatter_y<-as.vector(incr_cost)
scatter_x<-as.vector(incr_effect)
scatter<-as.data.frame(cbind(scatter_y,scatter_x))
scatter<-transform(scatter,slope=scatter_y/scatter_x)
slope_quantiles<-quantile(scatter$slope,c(0.025,0.975))
slope_mean<-mean(scatter$slope)
p <- ggplot(scatter, aes(scatter_x/1000, scatter_y/1000)) 
# add titles, add threshhold based on GDP per capita (replace slope), make points slightly transparent so we can see where they overlap
r <- p + geom_point(alpha = 0.2) +
  geom_abline(intercept=0, slope=slope_mean, colour = "blue")  +
  geom_abline(intercept=0, slope=slope_quantiles[1], colour = "blue", linetype="dotdash")  +
  geom_abline(intercept=0, slope=slope_quantiles[2], colour = "blue",  linetype="dotdash")  +
  scale_x_continuous(limits=c(0,25000), "Incremental effects (thousands of life-years of yaws averted)") +
  scale_y_continuous(limits=c(0,quantile(scatter_y/1000,c(0.99))), "Incremental costs (thousands US$)") + theme_tufte() 
r

# expected incremental benefit
wtp<-20
eib<-(wtp*incr_effect-incr_cost)
mean(eib)
quantile(eib,c(0.025,0.975))
hist(eib/1000000, freq=FALSE)

# expected value of perfect information (EVPI)
# opportunity loss
ol<--nmb
ol[ol<0]<-0
mean(ol)
# so what if we just look at the EVPPI (perfect partial information), focussing just on the uncertainty dimension 
nmb_mean<-rowMeans(nmb)
ol_part<--nmb_mean
ol_part[ol_part<0]<-0
mean(ol_part)
# or the variability dimension
nmb_mean<-colMeans(nmb)
ol_part<--nmb_mean
ol_part[ol_part<0]<-0
mean(ol_part)

# summary plots from the BCEA package (uncertainty and variability combined)
library(BCEA, lib="D:/Rlibrary")
# we first turn the 2D matrices into the BCEA vector format
fx<-cbind(as.vector(-apply(early_a_late_disc_collapse,c(2,3),sum)), as.vector(-apply(early_a_late_erad_disc_collapse,c(2,3),sum)))
# note here that the effects are negative (not QALYs)
costs<-cbind(0,as.vector(apply(total_usd_disc_collapse, c(2,3), sum)))
ints<-c("Baseline", "Eradication")
bcea<-bcea(e=fx,c=costs,ref=2,interventions=ints, Kmax=wtp)
# if matrix is too large, select 1000 rows randomly
# sample<-sample(nrow(fx),1000)
# bcea<-bcea(e=fx[sample,],c=costs[sample,],ref=1,interventions=ints, Kmax=30000)
print(bcea$ICER)
ceplane.plot(bcea,comparison=1,wtp=wtp)
eib.plot(bcea)
ceac.plot(bcea)
evi.plot(bcea)

############# TOOL BOX ##################

# impose correlation structure in PSA
# use cornode (Iman-Conover method) from mc2d to impose correlation structures
mat <- cbind(disc_c,disc_e)
cor(mat)
corr1 <- matrix(c(1, 0.99, 0.99, 1), ncol=2)
matc <- cornode(mat, target=corr1)
cor(matc)
# or use package copula
install.packages("copula", lib="D:/Rlibrary", dependencies=TRUE)
library(copula, lib="D:/Rlibrary")
rho<-0.5 # Spearman's rho
# For normal and T elliptical copulas, rho=sin(pi/2*tau) 
tau<-asin(rho)*2/pi # Kendall's tau
temp <- mvdc(normalCopula(rho), c("gamma", "norm"), list(list(shape = 1, rate = 2), list(mean = 0, sd =2))) 
dat <- rMvdc(1000, temp) 
cor(dat)
# The Frank copula is a symmetric Archimedean copula. The relationship between Kendall's tau rank correlation coefficient and the Frank copula parameter is given by (1-tau)/4
# The Clayton copula is an asymmetric Archimedean copula, exhibiting greater dependence in the negative tail than in the positive. The relationship between Kendall's tau and the Clayton copula parameter is given by 2*tau/(1-tau) 
# The Gumbel copula (a.k.a. Gumbel-Hougard copula) is an asymmetric Archimedean copula, exhibiting greater dependence in the positive tail than in the negative. The relationship is given by 1/(1-tau)
# http://www.vosesoftware.com/ModelRiskHelp/index.htm#Modeling_correlation/Copulas.htm

# import data from the Global Health Observatory
# you may need to run the following before downloading (depending on proxy settings)
setInternet2()
# WHO GHO
# http://apps.who.int/gho/data/view.main
whoData1<-read.csv(url("http://apps.who.int/gho/athena/data/data-verbose.csv?target=GHO/NTD_1&profile=verbose&filter=COUNTRY:*"))
whoData2<-read.csv(url("http://apps.who.int/gho/athena/data/data-verbose.csv?target=GHO/NTD_2&profile=verbose&filter=COUNTRY:*"))

# import data from WDI
# http://data.worldbank.org/data-catalog/world-development-indicators
# http://cran.r-project.org/web/packages/WDI/index.html
library(WDI, lib="D:/Rlibrary")
library(countrycode, lib="D:/Rlibrary")
# Use the WDIsearch function to get a list of indicators
indicatorMetaData <- WDIsearch("Population", field="name", short=FALSE)
# Define a list of countries for which to pull data
countries <- c("United States", "Britain", "Sweden", "Germany")
# Convert the country names to iso2c format used in the World Bank data
iso2cNames <- countrycode(countries, "country.name", "iso2c")
# Pull data for individual or all countries, single or multiple indicators
wdi_gdp2 <- WDI(country="all", indicator = "SP.POP.TOTL", start=2001, end=2050, extra=TRUE)
wdi_gdp3 <- WDI(country="all", indicator=c("NV.AGR.PCAP.KD.ZG","GDPPCKD"), start=2001, end=2011, extra=TRUE)
## Filter out the aggregates
subData <- subset(wdi_gdp3, !iso3c %in% NA) 

# get distribution parameters
# http://cran.r-project.org/web/packages/rriskDistributions/
test<-rnorm(10000,100,20)
hist(test)
# load GUI to test all distributions
fit.cont(test)
# get parameters for a specific distribution by maximum likelihood  
# including: "norm", "exp", "lnorm", "logis", "gamma", "weibull", "beta", "chisq", "t", "f", "cauchy", "gompertz"
a<-rriskFitdist.cont(test, "norm", method = c("mle"))
a$estimate                  
# get parameters for a specific distribution by method of moments
b<-rriskFitdist.cont(test, "norm", method = c("mme"))
b$estimate                  
a<-fit.perc(p=c(0.025,0.5,0.6,0.975),q=c(1,3.5,4,7))
a$fittedParams
# or if you think you know the distribution
# triangle, from three or more quantiles
b<-get.triang.par(p=c(0.025,0.5,0.975),q=c(1,3,5))
b
# pert, from four or more quantiles
c<-get.pert.par(p=c(0.025,0.5,0.6,0.975),q=c(1,3.5,4,7))
c
# beta, from two or more quantiles
d<-get.beta.par(p=c(0.025,0.975),q=c(0.2,0.8))
d
# gamma, from two or more quantiles
e<-get.gamma.par(p=c(0.025,0.975),q=c(1,5))
e
# log normal, from two or more
f<-get.lnorm.par(p=c(0.025,0.5,0.975),q=c(0.88,1.88,4.31))







