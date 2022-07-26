
# Code for: "Honeyguides and other wax-eaters: interspecific 
# competitors for beeswax stabilise rather than jeopardize a bird-human mutualism"

# Data stored in two files: "waxdata2022.csv" and "wax_td_format.csv"

library(tidyverse) # load packages
library(survival)
library(survminer)
library(hrbrthemes)
library(ggplot2)
library(ggsci)
library(lubridate)
library(forcats)
library(knitr)
library(dplyr)
library(plotrix)
library(lme4)
library (emmeans)
library(jtools)
library(suncalc)


working_directory <- "G:/Documents"
setwd(working_directory)
waxdata <- read.csv("waxdata2022.csv", sep=",", comment.char="#")

# ----- Clean and filter data ---- 
# remove two rows where visits to wax may have been affected by close placement of EFID readers

waxdata <- waxdata %>%
  filter(!(row_number() %in% 74:75)) 

waxdata <- waxdata %>%
  filter(!(trial == "11")) # to remove single crowned hornbill record for most analyses as it differs from other wax trials

# ----- Summary statistics for methods -----
# number of wax visits for each site (n=26 plus crowned hornbill observations)
# check number of species

waxdata %>%
  count(trial)

# number of different harvest sites (geographically unique locations) n = 26 sites
waxdata %>%
  count(lat,lon)

# number of trees we were guided by a honeyguide (guided: n= 22)
waxdata %>%
  count(trial, guided)

#count the number of non-zero days after harvest when wax was placed
data.frame(waxdata %>% 
             count(trial,time_from_harvesttoplaced,guided))

#average number of days + nights (24hr periods) camera-traps were set for
sum <- waxdata %>% # get the first feeding event
  group_by(trial) %>%
  filter(row_number()==1) 
mean(sum$trap_duration_days,na.rm = TRUE)
std.error(sum$trap_duration_days,na.rm = TRUE)
range(sum$trap_duration_days,na.rm = TRUE)
sum(sum$trap_duration_days,na.rm = TRUE)

# ---- Average number of species/site ---- 

sp_visits <- waxdata %>%
  filter (sp_com =='Greater honeyguide'|sp_com =='Scaly-throated honeyguide'|sp_com =='Striped bush squirrel'
          |sp_com =='Crowned hornbill'|sp_com =='African civet' |sp_com =='Honey badger'|sp_com =='Slender mongoose'
          |sp_com =='Yellow baboon'|sp_com =='Lesser honeyguide'|sp_com =='Mellers mongoose') %>%
  group_by(trial) %>%
  summarise(No_sp = n_distinct(sp_com))
sp_visits
mean <- mean(sp_visits$No_sp)
mean
median(sp_visits$No_sp)
sd(sp_visits$No_sp)
se <- sd(sp_visits$No_sp)/sqrt(length(sp_visits$No_sp))
se
range(sp_visits$No_sp)


# ----  Average feeding time for honeyguides --------

hg_mean <- waxdata %>%
  filter(sp_com =='Greater honeyguide') %>%
  filter (food_type == 'wax'|food_type == 'wax_brood') %>%
  filter(wax_visible == '1')
hg_mean
mean <-mean(hg_mean$duration)
mean/60 # to get result in minutes
sd(hg_mean$duration)
se <- sd(hg_mean$duration)/sqrt(length(hg_mean$duration))
length(hg_mean$duration)

squirrel_mean <- waxdata %>%
  filter(sp_com =='Striped bush squirrel') %>%
  filter (food_type == 'wax'|food_type == 'wax_brood') %>%
  filter(wax_visible == '1')
mean <-mean(squirrel_mean$duration)
mean/60 # to get result in minutes
sd(squirrel_mean$duration)
se <- sd(squirrel_mean$duration)/sqrt(length(squirrel_mean$duration))
mean
se

# number of squirrel visits per site
squirrel_count <- squirrel_mean %>% group_by(trial) %>% count()
mean <- mean(squirrel_count$n)
sd(squirrel_count$n)
se <- sd(squirrel_count$n)/sqrt(length(squirrel_count$n))
mean
se

# number of baboon visits per site
baboon_mean <- waxdata %>%
  filter(sp_com =='Yellow baboon') %>%
  filter (food_type == 'wax'|food_type == 'wax_brood') %>%
  filter(wax_visible == '1')
mean <-mean(baboon_mean$duration)
mean/60 # to get result in minutes
sd(baboon_mean$duration)
se <- sd(baboon_mean$duration)/sqrt(length(squirrel_mean$duration))
mean
se

# ---- Results section a ---- 
# ---- Hypothesis 1, prediction (i) ---- descriptive test: first arrivals at each site ----

waxdata_firsts <- waxdata %>%
  filter (food_type == 'wax'|food_type == 'wax_brood') 

waxdata_firsts %>%
  count(sp_com) # number of species which fed on wax

waxdata_first <- waxdata_firsts %>% # get the first feeding event
  group_by(trial) %>%
  filter(row_number()==1)
table_firsts <- waxdata_first %>% #table with a list of the first-arriving animal sp
  count(sp_com,wax_visible)
view(table_firsts)

#summarise final animal which depletes wax

waxdata_finalanimal <- waxdata_firsts %>% # get the first feeding event
  group_by(trial) %>%
  filter(row_number()==1) 
depletions <- waxdata_finalanimal %>% 
  count(final_animal)
table(depletions$final_animal)

4/7 # baboon
3/15 # mellers
2/6 # honeybadger
2/31 # civet
1/84 # minor 
2/173 # squirrel
6/633 # greater honeyguide
0 #scaly-throated honeyguide
0 #crowned hornbill 
0 #slender mongoose 

# number of times honeyguides feed during the first day (at a site that honeyguides arrived at): 

hg_visits <- waxdata %>%
 filter (food_type == 'wax'|food_type == 'wax_brood')

sixhrs <- hg_visits %>% # filter by any times that are prior to 21600 seconds (6 hrs)
  filter(timeafterplacement < 21600) %>%
  filter(sp_com == 'Greater honeyguide')
nrow(sixhrs)

twelvehrs <- hg_visits %>% # filter by any times that are prior to 43200 seconds (12 hrs)
  filter(timeafterplacement < 43200) %>%
  filter(sp_com == 'Greater honeyguide')
nrow(twelvehrs)

twentyfourhrs <- hg_visits %>% # filter by any times that are prior to 86400 seconds (24 hrs)
  filter(timeafterplacement < 86400) %>%
  filter(sp_com == 'Greater honeyguide')
nrow(twentyfourhrs)
356/633

# ---- Hypothesis 1, prediction (i) -- Cox PH model - do species arrive at significantly different times? ----

feedtimes <- waxdata %>% 
  filter(food_type == 'wax'|food_type == 'wax_brood') 

feedtimes$grp[feedtimes$sp_com == "Greater honeyguide"] <- "B" #group honey major competitors 
feedtimes$grp[feedtimes$sp_com == "Honey badger"] <- "A"
feedtimes$grp[feedtimes$sp_com == "Yellow baboon"] <- "A"
feedtimes$grp[feedtimes$sp_com == "Mellers mongoose"] <- "A"
feedtimes$grp[feedtimes$sp_com == "African civet"] <- "A" 
feedtimes$grp[feedtimes$sp_com == "Lesser honeyguide"] <- "A"
feedtimes$grp[feedtimes$sp_com == "Scaly-throated honeyguide"] <- "A"
feedtimes$grp[feedtimes$sp_com == "Striped bush squirrel"] <- "A"
feedtimes$grp[feedtimes$sp_com == "Slender mongoose"] <- "A"


feedtimes <- feedtimes %>% group_by(trial, sp_com) %>% slice(which.min(timeafterplacement), na.rm = TRUE)
feedtimes <- droplevels(feedtimes)
str(feedtimes)

cox.mod <- coxph(Surv(timeafterplacement/3600, wax_visible) ~ grp, data = feedtimes)
summary(cox.mod)
test.ph <- cox.zph(cox.mod)
test.ph
ggcoxzph(test.ph)


# ---- Figure 3a ----
# Probability of feeding times (for all for species with > 6 records). 

feedtimes <- waxdata %>% 
  #filter(food_type == 'wax'|food_type == 'wax_brood') %>%
  filter(sp_com =='Greater honeyguide'|sp_com =='Scaly-throated honeyguide'|
           sp_com =='Striped bush squirrel'|sp_com =='Lesser honeyguide'|
         sp_com =='Honey badger'|sp_com =='Mellers mongoose'|
           sp_com =='Yellow baboon'|sp_com =='African civet') 


feedtimes <- feedtimes %>% group_by(trial, sp_com) %>% slice(which.min(timeafterplacement), na.rm = TRUE)
feedtimes <- droplevels(feedtimes)
#str(feedtimes)

survplot <- survfit(Surv(timeafterplacement/3600, wax_visible) ~ sp_com, data = feedtimes)
survplot1 <-ggsurvplot(survplot, ggtheme=theme_classic2(),fun = 'event', 
                       xlab = "Time since wax availablity (h)",
                       ylab = "Probability of arrival at wax",
                       y.lim =c(0,1),
                       legend=c("right"), 
                       censor.size = c(5),
                       censor.shape = c(124),
                       axis.offset = FALSE,
                       #conf.int = TRUE,
                       #conf.int.style = c("step"),
                       size = 0.7, # line thickness
                       palette=c("#b00101","#0a68cd","#464646","#91bfdb","#e27728","#0a8f77","#464646","#fca70a"),
                       linetype = c("solid","solid","solid","solid","solid","solid","longdash","solid"),
                       censor =TRUE)

survplot1
# save as pdf - dimensions 5 x 10 inches


# ---- Hypothesis 1, prediction (ii) -- Are most important competitors diurnal? ----
visits <- waxdata 
# add in column whether visit were AM/PM
visits$hr <- hour(visits$timestart)
visits$min <- minute(visits$timestart)
visits$sec <- second(visits$timestart)
visits$AMPM <- "PM"
visits$AMPM[visits$hr > 5 & visits$hr < 17] <- "AM"
visits <- visits %>% group_by(trial) %>% slice(which.min(final_animal), na.rm = TRUE)
visits %>% select(final_animal,AMPM)


# ---- Hypothesis 1, prediction (iii) -- Do diurnal competitors consistently displace honeyguides from wax? ----

allvisits <- waxdata %>% #get only important diurnal competitors 
  filter(food_type == 'wax'|food_type == 'wax_brood') %>% # only wax-eating events
    filter (sp_com =='Greater honeyguide'|sp_com =='Scaly-throated honeyguide'|sp_com =='Striped bush squirrel'|sp_com =='Lesser honeyguide') 
allvisits %>%
  count(sp_com)

visits <- allvisits %>% 
  group_by(trial) %>%
  arrange(timestart_UNIX) %>%
  mutate(diff = timestart_UNIX - lag(timefinish_UNIX, default = NA)) %>%
  mutate(diffa = lead(timestart_UNIX) - timefinish_UNIX, default = NA)

visits <- visits %>% 
  mutate(simul = if_else(sp_com == "Greater honeyguide" & diff < 0 |sp_com == "Greater honeyguide" & diffa < 0, 'YES','NO')) # simultaneous feeding to squirrels, scaly-throated honeyguides or lesser honeyguides
visits <- visits %>% 
  mutate(befaf = if_else(sp_com == "Greater honeyguide" & diff <= 10 |sp_com == "Greater honeyguide" & diffa <= 10, 'YES','NO')) # before or after other species

#write.csv(visits, "G:/Documents/visits.csv", row.names = TRUE) # checked for accuracy 

# percentage of visits by diurnal competitors at wax (scaly-throated honeyguides, striped bush squirrels, yellow baboon which are
# simultaneous with greater honeyguides (identified through visual inspection of above dataset)
50/633 # n = 11 for lesser honeyguides, n= 17 for squirrels, n = 22 scaly-throated hg

# percentage of feeding visits cut short by either species 
1/633

# percentage of greater honeyguide visits which were immediately before or after either of above 3 species
24/633

# ---- Hypothesis 2, prediction (i) -- greater honeyguides are the first-arriving species after the wax has been exposed ----

# percentage of sites where greater honeyguides arrived first
waxdata_firsts <- waxdata %>%
  filter (food_type == 'wax'|food_type == 'wax_brood') 

waxdata_firsts <- waxdata_firsts %>% # get the first feeding event
  group_by(trial) %>%
  filter(row_number()==1)
table_firsts <- waxdata_firsts %>% #table with a list of the first-arriving animal sp
  count(sp_com,wax_visible)
table(table_firsts$sp_com)

18/26

# percentage of sites where greater honeyguides arrived last

waxdata_lasts <- waxdata_firsts %>% # get the first feeding event
  group_by(trial) %>%
  filter(row_number()==1) %>% 
  count(final_animal)
waxdata_lasts


# Probability of depletion by major wax competitors 

mwc_visits <- waxdata %>% #get only competitors sp
  filter(food_type == 'wax'|food_type == 'wax_brood') %>% # only wax-eating events
  filter (sp_com =='Honey badger'|sp_com =='Mellers mongoose'|sp_com =='African civet'| sp_com =='Yellow baboon') 
mwc_visits %>%
  count(sp_com)

9/28 #overall
2/6 #honeybadger
3/15 #mellers
4/7 #baboon


# ---- Hypothesis 2, prediction (ii) -- majority of greater honeyguide feeding events occur before those of major wax competitors ----

visits <- waxdata %>% 
  filter((sp_com == 'Greater honeyguide')) #get only greater honeyguides

visits <- visits %>% # get visits which are only wax eating for 'before other species arrive' and after, 
  # all visits after wax is complete (arrival factor level = 3)
  filter(!(sp_com == 'Greater honeyguide'& arrival == 'A' & food_type == 'not_eating')) %>%
  filter(!(sp_com == 'Greater honeyguide'& arrival == 'A' & food_type == 'eating_fooduncertain')) %>%
  filter(!(sp_com == 'Greater honeyguide'& arrival == 'A' & food_type == 'brood')) %>%
  filter(!(sp_com == 'Greater honeyguide'& arrival == 'B' & food_type == 'not_eating')) %>%
  filter(!(sp_com == 'Greater honeyguide'& arrival == 'B' & food_type == 'eating_fooduncertain')) %>%
  filter(!(sp_com == 'Greater honeyguide'& arrival == 'B' & food_type == 'brood')) %>%
  filter(!(sp_com == 'Greater honeyguide'& arrival == 'C' & food_type == 'eating_fooduncertain')) %>%
  filter(!(sp_com == 'Greater honeyguide'& arrival == 'C' & food_type == 'brood'))

visits1 <- visits %>% 
  group_by(trial, arrival,offset_days, .drop = FALSE) %>% 
  summarise(count = n(), duration = mean(duration, na.rm = TRUE))
visits1<-subset(visits1,duration!="NA")

h2p2 <- table(visits1$arrival)
h2p2

# Poisson GLMM
m1 <- glmer(count ~ arrival + offset(log(offset_days)) + (1|trial), family=poisson, data=visits1)
m1
summary(m1)

m2<-glmer(count~offset(log(offset_days))+(1|trial),data=visits1,family=poisson)
anova(m1,m2)

emm <- emmeans(m1, "arrival")
eff_size(emm, sigma = sigma(m1), edf = df.residual(m1))

# ---- Figure 4a -----
# Plot of feeding duration before, after other species and arrivals after wax is depleted

f1 <- effect_plot(m1, pred = arrival, interval = TRUE, plot.points = TRUE,set.offset=1,ylim=c(0,10))
pred.means <- f1$data$count;pred.means 
lower <- pred.means-f1$data$ymin 
upper <- f1$data$ymax-pred.means
lower
upper

#plot predicted means and error bars
library(gplots)
plotCI(x=c(1:3),
       y=pred.means,
       xlim=c(0.5,4),
       ylim=c(0,85),
       uiw=upper, gap = 0,#upper error interval width
       liw=lower,#lower error interval width
       pch=19,col="white", cex=0.7,lwd=1,las = 1,
       ylab="number of visits",xlab="phase",xaxt="n",
       yaxp=c(0,80,4), bty="l", cex.lab=1,cex.axis=1.)

#add plot x-values for different arrival factor levels, A=1,B=2,C=3
visits$arrival1 <-visits$arrival
visits1$arrival1<-c(1)
visits1$arrival1[visits1$arrival=="B"]<-c(2)
visits1$arrival1[visits1$arrival=="C"]<-c(3)
axis(1,at=c(1,2,3),labels=c("","",""),cex.axis=1)
points(jitter(visits1$arrival1,0.02, amount=0.1),visits1$count,pch=20,col="darkgrey",cex = 0.6)
#axis(side=2, at=seq(0, 85, by=20))
mtext(c("before \n major wax comp","after \n major wax comp","wax \n depleted"),at=c(1,2,3),side=1, line=1.5,cex=1)
max(visits1$count)

# save as pdf 3.5 x 3.5
# save as pdf 4.3x3"

# ---- Figure 4b -----
# Boxplot of feeding duration before, after other species and arrivals after wax is depleted 

visits2 <- subset(visits,duration!="NA")
visits2$arrival <- as.factor(visits2$arrival)
visits2$trial <- as.factor(visits2$trial)
visits2$duration <- as.numeric(visits2$duration)
visits2$duration <- visits2$duration/60
visits2$duration[which(is.nan(visits2$duration))] = NA
visits2$duration[which(visits2$duration==Inf)] = NA
visits2 <- subset(visits2,duration > 0)
sum(is.na(visits2$duration))
#str(visits2)

#Gamma GLMM

m3 <- glmer(duration ~ arrival + (1|trial), family=Gamma, data=visits2)
m4 <- glmer(duration ~ (1|trial), family=Gamma, data=visits2)
anova(m3,m4)
summary(m3)
isSingular(m3)

emm <- emmeans(m3, "arrival")
eff_size(emm, sigma = sigma(m1), edf = df.residual(m1))


mcheck<-function(obj,...){
  rs<-resid(obj)
  fv<-fitted(obj)
  par(mfrow=c(1,3))
  plot(fv,rs,xlab="FITTED VALUES",ylab="RESIDUALS")
  abline(h=0,lty=2,lwd=2)
  qqnorm(rs,xlab="NORMAL SCORES",ylab="ORDERED RESIDUALS",main="")
  qqline(rs,lty=2,lwd=2)
  hist(resid(obj),xlab="RESIDUALS",ylab="FREQUENCY",main="")
  par(mfrow=c(1,1))
  invisible(NULL)
}

mcheck(m3) #Check residuals
mcheck(m1) #Check residuals


f2 <- effect_plot(m3, pred = arrival,interval = TRUE, plot.points = TRUE)
pred.means2 <- f2$data$duration;pred.means2 
lower2 <- pred.means2-f2$data$ymax #effect_plot outputs coordinate errorbar tip, we need width
upper2 <- f2$data$ymin-pred.means2
max(visits2$duration)

P <- table(visits2$arrival)
P

#plot predicted means and error bars
plotCI(x=c(1:3),
       y=pred.means2,
       xlim=c(0.5,3.5),
       ylim=c(0,25),
       uiw=upper2, gap = 0,#upper error interval width
       liw=lower2,#lower error interval width
       pch=19,col="black", cex=1,lwd=1,las = 1,
       ylab="duration of feeding visits",xlab="phase",xaxt="n",
       yaxp=c(0,10,1), bty="l", cex.lab=1,cex.axis=1.)

visits2$arrival1 <-visits2$arrival
visits2$arrival1<-c(1)
visits2$arrival1[visits2$arrival=="B"]<-c(2)
visits2$arrival1[visits2$arrival=="C"]<-c(3)
points(jitter(visits2$arrival1,0.02,amount=0.1),visits2$duration,pch=20,col="darkgrey",cex=0.6)
axis(1,at=c(1,2,3),labels=c("","",""),cex.axis=1)
axis(side=2, at=seq(0, 25, by=5))
mtext(c("before \n other species","after \n other species","wax \n depleted"),at=c(1,2,3),side=1, line=1.5,cex=1)

#save as pdf: 4x3 inches

# ---- Hypothesis 2, prediction (iii) -- availability of wax should decline quickly, driven by heterospecific competitors ----

# ----- Figure 3b ----
# K-M Survival curve of the wax

waxsurv <- waxdata %>%
  distinct(waxdata$trial,.keep_all = TRUE) #filter to a single duration for each site
waxsurv$waxduration <- waxsurv$waxduration/3600 #convert sec to hr

waxsurv$grp <- "1" #with MWCs feeding
waxsurv$grp[waxsurv$trial=="1"]<-"2" #without MWCs feeding
waxsurv$grp[waxsurv$trial=="2"]<-"2"
waxsurv$grp[waxsurv$trial=="21"]<-"2"
waxsurv$grp[waxsurv$trial=="22"]<-"2"
waxsurv$grp[waxsurv$trial=="26"]<-"2"
waxsurv$grp[waxsurv$trial=="3"]<-"2"
waxsurv$grp[waxsurv$trial=="5"]<-"2"
waxsurv$grp[waxsurv$trial=="9"]<-"2"
waxsurv$grp[waxsurv$trial=="4"]<-"2"
waxsurv$grp[waxsurv$trial=="20"]<-"2"

waxsurvival <- Surv(time = waxsurv$waxduration, event = waxsurv$censored~1,data=waxsurv)

waxsurv <- waxsurv %>%
  distinct(waxsurv$trial,.keep_all = TRUE)
# ggplot (with censored data)
waxsurvplot <- survfit(Surv(time = waxsurv$waxduration, event = waxsurv$censored) ~ 1, data = waxsurv)
waxsurvplot1 <-ggsurvplot(waxsurvplot, ggtheme=theme_classic2(),
                          xlab = "Time since wax availablity (h)",
                          ylab = "Probability of Wax Survival",
                          y.lim =c(0,1),
                          legend= 'none', 
                          censor.size = c(3),
                          censor.shape = c(124),
                          axis.offset = FALSE,
                          conf.int = TRUE,
                          conf.int.style = c("step"), #confidence interval style. Allowed values include c("ribbon", "step").
                          size = 0.7, # line thickness
                          palette=c("#006D2C"),
                          linetype = c("solid"),
                          censor =TRUE)
#waxsurvplot1$plot + geom_vline(xintercept = c(c(0,24,48,72,96,120,144,168)), color = "lightgrey", size=0.1, alpha =0.3) #to add 
waxsurvplot1
#save as pdf - dimensions 4 x 8 inches

summary(survfit(Surv(time = waxsurv$waxduration, event = waxsurv$censored) ~ 1, data = waxsurv), times = 24)
summary(survfit(Surv(time = waxsurv$waxduration, event = waxsurv$censored) ~ 1, data = waxsurv), times = 48)
summary(survfit(Surv(time = waxsurv$waxduration, event = waxsurv$censored) ~ 1, data = waxsurv), times = 72)

# plotting sites in which MWC did and didn't feed
waxsurv <- waxsurv %>%
  distinct(waxsurv$trial,.keep_all = TRUE)
# ggplot (with censored data)
waxsurvplot <- survfit(Surv(time = waxsurv$waxduration, event = waxsurv$censored) ~ grp, data = waxsurv)
waxsurvplot1 <-ggsurvplot(waxsurvplot, ggtheme=theme_classic2(),
                          xlab = "Time since wax availablity (h)",
                          ylab = "Probability of Wax Survival",
                          y.lim =c(0,1),
                          legend= 'none', 
                          censor.size = c(3),
                          censor.shape = c(124),
                          axis.offset = FALSE,
                          conf.int = TRUE,
                          #conf.int.style = c("step"), #confidence interval style. Allowed values include c("ribbon", "step").
                          size = 0.7, # line thickness
                          palette=c("#e27728","#0a68cd"),
                          linetype = c("solid","solid"),
                          censor =TRUE)

#b00101","#0a68cd","#464646","#91bfdb","#e27728","#0a8f77","#464646","#fca70a"
#waxsurvplot1$plot + geom_vline(xintercept = c(c(0,24,48,72,96,120,144,168)), color = "lightgrey", size=0.1, alpha =0.3) #to add 
waxsurvplot1


#  ---- Cox PH model - do some major wax eater significantly affect the survival of wax?
# test this using time-dependant survival analysis

waxsurv <- waxdata %>% #get all waxeaters
  filter(food_type == 'wax'|food_type == 'wax_brood') # only wax-eating events
waxsurv$group1 <- "A" 
waxsurv$group1[waxsurv$sp_com == "African civet"] <- "A" #group honey major wax eaters
waxsurv$group1[waxsurv$sp_com == "Honey badger"] <- "B"
waxsurv$group1[waxsurv$sp_com == "Yellow baboon"] <- "B"
waxsurv$group1[waxsurv$sp_com == "Mellers mongoose"] <- "B"
waxsurv <- waxsurv %>% group_by(trial, sp_com) %>% slice(which.min(timeafterplacement), na.rm = TRUE)
tally <- waxsurv %>% group_by(trial,group1) %>% tally() # get count by group
tally

#write.csv(waxsurv, "G:/Documents/waxsurv.csv", row.names = TRUE)

# load data which has been reformatted for this analysis
waxtd <- read.csv("wax_td_format.csv", sep=",", comment.char="#")

newwaxtd <- tmerge(data1=waxtd[, 1:4],data2=waxtd[, 1:14],id=trial, tstop=futime)
newwaxtd <- tmerge(newwaxtd, waxtd, id=trial,
                 infect = event(etime1))
newwaxtd <- tmerge(newwaxtd, waxtd, id=trial,
                   infect = event(etime2))
newwaxtd <- tmerge(newwaxtd, waxtd, id=trial,
                   infect = event(etime3))
newwaxtd <- tmerge(newwaxtd, waxtd, id=trial,
                   infect = event(etime4))
newwaxtd <- tmerge(newwaxtd, waxtd, id=trial,
                   infect = event(etime5))
newwaxtd <- tmerge(newwaxtd, waxtd, id=trial,
                   infect = event(etime6))
newwaxtd <- tmerge(newwaxtd, waxtd, id=trial,
                   infect = event(etime7))
newwaxtd <- tmerge(newwaxtd, waxtd, id=trial,
                   infect = event(etime8))
newwaxtd <- tmerge(newwaxtd, waxtd, id=trial,
                   infect = event(etime9))
newwaxtd <- tmerge(newwaxtd, waxtd, id=trial,
                   infect = event(etime10))
attr(newwaxtd, "tcount") #check

# run Cox PH regression with a time-dependent covariate

m5 <- coxph(Surv(time = tstart, time2 = tstop, event = censored) ~ grp, 
  data = newwaxtd) %>% 
  gtsummary::tbl_regression(exp = TRUE)
m5

# wax availability with and without major wax competitors

mwc <- waxdata %>% # get the first feeding event
  group_by(trial) %>%
  filter(row_number()==1) 
mwc$grp <- mwc$grp
mwc$grp[mwc$final_animal == "NA"] <- "C" #group honey major competitiors 
mwc$grp[mwc$final_animal == "civet"] <- "A"
mwc$grp[mwc$final_animal == "baboon"] <- "B"
mwc$grp[mwc$final_animal == "honeybadger"] <- "B"
mwc$grp[mwc$final_animal == "indicator"] <- "A" 
mwc$grp[mwc$final_animal == "mellers"] <- "B"
mwc$grp[mwc$final_animal == "minor"] <- "A"
mwc$grp[mwc$final_animal == "squirrel"] <- "A"

means <- mwc %>% filter(grp=="A" | grp =="B") %>% group_by(grp) %>%
  summarise(mean=mean(waxduration))
means

185375/3600
174554/3600
(48.4/51.5)*100
100-93.9

# ---- Figure 2 ---- 

#get date range for all data and determine sunrise/sunset times to add to figure
suntimes <- waxdata %>%
  mutate(waxplacement = as.Date(waxplacement)) %>%
  group_by(trial) %>%
  summarise(Placement = min(waxplacement), .groups = 'drop')

getSunlightTimes(date = as.Date(c("2015-09-24","2015-10-25","2017-08-29","2017-10-15","2021-09-24","2021-10-7" )), lat = -12.131556, lon =  38.080503,
                 keep = c("sunrise", "sunset"), tz = "Africa/Johannesburg")

sun <-getSunlightTimes(date = as.Date(c("2015-09-24","2015-10-25","2017-08-29","2017-10-15","2021-09-24","2021-10-7" )), lat = -12.131556, lon =  38.080503,
                       keep = c("sunrise", "sunset"), tz = "Africa/Johannesburg")
median(sun$sunrise)
median(sun$sunset)

#select wax eating events    
waxeat <- waxdata %>%
  filter(food_type == 'wax'|food_type == 'wax_brood')

#time of feeding events and new column for times as fractions of an hour
waxeat$timestart_plot <- hour(waxeat$timestart) + minute(waxeat$timestart)/60 + second (waxeat$timestart)/3600
waxeat$timefinish_plot <- hour(waxeat$timefinish) + minute(waxeat$timefinish)/60 + second (waxeat$timefinish)/3600

#calculate wax availability times on first day for adding to plot and convert to fractions of an hour
waxavailability <- waxeat %>%
  group_by(trial) %>%
  distinct(waxplacement)
waxavailability$time <- hour(waxavailability$waxplacement) + 
  minute(waxavailability$waxplacement)/60 + second(waxavailability$waxplacement)/3600


#reorder the species for the plot
waxeat <- waxeat %>%
  mutate(sp_com = fct_relevel(sp_com, "Greater honeyguide", "Scaly-throated honeyguide", "Lesser honeyguide",
                              "Crowned hornbill","Striped bush squirrel","Slender mongoose","Yellow baboon",
                              "Mellers mongoose","African civet","Honey badger"))

# (b) plot feeding events over 24 hrs
ggplot(waxeat, aes(colour=sp_com)) + 
  geom_segment(aes(x=timestart_plot, xend=timefinish_plot, y=sp_com, yend=sp_com),
               size=7, alpha = 0.6, show.legend = FALSE) +
  xlab("Hour (UTC+2h)") +
  theme_classic()+
  scale_color_manual(values=c("darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen")) +
  geom_vline(xintercept = c(5.25,17.4), color = "black", size=0.6, linetype = "dashed") + 
  # geom_vline(xintercept = c(0,1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,18,19,20,21,22,23), color = "lightgrey", size=0.2, alpha =0.5) +
  scale_x_continuous(breaks= c(2,5,8,11,14,17,20,23,11.1,16.3,11.2,11.3,17.0,10.7,12.2,16.7,17.3,8.74,15.1,5.17,10.1,17.0,13.4,12.8,11.5,10.1,11.5), 
                     labels = (c("2","5","8","11","14","17","20","23","_","_","_","_","_","_","_","_","_","_","_","_","_","_","_","_","_","_","_")))

#save as pdf with 4 x 10 inch dimensions


# (a) histogram of visits by each hour of day
hist <- waxeat %>%
  ggplot(aes(x=timestart_plot)) +
  geom_histogram( binwidth=1, fill="grey", color="white", alpha= 1) +
  xlab("Time of day") +
  theme_classic() +
  scale_x_continuous(breaks= seq(0,24,3)) +
  theme(plot.title = element_text(size=15))
hist #save as pdf: 3 x 10 inches

#get n counts for right-hand column
n_count <- waxdata %>%
  filter (food_type == 'wax'|food_type == 'wax_brood') %>%
  count(sp_com)
n_count


# ---- ESM Figure S3. Matrix of species occurrence by wax site ----

library('scales')

waxdata <- waxdata %>%
  filter(!(trial == "11"))
waxeaters <- subset(waxdata, sp_com =='Greater honeyguide'|sp_com =='Scaly-throated honeyguide'|sp_com =='Lesser honeyguide'|sp_com =='Striped bush squirrel' |sp_com =='Slender mongoose' |
                      sp_com =='Yellow baboon' |sp_com =='Mellers mongoose'|sp_com =='African civet' |sp_com =='Honey badger')
mat <- data.frame(waxeaters$trial,waxeaters$sp_com)
mat <- as.data.frame.matrix(table(mat))
mat <- mat[,c('Greater honeyguide','Scaly-throated honeyguide','Lesser honeyguide','Striped bush squirrel',
              'Slender mongoose','Yellow baboon','Mellers mongoose','African civet','Honey badger')]
#mat <- mat[-c(12, 13, 15,20,22,23,26),]
mat <- replace(mat,mat>0,'1') #convert to presence/absence data
mat <- replace(mat,mat==0,'0')
mat <- mat %>% 
  as.data.frame() %>%
  rownames_to_column("site") %>%
  pivot_longer(-c(site), names_to = "species", values_to = "counts")%>%
  ggplot(aes(x=site, y=species, fill=(counts))) + 
  geom_raster() +
  scale_fill_manual(values=c("white", "black")) +
  geom_vline(xintercept = c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,
                            11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5,21.5,22.5,23.5,24.5,25.5,26.5), color = "white", size=0.6, linetype = "solid")+ #to get gridlines
  geom_hline(yintercept = c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,
                            11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5), color = "white", size=0.6, linetype = "solid")+
  theme_bw()
mat #save as pdf 3 x 9 inches

# looking at wax feeding only

waxdata <- waxdata %>%
  filter(!(trial == "11"))
waxdata <- waxdata %>% 
  filter(food_type == 'wax'|food_type == 'wax_brood') # only wax-eating events
waxeaters <- subset(waxdata, sp_com =='Greater honeyguide'|sp_com =='Scaly-throated honeyguide'|sp_com =='Lesser honeyguide'|sp_com =='Striped bush squirrel' |sp_com =='Slender mongoose' |
                      sp_com =='Yellow baboon' |sp_com =='Mellers mongoose'|sp_com =='African civet' |sp_com =='Honey badger')
mat <- data.frame(waxeaters$trial,waxeaters$sp_com)
mat <- as.data.frame.matrix(table(mat))
mat <- mat[,c('Greater honeyguide','Scaly-throated honeyguide','Lesser honeyguide','Striped bush squirrel',
              'Slender mongoose','Yellow baboon','Mellers mongoose','African civet','Honey badger')]
#mat <- mat[-c(12, 13, 15,20,22,23,26),]
mat <- replace(mat,mat>0,'1') #convert to presence/absence data
mat <- replace(mat,mat==0,'0')
mat <- mat %>% 
  as.data.frame() %>%
  rownames_to_column("site") %>%
  pivot_longer(-c(site), names_to = "species", values_to = "counts")%>%
  ggplot(aes(x=site, y=species, fill=(counts))) + 
  geom_raster() +
  scale_fill_manual(values=c("white", "black")) +
  geom_vline(xintercept = c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,
                            11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5,21.5,22.5,23.5,24.5,25.5,26.5), color = "white", size=0.6, linetype = "solid")+ #to get gridlines
  geom_hline(yintercept = c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,
                            11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5), color = "white", size=0.6, linetype = "solid")+
  theme_bw()
mat #save as pdf 3.6 x 8 inches

