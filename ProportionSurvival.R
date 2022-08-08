#Load Packages and Datasheets
{library(AER)
  library(ecotox)
  library(plyr)
  library(dplyr)
  library(ggplot2)
  library(Rmisc)
  library(RColorBrewer)
  library(emmeans)
  library(gam)
  library(mgcv)
  library(multcomp)
}

# Survival curves for Culex pipiens, only need to change the imported data sheet name
ds<-CulexSurvivalFinal
attach(ds)
ds$RH<-as.factor(ds$RH)
ds$BF<-as.factor(ds$BF)
ds
attach(ds)

#Split the dataset into appropriate groups (75% and 100% RH)
ds2<-dlply(ds,.(RH))
ds100<-ds2$'100'
ds75<-ds2$'75'
ds75
ds100
plot(PropSurv~Time,data=ds75)
plot(PropSurv~Time,data=ds100)

ds3<-dlply(ds75,.(BF))
ds75bf<-ds3$'Y'
ds75nbf<-ds3$'N'
ds75bf
ds75nbf
plot(PropSurv~Time,data=ds75bf)
plot(PropSurv~Time,data=ds75nbf)

ds4<-dlply(ds100,.(BF))
ds100bf<-ds4$'Y'
ds100nbf<-ds4$'N'
ds100bf
ds100nbf
plot(PropSurv~Time,data=ds100bf)
plot(PropSurv~Time,data=ds100nbf)

Totds<-rbind(ds75nbf,ds75bf,ds100nbf,ds100bf)
Totds
attach(Totds)

Totds2<-mutate(Totds,Mutated=paste(RH,BF,sep="_"))
attach(Totds2)
summary(Totds2)

CTotds<-Totds2 %>% mutate(Treatment = factor(Mutated, levels = c('75_N','75_Y','100_N','100_Y')))
attach(CTotds)
summary(CTotds)


#Start Here
#Geom_smooth uses stats::loess() for < 1,000 observations; otherwise mgcv::gam()

DS1<-ggplot(ds75nbf, aes(x=Time, y=PropSurv)) +
  geom_smooth(aes(x=Time, y=PropSurv),method=mgcv::gam,formula = y ~ s(x, bs = "cs"))
DS1

#Will generate the same formula, check 'ggplot' against 'plot' graph for verification
cstat1<-mgcv::gam(data=ds75nbf,formula = PropSurv ~ s(Time, bs = "cs"))
cstat1
plot(cstat1)

DS2<-ggplot(ds75bf, aes(x=Time, y=PropSurv)) +
  geom_smooth(aes(x=Time, y=PropSurv),method=mgcv::gam,formula = y ~ s(x, bs = "cs"))
DS2

cstat2<-mgcv::gam(data=ds75bf,formula = PropSurv ~ s(Time, bs = "cs"))
cstat2
plot(cstat2)

DS3<-ggplot(ds100nbf, aes(x=Time, y=PropSurv)) +
  geom_smooth(aes(x=Time, y=PropSurv),method=mgcv::gam,formula = y ~ s(x, bs = "cs"))
DS3

cstat3<-mgcv::gam(data=ds100nbf,formula = PropSurv ~ s(Time, bs = "cs"))
cstat3
plot(cstat3)

DS4<-ggplot(ds100bf, aes(x=Time, y=PropSurv)) +
  geom_smooth(aes(x=Time, y=PropSurv),method=mgcv::gam,formula = y ~ s(x, bs = "cs"))
DS4

cstat4<-mgcv::gam(data=ds100bf,formula = PropSurv ~ s(Time, bs = "cs"))
cstat4
plot(cstat4)


#Make a color palette (I have an affinity for grey) and create the combined figure
TreatColor = c('75_N'='grey30','75_Y'='grey50','100_N'='grey70','100_Y'='grey90')

DS5<-ggplot(CTotds,aes(x=Time, y=PropSurv, color=Treatment)) +
  geom_smooth(method=mgcv::gam,formula = y ~ s(x, bs = "cs")) +
  xlab("\nTime (hours)") +
  ylab("Proportion Surviving\n") +
  theme_classic() +
  theme(legend.position = c(0.75, 0.8)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(ylim=c(0,1.01),xlim=c(0,145)) +
  theme(axis.text.x=element_text(size=14,color="black"), axis.title=element_text(size=18,color="black")) +
  theme(axis.text.y=element_text(size=14,color="black"), axis.title=element_text(size=18,color="black")) +
  guides(guide_legend(title="Treatment", override.aes = list(size=.2))) +
  theme(panel.background=element_rect(fill="white",color="white")) +
  theme(axis.ticks=(element_line(color="black")),axis.line=(element_line(color="black"))) + 
  scale_color_manual(values=TreatColor) 
DS5
ggsave(filename="PropSurvSmooth.png",units=c("in"),width=9,height=7,plot=DS5)

#Identify points along the smoothed curve to compare

attach(CulexSurvivalFinal)
cs<-CulexSurvivalFinal

cs2<- mutate(cs, Treatment=paste(RH,BF, sep="_"))
cs2
attach(cs2)
as.numeric(Dead)
as.numeric(Total)
as.numeric(Time)
as.factor(Treatment)

m1 <- LC_probit((Dead / Total) ~ Time, p = c(25,50,75,99),
                weights = Total,
                data = cs2[cs2$Time != 0, ],
                subset = c(Treatment == "75_N"),
                log_x = FALSE)

m2 <- LC_probit(cbind(Dead, Alive) ~ Time, p = c(25,50,75,99),
                data = cs2[cs2$Time != 0, ],
                subset = c(Treatment == "75_N"),
                log_x = FALSE)

m1
m2

# dose-response curve can be plotted using 'ggplot2'

DN <- subset(cs2, Treatment %in% c("75_N"))
p1 <- ggplot(data = DN[DN$Time != 0, ],
             aes(x = Time, y = (Dead / Total))) +
  geom_point() +
  geom_smooth(method = "loess",
              #              method.args = list(family = quasibinomial(link = "probit")),
              aes(weight = Total), colour = "#FF0000", se = TRUE)
p1

# calculate LCs (25,50,75,99) for multiple treatments to determine differences between groups at those points
DN <- LC_probit((Dead / Total) ~ Time, p = c(25,50,75,99),
                weights = Total,
                conf_level = 0.95,
                data = cs2[cs2$Time != 0, ],
                subset = c(Treatment == "75_N"),
                log_x = FALSE)
DY <- LC_probit((Dead / Total) ~ Time, p = c(25,50,75,99),
                weights = Total,
                conf_level = 0.95,
                data = cs2[cs2$Time != 0, ],
                subset = c(Treatment == "75_Y"),
                log_x = FALSE)
NN <- LC_probit((Dead / Total) ~ Time, p = c(25,50,75,99),
                weights = Total,
                conf_level = 0.95,
                data = cs2[cs2$Time != 0, ],
                subset = c(Treatment == "100_N"),
                log_x = FALSE)
NY <- LC_probit((Dead / Total) ~ Time, p = c(25,50,75,99),
                weights = Total,
                conf_level = 0.95,
                data = cs2[cs2$Time != 0, ],
                subset = c(Treatment == "100_Y"),
                log_x = FALSE)

# group results together in a dataframe to plot with 'ggplot2'
results <- rbind(DN[, c(1, 3:8, 11)], DY[,c(1, 3:8, 11)],
                 NN[, c(1, 3:8, 11)], NY[, c(1, 3:8, 11)])

results$Treatment <- factor(c(rep("75_N", 4), rep("75_Y", 4),
                              rep("100_N", 4), rep("100_Y", 4)),
                            levels = c("75_N", "75_Y", "100_N", "100_Y"))

palettegrey<-c("#444444","#666666","#999999","#CCCCCC")

p2 <- ggplot(data = results, aes(x = Treatment, y = dose,
                                 group = factor(p), fill = factor(p))) +
  geom_col(position = position_dodge(width = 0.9), colour = "#000000") +
  geom_errorbar(aes(ymin = (LCL), ymax = (UCL)),
                size = 0.4, width = 0.06,
                position = position_dodge(width = 0.9)) +
  xlab("\nTreatment") +
  ylab("Time (Hours)\n") +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  theme_classic() +
  theme(axis.text.y=element_text(size=14,color="black"), axis.title=element_text(size=18,color="black")) +
  theme(axis.text=element_text(size=12,color="black"), axis.title=element_text(size=18,color="black")) +
  theme(panel.background=element_rect(fill="white",color="white"))+
  theme(legend.text=(element_text(size=10)))+
  theme(axis.text=(element_text(size=12,color="black")),axis.title=(element_text(size=14,color="black")))+
  theme(axis.ticks=(element_line(color="black")),axis.line=(element_line(color="black"))) +
  scale_fill_manual(values=c(palettegrey),name="Lethal Time",labels=c("25","50","75","99")) + 
  scale_color_manual(values=c(palettegrey),name="Lethal Time",labels=c("25","50","75","99")) 
p2
ggsave(filename="LDComparisons.png",units=c("in"),width=9,height=7,plot=p2)

# Source: Hlina, B.L., Birceanu, O., Robinson, C.S., Dhiyebi, H., Wilkie, M.P. In Reivew.
# Seasonal Variation in the Sensitivity of Invasive Sea Lampreys to the Lampricide 
# FM: Importance of Energy Reserves and Temperature. North American Journal of Fisheries Management

#Now for the stats (To determine the overall similarities/differences between GAMs)
#Make the models
#Culex for emmeans
CComp<-mgcv::gam(data=CTotds,
                  formula = PropSurv ~ s(Time, bs = "cs", by = Treatment),
                  na.action = na.gam.replace)

#check for normality?
CComp
plot(CComp)
plot(CComp)
plot(CComp)
plot(CComp)

#Compare the survival treatments
#Write and export all important data from here on out

#Culex
#Using emtrends looks at many means over the curves to compare entire trends between groups
Cestimatedtrends=emtrends(CComp, pairwise ~ "Treatment", var="Time")
Cestimatedtrends
write.csv(Cestimatedtrends$contrasts,"Culex_GAM_Trends.csv",row.names=TRUE)

#Compare the LDs between treatments
#View the data generated by ecotox (results) and use that data to compare LDs

results

#First compare the curves individually
#The following time points were determined via previous probit analyses (e.g. DN)
#These values are found in the 'dose' column
#These numbers can also be called in other manners, I just manually inserted them here
CSE1<-emmeans(cstat1, specs = pairwise ~ Time,
              at=list(Time=c(0,16.7,38.4,60,113)))
CSE1
write.csv(CSE1$contrasts,"C75N_Contrasts.csv",row.names=TRUE)
write.csv(CSE1$emmeans,"C75N_Emmeans.csv",row.names=TRUE)

CSE2<-emmeans(cstat2, specs = pairwise ~ Time,
              at=list(Time=c(0,35,51.7,68.3,109)))
CSE2
write.csv(CSE2$contrasts,"C75Y_Contrasts.csv",row.names=TRUE)
write.csv(CSE2$emmeans,"C75Y_Emmeans.csv",row.names=TRUE)

CSE3<-emmeans(cstat3, specs = pairwise ~ Time,
              at=list(Time=c(0,41.7,70.6,99.4,170)))
CSE3
write.csv(CSE3$contrasts,"C100N_Contrasts.csv",row.names=TRUE)
write.csv(CSE3$emmeans,"C100N_Emmeans.csv",row.names=TRUE)

CSE4<-emmeans(cstat4, specs = pairwise ~ Time,
              at=list(Time=c(0,65.7,89.1,113,170)))
CSE4
write.csv(CSE4$contrasts,"C100Y_Contrasts.csv",row.names=TRUE)
write.csv(CSE4$emmeans,"C100Y_Emmeans.csv",row.names=TRUE)

#Then compare all treatments and LDs, selecting the relevant comparisons from Excel later :(
#A more streamlined method can be employed here, I just saw the light at the end of the tunnel and ran for it.
CContrasts<-emmeans(CComp, specs = pairwise ~ Time + factor(Treatment),
        at=list(Time=c(0,16.7,38.4,60.0,113,35.0,51.7,68.3,109,41.7,70.6,99.4,170,65.7,89.1)))
CContrasts
write.csv(CContrasts$emmeans,"Culex_LT_Emmeans.csv",row.names=TRUE)
write.csv(CContrasts$contrasts,"Culex_LT_Contrasts.csv",row.names=TRUE)