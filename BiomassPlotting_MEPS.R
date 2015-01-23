## Author: Mitchell Lyons (mitchell.lyons@gmail.com)
## No licence, though please attribute me and publish any derivitives open source =)

library(Cairo)
library(ggplot2)
library(stringr)
library(gridExtra)
library(plyr)
library(mgcv)

# the plotting in this script requires that "BiomassModelling.R" has already been run

# model performance plots -------------------------------------------------
# Fig. 2

## obs vs. predicted plots
#CairoWin()
#CairoPNG(filename="strat.fit.val.png", width=600, height=600)
mixed = ggplot(bio, aes(x=T_AG_BM, y=linModel)) +
  geom_point(shape=19) +
  #geom_text(aes(label=row.names(bio)),hjust=0, vjust=0) +
  geom_smooth(method=lm, alpha=0.2) +
  geom_abline(intercept=0,slope=1, colour="red", linetype="twodash") +
  coord_cartesian(ylim=c(0,200), xlim=c(0,200)) +
  theme_classic() +
  xlab("Observed biomass (gDW.m^-2)") +
  ylab("Fitted biomass (gDW.m^-2)") +
  ggtitle("Mixed species linear model") +
  theme(legend.justification=c(0,1), legend.position=c(0,1), 
        legend.background = element_rect(fill="transparent")) +
  #add R/rms
  geom_text(data=NULL, x=100, y=25, 
            label = paste("Overall RMSE = ", round(sqrt(mean((bio$linModel-bio$T_AG_BM)^2)), 0)), hjust=0) +
  geom_text(data=NULL, x=175, y=10, label = "(a)", hjust=0) + 
  geom_vline(xintercept=c(25,50,75), linetype="dotted")

#dev.off()

# stratified model
#CairoWin()
#CairoPNG(filename="strat.fit.val.png", width=600, height=600)
dominant = ggplot(bio, aes(x=T_AG_BM, y=stratModel)) +
  geom_point(shape=19) +
  geom_smooth(method=lm, alpha=0.2) +
  geom_abline(intercept=0,slope=1, colour="red", linetype="twodash") +
  coord_cartesian(ylim=c(0,200), xlim=c(0,200)) +
  theme_classic() +
  xlab("Observed biomass (gDW.m^-2)") +
  ylab("Fitted biomass (gDW.m^-2)") +
  ggtitle("Dominant species biomass model") +
  theme(legend.justification=c(0,1), legend.position=c(0,1), 
        legend.background = element_rect(fill="transparent")) +
  #add R/rms
  geom_text(data=NULL, x=100, y=25, 
            label = paste("Overall RMSE =", round(sqrt(mean((bio$stratModel-bio$T_AG_BM)^2)), 0)), hjust=0) +
  geom_text(data=NULL, x=175, y=10, label = "(c)", hjust=0) + 
  geom_vline(xintercept=c(25,50,75), linetype="dotted")
#dev.off()


# summation model
#CairoWin()
#CairoPNG(filename="sum.fit.val.png", width=600, height=600)
component = ggplot(bio, aes(x=T_AG_BM, y=sumModel)) +
  geom_point(shape=19) +
  stat_smooth(method=lm, alpha=0.2) +
  geom_abline(intercept=0,slope=1, colour="red", linetype="twodash") +
  coord_cartesian(ylim=c(0,200), xlim=c(0,200)) +
  theme_classic() +
  xlab("Observed biomass (gDW.m^-2)") +
  ylab("Fitted biomass (gDW.m^-2)") +
  ggtitle("Species component biomass model") +
  theme(legend.justification=c(0,1), legend.position=c(0,1), 
        legend.background = element_rect(fill="transparent")) +
  #add R/rms
  geom_text(data=NULL, x=100, y=25, 
            label = paste("Overall RMSE =", round(sqrt(mean((bio$sumModel-bio$T_AG_BM)^2)), 0)), hjust=0) +
  geom_text(data=NULL, x=175, y=10, label = "(b)", hjust=0) + 
  geom_vline(xintercept=c(25,50,75), linetype="dotted")
#dev.off()

# in a single plot window
CairoPDF(file="../figures/figure2_R.pdf", width=15, height=5, bg="transparent")
#CairoWin()
grid.arrange(mixed, component, dominant, nrow=1)
dev.off()


# sample size simulation --------------------------------------------------
# Fig. 3, Fig. i

# number of simualtion runs
nBoot = 10000
# lower sample size limit
sample.min = 10
# total number of simulated points
sim.length = (nrow(bio)-(sample.min-1))*nBoot
# pre-allocate vectors
simN = numeric(sim.length)
sample.size = numeric(sim.length)
slope.perm = numeric(sim.length)
offset.perm = numeric(sim.length)
R2.perm = numeric(sim.length)
RMSE.perm = numeric(sim.length)
slope.boot = numeric(sim.length)
offset.boot = numeric(sim.length)
R2.boot = numeric(sim.length)
RMSE.boot = numeric(sim.length)
# fill ticker
pos = 0
for (sim in 1:nBoot) {                # simulation run
  for (i in sample.min:nrow(bio)){    # sample size, begining at n=sample.min
    pos = pos+1
    simN[pos] = sim
    sample.size[pos] = i
    # sample rows without replacement
    rows = sample(x=nrow(bio), size=i, replace=F)
    fit = lm(log(T_AG_BM) ~ T_SG, data=bio[rows,])
    slope.perm[pos] = coef(fit)[2]
    offset.perm[pos] = coef(fit)[1]
    R2.perm[pos] = summary(fit)$r.squared
    RMSE.perm[pos] = sqrt(mean((((exp(1)^predict.lm(fit))-(exp(1)^fit$model[,1]))^2)))
    # sample rows with replacement
    rows = sample(x=nrow(bio), size=i, replace=T)
    fit = lm(log(T_AG_BM) ~ T_SG, data=bio[rows,])
    slope.boot[pos] = coef(fit)[2]
    offset.boot[pos] = coef(fit)[1]
    R2.boot[pos] = summary(fit)$r.squared
    RMSE.boot[pos] = sqrt(mean((((exp(1)^predict.lm(fit))-(exp(1)^fit$model[,1]))^2)))
  }
}
sample.simulation = data.frame(simN, sample.size, 
                               slope.perm, offset.perm, R2.perm, RMSE.perm,
                               slope.boot, offset.boot, R2.boot, RMSE.boot)

## functions to do summary calculations
## mean function
MeanNumeric = function(x){
  ifelse(is.numeric(x), mean(x), x)
}
## 5% quantile function
quantile.025 = function(x){
  ifelse(is.numeric(x), quantile(x, 0.025), x)
}
## 95% quantile function
quantile.975 = function(x){
  ifelse(is.numeric(x), quantile(x, 0.975), x)
}


# claculate mean and CI's for simulation at each sample size
sample.simulation.means = ddply(sample.simulation[,-1], .(sample.size), colwise(MeanNumeric))
sample.simulation.025 = ddply(sample.simulation[,-1], .(sample.size), colwise(quantile.025))
sample.simulation.975 = ddply(sample.simulation[,-1], .(sample.size), colwise(quantile.975))

# plot it
size = 1

slope = ggplot(data=sample.simulation.means, aes(x=sample.size)) +
  geom_point(aes(y=slope.perm), size=size) +
  geom_errorbar(aes(ymin=sample.simulation.025$slope.boot,ymax=sample.simulation.975$slope.boot), colour="black", width=0.5) +
  geom_errorbar(aes(ymin=sample.simulation.025$slope.perm,ymax=sample.simulation.975$slope.perm), colour="red", width=0.5) +
  theme_classic()

offset = ggplot(data=sample.simulation.means, aes(x=sample.size)) +
  geom_point(aes(y=offset.perm), size=size) +
  geom_errorbar(aes(ymin=sample.simulation.025$offset.boot,ymax=sample.simulation.975$offset.boot), colour="black", width=0.5) +
  geom_errorbar(aes(ymin=sample.simulation.025$offset.perm,ymax=sample.simulation.975$offset.perm), colour="red", width=0.5) +
  theme_classic()

R2 = ggplot(data=sample.simulation.means, aes(x=sample.size)) +
  geom_point(aes(y=R2.perm), size=size) +
  geom_errorbar(aes(ymin=sample.simulation.025$R2.boot,ymax=sample.simulation.975$R2.boot), colour="black", width=0.5) +
  geom_errorbar(aes(ymin=sample.simulation.025$R2.perm,ymax=sample.simulation.975$R2.perm), colour="red", width=0.5) +
  theme_classic()

RMSE = ggplot(data=sample.simulation.means, aes(x=sample.size)) +
  geom_point(aes(y=RMSE.perm), size=size) +
  geom_errorbar(aes(ymin=sample.simulation.025$RMSE.boot,ymax=sample.simulation.975$RMSE.boot), colour="black", width=0.5) +
  geom_errorbar(aes(ymin=sample.simulation.025$RMSE.perm,ymax=sample.simulation.975$RMSE.perm), colour="red", width=0.5) +
  theme_classic()


#CairoWin()
CairoPDF(file="../figures/figure3_R.pdf", width=10, height=8, bg="transparent")
grid.arrange(slope, offset, R2, RMSE, ncol=2)
dev.off()


# photo transects panel plot ----------------------------------------------
# Fig. 4

# panel plot showing photo number histograms and autotrophic abline
binwidth = 5
yrange = c(0,1000)
xrange = c(0,115)
vs = 1.5 # vline thickness

# set plot colours - blue, red, green, purple, orange
values = c("AM"="#377eb8","CH"="#e41a1c","MA"="#4daf4a","MO"="#984ea3","WA"="#ff7f00")

# year to plot
photos.year = photos.j2012

Hs = ggplot(photos.year[photos.year$Hs_biomass>0&photos.year$corrSumModel!=0,], aes(x=Hs_biomass, fill=bank)) +
  scale_fill_manual(values=values) + geom_histogram(binwidth=binwidth) + guides(fill=FALSE) +
  theme_classic() + coord_cartesian(xlim=xrange) + 
  geom_vline(xintercept=41.8, linetype="dotted", size=vs) + geom_vline(xintercept=14.5, linetype="longdash", size=vs)

Cr = ggplot(photos.year[photos.year$Cr_biomass>0&photos.year$corrSumModel!=0,], aes(x=Cr_biomass, fill=bank)) +
  scale_fill_manual(values=values) + geom_histogram(binwidth=binwidth) + guides(fill=FALSE) +
  theme_classic() + coord_cartesian(xlim=xrange) + 
  geom_vline(xintercept=41.8, linetype="dotted", size=vs) + geom_vline(xintercept=67, size=vs)

ZmHu = ggplot(photos.year[photos.year$ZmHu_biomass>0&photos.year$corrSumModel!=0,], aes(x=ZmHu_biomass, fill=bank)) +
  scale_fill_manual(values=values) + geom_histogram(binwidth=binwidth) + guides(fill=FALSE) +
  theme_classic() + coord_cartesian(xlim=xrange) + 
  geom_vline(xintercept=41.8, linetype="dotted", size=vs) + geom_vline(xintercept=92.7, linetype="longdash", size=vs)

Ho = ggplot(photos.year[photos.year$Ho_biomass>0&photos.year$corrSumModel!=0,], aes(x=Ho_biomass, fill=bank)) +
  scale_fill_manual(values=values) + geom_histogram(binwidth=binwidth) + guides(fill=FALSE) +
  theme_classic() + coord_cartesian(xlim=xrange) + 
  geom_vline(xintercept=41.8, linetype="dotted", size=vs) + geom_vline(xintercept=14.5, size=vs)

Si = ggplot(photos.year[photos.year$Si_biomass>0&photos.year$corrSumModel!=0,], aes(x=Si_biomass, fill=bank)) +
  scale_fill_manual(values=values) + geom_histogram(binwidth=binwidth) + guides(fill=FALSE) +
  theme_classic() + coord_cartesian(xlim=xrange) + 
  geom_vline(xintercept=41.8, linetype="dotted", size=vs) + geom_vline(xintercept=32, size=vs)

total = ggplot(photos.year[photos.year$corrSumModel!=0,], aes(x=corrSumModel, fill=bank)) +
  scale_fill_manual(values=values) + geom_histogram(binwidth=binwidth) + guides(fill=FALSE) +
  theme_classic() + coord_cartesian(xlim=xrange) + 
  geom_vline(xintercept=41.8, linetype="dotted", size=vs)

# combine into single plot
CairoPDF(file="figure4_R.pdf", width=6, height=6, bg="transparent")
grid.arrange(ZmHu, Cr, Ho, Hs, Si, total, ncol=2)
dev.off()


# annual WA biomass plot --------------------------------------------------
# Fig. 5

## function to calculate bioamss dry weight (in tonnes) from DW.m^-2 and polygon area (for maps)
CalcBiomassDW.map = function(data){
  sum(data$biomass*data$Area, na.rm=TRUE)/1000000
}
# bank of interest
bank = "WA" 
# loop to get biomass totals for each year
biomass.map.ts = data.frame() #should pre-allocate here if you have many years or columns...
map.sets = c("j2004", "a2007", "2008", "d2009", "j2011", "f2012", "j2012", "f2013", "m2013")
for (i in 1:length(map.sets)){
  map = get(paste0("map.",map.sets[i]))
  biomass.map.ts[i,1] = paste0(map.sets[i],".WA")
  biomass.map.ts[i,2] = CalcBiomassDW.map(map[map$Bank==bank,])
  rm(map)
}
names(biomass.map.ts) = c("map.year", "map.biomass")

## function to calculate bioamss dry weight (in kg) from DW.m^-2 and photo area (for photos)
CalcBiomassDW.photo = function(data, PhotoArea){
  sum(data$corrSumModel*PhotoArea, na.rm=TRUE)/1000
}
# bank of interest
bank = "WA" 
# loop to get biomass totals for each year
biomass.photo.ts = data.frame() #should pre-allocate here if you have many years or columns...
photo.sets = c("j2004", "a2007", "j2011", "f2012", "j2012", "f2013", "m2013")
for (i in 1:length(photo.sets)){
  photo = get(paste0("photos.",photo.sets[i]))
  biomass.photo.ts[i,1] = paste0(photo.sets[i],".WA")
  biomass.photo.ts[i,2] = CalcBiomassDW.photo(photo[photo$bank==bank,], 1) #benthic photo area of 1 m^2
  rm(photo)
}
names(biomass.photo.ts) = c("photo.year", "photo.biomass")

# combine map and photo totals
biomass.ts = merge(biomass.map.ts, biomass.photo.ts, by.x="map.year", by.y="photo.year", all.x=TRUE)

# bar plot of map biomass DW
# order the map years
biomass.map.ts$map.year = as.factor(biomass.map.ts$map.year)
levels(biomass.map.ts$map.year) = map.sets

# plot it
CairoPDF(file="figure5_R.pdf", width=5, height=5, bg="transparent")
ggplot(data=biomass.map.ts, aes(x=factor(map.year), y=map.biomass)) + 
  geom_bar(stat="identity") + theme_classic() +
  ylab("biomass (tonnes DW)") + xlab("mapping month/year")
dev.off()

