library(Cairo)
library(ggplot2)
library(DAAG)
library(stringr)
library(gridExtra)
library(plyr)

# the plotting in this script requires that "BiomassModelling.R" has already been run

# model performance plots -------------------------------------------------
# Fig. 2

# linear regression models (Fig. 2b)
lin.fit.val = lm(bio$linModel ~ bio$T_AG_BM)
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
            label = paste("RMSE =", round(summary(lin.fit.val)$sigma, 1)), hjust=0) +
  geom_text(data=NULL, x=100, y=15, 
            label = paste("Pearson's R =", round(cor.test(bio$linModel, bio$T_AG_BM)$estimate, 2)), hjust=0) +
  geom_text(data=NULL, x=25, y=175, label = "(b)", hjust=0)

# stratified model (Fig. 2c)
strat.fit.val = lm(bio$stratModel ~ bio$T_AG_BM)
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
            label = paste("RMSE =", round(summary(strat.fit.val)$sigma, 2)), hjust=0) +
  geom_text(data=NULL, x=100, y=15, 
            label = paste("Pearson's R =", round(cor.test(bio$stratModel, bio$T_AG_BM)$estimate, 2)), hjust=0) +
  geom_text(data=NULL, x=25, y=175, label = "(c)", hjust=0)

# summation model (Fig. 2a)
sum.fit.val = lm(bio$sumModel ~ bio$T_AG_BM)
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
            label = paste("RMSE =", round(summary(sum.fit.val)$sigma, 1)), hjust=0) +
  geom_text(data=NULL, x=100, y=15, 
            label = paste("Pearson's R =", round(cor.test(bio$sumModel, bio$T_AG_BM)$estimate, 2)), hjust=0) +
  geom_text(data=NULL, x=25, y=175, label = "(a)", hjust=0)

# combine in single single plot
CairoPDF(file="figure2_R.pdf", width=15, height=5, bg="transparent")
grid.arrange(component, mixed, dominant, nrow=1)
dev.off()


# sample size simulation --------------------------------------------------
# Fig. 3, Fig. i

# number of simualtion runs
nBoot = 100
# lower sample size limit
sample.min = 10
# total number of simulated points
sim.length = (nrow(bio)-(sample.min-1))*nBoot
# pre-allocate vectors
simN = numeric(sim.length)
sample.size = numeric(sim.length)
slope = numeric(sim.length)
offset = numeric(sim.length)
R2 = numeric(sim.length)
RMSE = numeric(sim.length)
# fill ticker
pos = 0
for (sim in 1:nBoot) {                # simulation run
  for (i in sample.min:nrow(bio)){    # sample size, begining at n=sample.min
    pos = pos+1
    # sample rows according to simulated sample n and fit model
    rows = sample(nrow(bio), i)
    fit = summary(lm(log(T_AG_BM) ~ T_SG, data=bio[rows,]))
    # fill in info
    simN[pos] = sim
    sample.size[pos] = i
    slope[pos] = fit$coefficients[2]
    offset[pos] = fit$coefficients[1]
    R2[pos] = fit$r.squared
    RMSE[pos] = fit$sigma
  }
}
sample.simulation = data.frame(simN, sample.size, slope, offset, R2, RMSE)

# plot it
size = 0.5

slope = ggplot(data=sample.simulation, aes(x=sample.size, y=slope)) +
  geom_point(size=size) + theme_classic() + geom_density2d(colour="red")

offset = ggplot(data=sample.simulation, aes(x=sample.size, y=offset)) +
  geom_point(size=size) + theme_classic() + geom_density2d(colour="red")

R2 = ggplot(data=sample.simulation, aes(x=sample.size, y=R2)) +
  geom_point(size=size) + theme_classic() + geom_density2d(colour="red")

RMSE = ggplot(data=sample.simulation, aes(x=sample.size, y=RMSE)) +
  geom_point(size=size) + theme_classic() + geom_density2d(colour="red")

# combine into single plot
CairoPDF(file="figure3_R.pdf", width=10, height=10, bg="transparent")
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


# annual WA binned biomass plots ------------------------------------------
# Fig. 6

# create data frame to store biomass bin value counts
biomass.bins.labels = c("1-20","21-40","41-60","61-80","81-100","101-120","121-140", "140+")
biomass.bins.plot = data.frame(bin=character(0), percentage=numeric(0), year=character(0), data.source=character(0))

# apply to photos
## function to use to calculate the number of photos of a certain biomass bin
BiomassBinSum = function(data){
  length(data$biomass.bin)
}

photo.sets = c("j2004", "a2007", "j2011", "f2012", "j2012", "f2013", "m2013")
for (i in 1:length(photo.sets)){
  photo = get(paste0("photos.",photo.sets[i]))
  totals = ddply(photo, .(biomass.bin), BiomassBinSum)
  totals = totals[totals$biomass.bin %in% biomass.bins.labels,]
  totals$V1 = totals$V1/sum(totals$V1)
  names(totals) = c("bin","percentage")
  #fill in detail
  totals$year = rep(photo.sets[i], nrow(totals))
  totals$data.source = rep("photos", nrow(totals))
  # stack together
  biomass.bins.plot = rbind(biomass.bins.plot, totals)
}

# apply to maps (where photos exsist)
## function to use to calculate the area sum in a certain biomass bin
BiomassBinArea = function(data){
  sum(data$Area, na.rm=T)
}

map.sets = c("j2004", "a2007", "j2011", "f2012", "j2012", "f2013", "m2013")
for (i in 1:length(map.sets)){
  map = get(paste0("map.",map.sets[i]))
  totals = ddply(map, .(biomass.bin), BiomassBinArea)
  totals = totals[totals$biomass.bin %in% biomass.bins.labels,]
  totals$V1 = totals$V1/sum(totals$V1)
  names(totals) = c("bin","percentage")
  #fill in detail
  totals$year = rep(map.sets[i], nrow(totals))
  totals$data.source = rep("map", nrow(totals))
  # stack together
  biomass.bins.plot = rbind(biomass.bins.plot, totals)
}

# plot the bin data
biomass.bins.plot$year = as.factor(biomass.bins.plot$year)
levels(biomass.bins.plot$year) = map.sets

# create correlation plot between photos/map percentages
biomass.bin.corr = biomass.bins.plot
biomass.bin.corr$id = paste0(biomass.bin.corr$bin,".",biomass.bin.corr$year)
biomass.bin.corr = merge(biomass.bin.corr[biomass.bin.corr$data.source=="photos",], 
                         biomass.bin.corr[biomass.bin.corr$data.source=="map",],
                         by="id")

# plot it
CairoPDF(file="figure6_R.pdf", width=7, height=5, bg="transparent")
size = 3
ggplot(data=biomass.bin.corr, aes(x=percentage.x, y=percentage.y)) +
  geom_point(aes(fill=bin.x), size=size, shape=21) + theme_classic() +
  scale_fill_manual(name="Biomass: gDW/m2",
                      values=c("#de2d26","#fc9272","#edf8e9","#bae4b3","#74c476","#238b45","#000000")) +
  xlab("photo predicted % of landscape") + ylab("map predicted % of landscape") +
  geom_abline(intercept=0,slope=1, colour="red", linetype="twodash")
dev.off()
