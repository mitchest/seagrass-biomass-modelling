## Author: Mitchell Lyons (mitchell.lyons@gmail.com)
## No licence, though please attribute me and publish any derivitives open source =)

library(Cairo)
library(ggplot2)
library(stringr)
library(gridExtra)
library(plyr)
library(mgcv)

# set working directory to data and code location
setwd("C:/my/data/location")


# load data and prep for modelling ----------------------------------------

# load biomass core and coincident photo data
bio = read.csv("BiomassCoreData.csv", 
         header=TRUE, row.names=1, stringsAsFactors=FALSE)


Cth = 0.55 # define cover threshold for "dominant"

# calculate proportion of each species in biomass data
bio$pZmHu = bio$ZmHu_AG_BM/bio$T_AG_BM
bio$pHo = bio$Ho_AG_BM/bio$T_AG_BM
bio$pHs = bio$Hs_AG_BM/bio$T_AG_BM
bio$pCr = bio$Cr_AG_BM/bio$T_AG_BM
bio$pSi = bio$Si_AG_BM/bio$T_AG_BM

# calculate proportion of each species in percent cover data
bio$phZmHu = bio$SGZMHU/bio$T_SG
bio$phHo = bio$SGHO/bio$T_SG
bio$phHs = bio$SGHS/bio$T_SG
bio$phCr = bio$SGCR/bio$T_SG
bio$phSi = bio$SGSI/bio$T_SG

# create column for dominant species based on biomass and Cth threshold
dom.species.char = character(length(bio$T_SG))
for (i in 1:length(bio$T_SG)) {
  if (bio$pZmHu[i]>Cth) {
    dom.species.char[i] = "ZmHu"
  } else if (bio$pHo[i]>Cth) {
    dom.species.char[i] = "Ho"
  } else if (bio$pHs[i]>Cth) {
    dom.species.char[i] = "Hs"
  } else if (bio$pCr[i]>Cth) {
    dom.species.char[i] = "Cr"
  } else if (bio$pSi[i]>Cth) {
    dom.species.char[i] = "Si"
  } else {
    dom.species.char[i] = "mixed"
  }
}
bio$dom.biomass.species = as.factor(dom.species.char)

# create column for dominant species based on biomass and Cth threshold
dom.species.char = character(length(bio$T_SG))
for (i in 1:length(bio$T_SG)) {
  if (bio$phZmHu[i]>Cth) {
    dom.species.char[i] = "ZmHu"
  } else if (bio$phHo[i]>Cth) {
    dom.species.char[i] = "Ho"
  } else if (bio$phHs[i]>Cth) {
    dom.species.char[i] = "Hs"
  } else if (bio$phCr[i]>Cth) {
    dom.species.char[i] = "Cr"
  } else if (bio$phSi[i]>Cth) {
    dom.species.char[i] = "Si"
  } else {
    dom.species.char[i] = "mixed"
  }
}
bio$dom.photo.species = as.factor(dom.species.char)


# fitting total biomass ~ total cover -------------------------------------


# LM - one linear term
cover.lm.log = lm(log(T_AG_BM) ~ T_SG, data=bio)
cover.lm = lm(T_AG_BM ~ T_SG, data=bio)

PrintLM(cover.lm,10)
boot.lm(cover.lm,10,100)

PrintLM.log(cover.lm.log,10)
boot.lm.log(cover.lm.log,10,100)


# 2nd orer term
cover.lm2 = lm(log(T_AG_BM) ~ poly(T_SG,2), data=bio)
cover.lm2 = lm(log(T_AG_BM) ~ T_SG + I(T_SG^2), data=bio)
PrintLM.log(cover.lm2,10)
boot.lm.log(cover.lm2,10,100)

# 3rd order term
cover.lm2 = lm(log(T_AG_BM) ~ poly(T_SG,3), data=bio)
cover.lm3 = lm(log(T_AG_BM) ~ T_SG + I(T_SG^2) + I(T_SG^3), data=bio)
PrintLM.log(cover.lm3,10)
boot.lm.log(cover.lm3,10,100)

# GLM
cover.glm = glm(T_AG_BM ~ T_SG, family="Gamma", data=bio)
PrintGLM(cover.glm, 10)
boot.glm(cover.glm, 10, 100)

# GAM
cover.gam = gam(T_AG_BM ~ s(T_SG), family="Gamma", data=bio)
PrintGAM(cover.gam, 10)
boot.gam(cover.gam, 10, 100)

# write selected model predictions to bio data frame
bio$linModel = exp(1)^predict(cover.lm.log, bio)

# error at different biomass ranges
# cover.lm.log.25 = lm(log(T_AG_BM) ~ T_SG, data=bio[bio$T_AG_BM<25,])
# cover.lm.log.50 = lm(log(T_AG_BM) ~ T_SG, data=bio[bio$T_AG_BM<50,])
# cover.lm.log.75 = lm(log(T_AG_BM) ~ T_SG, data=bio[bio$T_AG_BM<75,])
# PrintLM.log(cover.lm.log.25,10)
# boot.lm.log(cover.lm.log.25,10,100)
# 
# PrintLM.log(cover.lm.log.50,10)
# boot.lm.log(cover.lm.log.50,10,100)
# 
# PrintLM.log(cover.lm.log.75,10)
# boot.lm.log(cover.lm.log.75,10,100)
# 
# PrintLM.log(cover.lm.log,10)

# varied k
for (i in 1:10) {
  print(i)
  PrintLM.log(cover.lm.log,10)
  boot.lm.log(cover.lm.log,10,100)
}

# dominant species model --------------------------------------------------

# # assume dominant species dominates total biomass
# ZmHu.dom.fit = lm(log(ZmHu_AG_BM) ~ T_SG, data=bio[bio$pZmHu>Cth,])
# Ho.dom.fit = lm(log(Ho_AG_BM) ~ T_SG, data=bio[bio$pHo>Cth,])
# Hs.dom.fit = lm(log(Hs_AG_BM) ~ T_SG, data=bio[bio$pHs>Cth,])
# Cr.dom.fit = lm(log(Cr_AG_BM) ~ T_SG, data=bio[bio$pCr>Cth,])
# Si.dom.fit = lm(log(Si_AG_BM) ~ T_SG, data=bio[bio$pSi>Cth,])

# assume total biomass follows "dominant type"
ZmHu.dom.fit = lm(log(T_AG_BM) ~ T_SG, data=bio[bio$pZmHu>Cth,])
Ho.dom.fit = lm(log(T_AG_BM) ~ T_SG, data=bio[bio$pHo>Cth,])
Hs.dom.fit = lm(log(T_AG_BM) ~ T_SG, data=bio[bio$pHs>Cth,])
Cr.dom.fit = lm(log(T_AG_BM) ~ T_SG, data=bio[bio$pCr>Cth,])
Si.dom.fit = lm(log(T_AG_BM) ~ T_SG, data=bio[bio$pSi>Cth,])

PrintLM.log(ZmHu.dom.fit,10)
boot.lm.log(ZmHu.dom.fit,10,100)

PrintLM.log(Ho.dom.fit,10)
boot.lm.log(Ho.dom.fit,10,100)

PrintLM.log(Hs.dom.fit,10)
boot.lm.log(Hs.dom.fit,10,100)

PrintLM.log(Cr.dom.fit,10)
boot.lm.log(Cr.dom.fit,10,100)

PrintLM.log(Si.dom.fit,10)
boot.lm.log(Si.dom.fit,10,100)

# create column of predicted values
strat.pred = numeric(length(bio$T_SG))
for (i in 1:length(bio$T_SG)) {
  if (bio$pZmHu[i]>Cth) {
    strat.pred[i] = coef(ZmHu.dom.fit)[2]*bio$T_SG[i] + coef(ZmHu.dom.fit)[1]
  } else if (bio$pHo[i]>Cth) {
    strat.pred[i] = coef(Ho.dom.fit)[2]*bio$T_SG[i] + coef(Ho.dom.fit)[1]
  } else if (bio$pHs[i]>Cth) {
    strat.pred[i] = coef(Hs.dom.fit)[2]*bio$T_SG[i] + coef(Hs.dom.fit)[1]
  } else if (bio$pCr[i]>Cth) {
    strat.pred[i] = coef(Cr.dom.fit)[2]*bio$T_SG[i] + coef(Cr.dom.fit)[1]
  } else if (bio$pSi[i]>Cth) {
    strat.pred[i] = coef(Si.dom.fit)[2]*bio$T_SG[i] + coef(Si.dom.fit)[1]
  } else {
    strat.pred[i] = coef(cover.lm.log)[2]*bio$T_SG[i] + coef(cover.lm.log)[1]
  }
}

bio$stratModel = exp(1)^strat.pred



# species component model -------------------------------------------------

# create individual species models
ZmHu.fit = lm(log(ZmHu_AG_BM) ~ SGZMHU, data=bio[bio$ZmHu_AG_BM > 0 & bio$SGZMHU > 0,])
Ho.fit = lm(log(Ho_AG_BM) ~ SGHO, data=bio[bio$Ho_AG_BM > 0 & bio$SGHO > 0,])
Hs.fit = lm(log(Hs_AG_BM) ~ SGHS, data=bio[bio$Hs_AG_BM > 0 & bio$SGHS > 0,])
Cr.fit = lm(log(Cr_AG_BM) ~ SGCR, data=bio[bio$Cr_AG_BM > 0 & bio$SGCR > 0,])
Si.fit = lm(log(Si_AG_BM) ~ SGSI, data=bio[bio$Si_AG_BM > 0 & bio$SGSI > 0,])

PrintLM.log(ZmHu.fit,10)
boot.lm.log(ZmHu.fit,10,100)

PrintLM.log(Ho.fit,10)
boot.lm.log(Ho.fit,10,100)

PrintLM.log(Hs.fit,10)
boot.lm.log(Hs.fit,10,100)

PrintLM.log(Cr.fit,10)
boot.lm.log(Cr.fit,10,100)

PrintLM.log(Si.fit,10)
boot.lm.log(Si.fit,10,100)

# create columns of predicted values for individual species, and a sum column
bio$mZmHu = exp(1)^predict.lm(object=ZmHu.fit, newdata=bio)
bio$mZmHu[bio$SGZMHU==0] = 0

bio$mHo = exp(1)^predict.lm(object=Ho.fit, newdata=bio)
bio$mHo[bio$SGHO==0] = 0

bio$mHs = exp(1)^predict.lm(object=Hs.fit, newdata=bio)
bio$mHs[bio$SGHS==0] = 0

bio$mCr = exp(1)^predict.lm(object=Cr.fit, newdata=bio)
bio$mCr[bio$SGCR==0] = 0

bio$mSi = exp(1)^predict.lm(object=Si.fit, newdata=bio)
bio$mSi[bio$SGSI==0] = 0

bio$sumModel = bio$mZmHu + bio$mHo + bio$mHs + bio$mCr + bio$mSi



# RMSE comparison ---------------------------------------------------------
sqrt(mean((bio$linModel[bio$T_AG_BM<25]-bio$T_AG_BM[bio$T_AG_BM<25])^2))
sqrt(mean((bio$linModel[bio$T_AG_BM<50]-bio$T_AG_BM[bio$T_AG_BM<50])^2))
sqrt(mean((bio$linModel[bio$T_AG_BM<75]-bio$T_AG_BM[bio$T_AG_BM<75])^2))
sqrt(mean((bio$linModel-bio$T_AG_BM)^2))

sqrt(mean((bio$stratModel[bio$T_AG_BM<25]-bio$T_AG_BM[bio$T_AG_BM<25])^2))
sqrt(mean((bio$stratModel[bio$T_AG_BM<50]-bio$T_AG_BM[bio$T_AG_BM<50])^2))
sqrt(mean((bio$stratModel[bio$T_AG_BM<75]-bio$T_AG_BM[bio$T_AG_BM<75])^2))
sqrt(mean((bio$stratModel-bio$T_AG_BM)^2))

sqrt(mean((bio$sumModel[bio$T_AG_BM<25]-bio$T_AG_BM[bio$T_AG_BM<25])^2))
sqrt(mean((bio$sumModel[bio$T_AG_BM<50]-bio$T_AG_BM[bio$T_AG_BM<50])^2))
sqrt(mean((bio$sumModel[bio$T_AG_BM<75]-bio$T_AG_BM[bio$T_AG_BM<75])^2))
sqrt(mean((bio$sumModel-bio$T_AG_BM)^2))



# application to photo transects ------------------------------------------

# load photo processed photo transect data
photo.sets = c("j2004", "a2007", "j2011", "f2012", "j2012", "f2013", "m2013")

for (i in 1:length(photo.sets)){
  photo.set = read.csv(paste0("photos_",photo.sets[i],".csv"), header=T, stringsAsFactors=FALSE)
  assign(paste0("photos.",photo.sets[i]), photo.set)
  rm(photo.set)
}

# explore the suitable range for biomass estimation, create max percentage cover threshold for estimation
for (i in c(2:7,16)) {
  print(names(bio)[i])
  print(max(bio[,i]))
}
# max %cover thresholds
max.Hu = 70
max.Ho = 30
max.Hs = 70

# apply component model and estimate biomass for each photo in the photo data sets
for (i in 1:length(photo.sets)){
  name = paste0("photos.",photo.sets[i])
  newPhoto = applyComponentModel(get(paste0("photos.",photo.sets[i])))
  assign(name, newPhoto)
  rm(name, newPhoto)
  write.csv(get(paste0("photos.",photo.sets[i])), 
            file=paste0("photos_",photo.sets[i],"_bioEstimates",".csv"), row.names=F)
}


# application to maps -----------------------------------------------------

# load photo processed photo transect data
map.sets = c("j2004", "a2007", "2008", "d2009", "j2011", "f2012", "j2012", "f2013", "m2013")
#map.sets = "m2013"

for (i in 1:length(map.sets)){
  map.set = read.csv(paste0("map_",map.sets[i],".csv"), header=T, stringsAsFactors=FALSE)
  assign(paste0("map.",map.sets[i]), map.set)
  rm(map.set)
}

# define max cover threshold for prediction
max.Hu = 70
max.Ho = 30
max.Hs = 70

# apply dominant model function to map data
for (i in 1:length(map.sets)){
  name = paste0("map.",map.sets[i])
  newmap = applyDominantModel(get(paste0("map.",map.sets[i])))
  assign(name, newmap)
  rm(name, newmap)
  write.csv(get(paste0("map.",map.sets[i])), 
            file=paste0("map_",map.sets[i],"_bioEstimates",".csv"), row.names=F)
}
