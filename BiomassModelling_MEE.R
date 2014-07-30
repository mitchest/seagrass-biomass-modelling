library(DAAG)
library(stringr)
library(plyr)

# set working directory to data and code location
setwd("C:/my/data/location")


# function definitions ----------------------------------------------------
# these can be moved to a separate source file if desired

## function for printing R2 and RMSE values for model fits
printFitCV = function(model.fit, kfold.fit){
  print(paste(length(kfold.fit$Predicted)," samples; ","rmse in gDW.m^-2:",sep=""))
  print(paste("overall rmse = ",exp(1)^summary(model.fit)$sigma,sep=""))
  print(paste("overall R^2 = ",summary(model.fit)$r.squared,sep=""))
  print(paste("k-fold rmse = ",exp(1)^sqrt(attr(kfold.fit, "ms"))),sep="")
  print(paste("predicted - CV-predicted rmse = ",
              abs((exp(1)^summary(model.fit)$sigma)-(exp(1)^sqrt(attr(kfold.fit, "ms")))),sep=""))
  print(paste("CV predicted Pearson's corr.  = ",cor.test(kfold.fit$cvpred, kfold.fit$Predicted)$estimate))
}

## function to estimate biomass in photo transect files
applyComponentModel = function(photos){
  # predict individual biomass values and sum
  photos$ZmHu_biomass = exp(1)^predict.lm(object=ZmHu.fit, newdata=photos)
  photos$ZmHu_biomass[photos$SGZMHU==0] = 0
  
  photos$Ho_biomass = exp(1)^predict.lm(object=Ho.fit, newdata=photos)
  photos$Ho_biomass[photos$SGHO==0] = 0
  
  photos$Hs_biomass = exp(1)^predict.lm(object=Hs.fit, newdata=photos)
  photos$Hs_biomass[photos$SGHS==0] = 0
  
  photos$Cr_biomass = exp(1)^predict.lm(object=Cr.fit, newdata=photos)
  photos$Cr_biomass[photos$SGCR==0] = 0
  
  photos$Si_biomass = exp(1)^predict.lm(object=Si.fit, newdata=photos)
  photos$Si_biomass[photos$SGSI==0] = 0
  
  photos$sumModel = photos$ZmHu_biomass + photos$Ho_biomass + photos$Hs_biomass + photos$Cr_biomass + photos$Si_biomass
  
  # remove values that are expected to be outside the confident prediction range
  photos$corrSumModel = photos$sumModel
  photos$corrSumModel[photos$SGHU>max.Hu] = 0
  photos$corrSumModel[photos$SGHO>max.Ho] = 0
  photos$corrSumModel[photos$SGHS>max.Hs] = 0
  
  # add binned biomass value for plotting
  photos$biomass.bin = cut(photos$corrSumModel, breaks=c(0,20,40,60,80,100,120,140,1000), 
                           labels=c("1-20","21-40","41-60","61-80","81-100","101-120","121-140", "140+"))
  
  return(photos)
}

## function to estimate biomass from cover/species label in map
applyDominantModel = function(map){
  # check cover/species values are correct
  print("Check cover labels are consistent")
  print(unique(map$Cover))
  print("Check species labels are consistent")
  print(unique(map$DomSpecies))
  # modify "1-10" class label to work with numeric conversion
  map$Cover[map$Cover=="1-10"] = "00-10"
  # create median value for estimation
  map$coverValue = as.numeric(substr(map$Cover,1,2))+5
  map$coverValueNA = ifelse(!is.na(map$coverValue), map$coverValue, 0)
  # conditionally apply dominant species model
  strat.pred = numeric(nrow(map))
  for (i in 1:nrow(map)){
    if (map$DomSpecies[i]=="Zm") {
      strat.pred[i] = coef(ZmHu.dom.fit)[2]*map$coverValue[i] + coef(ZmHu.dom.fit)[1]
    } else if (map$DomSpecies[i]=="Ho") {
      strat.pred[i] = coef(Ho.dom.fit)[2]*map$coverValue[i] + coef(Ho.dom.fit)[1]
      if (map$coverValueNA[i]>max.Ho) strat.pred[i] = coef(Ho.dom.fit)[2]*max.Ho + coef(Ho.dom.fit)[1]
    } else if (map$DomSpecies[i]=="Hs") {
      strat.pred[i] = coef(Hs.dom.fit)[2]*map$coverValue[i] + coef(Hs.dom.fit)[1]
      if (map$coverValueNA[i]>max.Hs) strat.pred[i] = coef(Hs.dom.fit)[2]*max.Hs + coef(Hs.dom.fit)[1]
    } else if (map$DomSpecies[i]=="Cr") {
      strat.pred[i] = coef(Cr.dom.fit)[2]*map$coverValue[i] + coef(Cr.dom.fit)[1]
    } else if (map$DomSpecies[i]=="Si") {
      strat.pred[i] = coef(Si.dom.fit)[2]*map$coverValue[i] + coef(Si.dom.fit)[1]
    } else {
      strat.pred[i] = NA
    } 
  }
  map$biomass = exp(1)^strat.pred
  # ensure area is a numeric value
  map$Area = as.numeric(map$Area)
  # add binned biomass value for plotting
  map$biomass.bin = cut(map$biomass, breaks=c(0,20,40,60,80,100,120,140,1000), 
                        labels=c("1-20","21-40","41-60","61-80","81-100","101-120","121-140", "140+"))
  
  return(map)
}


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


# linear regression for total bio ~ total cover ---------------------------

# one linear term
cover.fit1 = lm(log(T_AG_BM) ~ T_SG, data=bio)
cover.fit1.cv = CVlm(df=bio, cover.fit1, m=10)
printFitCV(cover.fit1, cover.fit1.cv)

# 2nd order term
cover.fit2 = lm(log(T_AG_BM) ~ T_SG + I(T_SG^2), data=bio)
cover.fit2.cv = CVlm(df=bio, cover.fit2, m=10)
printFitCV(cover.fit2, cover.fit2.cv)

# 3rd order term
cover.fit3 = lm(log(T_AG_BM) ~ T_SG + I(T_SG^2) + I(T_SG^3), data=bio)
cover.fit3.cv = CVlm(df=bio, cover.fit3, m=10)
printFitCV(cover.fit3, cover.fit3.cv)

# write predictions to bio data frame
bio$linModel = exp(1)^predict(cover.fit1, bio)
bio$linModelCV = exp(1)^cover.fit1.cv$cvpred
bio$poly2Model = exp(1)^predict(cover.fit2, bio)
bio$poly3Model = exp(1)^predict(cover.fit3, bio)


# dominant species model --------------------------------------------------

# fit dominant species models
ZmHu.dom.fit = lm(log(T_AG_BM) ~ T_SG, data=bio[bio$pZmHu>Cth,])
Ho.dom.fit = lm(log(T_AG_BM) ~ T_SG, data=bio[bio$pHo>Cth,])
Hs.dom.fit = lm(log(T_AG_BM) ~ T_SG, data=bio[bio$pHs>Cth,])
Cr.dom.fit = lm(log(T_AG_BM) ~ T_SG, data=bio[bio$pCr>Cth,])
Si.dom.fit = lm(log(T_AG_BM) ~ T_SG, data=bio[bio$pSi>Cth,])

# cross validated models
ZmHu.dom.fit.cv = CVlm(df=bio[bio$pZmHu>Cth,], ZmHu.dom.fit, m=10)
Ho.dom.fit.cv = CVlm(df=bio[bio$pHo>Cth,], Ho.dom.fit, m=10)
Hs.dom.fit.cv = CVlm(df=bio[bio$pHs>Cth,], Hs.dom.fit, m=10)
Cr.dom.fit.cv = CVlm(df=bio[bio$pCr>Cth,], Cr.dom.fit, m=10)
Si.dom.fit.cv = CVlm(df=bio[bio$pSi>Cth,], Si.dom.fit, m=10)

# check out the fit and CV model fit stats
printFitCV(ZmHu.dom.fit, ZmHu.dom.fit.cv)
printFitCV(Ho.dom.fit, Ho.dom.fit.cv)
printFitCV(Hs.dom.fit, Hs.dom.fit.cv)
printFitCV(Cr.dom.fit, Cr.dom.fit.cv)
printFitCV(Si.dom.fit, Si.dom.fit.cv)

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
    strat.pred[i] = coef(cover.fit1)[2]*bio$T_SG[i] + coef(cover.fit1)[1]
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

# cross validated models
ZmHu.fit.cv = CVlm(df=bio[bio$ZmHu_AG_BM > 0 & bio$SGZMHU > 0,], ZmHu.fit, m=10)
Ho.fit.cv = CVlm(df=bio[bio$Ho_AG_BM > 0 & bio$SGHO > 0,], Ho.fit, m=10)
Hs.fit.cv = CVlm(df=bio[bio$Hs_AG_BM > 0 & bio$SGHS > 0,], Hs.fit, m=10)
Cr.fit.cv = CVlm(df=bio[bio$Cr_AG_BM > 0 & bio$SGCR > 0,], Cr.fit, m=10)
Si.fit.cv = CVlm(df=bio[bio$Si_AG_BM > 0 & bio$SGSI > 0,], Si.fit, m=10)

# check out the fit and CV model fit stats
printFitCV(ZmHu.fit, ZmHu.fit.cv)
printFitCV(Ho.fit, Ho.fit.cv)
printFitCV(Hs.fit, Hs.fit.cv)
printFitCV(Cr.fit, Cr.fit.cv)
printFitCV(Si.fit, Si.fit.cv)


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
