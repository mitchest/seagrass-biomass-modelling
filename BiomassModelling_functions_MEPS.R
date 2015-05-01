## Author: Mitchell Lyons (mitchell.lyons@gmail.com)
## No licence, though please attribute me and publish any derivitives open source =)

# function definitions ----------------------------------------------------


# # function for printing R2 and RMSE values for model fits
# printFitCV = function(model.fit, kfold.fit){
#   print(paste(length(kfold.fit$Predicted)," samples; ","rmse in gDW.m^-2:",sep=""))
#   print(paste("overall rmse = ",exp(1)^summary(model.fit)$sigma,sep=""))
#   print(paste("overall R^2 = ",summary(model.fit)$r.squared,sep=""))
#   print(paste("k-fold rmse = ",exp(1)^sqrt(attr(kfold.fit, "ms"))),sep="")
#   print(paste("predicted - CV-predicted rmse = ",
#               abs((exp(1)^summary(model.fit)$sigma)-(exp(1)^sqrt(attr(kfold.fit, "ms")))),sep=""))
#   print(paste("CV predicted Pearson's corr.  = ",cor.test(kfold.fit$cvpred, kfold.fit$Predicted)$estimate))
# }

## function for a random K-fold CV on lm model fit
cv.lm = function(model.fit, K) {
  # get data used to fit overall model
  data = model.fit$model
  # pre-allocate vector for sotring results
  model.r2 = numeric(K)
  model.rmse = numeric(K)
  pred.rmse = numeric(K)
  # randomly shuffle data and create fold index
  data = data[sample(nrow(data)),]
  folds = cut(seq(1,nrow(data)),breaks=K,labels=FALSE)
  # split data and calculate desired stats
  for (i in 1:K) {
    # make test/model data
    test = which(folds==i,arr.ind=TRUE)
    test.data = data[test,]
    model.data = data[-test,]
    # fit models and store stats
    fold.fit = lm(formula(model.fit), data=model.data)
    model.r2[i] = summary(fold.fit)$r.squared
    model.rmse[i] = sqrt(mean(residuals.lm(fold.fit, type="response")^2))
    # test predicted values against test set
    predicted = predict.lm(fold.fit, newdata=test.data, type="response")
    pred.rmse[i] = sqrt(mean((predicted-test.data[,1])^2))
  }
  return(data.frame(model.r2,model.rmse,pred.rmse))
}

## function for a random K-fold CV on lm model fit - hacked to have log(y)
cv.lm.log = function(model.fit, K) {
  # get data used to fit overall model
  data = model.fit$model
  # hack y back to orginal
  data[,1] = exp(1)^data[,1]
  names(data)[1] = gsub(".*\\((.*)\\).*", "\\1", names(data)[1])
  # pre-allocate vector for sotring results
  model.r2 = numeric(K)
  model.rmse = numeric(K)
  pred.rmse = numeric(K)
  # randomly shuffle data and create fold index
  data = data[sample(nrow(data)),]
  folds = cut(seq(1,nrow(data)),breaks=K,labels=FALSE)
  # split data and calculate desired stats
  for (i in 1:K) {
    # make test/model data
    test = which(folds==i,arr.ind=TRUE)
    test.data = data[test,]
    model.data = data[-test,]
    # fit models and store stats
    fold.fit = lm(formula(model.fit), data=model.data)
    model.r2[i] = summary(fold.fit)$r.squared
    model.resid = exp(1)^predict.lm(fold.fit)-model.data
    model.rmse[i] = sqrt(mean(model.resid^2))
    # test predicted values against test set
    predicted = exp(1)^predict.lm(fold.fit, newdata=test.data)
    pred.rmse[i] = sqrt(mean((predicted-test.data[,1])^2))
  }
  return(data.frame(model.r2,model.rmse,pred.rmse))
}

## function for a random K-fold CV on glm model fit
cv.glm = function(model.fit, K) {
  # get data used to fit overall model
  data = model.fit$model
  # pre-allocate vector for sotring results
  model.rmse = numeric(K)
  pred.rmse = numeric(K)
  # randomly shuffle data and create fold index
  data = data[sample(nrow(data)),]
  folds = cut(seq(1,nrow(data)),breaks=K,labels=FALSE)
  # split data and calculate desired stats
  for (i in 1:K) {
    # make test/model data
    test = which(folds==i,arr.ind=TRUE)
    test.data = data[test,]
    model.data = data[-test,]
    # fit models and store stats
    fold.fit = glm(formula(model.fit), data=model.data, family=family(model.fit)$family)
    model.rmse[i] = sqrt(mean(residuals.glm(fold.fit, type="response")^2))
    # test predicted values against test set
    predicted = predict.glm(fold.fit, newdata=test.data, type="response")
    pred.rmse[i] = sqrt(mean((predicted-test.data[,1])^2))
  }
  return(data.frame(model.rmse,pred.rmse))
}

## function for a random K-fold CV on gam model fit
cv.gam = function(model.fit, K) {
  # get data used to fit overall model
  data = model.fit$model
  # pre-allocate vector for sotring results
  model.rmse = numeric(K)
  pred.rmse = numeric(K)
  # randomly shuffle data and create fold index
  data = data[sample(nrow(data)),]
  folds = cut(seq(1,nrow(data)),breaks=K,labels=FALSE)
  # split data and calculate desired stats
  for (i in 1:K) {
    # make test/model data
    test = which(folds==i,arr.ind=TRUE)
    test.data = data[test,]
    model.data = data[-test,]
    # fit models and store stats
    fold.fit = gam(formula(model.fit), data=model.data, family=family(model.fit)$family)
    model.rmse[i] = sqrt(mean(residuals.gam(fold.fit, type="response")^2))
    # test predicted values against test set
    predicted = predict.gam(fold.fit, newdata=test.data, type="response")
    pred.rmse[i] = sqrt(mean((predicted-test.data[,1])^2))
  }
  return(data.frame(model.rmse,pred.rmse))
}

## function for printing out lm and cv.lm statistics (need to specify K for k-fold CV)
PrintLM = function(model.fit, K){
  print(paste0("LM fit: ",Reduce(paste, deparse(formula(model.fit)))))
  print(paste0(nrow(model.fit$model)," samples; ","rmse in mg.m^-2:"))
  print(paste0("overall R^2 = ",summary(model.fit)$r.squared))
  print(paste0("overall rmse = ",sqrt(mean(residuals.lm(model.fit, type="response")^2))))
  cv = cv.lm(model.fit, K)
  print(paste0("mean k-fold model rmse = ", mean(cv$model.rmse)))
  print(paste0("mean k-fold prediction rmse = ", mean(cv$pred.rmse)))
  # make an observed vs. fitted plot
  plot = qplot(y=predict(model.fit,type="response"),x=model.fit$model[,1])
  plot = plot + theme_classic() +
    ggtitle(paste0("LM fit: ",Reduce(paste, deparse(formula(model.fit))))) +
    ylab("Fitted biomass") + xlab("Observed biomass") +
    geom_smooth(method=lm, alpha=0.2) +
    geom_abline(intercept=0,slope=1, colour="red", linetype="twodash") +
    coord_cartesian(ylim=c(0,200), xlim=c(0,150))
  return(plot)
}

## function for printing out lm and cv.lm statistics (need to specify K for k-fold CV) - with the log(y) cv.lm
PrintLM.log = function(model.fit, K){
  print(paste0("LM fit: ",Reduce(paste, deparse(formula(model.fit)))))
  print(paste0(nrow(model.fit$model)," samples; ","rmse in mg.m^-2:"))
  print(paste0("overall R^2 = ",summary(model.fit)$r.squared))
  residuals = exp(1)^predict(model.fit,type="response")-exp(1)^model.fit$model[,1]
  print(paste0("overall rmse = ",sqrt(mean(residuals^2))))
  cv = cv.lm.log(model.fit, K)
  print(paste0("mean k-fold model rmse = ", mean(cv$model.rmse)))
  print(paste0("mean k-fold prediction rmse = ", mean(cv$pred.rmse)))
  # make an observed vs. fitted plot
  plot = qplot(y=exp(1)^predict(model.fit,type="response"),x=exp(1)^model.fit$model[,1])
  plot = plot + theme_classic() +
    ggtitle(paste0("LM fit: ",Reduce(paste, deparse(formula(model.fit))))) +
    ylab("Fitted biomass") + xlab("Observed biomass") +
    geom_smooth(method=lm, alpha=0.2) +
    geom_abline(intercept=0,slope=1, colour="red", linetype="twodash") +
    coord_cartesian(ylim=c(0,200), xlim=c(0,150))
  return(plot)
}

## function for printing out glm stats (and cv.glm stats eventually...)
PrintGLM = function(model.fit, K) {
  print(paste0("GLM fit: ",Reduce(paste, deparse(formula(model.fit)))))
  print(paste0(nrow(model.fit$model)," samples; ","rmse in mg.m^-2:"))
  print(paste0("overall rmse = ",sqrt(mean(residuals.glm(model.fit, type="response")^2))))
  cv = cv.glm(model.fit, K)
  print(paste0("mean k-fold model rmse = ", mean(cv$model.rmse)))
  print(paste0("mean k-fold prediction rmse = ", mean(cv$pred.rmse)))
  # make an observed vs. fitted plot
  plot = qplot(y=predict(model.fit,type="response"),x=model.fit$model[,1])
  plot = plot + theme_classic() +
    ggtitle(paste0("GLM fit: ",Reduce(paste, deparse(formula(model.fit))))) +
    ylab("Fitted biomass") + xlab("Observed biomass") +
    geom_smooth(method=lm, alpha=0.2) +
    geom_abline(intercept=0,slope=1, colour="red", linetype="twodash") +
    coord_cartesian(ylim=c(0,200), xlim=c(0,150))
  return(plot)
}

## function for printing out gam stats (and cv.glm stats eventually...)
PrintGAM = function(model.fit, K) {
  print(paste0("GAM fit: ",Reduce(paste, deparse(formula(model.fit)))))
  print(paste0(nrow(model.fit$model)," samples; ","rmse in mg.m^-2:"))
  print(paste0("GAM deviance explained = ",summary(model.fit)$dev.expl))
  print(paste0("GAM adj. R^2 = ",summary(model.fit)$r.sq))
  print(paste0("overall rmse = ",sqrt(mean(residuals.gam(model.fit, type="response")^2))))
  cv = cv.gam(model.fit, K)
  print(paste0("mean k-fold model rmse = ", mean(cv$model.rmse)))
  print(paste0("mean k-fold prediction rmse = ", mean(cv$pred.rmse)))
  # make an observed vs. fitted plot
  plot = qplot(y=predict(model.fit,type="response"),x=model.fit$model[,1])
  plot = plot + theme_classic() +
    ggtitle(paste0("GAM fit: ",Reduce(paste, deparse(formula(model.fit))))) +
    ylab("Fitted biomass") + xlab("Observed biomass") +
    geom_smooth(method=lm, alpha=0.2) +
    geom_abline(intercept=0,slope=1, colour="red", linetype="twodash") +
    coord_cartesian(ylim=c(0,200), xlim=c(0,150))
  return(plot)
}

## function that calculates bootstrapped mean and 95% CI's for CV prediciton rmse for a given LM
boot.lm = function(model.fit, K, nBoot) {
  pred.rmse = numeric(K*nBoot)
  pos = 1
  for (i in 1:nBoot) {
    pred.rmse[pos:(pos+(K-1))] = cv.lm(model.fit,K)$pred.rmse
    pos = pos+K
  }
  print(paste0("mean prediction RMSE: ",mean(pred.rmse)))
  plot = plot(sort(pred.rmse),
              main=paste0("LM fit: ",Reduce(paste, deparse(formula(model.fit)))),
              ylab="bootstrapped prediction RMSE") +
    abline(h=mean(pred.rmse), col="red") +
    abline(h=quantile(pred.rmse,probs=0.95), lty=2)
  return(plot)
}

## function that calculates bootstrapped mean and 95% CI's for CV prediciton rmse for a given LM
boot.lm.log = function(model.fit, K, nBoot) {
  pred.rmse = numeric(K*nBoot)
  pos = 1
  for (i in 1:nBoot) {
    pred.rmse[pos:(pos+(K-1))] = cv.lm.log(model.fit,K)$pred.rmse
    pos = pos+K
  }
  print(paste0("mean prediction RMSE: ",mean(pred.rmse)))
  plot = plot(sort(pred.rmse),
              main=paste0("LM fit: ",Reduce(paste, deparse(formula(model.fit)))),
              ylab="bootstrapped prediction RMSE") +
    abline(h=mean(pred.rmse), col="red") +
    abline(h=quantile(pred.rmse,probs=0.95), lty=2)
  return(plot)
}

## function that calculates bootstrapped mean and 95% CI's for CV prediciton rmse for a given GLM
boot.glm = function(model.fit, K, nBoot) {
  pred.rmse = numeric(K*nBoot)
  pos = 1
  for (i in 1:nBoot) {
    pred.rmse[pos:(pos+(K-1))] = cv.glm(model.fit,K)$pred.rmse
    pos = pos+K
  }
  print(paste0("mean prediction RMSE: ",mean(pred.rmse)))
  plot = plot(sort(pred.rmse),
              main=paste0("GLM fit: ",Reduce(paste, deparse(formula(model.fit)))),
              ylab="bootstrapped prediction RMSE") +
    abline(h=mean(pred.rmse), col="red") +
    abline(h=quantile(pred.rmse,probs=0.95), lty=2)
  return(plot)
}

## function that calculates bootstrapped mean and 95% CI's for CV prediciton rmse for a given GAM
boot.gam = function(model.fit, K, nBoot) {
  pred.rmse = numeric(K*nBoot)
  pos = 1
  for (i in 1:nBoot) {
    pred.rmse[pos:(pos+(K-1))] = cv.gam(model.fit,K)$pred.rmse
    pos = pos+K
  }
  print(paste0("mean prediction RMSE: ",mean(pred.rmse)))
  plot = plot(sort(pred.rmse),
              main=paste0("GAM fit: ",Reduce(paste, deparse(formula(model.fit)))),
              ylab="bootstrapped prediction RMSE") +
    abline(h=mean(pred.rmse), col="red") +
    abline(h=quantile(pred.rmse,probs=0.95), lty=2)
  return(plot)
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