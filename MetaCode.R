##########################################################
#### section 1 - prepare environment ###
########################################

#clear environment
cat("\014")
rm(list = ls())
graphics.off()

## sets warnings
options(warn=1)

#setwd
setwd("DATA FILE PATHWAY HERE - DATA USED HERE 'MetaData'")
getwd()

## loads all the required packages
list.of.packages <- c("ggplot2","dplyr","meta","metafor","esc","MuMIn","glmulti","Hmisc")
lapply(list.of.packages, require, character.only = TRUE)

## install dmetar
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("MathiasHarrer/dmetar")

## prints the session info for diagnostics
print(sessionInfo())
rm(list.of.packages)

########################################
#### section 2 - load and check data ###
########################################

## Data version 9 ("NiSpp_v9")
## This includes data from 14 studies, 34 treatments
## all have a control (0N) and a treatment (XN)
df1 <- read.csv("Data_for_ms.csv", sep=",")

### change data types
df1$author<-as.factor(df1$author)
df1$es.id<-as.factor(df1$es.id)
df1$study_Nr<-as.factor(df1$study_Nr)
df1$application_cat<-as.factor(df1$application_cat)
#glimpse(df1)



#####################################################################
#### section 3 - calculate mean difference and sampling variances ###
#####################################################################

df1<-escalc(measure="MD", 
            m2i=c.mean_sr,
            m1i=tr.mean_sr,
            sd2i=c.sd,
            sd1i = tr.sd,
            n2i = c.n,
            n1i = tr.n,
            data=df1,
            append = TRUE)


####################################
#### section 4 - initial exploration
####################################

# Check all model variables for correlations

cor_data <- df1[, c(5:14)]
res<-rcorr(as.matrix(cor_data))
res



####################################################################
#### section 5 - fit full model and check for 2- vs 3- level fit ###
####################################################################
full.model <- rma.mv(yi = yi, 
                 V = vi, 
                 slab = author,
                 data = df1,
                 mods = 
                   ~ tr.nitrogen
                   +I(tr.nitrogen^2)
                 +defoliation_N
                 +c.mean_sr
                 +plot_size_m2
                 +experiment_years
                 +application_cat
                 +tr.nitrogen:defoliation_N
                 ,
                 random = ~ 1 | author/es.id, 
                 test = "t", 
                 method = "ML")

### See if a 2-level model is a better (or as good) a fit
l3.removed<-rma.mv(yi = yi, 
                      V = vi, 
                      slab = author,
                      data = df1,
                      mods = ~ tr.nitrogen
                      +defoliation_N
                      +c.mean_sr
                      +plot_size_m2
                      +experiment_years
                      +application_cat
                   +tr.nitrogen:defoliation_N
                      ,
                      random = ~ 1 | author/es.id, 
                      test = "t", 
                      method = "ML",
                      sigma2 = c(0,NA))
### Compare this with a full (3 level) model
anova(full.model,l3.removed)
## This shows the full (3-level) model is a significantly better fit

### Look at heterogeneity variance captured by each level in model


library(dmetar)
i2<-var.comp(model1)
summary(i2)
plot(i2)

##################################################
#### section 6 - model averaging using glmulti ###
##################################################

## Followed approach of Wolfgang Viechtbauer here:
## https://www.metafor-project.org/doku.php/

library(glmulti)

## Write function, specific for our rma.mv full model
rma.glmulti <- function(formula, data, ...)
  rma.mv(formula, vi, random= ~ 1|author/es.id, data=df1, method="ML", ...)

## Fit all possible models
res<-glmulti(yi~
               tr.nitrogen
             +I(tr.nitrogen^2)
             +defoliation_N
             +c.mean_sr
             +plot_size_m2
             +experiment_years
             +application_cat
             +tr.nitrogen:defoliation_N,
             data=df1,
             level=1,
             fitfunction = rma.glmulti,
             crit="aicc",
             confsetsize = 64)

## Summary of results
print(res)

## Inspect these top 3 models
top <- weightable(res)
top <- top[top$aicc <=min(top$aicc)+2,]
top

## Inspect the "best" model
summary(res@objects[[1]])

## Look at relative importance of predictors (<0.5 considered unimportant)
plot(res, type="s")

## Mulitmodel inference
eval(metafor:::.glmulti)
## Tidied output
mmi <- as.data.frame(coef(res, varweighting="Johnson"))
mmi <- data.frame(Estimate=mmi$Est, SE=sqrt(mmi$Uncond), 
                  Importance=mmi$Importance, row.names=row.names(mmi))
mmi$z <- mmi$Estimate / mmi$SE
mmi$p <- 2*pnorm(abs(mmi$z), lower.tail=FALSE)
names(mmi) <- c("Estimate", "Std. Error", "Importance", "z value", "Pr(>|z|)")
mmi$ci.lb <- mmi[[1]] - qnorm(.975) * mmi[[2]]
mmi$ci.ub <- mmi[[1]] + qnorm(.975) * mmi[[2]]
mmi <- mmi[order(mmi$Importance, decreasing=TRUE), c(1,2,4:7,3)]
round(mmi, 4)

## "Best" model includes 3 variables (all included in the top 3 candidate models too)
## Here we fit a reduced model including these variables

reduced.linear.model<- rma.mv(yi = yi, 
                         V = vi, 
                         slab = author,
                         data = df1,
                         mods = 
                           ~ tr.nitrogen
                         +defoliation_N
                         +c.mean_sr,
                         random = ~ 1 | author/es.id, 
                         test = "t", 
                         method = "ML")
reduced.nonlinear.model<- rma.mv(yi = yi, 
                              V = vi, 
                              slab = author,
                              data = df1,
                              mods = 
                                ~ tr.nitrogen
                              +I(tr.nitrogen^2)
                              +defoliation_N
                              +c.mean_sr,
                              random = ~ 1 | author/es.id, 
                              test = "t", 
                              method = "ML")
anova(reduced.linear.model,reduced.nonlinear.model)

##Final model 

final.model<- rma.mv(yi = yi, 
                              V = vi, 
                              slab = author,
                              data = df1,
                              mods = 
                                ~ tr.nitrogen
                              +defoliation_N
                              +c.mean_sr,
                              random = ~ 1 | author/es.id, 
                              test = "t", 
                              method = "REML")
summary(final.model)


## Exclude the Sammul et al study to see if the baseline 
## SR result holds without this (return to line 50 and repeat after exclusion)
#df1<-filter(df1, study_Nr != 210)


###ENDS###





