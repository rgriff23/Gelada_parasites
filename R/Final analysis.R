#################################################################################################################
# PREPARATIONS
#################################################################################################################

# Load packages
library(plyr) # For ddply
library(Deducer) # For G-test
library(survival) # For Cox model
library(simPH) # For plotting Cox model 
library(MASS) # For negative binomial glm
library(lmtest) # For likelihood ratio test

# Read data into R
lh_wide <- read.csv("~/Desktop/GitHub/Gelada_parasites/Data/lh_wide.csv", header=TRUE, stringsAsFactors=FALSE)
lh_long <- read.csv("~/Desktop/GitHub/Gelada_parasites/Data/lh_long.csv", header=TRUE, stringsAsFactors=FALSE)
lh_ibi <- read.csv("~/Desktop/GitHub/Gelada_parasites/Data/lh_ibi.csv", header=TRUE, stringsAsFactors=FALSE) 
samps <- read.csv("~/Desktop/GitHub/Gelada_parasites/Data/samps.csv", header=TRUE, stringsAsFactors=FALSE)

# Format dates
lh_wide$dob <- as.POSIXct(lh_wide$dob, format="%Y-%m-%d")
lh_wide$cyst <- as.POSIXct(lh_wide$cyst, format="%Y-%m-%d")
lh_wide$dod <- as.POSIXct(lh_wide$dod, format="%Y-%m-%d")
lh_wide$end <- as.POSIXct(lh_wide$end, format="%Y-%m-%d")

#################################################################################################################
# PERIOD PREVALENCE AND G-TESTS FOR LIFE HISTORY DATA (use lh_wide)
#################################################################################################################

############
## CLEAN  #
##########

# How many individuals to start?
nrow(lh_wide) # 621

# Exclude individuals who disappeared 
sum(lh_wide$dis) # 37 disappeared
lh_wide <- lh_wide[!lh_wide$dis==1,-which(names(lh_wide)=="dis")] # 584 remaining

# Exclude individuals who END before the age of four from data
drop <- lh_wide$name[(c(lh_wide$end - lh_wide$dob)/60/60/24/365) < 4]
length(drop) # 233 <4 years old
sum(!is.na(lh_wide[lh_wide$name %in% drop,"cyst"])) # 0 cysts among the dropped individuals
lh_wide <- lh_wide[!(lh_wide$name %in% drop),] # Drop 233 of 584 individuals, leaving 351 adults

###############
## SUMMARIZE #
#############

# Number of adults, males, and females in the long term data
nrow(lh_wide) # 351 individuals surviving past age 4 by the end of the study
sum(lh_wide$sex=="M") # 133 males
sum(lh_wide$sex=="F") # 218 females

# Number/proportion of adults, males, and females with cysts
tab0 <- table(lh_wide$sex, !is.na(lh_wide$cyst))
sum(tab0[,2])/sum(tab0)	# Overall period prevalence = 0.15
sum(tab0["M",2])/sum(tab0["M",])	# Male period prevalence = 0.09
sum(tab0["F",2])/sum(tab0["F",])	# Female period prevalence = 0.18

# Number/proportion of adults dying with/without cysts over study period
tab1 <- table(!is.na(lh_wide$dod), !is.na(lh_wide$cyst))
tab1[2,]/colSums(tab1)	# 15.7% without cysts die; 61.5% with cysts die
tab2 <- table(!is.na(lh_wide$dod[lh_wide$sex=="M"]), !is.na(lh_wide$cyst[lh_wide$sex=="M"]))
tab2[2,]/colSums(tab2)	# 7.4% without cysts die; 25% with cysts die
tab3 <- table(!is.na(lh_wide$dod[lh_wide$sex=="F"]), !is.na(lh_wide$cyst[lh_wide$sex=="F"]))
tab3[2,]/colSums(tab3) # 21.3% without cysts die; 72.5% with cysts die

#############
## ANALYZE #
###########

# G-test for independence between cysts and sex
likelihood.test(tab0) # Significant association (p < 0.05)

# G-test for difference between Guassa and SMNP sex-specific prevalences
likelihood.test(rbind(tab0["F",], c(68, 31))) # Significant association (p < 0.05)
likelihood.test(rbind(tab0["M",], c(49, 19))) # Highly significant association (p < 0.001)

# G-test for association between cysts and deaths over study period
likelihood.test(tab1) # Highly significant across combined sexes (p < 0.001)
likelihood.test(tab2) # Approaching significance in males (p = 0.08)
likelihood.test(tab3) # Highly significant in females  (p < 0.001)

###############
## VISUALIZE #
#############

# Barplot for cyst prevalence over course of study
bplot1 <- barplot(rev(tab0[,2]/rowSums(tab0)), col=c("dimgray", "lightgray"), names.arg=c("Male (n = 133)", "Female (n = 218)"), ylim=c(0, 1), ylab="% adults with cysts during study")
errs.m <- prop.test(tab0["M",2], rowSums(tab0)["M"])
errs.f <- prop.test(tab0["F",2], rowSums(tab0)["F"])
arrows(x0=c(bplot1), y0=c(errs.m$conf.int[1], errs.f$conf.int[1]), x1=c(bplot1), y1=c(errs.m$conf.int[2], errs.f$conf.int[2]), length=0.07, angle=90, code=3)


# Barplots for mortality among males/females with/without cysts
layout(matrix(1:2, 1, 2))
bplot2 <- barplot(rev(tab3[2,]/colSums(tab3)), col=c("dimgray", "lightgray"), ylim=c(0,1), xlab="Females", ylab="% adults who died during study", names.arg=c("Cyst (29/40)", "No cyst (38/178)"))
errs.c <- prop.test(tab3[2,2], colSums(tab3)[2])
errs.n <- prop.test(tab3[2,1], colSums(tab3)[1])
arrows(x0=c(bplot2), y0=c(errs.c$conf.int[1], errs.n$conf.int[1]), x1=c(bplot2), y1=c(errs.c$conf.int[2], errs.n$conf.int[2]), length=0.07, angle=90, code=3)
bplot3 <- barplot(rev(tab2[2,]/colSums(tab2)), col=c("dimgray", "lightgray"), ylim=c(0,1), xlab="Males", names.arg=c("Cyst (3/12)", "No cyst (9/121)"))
errs.c2 <- prop.test(tab2[2,2], colSums(tab2)[2])
errs.n2 <- prop.test(tab2[2,1], colSums(tab2)[1])
arrows(x0=c(bplot3), y0=c(errs.c2$conf.int[1], errs.n2$conf.int[1]), x1=c(bplot3), y1=c(errs.c2$conf.int[2], errs.n2$conf.int[2]), length=0.07, angle=90, code=3)

#################################################################################################################
# COX PROPORTIONAL HAZARDS MODEL FOR LIFE HISTORY DATA (use lh_long)
#################################################################################################################

#############
## ANALYZE #
###########

# Males and females
mod0 <- coxph(Surv(start, stop, death) ~ sex + cyst, data=lh_long)
summary(mod0)
cox.zph(mod0) # Proportional hazards assumption not met for cysts
plot(cox.zph(mod0))
# Add interaction between cyst and time, cyst_log
lh_long1 <- tvc(lh_long, b="cyst", tvar="stop", tfun="log")
mod1 <- coxph(Surv(start, stop, death) ~ sex + cyst + cyst_log, data=lh_long1)
summary(mod1)
cox.zph(mod1) # Check
simGG(coxsimtvc(obj=mod1, b="cyst", btvc="cyst_log", tfun="log", qi="Relative Hazard", Xj=c(1), from=1, to=300, by=1, nsim=100))

# Males only
mod.m <- coxph(Surv(start, stop, death) ~ cyst, data=lh_long[lh_long$sex=="M",])
summary(mod.m)
cox.zph(mod.m)  # Proportional hazards assumption nearly not met for cysts
plot(cox.zph(mod.m))
mod.m1 <- coxph(Surv(start, stop, death) ~ cyst + cyst_log, data=lh_long1[lh_long1$sex=="M",])
summary(mod.m1)
cox.zph(mod.m1) # Check
simGG(coxsimtvc(obj=mod.m1, b="cyst", btvc="cyst_log", tfun="log", qi="Relative Hazard", Xj=c(1), from=1, to=200, by=1, nsim=100))


# Females only
mod.f <- coxph(Surv(start, stop, death) ~ cyst, data=lh_long[lh_long$sex=="F",])
summary(mod.f)
cox.zph(mod.f)   # Proportional hazards assumption not met for cysts
plot(cox.zph(mod.f))
mod.f1 <- coxph(Surv(start, stop, death) ~ cyst + cyst_log, data=lh_long1[lh_long1$sex=="F",])
summary(mod.f1)
cox.zph(mod.f1)	# Check
simGG(coxsimtvc(obj=mod.f1, b="cyst", btvc="cyst_log", tfun="log", qi="Hazard Ratio", Xj=c(1), from=1, to=300, by=1, nsim=100))

###############
## VISUALIZE #
#############

# Separate curves for cyst/no cyst
cyst <- data.frame(cyst=c(0, 1))

# Survival curves for males and females
layout(matrix(1:2, 1, 2))
plot(survfit(mod.m1, newdata=cyst), col=c("dimgray", "lightgray"), xlab="Months", ylab="Proportion alive", main="Males")
legend("bottomleft", legend=c("No cyst", "Cyst"), fill=c("dimgray", "lightgray"))
plot(survfit(mod.f1, newdata=cyst), col=c("dimgray", "lightgray"), xlab="Months", ylab="", main="Females")

#################################################################################################################
# SUMMARIZING REPRODUCTION DATA (use lh_ibi)
#################################################################################################################

#############
## FEMALES #
###########

# Females only
lh_ibi.f <- lh_ibi[lh_ibi$SEX=="F",]

# Ignore "first birth" designation
lh_ibi.f[lh_ibi.f$EVENT %in% c("First Birth", "First Birth?"),"EVENT"] <- "Birth"

# Trim data to mothers
mothers <- ddply(lh_ibi.f, .(NAME), function(x) {sum(x$EVENT=="Birth")})
lh_ibi.f <- lh_ibi.f[lh_ibi.f $NAME %in% mothers$NAME,]

# Convert time data to POSIXct
lh_ibi.f$DATE <- as.POSIXct(lh_ibi.f$DATE, format="%Y-%m-%d")

# Look at life history of all 40 females with cysts who ever gave birth
t <- ddply(lh_ibi.f, .(NAME), function(x) {sum(x$EVENT=="Cyst")})
t <- t[t[,2]==1,]
length(unique(lh_ibi.f[lh_ibi.f $NAME %in% t$NAME,3:7]$NAME)) # 40 moms
lh_ibi.f[lh_ibi.f $NAME %in% t$NAME,3:7]

###########
## MALES #
#########

# Males only
lh_ibi.m <- lh_ibi[lh_ibi$SEX=="M",]

# Convert time data to POSIXct
lh_ibi.m$DATE <- as.POSIXct(lh_ibi.m$DATE, format="%Y-%m-%d")

# Look at life hitory of all males who ever had cysts
t2 <- ddply(lh_ibi.m, .(NAME), function(x) {sum(x$EVENT=="Cyst")})
t2 <- t2[t2[,2]==1,]
length(unique(lh_ibi.m[lh_ibi.m$NAME %in% t2$NAME,3:7]$NAME)) # 18 cyst males
lh_ibi.m[lh_ibi.m$NAME %in% t2$NAME,3:8]

#################################################################################################################
# G-TESTS AND GLMs FOR SAMPLE PREVALENCE (use samps)
#################################################################################################################


#################################################################################################################
# END
#################################################################################################################