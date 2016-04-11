#################################################################################################################
# PREPARATIONS
#################################################################################################################

# Load packages
library(plyr) # For ddply
library(Deducer) # For G-test
library(survival) # For Cox model
library(timereg) # For Aalen and Cox models
library(simPH) # For plotting Cox model 
library(MASS) # For negative binomial glm
library(lmtest) # For likelihood ratio test

# Read data into R
lh_wide <- read.csv("~/Desktop/GitHub/Gelada_parasites/Data/lh_wide.csv", header=TRUE, stringsAsFactors=FALSE)
lh_long <- read.csv("~/Desktop/GitHub/Gelada_parasites/Data/lh_long.csv", header=TRUE, stringsAsFactors=FALSE)
lh_long2 <- read.csv("~/Desktop/GitHub/Gelada_parasites/Data/lh_long2.csv", header=TRUE, stringsAsFactors=FALSE) # Starting at 4 years old
lh_ibi <- read.csv("~/Desktop/GitHub/Gelada_parasites/Data/lh_ibi.csv", header=TRUE, stringsAsFactors=FALSE) 
samps <- read.csv("~/Desktop/GitHub/Gelada_parasites/Data/samps.csv", header=TRUE, stringsAsFactors=FALSE)

# Format dates in 'wide' dataset
lh_wide$dob <- as.POSIXct(lh_wide$dob, format="%Y-%m-%d")
lh_wide$cyst <- as.POSIXct(lh_wide$cyst, format="%Y-%m-%d")
lh_wide$dod <- as.POSIXct(lh_wide$dod, format="%Y-%m-%d")
lh_wide$end <- as.POSIXct(lh_wide$end, format="%Y-%m-%d")

#################################################################################################################
# PERIOD PREVALENCE AND G-TESTS FOR LIFE HISTORY DATA (use lh_wide)
#################################################################################################################

############
## clean  #
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
## summarize #
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
## analyze #
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

#################################################
## barplots comparing prevalence and mortality #
###############################################

# Barplot for cyst prevalence over course of study
quartz()
bplot1 <- barplot(rev(tab0[,2]/rowSums(tab0)), col=c("dimgray", "lightgray"), names.arg=c("Male (n = 133)", "Female (n = 218)"), ylim=c(0, 1), ylab="% adults with cysts during study")
errs.m <- prop.test(tab0["M",2], rowSums(tab0)["M"])
errs.f <- prop.test(tab0["F",2], rowSums(tab0)["F"])
arrows(x0=c(bplot1), y0=c(errs.m$conf.int[1], errs.f$conf.int[1]), x1=c(bplot1), y1=c(errs.m$conf.int[2], errs.f$conf.int[2]), length=0.07, angle=90, code=3)

# Barplots for mortality among males/females with/without cysts
quartz()
layout(matrix(1:2, 1, 2))
bplot2 <- barplot(rev(tab3[2,]/colSums(tab3)), col=c("dimgray", "lightgray"), ylim=c(0,1), xlab="Females", ylab="% adults who died during study", names.arg=c("Cyst (29/40)", "No cyst (38/178)"))
errs.c <- prop.test(tab3[2,2], colSums(tab3)[2])
errs.n <- prop.test(tab3[2,1], colSums(tab3)[1])
arrows(x0=c(bplot2), y0=c(errs.c$conf.int[1], errs.n$conf.int[1]), x1=c(bplot2), y1=c(errs.c$conf.int[2], errs.n$conf.int[2]), length=0.07, angle=90, code=3)
bplot3 <- barplot(rev(tab2[2,]/colSums(tab2)), col=c("dimgray", "lightgray"), ylim=c(0,1), xlab="Males", names.arg=c("Cyst (3/12)", "No cyst (9/121)"))
errs.c2 <- prop.test(tab2[2,2], colSums(tab2)[2])
errs.n2 <- prop.test(tab2[2,1], colSums(tab2)[1])
arrows(x0=c(bplot3), y0=c(errs.c2$conf.int[1], errs.n2$conf.int[1]), x1=c(bplot3), y1=c(errs.c2$conf.int[2], errs.n2$conf.int[2]), length=0.07, angle=90, code=3)

#################################################
## line plot showing individual life histories #
###############################################

# Subset to individuals that ever had cysts
lp <- lh_wide[!is.na(lh_wide$cyst),]
nrow(lp) # 52 individuals

# Sort first by sex and then by length of observation time
lp <- lp[order(factor(lp$sex, levels=c("F", "M")), lp$end - lp$dob, decreasing=FALSE),]

# Create plot window
quartz()
par(mar=c(5, 6, 4, 2))
plot(0:(nrow(lp)+1) ~ seq(0, 25, length.out=length(0:(nrow(lp)+1))), data=lp, type="n", yaxt="n", ylab="", xlab="Time (years)")
ynames <- c(lp$name[lp$sex=="F"], "", lp$name[lp$sex=="M"])
axis(side=2, at=1:length(ynames), labels=ynames, las=2, cex.axis=0.5, tcl=0, mgp=c(3, 0.3, 0)) # y axis
text(x=c(-5, -5), y=c(53, 40), labels=c("Males", "Females"), xpd=TRUE) # Indicate male and female names

# Draw life history timelines
yloc <- (1:length(ynames))[-(sum(lp$sex=="F")+1)]
life <- c(lp[,"end"] - lp[,"dob"])/365
cyststart <- c(lp[,"cyst"] - lp[,"dob"])/365
for (i in 1:length(yloc)) {
	segments(x0=0, y0=yloc[i], x1=life[i], y1=yloc[i], col="darkgray") # Life line
	segments(x0=cyststart[i], y0=yloc[i], x1=life[i], y1=yloc[i], lwd=3) # Cyst line
	points(x=life[i], y=yloc[i], pch=21, cex=0.8, bg="white") # Censoring indicator
	if (!is.na(lp[i,"dod"])) {points(x=life[i], y=yloc[i], pch=4, cex=0.8)}
}

#################################################################################################################
# SURVIVAL ANALYSIS FOR LIFE HISTORY DATA (use lh_long)
#################################################################################################################

#############################################
## Non-parametric survival curve estimates #
###########################################

# Kaplan-Meier curves for cyst vs no-cyst (males only, females only, combined)
quartz()
layout(matrix(1:2, 1, 2))
plot(survfit(Surv(start, stop, death==1) ~ cyst, data=lh_long[lh_long2$sex=="M",], type="kaplan-meier"), col=c("lightgray", "black"), main="Males")
plot(survfit(Surv(start, stop, death==1) ~ cyst, data=lh_long[lh_long2$sex=="F",], type="kaplan-meier"), col=c("lightgray", "black"), main="Females")

# Fleming-Harrington curves for cyst vs no-cyst (males only, females only, combined)
quartz()
layout(matrix(1:2, 1, 2))
plot(survfit(Surv(start, stop, death==1) ~ cyst, data=lh_long[lh_long2$sex=="M",], type="fleming-harrington"), col=c("lightgray", "black"), main="Males")
plot(survfit(Surv(start, stop, death==1) ~ cyst, data=lh_long[lh_long2$sex=="F",], type="fleming-harrington"), col=c("lightgray", "black"), main="Females")

# Why does the Kaplan-Meier estimator produce such an abrupt crash in the survival curve?

####################################
## Cox proportional hazards model #
##################################

# Cox proportional hazards model (Model 1)
coxM <- coxph(Surv(start, stop, death) ~ cyst, data=lh_long[lh_long$sex=="M",])
coxF <- coxph(Surv(start, stop, death) ~ cyst, data=lh_long[lh_long$sex=="F",])
summary(coxM) # Cyst significant (exp(coef) = 2.63, p = 0.21)
summary(coxF) # Cyst significant (exp(coef) = 10.16, p < 0.001 ***)

# Plot estimated survival curves for cyst vs no-cyst
# It appears that the severity of the effect of cysts is underestimated compared to the non-parametric estimations
# I guess this makes sense since the proportional hazards model with estimate a hazard that is intermediate between
# the early severe effect and the later mild effect
quartz()
layout(matrix(1:2, 1, 2))
cyst <- data.frame(cyst=c(0, 1)) # Separate curves for cyst/no cyst
plot(survfit(coxM, newdata=cyst), col=c("dimgray", "lightgray"), xlab="Months", ylab="Proportion alive", main="Males")
legend("bottomleft", legend=c("No cyst", "Cyst"), fill=c("dimgray", "lightgray"))
plot(survfit(coxF, newdata=cyst), col=c("dimgray", "lightgray"), xlab="Months", ylab="", main="Females")

# Test proportional hazards assumption
zphM <- cox.zph(coxM) # Cyst not significant (p = 0.08.)
zphF <- cox.zph(coxF) # Cyst not significant (p < 0.01**)
summary(zphM)
summary(zphF)

# Plot Schoenfeld residuals vs. time plots (point should fall along line with slope=0 if proportional hazards is valid)
quartz()
layout(matrix(1:2, 1, 2))
plot(zphM) # Looks not good (hazard decreases through time)
plot(zphF) # Looks not good (hazard decreases through time)

##############################################################
## Cox proportional hazards model with time-varying effects #
############################################################

# Cox proportional hazards model (Model 2) with time-varying effect of cysts
coxM2 <- timecox(Surv(start, stop, death) ~ cyst, data=lh_long[lh_long$sex=="M",], residuals=1)
coxF2 <- timecox(Surv(start, stop, death) ~ cyst, data=lh_long[lh_long$sex=="F",], residuals=1)
summary(coxM2) # Significant varying effect of cysts (p < 0.01**)
summary(coxF2) # Significant varying effect of cysts (p < 0.001***)

# Top: Cumulative regression coefficients for males and females
# Bottom: Observed test process compared to 50 simulated test processes for males and females
quartz()
layout(matrix(1:4, 2, 2))
plot(coxM2, specific.comps=2, main="Males", xlab="") # Cyst has stronger effect early on
plot(coxM2, score=T, specific.comps=2, ylab="") # Constant hazard rate violated
plot(coxF2, specific.comps=2, main="Females", xlab="") # Cyst has stronger effect early on
plot(coxF2, score=T, specific.comps=2, ylab="") # Constant hazard rate violated
### Need to figure out how to change axis/main labels on plots ###

##############################################################
## Aalen's additive hazards model with time-varying effects #
############################################################

# Aalen's additive hazards model with time-varying effects (Model 3)
aalM <- aalen(Surv(start, stop, death) ~ cyst, data=lh_long[lh_long$sex=="M",], residuals=1)
aalF <- aalen(Surv(start, stop, death) ~ cyst, data=lh_long[lh_long$sex=="F",], residuals=1)
summary(aalM) # Significant varying effect of cysts (p < 0.001***)
summary(aalF) # Significant constant effect of cysts (p < 0.001***)

# Top: Cumulative regression coefficients for males and females
# Bottom: Observed test process compared to 50 simulated test processes for males and females
quartz()
layout(matrix(1:4, 2, 2))
plot(aalM, specific.comps=2, main="Males", xlab="") # Cyst has stronger effect early on
plot(aalM, score=T, specific.comps=2, ylab="") # Constant hazard coefficient violated
plot(aalF, specific.comps=2, main="Females", xlab="") # Cyst has a strong constant effect
plot(aalF, score=T, specific.comps=2, ylab="") # Constant hazard coefficient NOT violated

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