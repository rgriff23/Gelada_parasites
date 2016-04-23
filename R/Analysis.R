#################################################################################################################
# PREPARATIONS
#################################################################################################################

# Load packages
library(Deducer) # For G-test
library(survival) # For Cox model

# Read data into R
lh_wide <- read.csv("~/Desktop/GitHub/Gelada_parasites/Data/lh_wide.csv", header=TRUE, stringsAsFactors=FALSE)
lh_long <- read.csv("~/Desktop/GitHub/Gelada_parasites/Data/lh_long.csv", header=TRUE, stringsAsFactors=FALSE)

# Format dates in 'wide' dataset
lh_wide$dob <- as.POSIXct(lh_wide$dob, format="%Y-%m-%d")
lh_wide$cyst <- as.POSIXct(lh_wide$cyst, format="%Y-%m-%d")
lh_wide$dod <- as.POSIXct(lh_wide$dod, format="%Y-%m-%d")
lh_wide$end <- as.POSIXct(lh_wide$end, format="%Y-%m-%d")

# Exclude individuals who END before the age of four from data
drop <- lh_wide$name[(c(lh_wide$end - lh_wide$dob)/60/60/24/365) < 4]
length(drop) # 233 <4 years old
sum(!is.na(lh_wide[lh_wide$name %in% drop,"cyst"])) # 0 cysts among the dropped individuals
lh_wide <- lh_wide[!(lh_wide$name %in% drop),] 
lh_long <- lh_long[lh_long$name %in% lh_wide$name,]

#################################################################################################################
# G-TESTS FOR PREVALENCE AND MORTALITY (use lh_wide)
#################################################################################################################

#################
## prevalence  #
###############

# Number of adults, males, and females
nrow(lh_wide) # 388
nrow(lh_wide[lh_wide$sex=="M",]) # 170 males
nrow(lh_wide[lh_wide$sex=="F",]) # 218 females

# Prevalence of adults, males, and females with cysts
tab0 <- table(lh_wide$sex, !is.na(lh_wide$cyst))
sum(tab0[,2])/sum(tab0)	# Overall period prevalence = 0.15
sum(tab0["M",2])/sum(tab0["M",])	# Male period prevalence = 0.1
sum(tab0["F",2])/sum(tab0["F",])	# Female period prevalence = 0.18

# G-test for independence between cysts and sex
likelihood.test(tab0) # Significant association (p < 0.05)

# G-test for difference between Guassa and SMNP prevalences
likelihood.test(rbind(tab0["F",], c(68, 31))) # Significant association (p < 0.05)
likelihood.test(rbind(tab0["M",], c(49, 19))) # Highly significant association (p < 0.001)

###############
## mortality #
#############

# Exclude individuals who disappeared 
sum(!is.na(lh_wide$dod[lh_wide$sex=="M"])) # 13 males died
sum(is.na(lh_wide$dod[lh_wide$sex=="M"])) # 157 males did not die
sum(lh_wide$sex=="F" & lh_wide$dis==1)/sum(lh_wide$sex=="F") # 0 of 218 females disappeared
sum(lh_wide$sex=="M" & lh_wide$dis==1)/sum(lh_wide$sex=="M") # 36 of 170 (21.2% of males disappeared)
36/157 # 22.9% of censored males 'disappeared'
lh_wide <- lh_wide[!lh_wide$dis==1,-which(names(lh_wide)=="dis")]

# New number of adults, males, and females 
nrow(lh_wide) # 352 individuals
sum(lh_wide$sex=="M") # 134 males
sum(lh_wide$sex=="F") # 218 females

# Number/proportion of adults dying with/without cysts over study period
tab1 <- table(!is.na(lh_wide$dod), !is.na(lh_wide$cyst))
tab1[2,]/colSums(tab1)	# 16% without cysts die; 61.5% with cysts die
tab2 <- table(!is.na(lh_wide$dod[lh_wide$sex=="M"]), !is.na(lh_wide$cyst[lh_wide$sex=="M"]))
tab2[2,]/colSums(tab2)	# 8.2% without cysts die; 25% with cysts die
tab3 <- table(!is.na(lh_wide$dod[lh_wide$sex=="F"]), !is.na(lh_wide$cyst[lh_wide$sex=="F"]))
tab3[2,]/colSums(tab3) # 21.3% without cysts die; 72.5% with cysts die

# G-test for association between cysts and deaths over study period
likelihood.test(tab1) # Highly significant across combined sexes (p < 0.001)
likelihood.test(tab2) # Approaching significance in males (p = 0.1)
likelihood.test(tab3) # Highly significant in females  (p < 0.001)

################################################
## barplot comparing prevalence and mortality #
##############################################

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
life <- c(difftime(lp[,"end"], lp[,"dob"], unit="days"))/365.25
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

####################################
## Cox proportional hazards model #
##################################

# Male data
mdat <- lh_long[lh_long$sex=="M",]

# Female data (drop over 24 yrs old)
drop <- lh_long[lh_long$sex=="F" & lh_long$stop >= 24,"name"]
fdat <- lh_long[lh_long$sex=="F" & !(lh_long$name %in% drop),]

# Cox proportional hazards model (Model 1)
coxM <- coxph(Surv(start, stop, death) ~ cyst, data=mdat)
coxF <- coxph(Surv(start, stop, death) ~ cyst, data=fdat)
summary(coxM) # Cyst not significant (exp(coef) = 2.52, p = 0.23)
summary(coxF) # Cyst significant (exp(coef) = 10.16, p < 0.001 ***)

# Test proportional hazards assumption (model 1)
cox.zph(coxM) # Violated!
cox.zph(coxF) # Violated!

# Residual plots showing the nature of the proportional hazard violation
layout(matrix(1:2, 1,2))
plot(cox.zph(coxM))
plot(cox.zph(coxF))

# Cox proportional hazards with a time transform (Model 2)
coxM2 <- coxph(Surv(start, stop, death) ~ cyst + tt(cyst), data=mdat, tt=function(x, t, ...) x*t)
coxF2 <- coxph(Surv(start, stop, death) ~ cyst + tt(cyst), data=fdat, tt=function(x, t, ...) x*t)
summary(coxM2)
summary(coxF2)

# Test proportional hazards assumption (model 2)
cox.zph(coxM2) # Check
cox.zph(coxF2) # Check

############################
## Plot log hazard ratios #
##########################

# Hazard function? Why are coefficients so large???
layout(matrix(1:2,1,2))
s <- seq(0, 25, length.out=100)
ym <- function (x) {exp(coef(coxM2)[1] + coef(coxM2)[2]*x)}
yf <- function (x) {exp(coef(coxF2)[1] + coef(coxF2)[2]*x)}
ym.l <- function (x) {exp(log(6.42) + log(0.47)*x)}
ym.u <- function (x) {exp(log(21547.82) + log(0.95)*x)}
yf.l <- function (x) {exp(log(80.8) + log(0.73)*x)}
yf.u <- function (x) {exp(log(1466.88) + log(0.87)*x)}
plot(ym, ylim=c(0,22000), xlim=c(0, 25), xlab="Time", ylab="Hazard ratio", main="Males")
lines(s, ym.l(s), lty=2)
lines(s, ym.u(s), lty=2)
plot(yf, ylim=c(0,1500), xlim=c(0, 25), xlab="Time", ylab="", main="Females")
lines(s, yf.l(s), lty=2)
lines(s, yf.u(s), lty=2)
# log
layout(matrix(1:2,1,2))
s <- seq(0, 25, length.out=100)
ym <- function (x) {(coef(coxM2)[1] + coef(coxM2)[2]*x)}
yf <- function (x) {(coef(coxF2)[1] + coef(coxF2)[2]*x)}
ym.l <- function (x) {(log(6.42) + log(0.47)*x)}
ym.u <- function (x) {(log(21547.82) + log(0.95)*x)}
yf.l <- function (x) {(log(80.8) + log(0.73)*x)}
yf.u <- function (x) {(log(1466.88) + log(0.87)*x)}
plot(ym, ylim=c(-10,10), xlim=c(0, 15), xlab="Time", ylab="Log hazard ratio", main="Males", lwd=2)
lines(s, ym.l(s), lty=2)
lines(s, ym.u(s), lty=2)
abline(h=0, col="gray", lty=3)
plot(yf, ylim=c(-10,10), xlim=c(0, 25), xlab="Time", ylab="", main="Females", lwd=2)
lines(s, yf.l(s), lty=2)
lines(s, yf.u(s), lty=2)
abline(h=0, col="gray", lty=3)

#################################################################################################################
# END
#################################################################################################################

