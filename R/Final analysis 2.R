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
lh_ibi <- read.csv("~/Desktop/GitHub/Gelada_parasites/Data/lh_ibi.csv", header=TRUE, stringsAsFactors=FALSE) 
samps <- read.csv("~/Desktop/GitHub/Gelada_parasites/Data/samps.csv", header=TRUE, stringsAsFactors=FALSE)

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
# PERIOD PREVALENCE AND G-TESTS FOR LIFE HISTORY DATA (use lh_wide)
#################################################################################################################

############
## clean  #
##########

# How many individuals to start?
nrow(lh_wide) # 388

# Exclude individuals who disappeared 
sum(!is.na(lh_wide$dod[lh_wide$sex=="M"])) # 12 males died
sum(is.na(lh_wide$dod[lh_wide$sex=="M"])) # 158 males did not die
sum(lh_wide$sex=="F" & lh_wide$dis==1)/sum(lh_wide$sex=="F") # 0 of 218 females disappeared
sum(lh_wide$sex=="M" & lh_wide$dis==1)/sum(lh_wide$sex=="M") # 37 of 170 (21.8% of males disappeared)
37/158 # 23.4% of censored males 'disappeared'
lh_wide <- lh_wide[!lh_wide$dis==1,-which(names(lh_wide)=="dis")]
nrow(lh_wide) # 351 remaining

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
life <- c(difftime(lp[,"end"], lp[,"dob"], unit="days"))/365.25
cyststart <- c(lp[,"cyst"] - lp[,"dob"])/365
for (i in 1:length(yloc)) {
	linecol <- ifelse(!is.na(lp[i,"dod"]), "black", "blue") 
	segments(x0=0, y0=yloc[i], x1=life[i], y1=yloc[i], col="darkgray") # Life line
	segments(x0=cyststart[i], y0=yloc[i], x1=life[i], y1=yloc[i], lwd=3, col=linecol) # Cyst line
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
plot(survfit(Surv(start, stop, death==1) ~ cyst, data=lh_long[lh_long$sex=="M",], type="kaplan-meier"), col=c("lightgray", "black"), main="Males")
plot(survfit(Surv(start, stop, death==1) ~ cyst, data=lh_long[lh_long$sex=="F",], type="kaplan-meier"), col=c("lightgray", "black"), main="Females")

# Fleming-Harrington curves for cyst vs no-cyst (males only, females only, combined)
quartz()
layout(matrix(1:2, 1, 2))
plot(survfit(Surv(start, stop, death==1) ~ cyst, data=lh_long[lh_long$sex=="M",], type="fleming-harrington"), col=c("lightgray", "black"), main="Males")
plot(survfit(Surv(start, stop, death==1) ~ cyst, data=lh_long[lh_long$sex=="F",], type="fleming-harrington"), col=c("lightgray", "black"), main="Females")

####################################
## Cox proportional hazards model #
##################################

# Drop over 24 yrs old
drop <- lh_long[lh_long$sex=="F" & lh_long$stop >= 24,"name"]
fdat <- lh_long[lh_long$sex=="F" & !(lh_long$name %in% drop),]

# Cox proportional hazards model (Model 1)
coxM <- coxph(Surv(start, stop, death) ~ cyst, data=lh_long[lh_long$sex=="M",])
coxF <- coxph(Surv(start, stop, death) ~ cyst, data=fdat)
summary(coxM) # Cyst not significant (exp(coef) = 2.52, p = 0.23)
summary(coxF) # Cyst significant (exp(coef) = 10.16, p < 0.001 ***)
cox.zph(coxM)
cox.zph(coxF)
layout(matrix(1:2, 1,2))
plot(cox.zph(coxM))
plot(cox.zph(coxF))

# Cox proportional hazards with a time transform term (Model 1b)
coxMb <- coxph(Surv(start, stop, death) ~ cyst + tt(cyst), data=lh_long[lh_long$sex=="M",], tt=function(x, t, ...) x*t)
coxFb <- coxph(Surv(start, stop, death) ~ cyst + tt(cyst), data=fdat, tt=function(x, t, ...) x*t)
summary(coxMb)
summary(coxFb)
cox.zph(coxMb)
cox.zph(coxFb)
plot(cox.zph(coxMb))
plot(cox.zph(coxFb))

# Hazard function? Why are coefficients so large???
layout(matrix(1:2,1,2))
s <- seq(0, 25, length.out=100)
ym <- function (x) {exp(coef(coxMb)[1] + coef(coxMb)[2]*x)}
yf <- function (x) {exp(coef(coxFb)[1] + coef(coxFb)[2]*x)}
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
ym <- function (x) {(coef(coxMb)[1] + coef(coxMb)[2]*x)}
yf <- function (x) {(coef(coxFb)[1] + coef(coxFb)[2]*x)}
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


# Plot estimated survival curves for cyst vs no-cyst
# It appears that the severity of the effect of cysts is underestimated compared to the non-parametric estimations
# I guess this makes sense since the proportional hazards model will estimate a hazard that is intermediate between
# the early severe effect and the later mild effect
quartz()
layout(matrix(1:2, 1, 2))
new <- data.frame(name=c("A", "A", "B"), start=c(0, 15, 0), stop=c(15, 25, 25), death=c(0, 1, 1), cyst=c(1, 1, 0))
plot(survfit(coxM, newdata=new, id=c("A", "A", "B")), col=c("dimgray", "lightgray"), xlab="Years", ylab="Proportion alive", main="Males", conf.int=T, xlim=c(0, 18))
legend("bottomleft", legend=c("No cyst", "Cyst"), fill=c("dimgray", "lightgray"))
plot(survfit(coxF, newdata=new, id=c("A", "A", "B")), col=c("dimgray", "lightgray"), xlab="Years", ylab="", main="Females", conf.int=T, xlim=c(0, 25))

##############################################################
## Cox proportional hazards model with time-varying effects #
############################################################

# Cox proportional hazards model (Model 2) with time-varying effect of cysts
coxM2 <- timecox(Surv(start, stop, death) ~ cyst, data=lh_long[lh_long$sex=="M",], id=lh_long[lh_long$sex=="M","name"], residuals=1, n.sim=100)
coxF2 <- timecox(Surv(start, stop, death) ~ cyst, data=lh_long[lh_long$sex=="F",], residuals=1, n.sim=100)
summary(coxM2) # Significant varying effect of cysts (p < 0.01**)
summary(coxF2) # Significant varying effect of cysts (p < 0.001***)

# Top: Cumulative regression coefficients for males and females
# Bottom: Observed test process compared to 50 simulated test processes for males and females
quartz()
layout(matrix(1:4, 2, 2))
plot(coxM2, specific.comps=2, mains=F, xlab="") # Cyst has stronger effect early on
mtext("Males", font=2, line=1.5)
plot(coxM2, score=T, ylab="Test process", specific.comps=2) # Constant hazard rate violated
plot(coxF2, specific.comps=2, mains=F,xlab="") # Cyst has stronger effect early on
mtext("Females", font=2, line=1.5)
plot(coxF2, score=T, ylab="Test process", specific.comps=2, ylab="") # Constant hazard rate violated
### Need to figure out how to change axis/main labels on plots ###

## Just cumulative coefficients for India's talk
quartz()
layout(matrix(1:2, 1, 2))
plot(coxM2, specific.comps=2, mains=F, xlab="Time", start.time=4, stop.time=16, ylab="Cumulative coefficient") # Cyst has stronger effect early on
mtext("Males", font=2, line=1.5)
plot(coxF2, specific.comps=2, mains=F,xlab="Time", start.time=4, stop.time=24, ylab="") # Cyst has stronger effect early on
mtext("Females", font=2, line=1.5)

# cumulative residuals for M
temp <- lh_long[lh_long$sex=="M",]
m <- as.matrix(temp$cyst)
colnames(m) <- "Cyst"
rownames(m) <- rownames(temp)
test <- cum.residuals(coxM2, data=temp, cum.resid=0, modelmatrix=m)
quartz()
plot(test, score=1, conf.band=T)
# cumulative residuals for F
temp <- lh_long[lh_long$sex=="F",]
m <- as.matrix(temp$cyst)
colnames(m) <- "Cyst"
rownames(m) <- rownames(temp)
test <- cum.residuals(coxM2, data=temp, cum.resid=0, modelmatrix=m)
quartz()
plot(test, score=1, conf.band=T)

##############################################################
## Aalen's additive hazards model with time-varying effects #
############################################################

# Aalen's additive hazards model with time-varying effects (Model 3)
aalM <- aalen(Surv(start, stop, death) ~ cyst, data=lh_long[lh_long$sex=="M",], residuals=1, resample.iid=1, n.sim=100)
aalF <- aalen(Surv(start, stop, death) ~ cyst, data=lh_long[lh_long$sex=="F",], residuals=1, resample.iid=1, n.sim=100)
summary(aalM) # Significant varying effect of cysts (p < 0.001***)
summary(aalF) # Significant constant effect of cysts (p < 0.001***)

# Top: Cumulative regression coefficients for males and females
# Bottom: Observed test process compared to 50 simulated test processes for males and females
quartz()
layout(matrix(1:4, 2, 2))
plot(aalM, specific.comps=2, main="Males", xlab="") # Cyst has stronger effect early on
plot(aalM, score=T, ylab="Test process", specific.comps=2, ylab="") # Constant hazard coefficient violated
plot(aalF, specific.comps=2, main="Females", xlab="") # Cyst has a strong constant effect
plot(aalF, score=T, ylab="Test process", specific.comps=2, ylab="") # Constant hazard coefficient NOT violated

# cumulative residuals for M
temp <- lh_long[lh_long$sex=="M",]
m <- as.matrix(temp$cyst)
colnames(m) <- "Cyst"
rownames(m) <- rownames(temp)
test <- cum.residuals(aalM, data=temp, cum.resid=0, modelmatrix=m)
quartz()
layout(matrix(1:2, 1, 2))
plot(test, score=1)
temp <- lh_long[lh_long$sex=="F",]
m <- as.matrix(temp$cyst)
colnames(m) <- "Cyst"
rownames(m) <- rownames(temp)
test <- cum.residuals(aalF, data=temp, cum.resid=0, modelmatrix=m)
plot(test, score=1)

##############################################################
## Survival curves comparing:  #
## 	1.	Fleming-Harrington estimator
##	2.	Cox proportional hazards model
##	3.	Time-varying Cox proportional hazards model
##	4.	Time-varying Aalen's additive hazards model
############################################################

#################################################################################################################
# SURVIVAL ANALYSIS FOR ONLY FEMALES WITH CYSTS (use modified lh_long)
#################################################################################################################

# Create modified data frame
lh_cyst <- lh_long[lh_long$cyst==T & lh_long$sex=="F",]
time <- lh_cyst$stop - lh_cyst$start
names(lh_cyst)[names(lh_cyst)=="start"] <- "age"
lh_cyst <- cbind(lh_cyst, time)

# Fit cox ph model
coxC <- coxph(Surv(time, death) ~ age, data=lh_cyst)
summary(coxC)

# Test proportional hazards assumption
zphC <- cox.zph(coxC)
zphC
plot(zphC)

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
# END
#################################################################################################################