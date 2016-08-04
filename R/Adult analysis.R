#################################################################################################################
# PREPARATIONS
#################################################################################################################

# Load packages
library(Deducer) # For G-test
library(survival) # For survival models

# Read data into R
adult_wide <- read.csv("~/Desktop/GitHub/Gelada_parasites/Data/adult_wide.csv", header=TRUE, stringsAsFactors=FALSE)
adult_long <- read.csv("~/Desktop/GitHub/Gelada_parasites/Data/adult_long.csv", header=TRUE, stringsAsFactors=FALSE)

###############################################################
## Line plot showing individual life histories in adult_wide #
#############################################################

# Subset to individuals that ever had cysts
lp <- adult_wide[!is.na(adult_wide$cyst),]
lp <- lp[-which(lp$name=="Delia"),] # Sole individual who developed a cyst and recovered
nrow(lp) # 54 individuals
nrow(lp[lp$sex=="M",]) # 16 males
nrow(lp[lp$sex=="F",]) # 38 females

# Sort first by sex and then by length of observation time
lp <- lp[order(factor(lp$sex, levels=c("F", "M")), lp$stop, decreasing=FALSE),]

# Create plot window
quartz()
par(mar=c(5, 6, 4, 2))
plot(0:(nrow(lp)+1) ~ seq(0, 25, length.out=length(0:(nrow(lp)+1))), data=lp, type="n", yaxt="n", ylab="", xlab="Time (years)")
ynames <- c(toupper(substr(lp$name[lp$sex=="F"], 1, 3)), "", toupper(substr(lp$name[lp$sex=="M"], 1, 3)))
axis(side=2, at=1:length(ynames), labels=ynames, las=2, cex.axis=0.5, tcl=0, mgp=c(3, 0.3, 0)) # y axis
text(x=c(-6, -6), y=c(57, 39), labels=c("Males", "Females"), xpd=TRUE) # Indicate male and female names

# Draw life history timelines
yloc <- (1:length(ynames))[-(sum(lp$sex=="F")+1)]
for (i in 1:length(yloc)) {
	ifelse (lp[i,"cyst"] < lp[i,"start"], cyststart <- lp[i,"start"], cyststart <- lp[i,"cyst"])
	segments(x0=lp[i,"start"], y0=yloc[i], x1=lp[i,"stop"], y1=yloc[i], col="darkgray") # Life line
	segments(x0=lp[i,"cyst"], y0=yloc[i], x1=lp[i,"stop"], y1=yloc[i], lwd=3) # Cyst line
	points(x=lp[i,"stop"], y=yloc[i], pch=21, cex=0.8, bg="white") # Censoring indicator
	if (!is.na(lp[i,"dod"])) {points(x=lp[i,"stop"], y=yloc[i], pch=4, cex=0.8)}
}

#################################################################################################################
# G-TESTS FOR PREVALENCE AND MORTALITY (use adult_wide)
#################################################################################################################

# Number of adults, males, and females
nrow(adult_wide) # 387
nrow(adult_wide[adult_wide$sex=="M",]) # 170 males
nrow(adult_wide[adult_wide$sex=="F",]) # 217 females

# Limit to last 6.5 years of study
adult_wide <- adult_wide[adult_wide$drop == 0,]

# New number of adults, males, and females (dropped 36: 9 males and 27 females)
nrow(adult_wide) # 351
nrow(adult_wide[adult_wide$sex=="M",]) # 161 males
nrow(adult_wide[adult_wide$sex=="F",]) # 190 females

#################
## Prevalence  #
###############

# Prevalence of adults, males, and females with cysts
tab0 <- table(adult_wide$sex, !is.na(adult_wide$cyst))
sum(tab0[,2])/sum(tab0)	# Overall period prevalence = 0.13
sum(tab0["M",2])/sum(tab0["M",]) 	# Male period prevalence = 0.09
sum(tab0["F",2])/sum(tab0["F",]) 	# Female period prevalence = 0.16

# G-test for independence between Guassa vs SMNP male prevalence
likelihood.test(rbind(tab0["M",], c(49, 19))) # G = 13.16, p < 0.001***

# G-test for independence between Guassa vs SMNP female prevalence
likelihood.test(rbind(tab0["F",], c(68, 31))) # G = 8.38, p < 0.01**

# G-test for independence between Guassa and SMNP adult versus immature prevalences
# (note these numbers were compiled from a data sheet that is not included here)
likelihood.test(rbind(c(174, 7),c(405, 2))) # G = 8.58, p < 0.01**

# G-test for independence between cysts and sex at SMNP
likelihood.test(tab0) # G = 4.66, p < 0.05*

# G-test for independence between cysts and age at SMNP
likelihood.test(rbind(c(306, 45),c(405, 2))) # G = 58.3, p < 0.001***

# Test for a unit-level effect
unit_test <- adult_wide[! adult_wide$unit %in% c("T?", "unknown"),]
length(unique(unit_test$unit)) # 27 units
likelihood.test(unit_test$unit, is.na(unit_test$cyst)) # G = 24.91, p > 0.05

###############
## Mortality #
#############

# Drop individuals who disappeared during the study
sum(!is.na(adult_wide$dod[adult_wide$sex=="M"])) # 12 males died
sum(is.na(adult_wide$dod[adult_wide$sex=="M"])) # 149 males did not die
sum(adult_wide$sex=="F" & adult_wide$dis==1)/sum(adult_wide$sex=="F") # 0 of 190 females disappeared
sum(adult_wide$sex=="M" & adult_wide$dis==1)/sum(adult_wide$sex=="M") # 31 of 161 (21.2% of males disappeared)
adult_wide2 <- adult_wide[!adult_wide$dis==1,-which(names(adult_wide)=="dis")]

# New number of adults, males, and females 
nrow(adult_wide2) # 320 individuals
sum(adult_wide2$sex=="M") # 130 males
sum(adult_wide2$sex=="F") # 190 females

# Number/proportion of adults dying with/without cysts over study period
tab1 <- table(!is.na(adult_wide2$dod), !is.na(adult_wide2$cyst))
tab1[2,]/colSums(tab1)	# 15% without cysts die; 59% with cysts die
tab2 <- table(!is.na(adult_wide2$dod[adult_wide2$sex=="M"]), !is.na(adult_wide2$cyst[adult_wide2$sex=="M"]))
tab2[2,]/colSums(tab2)	# 8% without cysts die; 30% with cysts die
tab3 <- table(!is.na(adult_wide2$dod[adult_wide2$sex=="F"]), !is.na(adult_wide2$cyst[adult_wide2$sex=="F"]))
tab3[2,]/colSums(tab3) # 20% without cysts die; 68% with cysts die

# G-test for association between cysts and deaths over study period
likelihood.test(tab1) # Adults: G = 34.47, p < 0.001***
likelihood.test(tab2) # Males: G = 3.89, p < 0.05*
likelihood.test(tab3) # Females: G = 26.28, p < 0.001***

# G-test for difference in cyst-associated mortality between SMNP and Guassa
likelihood.test(rbind(c(3, 10), c(9, 12))) # Males: G = 1.42, p > 0.05
likelihood.test(rbind(c(21, 31), c(9, 22))) # Females: G = 1.1, p > 0.05

#################################
## Barplot comparing mortality #
###############################

# Barplots for mortality among males/females with/without cysts
quartz()
layout(matrix(1:2, 1, 2))
bplotF <- barplot(rev(tab2[2,]/colSums(tab2)), col=c("dimgray", "lightgray"), ylim=c(0,1), xlab="", main="Males", names.arg=c("Cyst (3/10)", "No cyst (9/120)"), ylab="% adults that died during study")
errs.c2 <- prop.test(tab2[2,2], colSums(tab2)[2])
errs.n2 <- prop.test(tab2[2,1], colSums(tab2)[1])
arrows(x0=c(bplotF), y0=c(errs.c2$conf.int[1], errs.n2$conf.int[1]), x1=c(bplotF), y1=c(errs.c2$conf.int[2], errs.n2$conf.int[2]), length=0.07, angle=90, code=3)
mtext("A", line=2, cex=1.5, adj=0)
bplotM <- barplot(rev(tab3[2,]/colSums(tab3)), col=c("dimgray", "lightgray"), ylim=c(0,1), xlab="", ylab="", main="Females", names.arg=c("Cyst (21/31)", "No cyst (32/159)"))
errs.c <- prop.test(tab3[2,2], colSums(tab3)[2])
errs.n <- prop.test(tab3[2,1], colSums(tab3)[1])
arrows(x0=c(bplotM), y0=c(errs.c$conf.int[1], errs.n$conf.int[1]), x1=c(bplotM), y1=c(errs.c$conf.int[2], errs.n$conf.int[2]), length=0.07, angle=90, code=3)
mtext("B", line=2, cex=1.5, adj=0)

#################################################################################################################
# SURVIVAL ANALYSIS FOR LIFE HISTORY DATA (use adult_long)
#################################################################################################################

###############
## Data prep #
#############

# Male data
mdat <- adult_long[adult_long$sex=="M",]
fdat <- adult_long[adult_long$sex=="F",]

# Number censored
nmales <- length(unique(mdat$name))
nfemales <- length(unique(fdat$name))
(nmales - sum(mdat$death))/nmales # 92.4% censored
(nfemales - sum(fdat$death))/nfemales # 72.2% censored

# Number with cysts
sum(fdat$cyst) # 38
sum(mdat$cyst) # 16

################
## Cox models #
##############

# Cox proportional hazards model (Model 1)
coxM <- coxph(Surv(start, stop, death) ~ cyst, data=mdat, x=TRUE)
coxF <- coxph(Surv(start, stop, death) ~ cyst, data=fdat, x=TRUE)
summary(coxM) # Cyst not significant (exp(coef) = 1.81, p = 0.41)
summary(coxF) # Cyst significant (exp(coef) = 5.77, p < 0.001 ***)

# Cox-Snell residuals
rr <- fdat$death - coxF$residual
fit <- survfit(Surv(rr, fdat$death)~1)
s <- fit$surv
t <- fit$time
plot(s ~ t)
rr <- mdat$death - coxM$residual
fit <- survfit(Surv(rr, mdat$death)~1)
s <- fit$surv
t <- fit$time
plot(s ~ t)

# Test proportional hazards assumption (model 1)
cox.zph(coxM, transform="identity") # Check
cox.zph(coxF, transform="identity") # p < 0.001 ***, Violated!

# Cox proportional hazards with a time transform (Model 2)
coxM2 <- coxph(Surv(start, stop, death) ~ cyst + tt(cyst), data=mdat, tt=function(x, t, ...) x*t, x=TRUE)
coxF2 <- coxph(Surv(start, stop, death) ~ cyst + tt(cyst), data=fdat, tt=function(x, t, ...) x*t, x=TRUE)
summary(coxM2) # Cyst approaching significance, p = 0.06
summary(coxF2) # Cyst and time transform significant, p < 0.001 ***

######################################
## Plot scaled Schoenfeld residuals #
####################################

# Residual plots showing the nature of the proportional hazard violation
layout(matrix(1:2, 1, 2, byrow=T))
plot(cox.zph(coxM, transform="identity"), xlab="Time (years)", main="Males")
mtext("A", line=2, cex=1.5, adj=0)
plot(cox.zph(coxF, transform="identity"), xlab="Time (years)", main="Females")
mtext("B", line=2, cex=1.5, adj=0)

############################
## Plot log hazard ratios #
##########################

layout(matrix(1:2,1,2))
s <- seq(0, 25, length.out=100)
ym <- function (x) {(coef(coxM2)[1] + coef(coxM2)[2]*x)}
yf <- function (x) {(coef(coxF2)[1] + coef(coxF2)[2]*x)}
ym.l <- function (x) {log(summary(coxM2)$conf.int["cyst","lower .95"]) + log(summary(coxM2)$conf.int["tt(cyst)","lower .95"])*x}
ym.u <- function (x) {log(summary(coxM2)$conf.int["cyst","upper .95"]) + log(summary(coxM2)$conf.int["tt(cyst)","upper .95"])*x}
yf.l <- function (x) {log(summary(coxF2)$conf.int["cyst","lower .95"]) + log(summary(coxF2)$conf.int["tt(cyst)","lower .95"])*x}
yf.u <- function (x) {log(summary(coxF2)$conf.int["cyst","upper .95"]) + log(summary(coxF2)$conf.int["tt(cyst)","upper .95"])*x}
plot(ym, ylim=c(-10,10), xlim=c(4, 15), ylab="Log hazard ratio", lwd=2.5, xlab="Time (years)", main="Males")
lines(s, ym.l(s))
lines(s, ym.u(s))
abline(h=coef(coxM)[1], lty=2)
abline(h=0, col="gray", lty=3)
mtext("A", line=2, cex=1.5, adj=0)
plot(yf, ylim=c(-10,10), xlim=c(4, 24), ylab="", lwd=2.5, xlab="Time (years)", main="Females")
lines(s, yf.l(s))
lines(s, yf.u(s))
abline(h=0, col="gray", lty=3)
abline(h=coef(coxF)[1], lty=2)
mtext("B", line=2, cex=1.5, adj=0)

#################################################################################################################
# END
#################################################################################################################
