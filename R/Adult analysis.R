#################################################################################################################
# PREPARATIONS
#################################################################################################################

# Load packages
library(Deducer) # For G-test
library(survival) # For survival models

# Read data into R
lh_wide <- read.csv("~/Desktop/GitHub/Gelada_parasites/Data/lh_wide.csv", header=TRUE, stringsAsFactors=FALSE)
lh_long <- read.csv("~/Desktop/GitHub/Gelada_parasites/Data/lh_long.csv", header=TRUE, stringsAsFactors=FALSE)

#################################################################################################################
# G-TESTS FOR PREVALENCE AND MORTALITY (use lh_wide)
#################################################################################################################

#################
## Prevalence  #
###############

# Number of adults, males, and females
nrow(lh_wide) # 387
nrow(lh_wide[lh_wide$sex=="M",]) # 170 males
nrow(lh_wide[lh_wide$sex=="F",]) # 217 females

# Prevalence of adults, males, and females with cysts
tab0 <- table(lh_wide$sex, !is.na(lh_wide$cyst))
sum(tab0[,2])/sum(tab0)	# Overall period prevalence = 0.14
sum(tab0["M",2])/sum(tab0["M",])	# Male period prevalence = 0.1
sum(tab0["F",2])/sum(tab0["F",])	# Female period prevalence = 0.18

# G-test for difference between Guassa and SMNP prevalences
likelihood.test(rbind(tab0["F",], c(68, 31))) # G = 6.74, p < 0.001
likelihood.test(rbind(tab0["M",], c(49, 19))) # G = 11.15, p < 0.001)

# G-test for independence between cysts and sex
likelihood.test(tab0) # G = 5.05, p < 0.05

###############
## Mortality #
#############

# Exclude individuals who disappeared 
sum(!is.na(lh_wide$dod[lh_wide$sex=="M"])) # 13 males died
sum(is.na(lh_wide$dod[lh_wide$sex=="M"])) # 157 males did not die
sum(lh_wide$sex=="F" & lh_wide$dis==1)/sum(lh_wide$sex=="F") # 0 of 216 females disappeared
sum(lh_wide$sex=="M" & lh_wide$dis==1)/sum(lh_wide$sex=="M") # 36 of 170 (21.2% of males disappeared)
lh_wide2 <- lh_wide[!lh_wide$dis==1,-which(names(lh_wide)=="dis")]

# New number of adults, males, and females 
nrow(lh_wide2) # 351 individuals
sum(lh_wide2$sex=="M") # 134 males
sum(lh_wide2$sex=="F") # 217 females

# Number/proportion of adults dying with/without cysts over study period
tab1 <- table(!is.na(lh_wide2$dod), !is.na(lh_wide2$cyst))
tab1[2,]/colSums(tab1)	# 15.7% without cysts die; 60.8% with cysts die
tab2 <- table(!is.na(lh_wide2$dod[lh_wide2$sex=="M"]), !is.na(lh_wide2$cyst[lh_wide2$sex=="M"]))
tab2[2,]/colSums(tab2)	# 8.2% without cysts die; 25% with cysts die
tab3 <- table(!is.na(lh_wide2$dod[lh_wide2$sex=="F"]), !is.na(lh_wide2$cyst[lh_wide2$sex=="F"]))
tab3[2,]/colSums(tab3) # 21.3% without cysts die; 71.8% with cysts die

# G-test for association between cysts and deaths over study period
likelihood.test(tab1) # Highly significant across combined sexes (p < 0.001)
likelihood.test(tab2) # Approaching significance in males (p = 0.1)
likelihood.test(tab3) # Highly significant in females  (p < 0.001)

################################################
## Barplot comparing prevalence and mortality #
##############################################

# Barplots for mortality among males/females with/without cysts
quartz()
layout(matrix(1:2, 1, 2))
bplotF <- barplot(rev(tab2[2,]/colSums(tab2)), col=c("dimgray", "lightgray"), ylim=c(0,1), xlab="Males", names.arg=c("Cyst (3/12)", "No cyst (10/122)"), ylab="% adults who died during study")
errs.c2 <- prop.test(tab2[2,2], colSums(tab2)[2])
errs.n2 <- prop.test(tab2[2,1], colSums(tab2)[1])
arrows(x0=c(bplotF), y0=c(errs.c2$conf.int[1], errs.n2$conf.int[1]), x1=c(bplotF), y1=c(errs.c2$conf.int[2], errs.n2$conf.int[2]), length=0.07, angle=90, code=3)
mtext("A", line=2, cex=1.5, adj=0)
bplotM <- barplot(rev(tab3[2,]/colSums(tab3)), col=c("dimgray", "lightgray"), ylim=c(0,1), xlab="Females", ylab="", names.arg=c("Cyst (28/39)", "No cyst (38/178)"))
errs.c <- prop.test(tab3[2,2], colSums(tab3)[2])
errs.n <- prop.test(tab3[2,1], colSums(tab3)[1])
arrows(x0=c(bplotM), y0=c(errs.c$conf.int[1], errs.n$conf.int[1]), x1=c(bplotM), y1=c(errs.c$conf.int[2], errs.n$conf.int[2]), length=0.07, angle=90, code=3)
mtext("B", line=2, cex=1.5, adj=0)

#################################################
## Line plot showing individual life histories #
###############################################

# Subset to individuals that ever had cysts
lp <- lh_wide[!is.na(lh_wide$cyst),]
lp <- lp[-which(lp$name=="Delia"),] 
nrow(lp) # 56 individuals

# Sort first by sex and then by length of observation time
lp <- lp[order(factor(lp$sex, levels=c("F", "M")), lp$stop, decreasing=FALSE),]

# Create plot window
quartz()
par(mar=c(5, 6, 4, 2))
plot(0:(nrow(lp)+1) ~ seq(0, 25, length.out=length(0:(nrow(lp)+1))), data=lp, type="n", yaxt="n", ylab="", xlab="Time (years)")
ynames <- c(lp$name[lp$sex=="F"], "", lp$name[lp$sex=="M"])
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
# SURVIVAL ANALYSIS FOR LIFE HISTORY DATA (use lh_long)
#################################################################################################################

###############
## Data prep #
#############

# Male data
mdat <- lh_long[lh_long$sex=="M",]
fdat <- lh_long[lh_long$sex=="F",]

# Number censored
nmales <- length(unique(mdat$name))
nfemales <- length(unique(fdat$name))
(nmales - sum(mdat$death))/nmales # 92.4% censored
(nfemales - sum(fdat$death))/nfemales # 72.2% censored

# Number with cysts
sum(fdat$cyst) # 38
sum(mdat$cyst) # 17

################
## Cox models #
##############

# Cox proportional hazards model (Model 1)
coxM <- coxph(Surv(start, stop, death) ~ cyst, data=mdat)
coxF <- coxph(Surv(start, stop, death) ~ cyst, data=fdat)
summary(coxM) # Cyst not significant (exp(coef) = 1.81, p = 0.41)
summary(coxF) # Cyst significant (exp(coef) = 5.77, p < 0.001 ***)

# Test proportional hazards assumption (model 1)
cox.zph(coxM, transform="identity") # Check
cox.zph(coxF, transform="identity") # p < 0.001 ***, Violated!

# Cox proportional hazards with a time transform (Model 2)
coxM2 <- coxph(Surv(start, stop, death) ~ cyst + tt(cyst), data=mdat, tt=function(x, t, ...) x*t)
coxF2 <- coxph(Surv(start, stop, death) ~ cyst + tt(cyst), data=fdat, tt=function(x, t, ...) x*t)
summary(coxM2) # Cyst approaching significance, p = 0.06
summary(coxF2) # Cyst and time transform significant, p < 0.001 ***

######################################
## Plot scaled Schoenfeld residuals #
####################################

# Residual plots showing the nature of the proportional hazard violation
layout(matrix(1:2, 1, 2, byrow=T))
plot(cox.zph(coxM, transform="identity"))
mtext("A", line=2, cex=1.5, adj=0)
plot(cox.zph(coxF, transform="identity"))
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
plot(ym, ylim=c(-10,10), xlim=c(4, 15), xlab="Time", ylab="Log hazard ratio", lwd=2.5)
lines(s, ym.l(s))
lines(s, ym.u(s))
abline(h=coef(coxM)[1], lty=2)
abline(h=0, col="gray", lty=3)
mtext("A", line=2, cex=1.5, adj=0)
plot(yf, ylim=c(-10,10), xlim=c(4, 24), xlab="Time", ylab="", lwd=2.5)
lines(s, yf.l(s))
lines(s, yf.u(s))
abline(h=0, col="gray", lty=3)
abline(h=coef(coxF)[1], lty=2)
mtext("B", line=2, cex=1.5, adj=0)

#################################################################################################################
# END
#################################################################################################################

