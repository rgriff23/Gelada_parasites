################
# PREPARATIONS #
################

# Load packages
library(survival) 
library(plyr)

# Read data into R
adult <- read.csv("https://raw.githubusercontent.com/rgriff23/Gelada_parasites/master/Data/adult.csv", header=TRUE, stringsAsFactors=FALSE)

#################
# SUMMARY STATS #
#################

# Format data so there is one row per individual
adult1 <- ddply(adult, .(name), function (x) {data.frame(name=x$name[1], sex=x$sex[1], cyst=sum(x$cyst))})

# Number of adults, males, females, and infants
nrow(adult1) # 386
nrow(adult1[adult1$sex=="M",]) # 170 males
nrow(adult1[adult1$sex=="F",]) # 216 females
nrow(infant)

# Number of cysts among adults, males, and females
sum(adult1[,"cyst"]) # 54
sum(adult1[adult1$sex=="M","cyst"]) # 16
sum(adult1[adult1$sex=="F","cyst"]) # 38

# Prevalence of cysts among adults, males, and females
sum(adult1[,"cyst"])/nrow(adult1) # 13.99%
sum(adult1[adult1$sex=="M","cyst"])/nrow(adult1[adult1$sex=="M",]) # 9.41%
sum(adult1[adult1$sex=="F","cyst"])/nrow(adult1[adult1$sex=="F",]) # 17.59%

#################
# COX PH MODELS #
#################

# Subset data for males and females
mdat <- adult[adult$sex=="M",]
fdat <- adult[adult$sex=="F",]

# Number of censored individuals
nmales <- length(unique(mdat$name))
nfemales <- length(unique(fdat$name))
(nmales - sum(mdat$death))/nmales # 92.4% censored
(nfemales - sum(fdat$death))/nfemales # 72.2% censored

# Cox proportional hazards model (Model 1)
coxM <- coxph(Surv(start, stop, death) ~ cyst, data=mdat, x=TRUE)
coxF <- coxph(Surv(start, stop, death) ~ cyst, data=fdat, x=TRUE)
summary(coxM) # Cyst not significant (exp(coef) = 1.81, p = 0.41)
summary(coxF) # Cyst significant (exp(coef) = 5.77, p < 0.001 ***)

# Test proportional hazards assumption (Model 1)
cox.zph(coxM, transform="identity") # p = 0.18
cox.zph(coxF, transform="identity") # p < 0.001 ***, Violated!

# Cox proportional hazards with a time transform (Model 2)
coxM2 <- coxph(Surv(start, stop, death) ~ cyst + tt(cyst), data=mdat, tt=function(x, t, ...) x*t, x=TRUE)
coxF2 <- coxph(Surv(start, stop, death) ~ cyst + tt(cyst), data=fdat, tt=function(x, t, ...) x*t, x=TRUE)
summary(coxM2) # Cyst approaching significance, p = 0.06
summary(coxF2) # Cyst and time transform significant, p < 0.001 ***

# Schoenfeld residual plots showing the nature of the proportional hazard violation (Figure 3)
layout(matrix(1:2, 1, 2, byrow=T))
plot(cox.zph(coxM, transform="identity"), xlab="Age (years)", main="Males")
mtext("A", line=2, cex=1.5, adj=0)
plot(cox.zph(coxF, transform="identity"), xlab="Age (years)", main="Females")
mtext("B", line=2, cex=1.5, adj=0)

# Plot log hazard ratios (Figure 4)
layout(matrix(1:2,1,2))
s <- seq(0, 25, length.out=100)
ym <- function (x) {(coef(coxM2)[1] + coef(coxM2)[2]*x)}
yf <- function (x) {(coef(coxF2)[1] + coef(coxF2)[2]*x)}
ym.l <- function (x) {log(summary(coxM2)$conf.int["cyst","lower .95"]) + log(summary(coxM2)$conf.int["tt(cyst)","lower .95"])*x}
ym.u <- function (x) {log(summary(coxM2)$conf.int["cyst","upper .95"]) + log(summary(coxM2)$conf.int["tt(cyst)","upper .95"])*x}
yf.l <- function (x) {log(summary(coxF2)$conf.int["cyst","lower .95"]) + log(summary(coxF2)$conf.int["tt(cyst)","lower .95"])*x}
yf.u <- function (x) {log(summary(coxF2)$conf.int["cyst","upper .95"]) + log(summary(coxF2)$conf.int["tt(cyst)","upper .95"])*x}
plot(ym, ylim=c(-10,10), xlim=c(4, 15), ylab="Log hazard ratio", lwd=2.5, xlab="Age (years)", main="Males", col= "turquoise4")
lines(s, ym.l(s), col= "turquoise")
lines(s, ym.u(s), col= "turquoise")
abline(h=coef(coxM)[1], lty=2)
abline(h=0, col="gray", lty=3)
mtext("A", line=2, cex=1.5, adj=0)
plot(yf, ylim=c(-10,10), xlim=c(4, 24), ylab="", lwd=2.5, xlab="Age (years)", main="Females", col= "turquoise4")
lines(s, yf.l(s), col= "turquoise")
lines(s, yf.u(s), col= "turquoise")
abline(h=0, col="gray", lty=3)
abline(h=coef(coxF)[1], lty=2)
mtext("B", line=2, cex=1.5, adj=0)

#######
# END #
#######
