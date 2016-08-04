#################################################################################################################
# PREPARATIONS
#################################################################################################################

# Load package for survival analysis
library(survival) 

# Read data into R
inf_wide <- read.csv("~/Desktop/GitHub/Gelada_parasites/Data/inf_wide.csv", header=TRUE, stringsAsFactors=FALSE)
inf_long <- read.csv("~/Desktop/GitHub/Gelada_parasites/Data/inf_long.csv", header=TRUE, stringsAsFactors=FALSE)

# Summary stats for full data
nrow(inf_wide) # 375 infants
nrow(inf_wide[inf_wide$cyst==0,]) # 339 infants born to
length(unique(inf_wide[inf_wide$cyst==0,"mom"])) # 154 non-cyst mothers
nrow(inf_wide[inf_wide$cyst==1,]) # 36 infants born to 
length(unique(inf_wide[inf_wide$cyst==1,"mom"])) # 27 cyst mothers
median(table(inf_wide$mom)) # median number of births per mother is 2
max(table(inf_wide$mom)) # maximum number of births per mother is 7

# Summary stats for data excluding infants who died with moms
inf_wide2 <- inf_wide[inf_wide$momdied==0,]
nrow(inf_wide2[inf_wide2$cyst==0,]) # 325 infants born to
length(unique(inf_wide2[inf_wide2$cyst==0,"mom"])) # 149 non-cyst mothers
nrow(inf_wide2[inf_wide2$cyst==1,]) # 149 infants born to 
length(unique(inf_wide2[inf_wide2$cyst==1,"mom"])) # 17 cyst mothers


#################################################################################################################
# KAPLAN-MEIER ANALYSIS (use inf_wide)
#################################################################################################################

# Test for significant difference in KM survival estimates
survdiff(Surv(stop, death) ~ cyst, data=inf_wide) # Difference is significant

# Test for significant difference in KM survival estimates, dropping cyst moms that died
survdiff(Surv(stop, death) ~ cyst, data=infant.wide[inf_wide$momdied==0,]) # Difference not significant

# Plot KM survival curves
layout(matrix(1:2, 1, 2))
plot(survfit(Surv(start, stop/365.25, death) ~ cyst, data=inf_wide), xlab="Time (years)", ylab="Proportion of infants alive", conf.int=T, col=c("skyblue", "darkblue"), mark.time=TRUE, lwd=c(1.5, 2))
legend("bottomleft", legend=c("Non-cyst mother", "Cyst mother"), col=c("skyblue", "darkblue"), lwd=2, cex=0.6)
mtext("A", line=2, cex=1.5, adj=0)
plot(survfit(Surv(start, stop/365.25, death) ~ cyst, data=inf_wide[inf_wide$momdied==0,]), xlab="Time (years)", conf.int=T, col=c("skyblue", "darkblue"), mark.time=TRUE, lwd=c(1.5, 2))
mtext("B", line=2, cex=1.5, adj=0)

#################################################################################################################
# COX PH MODELS (use inf_long)
#################################################################################################################

# Regular Cox model 
coxI0 <- coxph(Surv(start, stop, death) ~ cyst, data=inf_long)
summary(coxI0) # p < 0.001, HR = 2.49
cox.zph(coxI0) # not significant

# Regular Cox model excluding mom-infant deaths
momalive <- inf_wide[inf_wide$momdied==0,"name"]
inf_long2 <- inf_long[inf_long$name%in%momalive,]
coxI0b <- coxph(Surv(start, stop, death) ~ cyst, data=inf_long2)
summary(coxI0b) # p < 0.001, HR = 2.49
cox.zph(coxI0b) # not significant

# Regular Cox model with clustering by mom ID
coxI0c <- coxph(Surv(start, stop, death) ~ frailty(mom) + cyst, data=inf_long)
summary(coxI0c) # p < 0.001, HR = 2.49
cox.zph(coxI0c) # not significant

# Regular Cox model controlling for takeover
coxI1 <- coxph(Surv(start, stop, death) ~ cyst + takeover, data= inf_long)
summary(coxI1) # cyst p < 0.01, HR = 2.38; takeover p < 0.001, HR = 4.76
cox.zph(coxI1) # significant time varying effect of takeover
plot(cox.zph(coxI1)) # effect of takeover declines with infant age

# Extended Cox model controlling for takeover with a time varying effect
coxI2 <- coxph(Surv(start, stop, death) ~ cyst + takeover + tt(takeover), data= inf_long, tt=function(x, t, ...) x*t)
summary(coxI2) # cyst p < 0.001, HR = 2.46; takeover p < 0.001, HR = 12.92; tt(takeover) p < 0.01, HR = 0.997

# Extended Cox model controlling for takeover with a time varying effect and clustering by mom ID
coxI3 <- coxph(Surv(start, stop, death) ~ cluster(mom) + cyst + takeover + tt(takeover), data= inf_long, tt=function(x, t, ...) x*t)
summary(coxI3) # cyst p < 0.01, HR = 2.46; takeover p < 0.001, HR = 12.92; tt(takeover) p < 0.01, HR = 0.997

#################################################################################################################
# END
#################################################################################################################