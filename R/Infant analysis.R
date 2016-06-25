#################################################################################################################
# PREPARATIONS
#################################################################################################################

# Load package for survival analysis
library(survival) 

# Read data into R
inf_wide <- read.csv("~/Desktop/GitHub/Gelada_parasites/Data/inf_wide.csv", header=TRUE, stringsAsFactors=FALSE)
inf_long <- read.csv("~/Desktop/GitHub/Gelada_parasites/Data/inf_long.csv", header=TRUE, stringsAsFactors=FALSE)

# Summary stats
nrow(inf_wide) # 375 infants
table(inf_wide $sex) # 203 females, 172 males
table(inf_wide $cyst) # 339 to non-cyst moms, 36 to cyst moms
sum(inf_long$takeover) # 107 experienced takeovers, 268 did not
table(inf_wide$death) # 288 were censored, 87 died
length(unique(inf_wide$mom)) # 170 unique mothers
median(table(inf_wide$mom)) # median number of births per mother is 2

#################################################################################################################
# KAPLAN-MEIER ANALYSIS (use inf_wide)
#################################################################################################################

# Test for significant difference in KM survival estimates
survdiff(Surv(stop, death) ~ cyst, data=inf_wide) # Difference is significant

# Test for significant difference in KM survival estimates, dropping cyst moms that died
survdiff(Surv(stop, death) ~ cyst, data=infant.wide[inf_wide$momdied==0,]) # Difference not significant

# Plot KM survival curves
layout(matrix(1:2, 1, 2))
plot(survfit(Surv(start, stop/365.25, death) ~ cyst, data=inf_wide), xlab="Time (years)", ylab="Proportion alive", conf.int=T, col=c("skyblue", "darkblue"), mark.time=TRUE, lwd=c(1.5, 2))
legend("bottomleft", legend=c("Non-cyst", "Cyst"), col=c("skyblue", "darkblue"), lwd=2)
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
momalive <- infant.wide[inf_wide$momdied==0,"name"]
coxI0b <- coxph(Surv(start, stop, death) ~ cyst, data=inf_long[inf_long$name%in%momalive,])
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