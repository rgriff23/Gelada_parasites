################
# PREPARATIONS #
################

# Load packages
library(survival) 

# Read data into R
infant <- read.csv("~/Desktop/GitHub/Gelada_parasites/Data/infant.csv", header=TRUE, stringsAsFactors=FALSE)

#################
# SUMMARY STATS #
#################

# Number of infants and mothers with/without cysts (full data)
nrow(infant) # 375 infants
nrow(infant[infant$cyst==0,]) # 339 infants born to
length(unique(infant[infant$cyst==0,"mom"])) # 154 non-cyst mothers
nrow(infant[infant$cyst==1,]) # 36 infants born to 
length(unique(infant[infant$cyst==1,"mom"])) # 27 cyst mothers
median(table(infant$mom)) # median number of births per mother is 2
max(table(infant$mom)) # maximum number of births per mother is 7

# Number of infants and mothers with/without cysts (exclude infants that died with mothers)
infant2 <- infant[infant$momdied==0,]
nrow(infant2[infant2$cyst==0,]) # 325 infants born to non-cyst mothers
length(unique(infant2[infant2$cyst==0,"mom"])) # 149 non-cyst mothers
nrow(infant2[infant2$cyst==1,]) # 25 infants born to cyst mothers
length(unique(infant2[infant2$cyst==1,"mom"])) # 17 cyst mothers

#########################
# KAPLAN-MEIER ANALYSIS #
#########################

# Test for significant difference in KM survival estimates
survdiff(Surv(stop, death) ~ cyst, data= infant) # p < 0.001***

# Test for significant difference in KM survival estimates, dropping cyst moms that died
survdiff(Surv(stop, death) ~ cyst, data= infant[infant$momdied==0,]) # p = 0.465

# Plot KM survival curves
layout(matrix(1:2, 1, 2))
plot(survfit(Surv(start, stop/365.25, death) ~ cyst, data=infant), xlab="Age (years)", ylab="Proportion of infants alive", conf.int=T, col=c("skyblue", "darkblue"), mark.time=TRUE, lwd=c(1.5, 2))
legend("bottomleft", title="Mother status", legend=c("No cyst", "Cyst"), col=c("skyblue", "darkblue"), lwd=2, text.font=3)
mtext("A", line=2, cex=1.5, adj=0)
plot(survfit(Surv(start, stop/365.25, death) ~ cyst, data=infant[infant$momdied==0,]), xlab="Age (years)", conf.int=T, col=c("skyblue", "darkblue"), mark.time=TRUE, lwd=c(1.5, 2))
mtext("B", line=2, cex=1.5, adj=0)

#################
# COX PH MODELS #
#################

# Cox model with full data
cox1 <- coxph(Surv(start, stop, death) ~ cyst, data=infant)
summary(cox1) # p < 0.001, HR = 2.49
cox.zph(cox1) # p = 0.23

# Cox model excluding mom-infant deaths
cox2 <- coxph(Surv(start, stop, death) ~ cyst, data=infant2)
summary(cox2) # p = 0.467, HR = 1.37
cox.zph(cox2) # p = 0.471

# Cox model with clustering by mom ID
cox3 <- coxph(Surv(start, stop, death) ~ frailty(mom) + cyst, data=infant)
summary(cox3) # p < 0.001, HR = 2.49
cox.zph(cox3) # p = 0.23

#######
# END #
#######