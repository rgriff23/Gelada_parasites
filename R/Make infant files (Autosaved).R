#################################################################################################################
# PREPARATIONS
#################################################################################################################

# Load packages
library(survival)
library(plyr)
library(Deducer)

# Read data into R
data <- read.csv("~/Desktop/GitHub/Gelada_parasites/Data/Infant life history.csv", header=TRUE, stringsAsFactors=FALSE)

# Trim unwanted columns
data <- data[,c("NAME", "SEX", "DOB", "EndDate", "Dead..1..or.Gone.ALIVE..0.", "MomCyst", "TAKEOVER.DATE.1", "MOM", "Died.W.Mom")]
names(data) <- c("name", "sex", "dob", "end", "death", "cyst", "takeover", "mom", "momdied")
data$momdied <- ifelse(data$momdied=="1", 1, 0)
data$momdied[is.na(data$momdied)] <- 0

# Convert time data to POSIXct
data$dob <- as.POSIXct(data$dob, format="%m/%d/%y")
data$end <- as.POSIXct(data$end, format="%m/%d/%y")
data$takeover <- as.POSIXct(data$takeover, format="%m/%d/%y")

# Trim unwanted rows
data <- data[data$sex %in% c("M", "F"),]
data <- data[data$mom != "UN",]
data <- data[data$cyst != "UN",]
data <- data[!is.na(data$dob),]
data <- data[!is.na(data$end),]

#################################################################################################################
# CHECKING DATA
#################################################################################################################

# Data is not NA
sum(is.na(data$mom))
sum(is.na(data$cyst))
sum(is.na(data$dob))
sum(is.na(data$end))

# Time data makes sense (dob < takeover < end)
sum((data$takeover - data$dob) < 0, na.rm=T)
sum((data$end - data$dob) < 0, na.rm=T)
sum((data$end - data$takeover) < 0, na.rm=T)

# Summary stats
nrow(data) # 375 infants
table(data$sex) # 203 females, 172 males
table(data$cyst) # 339 to non-cyst moms, 36 to cyst moms
table(is.na(data$takeover)) # 137 experienced takeovers, 238 did not
table(data$death) # 105 were censored, 270 died
length(unique(data$mom)) # 170 unique mothers
median(table(data$mom)) # median number of births per mother is 2

#################################################################################################################
# CREATE WIDE-FORM DATA FOR KAPLAN-MEIER ANALYSIS
#################################################################################################################

# Number of years for right truncation
rt <- 3

# Create wide form data
infant.wide <- data.frame()
for (i in 1:nrow(data)) {
	stop <- difftime(data[i,"end"], data[i,"dob"], units="days")
	if (stop > 365.25*rt) {stop <- 365.25*rt}
	if (stop < 365.25*rt & data[i,"death"]==1) {death <- 1} else {death <- 0}
	new <- data.frame(name=data[i,"name"], sex=data[i,"sex"], start=0, stop=stop, cyst=data[i,"cyst"], death=death, mom=data[i,"mom"], momdied=data[i,"momdied"])
	infant.wide <- rbind(infant.wide, new)
}

# Plot histograms and median lifespans of infants of cyst and non-cyst moms
quartz()
layout(matrix(1:2,1,2))
cyst.inf <- ddply(infant.wide[infant.wide$cyst==1,], .(name), function(x) {tail(x$stop,1)})[,2]
ncyst.inf <- ddply(infant.wide[infant.wide$cyst==0,], .(name), function(x) {tail(x$stop,1)})[,2]
hist(cyst.inf/365.25, main="Distribution of infant lives (cyst)", xlab="Life in years")
abline(v=median(cyst.inf)/365.25, col="red")
hist(ncyst.inf/365.25, main="Distribution of infant lives (no cyst)", xlab="Life in years")
abline(v=median(ncyst.inf)/365.25, col="red")

# Write file
write.csv(infant.wide, file="~/Desktop/GitHub/Gelada_parasites/Data/inf_wide.csv", row.names=FALSE)

# Plot KM survival curves
quartz()
plot(survfit(Surv(start, stop/365.25, death) ~ cyst, data=infant.wide))

# Test for significant difference in KM survival estimates
survdiff(Surv(stop, death) ~ cyst, data=infant.wide) # Difference is significant

# Test for significant difference in KM survival estimates, dropping cyst moms that died
survdiff(Surv(stop, death) ~ cyst, data=infant.wide[infant.wide$momdied==0,]) # Difference not significant

#################################################################################################################
# CREATE LONG-FORM DATA FOR COX MODELS
#################################################################################################################

# Number of years for right truncation
rt <- 3

# Create long form data
infant.long <- data.frame()
for (i in 1:nrow(data)) {
	name <- data[i,"name"]
	sex <- data[i,"sex"]
	start <- 0
	stop <- difftime(data[i,"end"], data[i,"dob"], units="days")
	if (stop > 365.25*rt) {stop <- 365.25*rt}
	if (stop < 365.25*rt & data[i,"death"]==1) {death <- 1} else {death <- 0}
	cyst <- data[i,"cyst"]
	mom <- data[i,"mom"]
	takeover <- 0
	if (!is.na(data[i,"takeover"])) {
		takeovertime <- difftime(data[i,"takeover"], data[i,"dob"], units="days")
		if (takeovertime < 365.25*rt) {
			name <- c(name, name)
			sex <- c(sex, sex)
			start <- c(start, takeovertime)
			stop <- c(takeovertime, stop)
			death <- c(0, death)
			cyst <- c(cyst, cyst)
			mom <- c(mom, mom)
			takeover <- c(0, 1)	
		}
	}
	new <- data.frame(name=name, sex=sex, start=start, stop=stop, cyst=cyst, death=death, mom=mom, takeover=takeover)
	infant.long <- rbind(infant.long, new)
}

# Write file
write.csv(infant.long, file="~/Desktop/GitHub/Gelada_parasites/Data/inf_long.csv", row.names=FALSE)

# Regular Cox model 
coxI0 <- coxph(Surv(start, stop, death) ~ cyst, data=infant.long)
summary(coxI0) # p < 0.001, HR = 2.49
cox.zph(coxI0) # not significant

# Regular Cox model controlling for takeover
coxI1 <- coxph(Surv(start, stop, death) ~ cyst + takeover, data=infant.long)
summary(coxI1) # cyst p < 0.01, HR = 2.38; takeover p < 0.001, HR = 4.76
cox.zph(coxI1) # significant time varying effect of takeover
plot(cox.zph(coxI1)) # effect of takeover declines with infant age

# Extended Cox model controlling for takeover with a time varying effect
coxI2 <- coxph(Surv(start, stop, death) ~ cyst + takeover + tt(takeover), data=infant.long, tt=function(x, t, ...) x*t)
summary(coxI2) # cyst p < 0.001, HR = 2.46; takeover p < 0.001, HR = 12.92; tt(takeover) p < 0.01, HR = 0.997

# Extended Cox model controlling for takeover with a time varying effect and clustering by mom ID
coxI3 <- coxph(Surv(start, stop, death) ~ cluster(mom) + cyst + takeover + tt(takeover), data=infant.long, tt=function(x, t, ...) x*t)
summary(coxI3) cyst p < 0.01, HR = 2.46; takeover p < 0.001, HR = 12.92; tt(takeover) p < 0.01, HR = 0.997

#################################################################################################################
# END
#################################################################################################################

