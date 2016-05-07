#################################################################################################################
# PREPARATIONS
#################################################################################################################

# Load packages
library(survival)
library(plyr)
library(Deducer)

# Read data into R
data <- read.csv("~/Desktop/GitHub/Gelada_parasites/Data/survival_050416.csv", header=TRUE, stringsAsFactors=FALSE)

# Convert time data to POSIXct
data$DATE <- as.POSIXct(data$DATE, format="%m/%d/%y")
hist(data$DATE, breaks=100, col="gray")

# Convert missing sex data to NA
data[data$SEX %in% c("?", ""),"SEX"] <- NA

#################################################################################################################
# CHECKING DATA
#################################################################################################################

# Check 1 DOB and 1 START and 1 END for each individual?
test <- ddply(data, .(NAME), function(x) {
	data.frame(
	sum(x$EVENT=="DOB"),
	sum(x$EVENT=="START"),
	sum(x$EVENT=="END")
	)
	})
names(test) <- c("NAME", "DOB", "START","END")
drop1 <- test[test[,2]!=1 | test[,3]!=1 | test[,4]!=1,"NAME"]	 # drop
drop1
data <- data[!(data$NAME %in% drop1),]

# Check DOB and START and END are not NA
test1 <- data[data$EVENT=="DOB", c("NAME", "DATE")]
test1 <- test1[is.na(test1[,2]),] # eliminate these individuals from analysis
test2 <- data[data$EVENT=="START", c("NAME", "DATE")]
test2 <- test2[is.na(test2[,2]),] # eliminate these individuals from analysis
test3 <- data[data$EVENT=="END", c("NAME", "DATE")]
test3 <- test3[is.na(test3[,2]),] # eliminate these individuals from analysis
nrow(test1)	# 10 missing DOB
nrow(test2)	# 4 missing START
nrow(test3)	# 58 missing END

# Drop data with DOB=NA and START=NA and END=NA
drop2 <- c(test1$NAME, test2$NAME, test3$NAME)
data <- data[!(data$NAME %in% drop2),]

# Confirm DOB, START, Cyst, Cyst End, and DOD are not NA
sum(is.na(data[data$EVENT=="DOB","DATE"]))	# Check
sum(is.na(data[data$EVENT=="START","DATE"]))	# Check
sum(is.na(data[data$EVENT=="Cyst","DATE"]))	# Check
sum(is.na(data[data$EVENT=="Cyst End","DATE"]))	# Check
sum(is.na(data[data$EVENT=="DOD","DATE"]))	# Check

# Check that order of DOB, START, Cyst, Cyst End, and END make sense
test3 <- ddply(data, .(NAME), function(x) {
	dob <- start <- cyst <- end <- c()
	dob <- x[x$EVENT=="DOB", "DATE"]
	start <- x[x$EVENT=="START", "DATE"]
	if ("Cyst" %in% x$EVENT) {
		cyst <- x[x$EVENT=="Cyst", "DATE"]
	} else {cyst <- as.POSIXct(NA)}
	if ("Cyst End" %in% x$EVENT) {
		cystend <- x[x$EVENT=="Cyst End", "DATE"]
	} else {cystend <- as.POSIXct(NA)}
	end <- x[x$EVENT=="END", "DATE"]
	data.frame(dob, start, cyst, cystend, end)
	})
sum((test3$start - test3$dob)<0, na.rm=T)	# Check
sum((test3$cyst - test3$dob)<0, na.rm=T)	# Check
sum((test3$cystend - test3$dob)<0, na.rm=T)	# Check
sum((test3$cystend - test3$cyst)<0, na.rm=T)	# Check
sum((test3$end - test3$cyst)<0, na.rm=T)	# Check
sum((test3$end - test3$cystend)<0, na.rm=T)	# Check
sum((test3$end - test3$dob)<0, na.rm=T)	# Check
sum((test3$end - test3$start)<0, na.rm=T)	# Check

#################################################################################################################
# CREATE WIDE-FORM DATA FOR ESTIMATING PERIOD PREVALENCE AND PERFORMING G-TESTS
#################################################################################################################

# Longevity data, wide format (cyst is not time varying)
longevity.wide <- ddply(data, .(NAME), function(x) {
	if ("Cyst" %in% x$EVENT) {
		cyst <- x[x$EVENT=="Cyst","DATE"]
		} else {cyst <- as.POSIXct(NA)}
	if ("DOD" %in% x$EVENT) {
		dod <- x[x$EVENT=="DOD","DATE"]
		} else {dod <- as.POSIXct(NA)}
	if ("DIS" %in% x$EVENT) { 
		dis <- 1
	} else {dis <- 0}
	data.frame(
	sex <- x$SEX[1],
	dob <- x[x$EVENT=="DOB","DATE"],
	start <- x[x$EVENT=="START","DATE"],
	cyst <- cyst,
	dod <- dod,
	end <- x[x$EVENT=="END","DATE"],
	mom <- x$MOM[1],
	dis <- dis
	)
})
names(longevity.wide) <- c("name", "sex", "dob", "start", "cyst", "dod", "end", "mom", "dis")
longevity.wide$sex <- as.character(longevity.wide$sex)

# Check histogram of longevity
hist(c(longevity.wide$dod - longevity.wide$dob)/60/60/24/365, main="Longevity in years", xlab="Longevity", breaks=44, col="gray")

# Write file
write.csv(longevity.wide, file="~/Desktop/GitHub/Gelada_parasites/Data/lh_wide.csv", row.names=FALSE)

#################################################################################################################
# CREATE LONG-FORM DATA FOR SURVIVAL ANALYSIS
#################################################################################################################

# Longevity data, long format
names <- unique(data$NAME)
longevity.long <- data.frame()
for (i in 1:length(names)) {
	x <- data[data$NAME==names[i],]
	start <- difftime(x[x$EVENT=="START","DATE"], x[x$EVENT=="DOB","DATE"], units="days")
	stop <- difftime(x[x$EVENT=="END","DATE"], x[x$EVENT=="DOB","DATE"], units="days")
	if ("Cyst" %in% x$EVENT) {
		cystdate <- difftime(x[x$EVENT=="Cyst","DATE"], x[x$EVENT=="DOB","DATE"], units="days")
		if (cystdate > start) {
			start <- c(start, cystdate)
			if ("Cyst End" %in% x$EVENT) {
				cystenddate <- difftime(x[x$EVENT=="Cyst End","DATE"], x[x$EVENT=="DOB","DATE"], units="days")
				start <- c(start, cystenddate)
				cyst <- c(0, 1, 0)
				stop <- c(cystdate, cystenddate, stop)
				ifelse("DOD" %in% x$EVENT, death <- c(0,0,1), death <- c(0,0,0))
				} else {
					cyst <- c(0, 1)
					stop <- c(cystdate, stop)
					ifelse("DOD" %in% x$EVENT, death <- c(0,1), death <- c(0,0))
					}
			} else {
				cyst <- 1
				ifelse("DOD" %in% x$EVENT, death <- 1, death <- 0)
			}
		} else {
			cyst <- 0
			ifelse("DOD" %in% x$EVENT, death <- 1, death <- 0)
		}
	name <- rep(x$NAME[1], length(start))
	sex <- rep(x$SEX[1], length(start))	
	new <- data.frame(name=name, sex=sex, start=start, stop=stop, cyst=cyst, death=death)
	longevity.long <- rbind(longevity.long, new)
}

# Translate days into years
longevity.long$start <- c(longevity.long$start/365.25)
longevity.long$stop <- c(longevity.long$stop/365.25)

# Check that stop > start and remove problematic case
drop3 <- longevity.long[which(longevity.long $stop <= longevity.long $start),"name"] # Tripod is recorded with cyst at the moment of death, don't trust this
longevity.long <- longevity.long[!(longevity.long$name %in% drop3),]

# Write file
write.csv(longevity.long, file="~/Desktop/GitHub/Gelada_parasites/Data/lh_long.csv", row.names=FALSE)

#################################################################################################################
# END
#################################################################################################################
