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
drop1 # Patsy and Princess
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
	stop <- as.numeric(difftime(x[x$EVENT=="END","DATE"], x[x$EVENT=="DOB","DATE"], units="days")/365.25)
	if ("Cyst" %in% x$EVENT) {
		cyst <- as.numeric(difftime(x[x$EVENT=="Cyst","DATE"], x[x$EVENT=="DOB","DATE"], units="days"))/365.25
		if (cyst > stop) {cyst <- NA}
		} else {cyst <- NA}
	if ("DOD" %in% x$EVENT) {
		dod <- as.numeric(difftime(x[x$EVENT=="DOD","DATE"], x[x$EVENT=="DOB","DATE"], units="days"))/365.25
		} else {dod <- NA}	
	data.frame(
	sex <- x$SEX[1],
	start <- as.numeric(difftime(x[x$EVENT=="START","DATE"], x[x$EVENT=="DOB","DATE"], units="days"))/365.25,
	cyst <- cyst,
	dod <- dod,
	stop <- stop,
	mom <- x$MOM[1],
	dis <- ifelse("DIS" %in% x$EVENT, 1, 0),
	drop <- ifelse(x[x$EVENT=="END","DATE"] < as.POSIXct("07/01/08", format="%m/%d/%y"), 1, 0)
	)
})
names(longevity.wide) <- c("name", "sex", "start", "cyst", "dod", "stop", "mom", "dis", "drop")

# Check histogram of longevity
hist(longevity.wide$dod, main="Longevity in years", xlab="Longevity", col="gray")

# Left truncate at age 4
longevity.wide <- longevity.wide[longevity.wide$stop > 4,] # Drop those that didn't make it to 4
longevity.wide$start[longevity.wide$start < 4] <- 4 # Set pre-4 start dates to 4

# Write file
write.csv(longevity.wide, file="~/Desktop/GitHub/Gelada_parasites/Data/adult_wide.csv", row.names=FALSE)

#################################################################################################################
# CREATE LONG-FORM DATA FOR SURVIVAL ANALYSIS
#################################################################################################################

# Longevity data, long format
names <- unique(longevity.wide$name)
longevity.long <- data.frame()
for (i in 1:length(names)) {
	x <- data[data$NAME==names[i],]
	start <- as.numeric(difftime(x[x$EVENT=="START","DATE"], x[x$EVENT=="DOB","DATE"], units="days")/365.25)
	stop <- as.numeric(difftime(x[x$EVENT=="END","DATE"], x[x$EVENT=="DOB","DATE"], units="days")/365.25)
	if ("Cyst" %in% x$EVENT) {
		cystdate <- as.numeric(difftime(x[x$EVENT=="Cyst","DATE"], x[x$EVENT=="DOB","DATE"], units="days")/365.25)
		if (cystdate > start) {
			start <- c(start, cystdate)
			if ("Cyst End" %in% x$EVENT) {
				cystenddate <- as.numeric(difftime(x[x$EVENT=="Cyst End","DATE"], x[x$EVENT=="DOB","DATE"], units="days")/365.25)
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

# Left truncate at age 4
longevity.long$start[longevity.long$start < 4] <- 4 # Set pre-4 start dates to 4

# Right truncate at age 24
longevity.long$death[longevity.long$stop > 24] <- 0
longevity.long$stop[longevity.long$stop > 24] <- 24

# Check that stop > start and remove problematic case
longevity.long <- longevity.long[-which(longevity.long$stop <= longevity.long$start),] 

# Write file
write.csv(longevity.long, file="~/Desktop/GitHub/Gelada_parasites/Data/adult_long.csv", row.names=FALSE)

#################################################################################################################
# END
#################################################################################################################
