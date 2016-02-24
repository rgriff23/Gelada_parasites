#########################################################################
# PREPARATIONS
#########################################################################

# Load plyr package for ddply and MASS for neg binom GLM
library(plyr)
library(MASS)
library(lmtest)

# Read data into R
data <- read.csv("~/Dropbox/Geladas/data_021116.csv", header=TRUE, stringsAsFactors=FALSE)

# Convert time data to POSIXct
data$Smp.date <- as.POSIXct(data$Smp.date, format="%m/%d/%Y")
data$DOB <- as.POSIXct(data$DOB, format="%m/%d/%Y")

# Convert POSITIVE/negative to 1/0
data$Results[data$Results=="POSITIVE"] <- 1
data$Results[data$Results=="negative"] <- 0

# Convert cyst/no cyst to 1/0
data$Cyst <- ifelse(data$Cyst=="Y", 1, 0)

# Subset to SK monkeys with recorded birthdates and non-borderline results
data.sk <- data[data$Site=="SK" & !is.na(data$DOB) & data$Results %in% c(0,1),]

#########################################################################
# CHECKING AND SUBSETTING DATA
#########################################################################

# Checking data for SK
x <- ddply(data.sk, .(Name), function(x) {
  data.frame(
    Smps=length(unique(x$Smp.date)),
    Sex=x$Sex[1], 
    DOB=x$DOB[1])
}) # all indivs in SK
#nrow(x) # 191 individuals
#table(x$Sex) # 107 females, 84 males
#hist(x$DOB, breaks=10) # distribution of birthdates

# SK subsetted to indivs with >1 dates sampled
x2 <- x[x$Smps > 1,] 
#nrow(x2) # 107 individuals
#table(x2$Sex) # 57 females, 50 males
#hist(x2$DOB, breaks=10) # distribution of birthdates

# Probability that SK indivs were positive at some point
x3 <- ddply(data.sk, .(Name), function(x) {
  x <- x[!duplicated(x$Smp.date),]
  pos <- ifelse (1 %in% x$Results, 1, 0)
  numpos <- sum(x$Results=="1")
  data.frame(
    Pos=pos,
    NumPos=numpos,
    Cyst=unique(x$Cyst),
    SmpDates=length(x$Smp.date), 
    SmpSpan=c(x$Smp.date[nrow(x)] - x$Smp.date[1]),
    Sex=x$Sex[1], 
    Age=x$Age[1],
    DOB=x$DOB[1])
})
#table(x3$Pos) # 142 negative, 49 positive at least once

# Probability of a negative sample after a positive sample
# Exclude individuals who were sampled only once
# Or who tested negative until their last sample
x4 <- ddply(data.sk[data.sk$Name %in% x2$Name,], .(Name), function(x) {
  x <- x[!duplicated(x$Smp.date),]
  x <- x[order(x$Smp.date),]
  keep <- ifelse (1 %in% x$Results[-nrow(x)], 1, 0)
  trans <- ifelse (grepl("10", paste(x$Results, collapse="")), 1, 0)
  data.frame(Keep=keep,
             Trans=trans,
             Cyst=unique(x$Cyst),
             SmpDates=length(unique(x$Smp.date)), 
             SmpSpan=c(x$Smp.date[nrow(x)] - x$Smp.date[1]),
             Sex=x$Sex[1], 
             Age=x$Age[1],
             DOB=x$DOB[1])
})
# nrow(x4) # 107 individuals

#########################################################################
# VISUALIZE OD VALUE DISTRIBUTIONS FOR CYST AND NON-CYST INDIVS
#########################################################################

# Define bar heights and combine into matrix for barplot
breaks <- seq(0, 6, by=0.1)
cyst0 <- hist(log(data.sk$ODV[data.sk$Cyst==0] + 16), breaks=breaks, plot=F)
cyst1 <- hist(log(data.sk$ODV[data.sk$Cyst==1] + 16), breaks=breaks, plot=F)
cyst01 <- rbind(cyst0$counts, cyst1$counts)

# Stacked barplot
quartz()
barplot(cyst01, space=0, legend.text=c("No Cyst", "Cyst"), col=c("gray28", "gray87"), xlab="log(OD Value)", ylab="Count (stacked)", ylim=c(0,50))
axis(side=1, at=(0:6)*10, labels=0:6)
abline(v=log(25.56 - 5 + 16)*10, col="red", lwd=2)
abline(v=log(25.56 + 5 + 16)*10, col="red", lwd=2)

# Which individual has a cyst and lies below the cutoff
data.sk$Name[data.sk$Cyst==1 & data.sk$ODV<37.5] # Mary

#########################################################################
# VISUALIZE SAMPLE RESULTS AND DATES FOR ALL SK INDIVIDUALS
#########################################################################

# Create data for figures
fig1 <- data.sk[data.sk$Name %in% x2$Name,]
fig1m <- fig1[fig1$Sex == "M",]
fig1f <- fig1[fig1$Sex == "F",]

# Names ordered by age
m.names <- ddply(fig1m, .(Name), function(x){
  x <- x[order(x$Smp.date),]
  data.frame(x$DOB[1], x$Age[1])
})
f.names <- ddply(fig1f, .(Name), function(x){
  x <- x[order(x$Smp.date),]
  data.frame(x$DOB[1], x$Age[1])
})
m.names <- m.names[order(m.names[,2], decreasing=T),]
m.names[,3] <- as.character(m.names[,3])
f.names <- f.names[order(f.names[,2], decreasing=T),]
f.names[,3] <- as.character(f.names[,3])
m.names[m.names[,3]=="I",3] <- paste(m.names[m.names[,3]=="I",3], 1:sum(m.names[,3]=="I"), sep="")
m.names[m.names[,3]=="J",3] <- paste(m.names[m.names[,3]=="J",3], 1:sum(m.names[,3]=="J"), sep="")
m.names[m.names[,3]=="SA",3] <- paste(m.names[m.names[,3]=="SA",3], 1:sum(m.names[,3]=="SA"), sep="")
m.names[m.names[,3]=="A",3] <- paste(m.names[m.names[,3]=="A",3], 1:sum(m.names[,3]=="A"), sep="")
f.names[f.names[,3]=="I",3] <- paste(f.names[f.names[,3]=="I",3], 1:sum(f.names[,3]=="I"), sep="")
f.names[f.names[,3]=="J",3] <- paste(f.names[f.names[,3]=="J",3], 1:sum(f.names[,3]=="J"), sep="")
f.names[f.names[,3]=="SA",3] <- paste(f.names[f.names[,3]=="SA",3], 1:sum(f.names[,3]=="SA"), sep="")
f.names[f.names[,3]=="A",3] <- paste(f.names[f.names[,3]=="A",3], 1:sum(f.names[,3]=="A"), sep="")

# Two panel plot
quartz()
layout(matrix(1:2, 1, 2))

# Males
plot(c(1:10)~c(1:10), ylim=c(0,nrow(m.names)), xlim=c(min(fig1m$Smp.date), max(fig1m$Smp.date)), type="n", xaxt="n", yaxt="n", xlab="Sample date", ylab="Individual", main="Males (n = 50)")
ddply(fig1[fig1$Sex=="M",], .(Name), function(x) {
  x <- x[!duplicated(x$Smp.date),]
  x <- x[order(x$Smp.date),]
  cols <- ifelse(x$Results==1, "black", "white")
  pos <- which(m.names[,1]==x$Name[1])
  #points(rep(pos, nrow(x)) ~ Smp.date, data=x, type="l", col="gray")
  abline(h=pos, col="gray", lty=3)
  points(rep(pos, nrow(x)) ~ Smp.date, data=x, pch=21, bg=cols, cex=0.5)
})
axis(side=2, at=1:nrow(m.names), labels=m.names[,3], las=2, cex.axis=0.3, tcl=-0.2, mgp=c(3,0.3,0))
xdates <- c("Aug 2014", "Sep 2014", "Oct 2014", "Apr 2015", "May 2015", "June 2015")
xpos <- as.POSIXct(c("2014/08/01", "2014/09/01", "2014/10/01", "2015/04/01", "2015/05/01", "2015/06/01"), format="%Y/%m/%d")
axis(side=1, at=xpos, labels=FALSE)
text(x=xpos, y=rep(-7, length(xpos)), labels=xdates, cex=0.75, srt=45, pos=1, xpd=T)

# Add legend
legend(1402047062, 63.27607, legend=c("Positive", "Negative"), pch=21, pt.bg=c("black", "white"), xpd=TRUE, cex=0.75)

# Females
plot(c(1:10)~c(1:10), ylim=c(0,nrow(f.names)), xlim=c(min(fig1f$Smp.date), max(fig1f$Smp.date)), type="n", xaxt="n", yaxt="n", xlab="Sample date", ylab="", main="Females (n = 57)")
ddply(fig1[fig1$Sex=="F",], .(Name), function(x) {
  x <- x[!duplicated(x$Smp.date),]
  x <- x[order(x$Smp.date),]
  cols <- ifelse(x$Results==1, "black", "white")
  pos <- which(f.names[,1]==x$Name[1])
  #points(rep(pos, nrow(x)) ~ Smp.date, data=x, type="l", col="gray")
  abline(h=pos, col="gray", lty=3)
  points(rep(pos, nrow(x)) ~ Smp.date, data=x, pch=21, bg=cols, cex=0.5)
})
axis(side=2, at=1:nrow(f.names), labels=f.names[,3], las=2, cex.axis=0.3, tcl=-0.2, mgp=c(3,0.3,0))
axis(side=1, at=xpos, labels=FALSE)
text(x=xpos, y=rep(-7, length(xpos)), labels=xdates, cex=0.75, srt=45, pos=1, xpd=T)

#########################################################################
# VISUALIZE INFECTION/RECOVERY RATE BY AGE AND SEX
#########################################################################

# Bar heights for infection
indivs <- table(x3$Sex, x3$Age)[,c(2,3,4,1)]
indivs <- cbind(indivs, rowSums(indivs))
positive <- table(x3[x3$Pos==1,"Sex"], x3[x3$Pos==1,"Age"])[,c(2,3,4,1)]
positive <- cbind(positive, rowSums(positive))
prop.pos <- positive/indivs
colnames(prop.pos) <- colnames(positive) <- colnames(indivs) <- c("Infant", "Juvenile", "Subadult", "Adult", "Total")
err.upper <- c()
err.lower <- c()
for (i in 1:10) {
  test <- prop.test(positive[i], indivs[i])
  err.upper[i] <- test$conf.int[2]
  err.lower[i] <- test$conf.int[1]
}

# Barplot with error bars for proportion with at least one positive
quartz()
pp.bplot <- barplot(prop.pos, ylim=c(0, 1), legend.text=c("Female", "Male"), axis.lty=1, beside=T, ylab="Proportion infected", xlab="Age category")
arrows(x0=c(pp.bplot), y0=err.lower, x1=c(pp.bplot), y1=err.upper, length=0.07, angle=90, code=3)

#########################################################################
# MODELS
#########################################################################

# Probability of having a cyst (n=191)
mod0 <- glm(Cyst ~ Sex + DOB, data=x3, family="binomial")
summary(mod0) # DOB p<0.01

# Probability of testing positive in at least one sample (n=191)
mod1 <- glm(Pos ~ Sex + DOB + SmpDates, data=x3, family="binomial")
summary(mod1) # DOB p<0.05, SmpDates p<0.01
mod2 <- glm(Pos ~ Sex + DOB + SmpSpan, data=x3, family="binomial")
summary(mod2) # DOB<0.01, SmpSpan p<0.1

# Number of positive samples (poisson model, n=191)
mod3.0 <- glm(NumPos ~ Sex + DOB + SmpDates, data=x3, family="poisson")
summary(mod3.0) # SexM p<0.1, DOB p<0.001, SmpDates p<0.001
mod4.0 <- glm(NumPos ~ Sex + DOB + SmpSpan, data=x3, family="poisson")
summary(mod4.0) # SexM p<0.01, DOB p<0.001, SmpSpan p<0.05

# Number of positive samples (neg binomial model, n=191)
mod3.1 <- glm.nb(NumPos ~ Sex + DOB + SmpDates, data=x3)
summary(mod3.1) # DOB p<0.001, SmpDates p<0.001
mod4.1 <- glm.nb(NumPos ~ Sex + DOB + SmpSpan, data=x3)
summary(mod4.1) # SexM p<0.05, DOB p<0.001, SmpSpan p<0.05

# Compare poisson and negative binomial
lrtest(mod3.1, mod3.0)
lrtest(mod4.1, mod4.0)

# Probability of testing negative after testing positive (n=32)
mod5 <- glm(Trans ~ Sex + DOB + SmpDates, data=x4[x4$Keep==1,], family="binomial")
summary(mod5) # nothing significant
mod6 <- glm(Trans ~ Sex + DOB + SmpSpan, data=x4[x4$Keep==1,], family="binomial")
summary(mod6) # SmpSpan p<0.1





