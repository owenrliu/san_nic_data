#Numerical Ecology With R, Chapter 2
## Data Exploration

#### Libraries and Data Import ####
library(ade4)
library(vegan)
library(gclus)
library(cluster)
library(FD)
library(RColorBrewer)
library(labdsv)

data(doubs)
env<-doubs$env
spa<-doubs$xy
spe<-doubs$fish
key <- doubs$species

#### 2.2.2 Species Data: First Contact ####
#minimum and maximum of baundance valuesi n the whole data set range (spe)
# coutn cases for each abundance class
ab<-table(unlist(spe))
ab

#Barplot of the distribution, all species confounded
barplot(ab, las=1, xlab="Abundance Class",ylab="Frequency",col=gray(5:0/5))
# Number of absences
sum(spe==0)
#Proportion of zeros in the community data set
sum(spe==0)/(nrow(spe)*ncol(spe))

#### 2.2.3 species data: a closer look ####
# Map of the positions of the sites
#***********************************

#Create an empty frame (proportional axes 1:1, with titles)
#Geographic coordinates x and y from the spa data frame
plot(spa,asp=1,type='n',main="Site Locations",xlab='x coordinate (km)',ylab='y coordinate (km)')

#add a blue line conencting the sites (representing the river)
lines(spa,col="light blue")

#add site labels
text(spa,row.names(spa),cex=0.8,col='red')

#Add text blocks
text(50,10,"Upstream",cex=1.2,col="red")
text(30,120, "Downstream", cex=1.2, col="red")

# Maps of some fish species
#**************************

# Divide the plot window in  4 frames, 2 per row
par(mfrow=c(2,2))
plot(spa,asp=1,col='brown',cex=spe$Satr, main="Brown Trout",xlab="x coord (km)",
     ylab="y coord (km)")
lines(spa,col='light blue')

plot(spa,asp=1,col='brown',cex=spe$Thth, main="Grayling",xlab="x coord (km)",
     ylab="y coord (km)")
lines(spa,col='light blue')

plot(spa,asp=1,col='brown',cex=spe$Baba, main="Barbel",xlab="x coord (km)",
     ylab="y coord (km)")
lines(spa,col='light blue')

plot(spa,asp=1,col='brown',cex=spe$Abbr, main="Common Bream",xlab="x coord (km)",
     ylab="y coord (km)")
lines(spa,col='light blue')

# Compare species: number of occurences
#**************************************
# At how many sites does each species occur?

# Compute the number of sites where each species is pr3esent
# to sum by columns, the second argument of apply(), MARGIN, is set to 2
spe.pres <- apply(spe>0,MARGIN=2,sum)
# sort in increasing order
sort(spe.pres)
# Compute percentage frequencies
spe.relf <- 100*spe.pres/nrow(spe)
round(sort(spe.relf),1) # Round to 1 digit

#Plot the histograms
par(mfrow=c(1,2))

hist(spe.pres, main="Species Occurences",right=FALSE, las=1, xlab="Number of occurences",
     ylab="Number of species",breaks=seq(0,30,by=5),col="bisque")

hist(spe.relf,main="Species Relative Freq", right=F,las=1,xlab="Frequency of occurrences(%)",
     ylab="Number of species",breaks=seq(0,100,by=10), col="bisque")

# Compare sites: species richness
# *******************************

# Compute the number of species at each site
# To sum by rows, the second argument of apply(), MARGIN, is set to 1
sit.pres <- apply(spe>0,1,sum)
# Sort the results in increasing order
sort(sit.pres)
par(mfrow=c(1,2))

# Plot species richness vs. position of the sites along the river
plot(sit.pres,type='s',las=1,col="gray",main="Species Richness vs. \n Upstream-Downstream Gradient",
     xlab="Positions of sites along the river", ylab="Species richness")
text(sit.pres,row.names(spe),cex=0.8,col='red')

# Use geographic coordinates to plot a bubble map
plot(spa,asp=1,main="Map of Species Richness", pch=21,col='white',bg="brown",cex=5*sit.pres/max(sit.pres),
     xlab='x coord (km)',ylab='y coord (km)')
lines(spa,col='light blue')

# Finally, we can look at diversity indices
# Compute diversity indices
# *************************

N0 <- rowSums(spe>0) # Species richness
H <- diversity(spe)  # Shannon entropy
N1 <- exp(H)         # Shannon diversity number
N2 <- diversity(spe,index='inv') # Simpson diversity number
J <- H/log(N0)       # Pielou evenness
E1 <- N1/N0          # Shannon evenness (Hill's ratio)
E2 <- N2/N0          # Simpson evenness (Hill's ratio)
div <- data.frame(N0,H, N1, N2, E1, E2, J)
div

#### 2.2.4 Species Data Transformation ####
# Data transformation and standardization
# *****************************************

# Simple transformations

# Partial view of the raw data (abundance codes)
spe[1:5, 2:4]
# Transform abundances to presence-absence (1-0)
spe.pa <- decostand (spe, method='pa')
spe.pa[1:5,2:4]

# Species profiles: 2 methods; presence-absence or abundance data
# *******************************

# Scale abundances by dividing them by the max value for each species
spe.scal <- decostand(spe, 'max')
spe.scal[1:5,2:4]
# Display the maximum by column
apply(spe.scal, MARGIN=2, max)

# Scale abundances by dividing them by the species totals (relative abundance per species)
spe.relsp <- decostand(spe, "total", MARGIN=2)
spe.relsp[1:5,2:4]
# Display the sum by column
apply(spe.relsp, MARGIN=2, sum)

#Site profiles: 3 methods; p-a or abundance data
#*******************************

# Scale abundances by dividing them by the site totals
# (relative abundances, or relative frequencies, per site)
spe.rel <- decostand(spe, "total", MARGIN=1)
spe.rel[1:5,2:4]
# Display the sum of rows to determine if the scaling worked properly
apply(spe.rel,MARGIN=1, sum)

# Give a length of 1 to each row vector (Euclidean norm)
spe.norm <- decostand(spe,"normalize")
spe.norm[1:5,2:4]
# Verify the norm of row vectors
normf <- function(x) sqrt(x%*%x)
apply(spe.norm,1,normf)
# The scaling above is called the 'chord transformation'. Useful for later analyses!

# Compute relative frequencies by rows (site profiles), then square root
spe.hel <- decostand(spe, 'hellinger')
spe.hel[1:5,2:4]
# Check the norm of row vectors
apply(spe.hel,MARGIN=1,function(x)sqrt(sum(x^2)))

#Standardization of both species and sites (double profiles)
# ******************************************

# Chi-square transformation
spe.chi <- decostand(spe,"chi.square")
spe.chi[1:5,2:4]
# Check what happened to site 8 where no species were found
spe.chi[7:9,]

# Note: the Euclidean distance function applied to chi-square-transformed data produces
# a chi-quare distance matrix (see Chapter 3)

# Wisconsin standardization: abundances are first ranged by
# species maxima and then by site totals

spe.wis <- wisconsin(spe)
spe.wis[1:5,2:4]

# Boxplots of transformed abundances of a commoin species (stone loach)
# *************************************

par(mfrow=c(2,2))

boxplot(spe$Neba, sqrt(spe$Neba),log1p(spe$Neba), las=1, main= "Simple transformation", 
        names=c('raw data',"sqrt","log"), col="bisque")

boxplot(spe.scal$Neba, spe.relsp$Neba, las=1, main="Standardization by species", names=c("max","total"), col="lightgreen")

boxplot(spe.hel$Neba,spe.rel$Neba,spe.norm$Neba,las=1, main="Standardization by sites",
        names=c("Hellinger","total","norm"), col="lightblue")

boxplot(spe.chi$Neba,spe.wis$Neba,las=1,main="Double Standardization",names=c("Chi-square","Wisconsin"),col='orange')

# Plot profiles along the upstream-downstream gradient
# ****************************************************

par(mfrow=c(2,2))
plot(env$dfs, spe$Satr, type='l', col=4, main="Raw data", xlab="Distance from the source [km]", ylab="Raw abundance code")
lines(env$dfs, spe$Thth, col=3)
lines(env$dfs, spe$Baba, col="orange")
lines(env$dfs, spe$Abbr, col=2)
lines(env$dfs, spe$Neba, col=1, lty='dotted')

plot(env$dfs, spe.scal$Satr, type='l', col=4, main="Species profiles (max)", xlab="Distance from the source [km]", 
     ylab="Standardized abundance")
lines(env$dfs, spe.scal$Thth, col=3)
lines(env$dfs, spe.scal$Baba, col="orange")
lines(env$dfs, spe.scal$Abbr, col=2)
lines(env$dfs, spe.scal$Neba, col=1, lty='dotted')

plot(env$dfs, spe.hel$Satr, type='l', col=4, main="Site profiles (Hellinger)", xlab="Distance from the source [km]", 
     ylab="Standardized abundance")
lines(env$dfs, spe.hel$Thth, col=3)
lines(env$dfs, spe.hel$Baba, col="orange")
lines(env$dfs, spe.hel$Abbr, col=2)
lines(env$dfs, spe.hel$Neba, col=1, lty='dotted')

plot(env$dfs, spe.chi$Satr, type='l', col=4, main="Double profiles (Chi-square)", xlab="Distance from the source [km]", 
     ylab="Standardized abundance")
lines(env$dfs, spe.chi$Thth, col=3)
lines(env$dfs, spe.chi$Baba, col="orange")
lines(env$dfs, spe.chi$Abbr, col=2)
lines(env$dfs, spe.chi$Neba, col=1, lty='dotted')

#### 2.2.5 Environmental Data ####
# **********************************

#Bubble maps of some environmental variables
# ******************************************

par(mfrow=c(2,2))

plot(spa, asp=1, main="Altitude", pch=21, col="white", bg="red", cex=5*env$alt/max(env$alt), xlab="x", ylab="y")
lines(spa, col='light blue')

plot(spa, asp=1, main="Discharge", pch=21, col="white", bg="blue", cex=5*env$flo/max(env$flo), xlab="x", ylab="y")
lines(spa, col='light blue')

plot(spa, asp=1, main="Oxygen", pch=21, col="white", bg="green3", cex=5*env$oxy/max(env$oxy), xlab="x", ylab="y")
lines(spa, col='light blue')

plot(spa, asp=1, main="Nitrate", pch=21, col="white", bg="brown", cex=5*env$nit/max(env$nit), xlab="x", ylab="y")
lines(spa, col='light blue')

## Variation in some descriptors along the stream
# Line plots
# ***************

par(mfrow=C(2,2))
plot(env$dfs, env$alt, type='l', xlab="Distance from the source (km)", ylab="Altitude (m)", col="red", main="Altitude")
plot(env$dfs, env$flo, type="l",xlab="Distance from the source (km)",ylab="Discharge (m3/s)", col="blue", main="Discharge")
plot(env$dfs, env$oxy, type='l', xlab="Distance from the source (km)", ylab="Oxygen (mg/L)", col="green3", main="Oxygen")
plot(env$dfs, env$nit, type="l",xlab="Distance from the source (km)",ylab="Discharge (m3/s)", col="brown", main="Nitrate")

# Scatter plots for all pairs of environmental variables
# ****************************************************

# Load additional functions from an R script
#### panelutils.R ####
#
# License: GPL-2 
# Author: Francois Gillet, 23 August 2012
#
## Put Pearson, Spearman or Kendall correlations on the upper panel
panel.cor <- function(x, y, method="pearson", digits=3, cex.cor=1.2, no.col=FALSE)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, method=method)
  ra <- cor.test(x, y, method=method)$p.value
  txt <- round(r, digits)
  prefix <- ""
  if(ra <= 0.1) prefix <- "."
  if(ra <= 0.05) prefix <- "*"
  if(ra <= 0.01) prefix <- "**"
  if(ra <= 0.001) prefix <- "***"
  if(no.col)
  {
    color <- 1
    if(r < 0) { if(ra <= 0.001) sig <- 4 else sig <- 3 }
    else { if(ra <= 0.001) sig <- 2 else sig <- 1 }
  }
  else
  {
    sig <- 1
    if(ra <= 0.001) sig <- 2
    color <- 2
    if(r < 0) color <- 4
  }
  txt <- paste(txt, prefix, sep="\n")
  text(0.5, 0.5, txt, cex = cex.cor, font=sig, col=color)
}


## Put histograms on the diagonal
panel.hist <- function(x, no.col=FALSE, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  his <- hist(x, plot=FALSE)
  breaks <- his$breaks; nB <- length(breaks)
  y <- his$counts
  y <- y/max(y)
  if(no.col) rect(breaks[-nB], 0, breaks[-1], y, col="gray", ...)
  else rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}


## Add black lowess curves to scatter plots
panel.smoothb <- function (x, y, col=par("col"), bg=NA, pch=par("pch"), 
                           cex=1, col.smooth="black", span=2/3, iter=3, ...) 
{
  points(x, y, pch=pch, col=col, bg=bg, cex=cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f=span, iter=iter), col=col.smooth, ...)
}


#Usage:
#pairs(num.mat, lower.panel=panel.smooth, upper.panel=panel.cor, diag.panel=panel.hist)
#pairs(num.mat, lower.panel=panel.smooth, upper.panel=panel.cor, method="kendall")

# Bivariate plots with histograms on the diagonal and smooth fitted curves
op <- par(mfrow=c(1,1), pty="s")
pairs(env, panel=panel.smooth, diag.panel = panel.hist, main= "Bivariate Plots and Histrograms and Smooth Curves")
par(op)

# Simple transformation of an environmental variable
# ************************************************

range(env$slo)
# Log-transformation of the variable 'phope' (y=ln(x))
# Compare histrograms and boxplots of raw and transformed values
par(mfrow=c(2,2))
hist(env$pho, col="bisque", right=F)
hist(log(env$pho), col="light green", right=F, main="Histogram of ln(env$pho)")
boxplot(env$pho, col="bisque", main="Boxplot of env$pho", ylab="env$pho")
boxplot(log(env$pho), col="light green", main="Boxplot of ln(env$pho)", ylab="log(env$pho")
