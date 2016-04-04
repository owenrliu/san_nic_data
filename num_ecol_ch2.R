#Numerical Ecology With R, Chapter 2
## Data Exploration

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

#### 2.2.2 Speciesw Data: First Contact ####
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
?decostand
