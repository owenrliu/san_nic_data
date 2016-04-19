## NUMERICAL ECOLOGY CHAPTER 5: UNCONSTRAINED ORDINATION ##

## While cluster analysis looks for discontinuities in a dataset, ordination extracts the main trends
## in the form of continuous axes.
## Objectives
## Learn how to choose among various ordination technique (PCA, CA, PCoA, NMDS), compute and interpret them
## Apply these techniques to the Doubs data
## Overlay the result of a cluster analysis on an ordination diagram to improve interpretation of clustering
## Interpret the sturcture in the species data using the env variables
## Write you own PCA function


#### helper functions ####
# Plot eigenvalues and percentages of variation of an ordination object
# Kaiser rule and broken stick model
# Usage:
# evplot(ev)
# where ev is a vector of eigenvalues

# License: GPL-2 
# Author: Francois Gillet, 25 August 2012

evplot <- function(ev)
{
  # Broken stick model (MacArthur 1957)
  n <- length(ev)
  bsm <- data.frame(j=seq(1:n), p=0)
  bsm$p[1] <- 1/n
  for (i in 2:n) bsm$p[i] <- bsm$p[i-1] + (1/(n + 1 - i))
  bsm$p <- 100*bsm$p/n
  # Plot eigenvalues and % of variation for each axis
  op <- par(mfrow=c(2,1))
  barplot(ev, main="Eigenvalues", col="bisque", las=2)
  abline(h=mean(ev), col="red")
  legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
  barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE, 
          main="% variation", col=c("bisque",2), las=2)
  legend("topright", c("% eigenvalue", "Broken stick model"), 
         pch=15, col=c("bisque",2), bty="n")
  par(op)
}


"cleanplot.pca" <- function(res.pca, ax1=1, ax2=2, point=FALSE, 
                            ahead=0.07, cex=0.7) 
{
  # A function to draw two biplots (scaling 1 and scaling 2) from an object 
  # of class "rda" (PCA or RDA result from vegan's rda() function)
  #
  # License: GPL-2 
  # Authors: Francois Gillet & Daniel Borcard, 24 August 2012
  
  require("vegan")
  
  par(mfrow=c(1,2))
  p <- length(res.pca$CA$eig)
  
  # Scaling 1: "species" scores scaled to relative eigenvalues
  sit.sc1 <- scores(res.pca, display="wa", scaling=1, choices=c(1:p))
  spe.sc1 <- scores(res.pca, display="sp", scaling=1, choices=c(1:p))
  plot(res.pca, choices=c(ax1, ax2), display=c("wa", "sp"), type="n", 
       main="PCA - scaling 1", scaling=1)
  if (point)
  {
    points(sit.sc1[,ax1], sit.sc1[,ax2], pch=20)
    text(res.pca, display="wa", choices=c(ax1, ax2), cex=cex, pos=3, scaling=1)
  }
  else
  {
    text(res.pca, display="wa", choices=c(ax1, ax2), cex=cex, scaling=1)
  }
  text(res.pca, display="sp", choices=c(ax1, ax2), cex=cex, pos=4, 
       col="red", scaling=1)
  arrows(0, 0, spe.sc1[,ax1], spe.sc1[,ax2], length=ahead, angle=20, col="red")
  pcacircle(res.pca)
  
  # Scaling 2: site scores scaled to relative eigenvalues
  sit.sc2 <- scores(res.pca, display="wa", choices=c(1:p))
  spe.sc2 <- scores(res.pca, display="sp", choices=c(1:p))
  plot(res.pca, choices=c(ax1,ax2), display=c("wa","sp"), type="n", 
       main="PCA - scaling 2")
  if (point) {
    points(sit.sc2[,ax1], sit.sc2[,ax2], pch=20)
    text(res.pca, display="wa", choices=c(ax1 ,ax2), cex=cex, pos=3)
  }
  else
  {
    text(res.pca, display="wa", choices=c(ax1, ax2), cex=cex)
  }
  text(res.pca, display="sp", choices=c(ax1, ax2), cex=cex, pos=4, col="red")
  arrows(0, 0, spe.sc2[,ax1], spe.sc2[,ax2], length=ahead, angle=20, col="red")
}



"pcacircle" <- function (pca) 
{
  # Draws a circle of equilibrium contribution on a PCA plot 
  # generated from a vegan analysis.
  # vegan uses special constants for its outputs, hence 
  # the 'const' value below.
  
  eigenv <- pca$CA$eig
  p <- length(eigenv)
  n <- nrow(pca$CA$u)
  tot <- sum(eigenv)
  const <- ((n - 1) * tot)^0.25
  radius <- (2/p)^0.5
  radius <- radius * const
  symbols(0, 0, circles=radius, inches=FALSE, add=TRUE, fg=2)
}

#### 5.2 Ordination Overview ####
## 5.2.1 Multidimensional space
## variables = dimensions, which are often more than 2
# The aim of ordination methods, therefore, is to represent the data along a reduced number of
# orthogonal axes, constructed in such a way that they represent, in decreasing order, the main
# trends of the data.
# The methods in this chapter are descriptive: i.e. no statistical test is provided to assess 
# the significance of the structures detected. For that, see Ch. 6.

## 5.2.2 Ordination in Reduced Space
# Most methods are based on the extraction of the eigenvectors of an association matrix
# The methods in this chapter:
#  * Principal component analysis: works on raw, quantitative data; preserves Euclidean distances among sites
#  * Correspondence analysis: works on frequency-like data; preserves chi-sq distance among rows or columns
#  * Principal coordinate analysis: devoted to ordinatino of distance matrices instead of site-by-variables tables
#  * Non-metric multidimensional scaling: not an eigenvector-based method. NMDS tries to represent the set of
#    objects along a predetermined number of aces while perserving the ordering relationships among them

#### 5.3 Principal Component Analysis ####
# remember- chooses first axis of the scatter that has the greatest dimension (variance), followed by the next,
# orthogonal axis.  PCA sort of carries out a 'rotation' of the orginial system of axes, such that successive
# new axes are orthogonal to one another.
# Works on a dispersion matrix S.

## 5.3.2 PCA on environmental variables using rda()
# since the bariables are expressed in different measurement scales, we compute a PCA on the correlation matrix
# (covariances of standardized variables)
# packages
library(ade4)
library(vegan)
library(gclus)
library(ape)


#data
data(doubs)
spe <- doubs$fish
env <- doubs$env
spa <- doubs$xy
spe <- spe[-8,]
env <- env[-8,]
spa <- spa[-8,]

# PCA on the full dataset
# ************************
env.pca<-rda(env,scale=TRUE) # scale argument calls for a standardization of the variables
env.pca
summary(env.pca) # Default scaling = 2 (see text for the difference between scaling types)
summary(env.pca, scaling=1)

# Examine and plot partial results form PCA output
# *************************************************
?cca.object
# Eigenvalues
(ev<-env.pca$CA$eig)
# Apply Kaiser-Guttman criterion to select axes
ev[ev>mean(ev)]

# Broken stick model
n <- length(ev)
bsm <- data.frame(j=seq(1:n), p=0)
bsm$p[1]<-1/n
for(i in 2:n) bsm$p[i] <- bsm$p[i-1] +(1/(n+1-i))
bsm$p <- 100*bsm$p/n
bsm

# Plot eigenvalues of % of variance for each axis
par(mfrow=c(2,1))
barplot(ev, main="Eigenvalues", col="bisque",las=2)
abline(h=mean(ev),col="red")
legend("topright", "Average eigenvalue",lwd=1,col=2,bty="n")
barplot(t(cbind(100*ev/sum(ev),bsm$p[n:1])),beside=TRUE,main="% variance", col=c("bisque",2),las=2)
legend("topright",c("% eigenvalue","Broken stick model"),pch=15,col=c("bisque",2),bty="n")

# Two PCA biplots: scaling 1 and 2
# ***********************************
# plots using biplot.rda()
par(mfrow=c(1,2))
biplot(env.pca, scaling=1, main="PCA - scaling 1")
biplot(env.pca,scaling=2,main="PCA - scaling 2")

# Plots using cleanplot.pca
# A rectangular graphic window is needed for the two plots
# with points and arrowheads
cleanplot.pca(env.pca, point=TRUE)
# With site labels only
cleanplot.pca(env.pca)
cleanplot.pca(env.pca,ahead = 0) # ... and without arrowheads

# lots of interpretation in text

## 5.3.2.4 Combining clustering and ordination results
# ***************************************************
# clustering the objects using the environmental data: Euclidean distance after
# standardizing the variables, followed by Ward clustering

env.w <- hclust(dist(scale(env)), "ward")

# Cut the dendrogram to 4 groups
gr<- cutree(env.w,k=4)
grl <- levels(factor(gr))

# Get the site scores, scaling 1
sit.scl <- scores(env.pca, display="wa", scaling=1) # "wa" is site scores

# Plot the sites with cluster symbols and colours (scaling 1)
p <- plot(env.pca,display="wa", scaling=1, type="n", main="PCA correlation + clusters")
abline(v=0, lty="dotted")
abline(h=0, lty="dotted")
for (i in 1:length(grl)) points(sit.scl[gr==i,], pch=(14+i), cex=2, col=i+1)
text(sit.scl, row.names(env), cex=0.7, pos=3)

# Add the dendrograms
ordicluster(p, env.w, col="dark grey")
legend(locator(1), paste("Group",c(1:length(grl))), pch=14+c(1:length(grl)), col=1+c(1:length(grl)), pt.cex=2)

## 5.3.3 PCA on Transformed Species Data
# PCA not naturally adpated to the analysis of species abundance data, but data can be pre-transformed

# PCA on the fish data
# **********************

# Hellinger pre-transformation
spe.h <- decostand(spe, "hellinger")
spe.h.pca <- rda(spe.h)
spe.h.pca

# Plot eigenvalues and % of variance for each axis
ev <- spe.h.pca$CA$eig
evplot(ev)

# PCA biplots
cleanplot.pca(spe.h.pca, ahead=0)
spe.pca <- rda(spe)
cleanplot.pca(spe.pca, ahead=0)

## Conditions for appropriate application of PCA -- see text

## 5.3.5 PCA using function PCA()
# PCA on the env data using PCA and biplot.PCA - see text (requires unique functions)

#### 5.4 Correspondence Analysis ####
# see text for introduction; but appropriate for raw abundance data, but the data
# must be frequency or frequency-like, dimensionally homogenous, and non-negative

## 5.4.2 CA using function cca() in vegan
# CA of the raw species dataset (original species abundances)
# *****************************************************

# Compute CA
spe.ca <- cca(spe)
spe.ca
summary(spe.ca)
summary(spe.ca, scaling=1)
(ev2 <- spe.ca$CA$eig)
evplot(ev2)

# CA biplots
# *************
par(mfrow=c(1,2))

# Scaling 1: sites are centroids of species
plot(spe.ca, scaling=1, main="CA fish abundances - biplot scaling 1")
# Scaling 2 (defaul): species are centroids of sites
plot(spe.ca, main="Ca fish abundances - biplot scaling 2")

## 5.4.2.2 Post Hoc explanation of axes using environmental variables
# There are ways to incorporate explanatory variables directly in ordination, but we can also
# try to interpret a simple ordination by means of external variables
# a posteriori projection of env var in a CA
# the last plot produced (CA scaling 2) must be active
spe.ca.env <-envfit(spe.ca, env)
plot(spe.ca.env)

## 5.4.2.3 Reordering the Data Table on the Basis of an Ordination Axis

# Species data table ordered after the CA result
# **************************************************

vegemite(spe,spe.ca)

## 5.4.4 Arch Effect and Detrended Correspondence Analysis (DCA)
# lots of problems- see text

#### 5.5 Principal Coordinate Analysis ####

### 5.5.1 Introduction
## can use any distance measure with PCoA: not just Euclidean or chi-sq (as in two previous methods)
## 5.5.2 
# PCoA on a Bray-Curtis dissimilarity matrix of fish species
# **********************************************************
spe.bray <- vegdist(spe)
spe.b.pcoa <- cmdscale(spe.bray, k=(nrow(spe)-1), eig=TRUE)

# Plot of the sites and weighted average projection of species
ordiplot(scores(spe.b.pcoa)[,c(1,2)], type="t", main= "PCoA with species")
abline(h=0, lty=3)
abline(v=0, lty=3)
# Add species
spe.wa <- wascores(spe.b.pcoa$points[,1:2], spe)
text(spe.wa, rownames(spe.wa), cex=0.7, col="red")

## 5.5.3, Application to the data using pcoa()
# PCoA and projection of species vectors using function pcoa()
# ***********************************************************
spe.h.pcoa <- pcoa(dist(spe.h))
 
# Biplots
par(mfrow=c(1,2))
# First biplot: Hellinger-transformed species data
biplot.pcoa(spe.h.pcoa, spe.h, dir.axis2=-1)
abline(h=0, lty=3)
abline(v=0,lty=3)
#Second biplot: standardized Hellinger-transformed species data
spe.std <- apply(spe.h,2,scale)
biplot.pcoa(spe.h.pcoa, spe.std, dir.axis2=-1)
abline(h=0, lty=3)
abline(v=0,lty=3)
 
## Note: PCoA should actually be reserved to situations where no Euclidean measure is available or selected
# Comparison of PCoA results with Euclidean and non-Euclidean dissimilarity matrices
# **************************************************************
 
#PCoA on a Hellinger distance matrix
is.euclid(dist(spe.h))
summary(spe.h.pcoa)
spe.h.pcoa$values
 
# PCoA on a Bray-Curtis dissimilarity matrix
is.euclid(spe.bray)
spe.bray.pcoa <- pcoa(spe.bray)
spe.bray.pcoa$values
 
# PCoA on the square root of a Bray-Curtis dissimilarity matrix
is.euclid(sqrt(spe.bray))
spe.braysq.pcoa <- pcoa(sqrt(spe.bray))
spe.braysq.pcoa$values

# PCoA on a Bray-Curtis dissimilarity matrix with Lingoes correction
spe.brayl.pcoa <- pcoa(spe.bray, correction="lingoes")
spe.brayl.pcoa$values

# PCoA on a Bray-Curtis dissimilarity matrix with Cailliez correction
spe.brayc.pcoa <- pcoa(spe.bray, correction="cailliez")
spe.brayc.pcoa$values

#### 5.6 Nonmetric Mltidimensional Scaling ####

# if the priority is not to preserve the exact distances, but rather to represent well the ordering 
# relationships among objects in a small and specified number of axes, NMDS may be the solution
# Can use any distance matrix
# Can cope with missing distances
# Not an eigenvalue technique
# see text for generalized procedure

## 5.6.2 Application to the fish data 
# NMDS applied to the fish species - Bray-Curtis distance matrix
# ****************************************************************

spe.nmds <- metaMDS(spe,distance="bray")
spe.nmds
spe.nmds$stress
plot(spe.nmds,type="t",main=paste("NMDS/Bray - Stress = ", round(spe.nmds$stress,3)))
# Goodness-of-fit can be measured with a regression of NMDS distances on the original ones
# Shepard plot and g.o.f
# ************************
par(mfrow=c(1,2))
stressplot(spe.nmds,main="Shephard Plot")
gof=goodness(spe.nmds)
plot(spe.nmds,type="t",main="Goodness of fit")
points(spe.nmds,display="sites",cex=gof*100)

# Add colours from a clustering results to an NMDS plot
# ****************************************************

# Ward clustering of Bray-Curtis dissimilarity matrix
# and extraction of four groups
spe.bray.ward <- hclust(spe.bray, "ward")
spe.bw.groups <- cutree(spe.bray.ward,k=4)
grp.lev <- levels(factor(spe.bw.groups))

# Combination with NMDS
sit.sc <- scores(spe.nmds)
p <- ordiplot(sit.sc, type="n", main="NMDS/Bray + clusters Ward/Bray")
for(i in 1:length(grp.lev)) points(sit.sc[spe.bw.groups==i,],pch=(14+i), cex=2, col=i+1)
text(sit.sc,row.names(spe),pos=4,cex=0.7)

#Add the dendrogram
ordicluster(p, spe.bray.ward,col="dark grey")
legend(locator(1), paste("Group",c(1:length(grp.lev))),pch=14+c(1:length(grp.lev)), col=1+c(1:length(grp.lev)), pt.cex=2)

#### 5.7 Handwritten Ordination Function ####
# See text for explanation and code...