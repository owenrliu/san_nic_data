library(ade4)
library(vegan)
library(gclus)
library(cluster)
library(FD)
library(RColorBrewer)
library(labdsv)



data(doubs)
spe <- doubs$fish
spe <- spe[-8,]

#### Q MODE ####

#Bray-Curtis dissimilarity
spe.db <- vegdist(spe)

#Bray-Curtis (log trans)
spe.dbln <- vegdist(log1p(spe))

#chord distance matrix
spe.norm <- decostand(spe,'nor')
spe.dc <- dist(spe.norm)

#Hellinger distance matrix
spe.hel <- decostand(spe,'hel')
spe.dh <- dist(spe.hel)

# coldiss()
# Color plots of a dissimilarity matrix, without and with ordering
#
# License: GPL-2 
# Author: Francois Gillet, 23 August 2012
#

"coldiss" <- function(D, nc = 4, byrank = TRUE, diag = FALSE)
{
  require(gclus)
  
  if (max(D)>1) D <- D/max(D)
  
  if (byrank) {
    spe.color <- dmat.color(1-D, cm.colors(nc))
  }
  else {
    spe.color <- dmat.color(1-D, byrank=FALSE, cm.colors(nc))
  }
  
  spe.o <- order.single(1-D)
  speo.color <- spe.color[spe.o, spe.o]
  
  op <- par(mfrow=c(1,2), pty="s")
  
  if (diag) {
    plotcolors(spe.color, rlabels=attributes(D)$Labels, 
               main="Dissimilarity Matrix", 
               dlabels=attributes(D)$Labels)
    plotcolors(speo.color, rlabels=attributes(D)$Labels[spe.o], 
               main="Ordered Dissimilarity Matrix", 
               dlabels=attributes(D)$Labels[spe.o])
  }
  else {
    plotcolors(spe.color, rlabels=attributes(D)$Labels, 
               main="Dissimilarity Matrix")
    plotcolors(speo.color, rlabels=attributes(D)$Labels[spe.o], 
               main="Ordered Dissimilarity Matrix")
  }
  
  par(op)
}

# Usage:
# coldiss(D = dissimilarity.matrix, nc = 4, byrank = TRUE, diag = FALSE)
# If D is not a dissimilarity matrix (max(D) > 1), then D is divided by max(D)
# nc 							number of colours (classes)
# byrank= TRUE		equal-sized classes
# byrank= FALSE		equal-length intervals
# diag = TRUE			print object labels also on the diagonal

# Example:
# coldiss(spe.dj, nc=9, byrank=F, diag=T)

# Function hcoplot()
# Reorder and plot dendrogram with colors for groups and legend
#
# Usage:
# hcoplot(tree = hclust.object, diss = dissimilarity.matrix, k = nb.clusters, 
#	title = paste("Reordered dendrogram from",deparse(tree$call),sep="\n"))
#
# License: GPL-2 
# Author: Francois Gillet, 23 August 2012

"hcoplot" <- function(tree, diss, k, 
                      title=paste("Reordered dendrogram from", deparse(tree$call), sep="\n"))
{
  require(gclus)
  gr <- cutree(tree, k=k)
  tor <- reorder.hclust(tree, diss)
  plot(tor, hang=-1, xlab=paste(length(gr),"sites"), sub=paste(k,"clusters"), 
       main=title)
  so <- gr[tor$order]
  gro <- numeric(k)
  for (i in 1:k)
  {
    gro[i] <- so[1]
    if (i<k) so <- so[so!=gro[i]]
  }
  rect.hclust(tor, k=k, border=gro+1, cluster=gr)
  legend("topright", paste("Cluster",1:k), pch=22, col=2:(k+1), bty="n")
}

#B-C untrans
coldiss(spe.db, byrank=FALSE,diag=TRUE)

#B-C logtrans
coldiss(spe.dbln, byrank=FALSE,diag=TRUE)

#Chord distance
coldiss(spe.dc,byrank=FALSE,diag=TRUE)

#Hellinger
coldiss(spe.dh,byrank=F,diag=T)

#### R MODE ####
spe.t <- t(spe)
spe.t.chi <- decostand(spe.t,'chi.square')
spe.t.D16 <- dist(spe.t.chi)
coldiss(spe.t.D16,diag=TRUE)

####Chapter 4 Cluster Analysis ####
data(doubs)
spe <- doubs$fish
env <- doubs$env
spa <- doubs$xy
spe <- spe[-8,]
env <- env[-8,]
spa <- spa[-8,]

# Chord distance
spe.norm <- decostand(spe,"normalize")
spe.ch <- vegdist(spe.norm,"euc")

#Single Linkage Agglomerative Clustering
spe.ch.single <- hclust(spe.ch, method='single')
plot(spe.ch.single)

#Complete Linkage Agglomerative Clustering
spe.ch.complete <- hclust(spe.ch, method='complete')
plot(spe.ch.complete)

#UPGMA (average distance) agglomerative clustering
spe.ch.UPGMA <- hclust(spe.ch,method='average')
plot(spe.ch.UPGMA)

#centroid clustering
spe.ch.centroid <- hclust(spe.ch, method='centroid')
plot(spe.ch.centroid)

#Ward's minimum variance clustering
spe.ch.ward <- hclust(spe.ch, method='ward')
plot(spe.ch.ward)
spe.ch.ward$height <- sqrt(spe.ch.ward$height)
plot(spe.ch.ward)

## Cophenetic correlation for the above ##
# single linkage clustering
spe.ch.single.coph <- cophenetic(spe.ch.single)
cor(spe.ch,spe.ch.single.coph)
#complete linkage clustering
spe.ch.comp.coph <- cophenetic(spe.ch.complete)
cor(spe.ch,spe.ch.comp.coph)
#average clustering
spe.ch.UPGMA.coph <- cophenetic(spe.ch.UPGMA)
cor(spe.ch,spe.ch.UPGMA.coph)
#ward clustering
spe.ch.ward.coph <- cophenetic(spe.ch.ward)
cor(spe.ch,spe.ch.ward.coph)
#Shepard-like diagrams to illustrate the relationship between a distance matrix and a set of cophenetic matrices
par(mfrow=c(2,2))
plot(spe.ch,spe.ch.single.coph,xlab='Chord Distance',ylab='Cophenetic Distance',asp=1,xlim=c(0,sqrt(2)),ylim=c(0,sqrt(2)),
     main=c('Single linkage', paste('Cophenetic correlation',round(cor(spe.ch,spe.ch.single.coph),3))))
abline(0,1)
lines(lowess(spe.ch,spe.ch.single.coph),col='red')

plot(spe.ch,spe.ch.comp.coph,xlab='Chord Distance',ylab='Cophenetic Distance',asp=1,xlim=c(0,sqrt(2)),ylim=c(0,sqrt(2)),
     main=c('Complete linkage', paste('Cophenetic correlation',round(cor(spe.ch,spe.ch.comp.coph),3))))
abline(0,1)
lines(lowess(spe.ch,spe.ch.comp.coph),col='red')

plot(spe.ch,spe.ch.UPGMA.coph,xlab='Chord Distance',ylab='Cophenetic Distance',asp=1,xlim=c(0,sqrt(2)),ylim=c(0,sqrt(2)),
     main=c('UPGMA', paste('Cophenetic correlation',round(cor(spe.ch,spe.ch.UPGMA.coph),3))))
abline(0,1)
lines(lowess(spe.ch,spe.ch.UPGMA.coph),col='red')

plot(spe.ch,spe.ch.ward.coph,xlab='Chord Distance',ylab='Cophenetic Distance',asp=1,xlim=c(0,sqrt(2)),ylim=c(0,sqrt(2)),
     main=c('Ward clustering', paste('Cophenetic correlation',round(cor(spe.ch,spe.ch.ward.coph),3))))
abline(0,1)
lines(lowess(spe.ch,spe.ch.ward.coph),col='red')

#Gower distance (another choice for determining best clustering)
gow.dist.single <- sum((spe.ch-spe.ch.single.coph)^2)
gow.dist.comp <- sum((spe.ch-spe.ch.comp.coph)^2)
gow.dist.UPGMA <- sum((spe.ch-spe.ch.UPGMA.coph)^2)
gow.dist.ward <- sum((spe.ch-spe.ch.ward.coph)^2)
gow.dist.single
gow.dist.comp
gow.dist.UPGMA
gow.dist.ward

#Graphs of fusion level values (where to cut the dendrogram?)
par(mfrow=c(2,2))
#Plot the fusion level values of the single linkage clustering
summary(spe.ch.single)
plot(spe.ch.single$height,nrow(spe):2,type='S',main='Fusion levels - Chord - Single',ylab='k (number of clusters)',
     xlab='h (node height)',col='grey')
text(spe.ch.single$height,nrow(spe):2,nrow(spe):2, col='red',cex=0.8)

#Plot the fusion level values of the complete linkage clustering
summary(spe.ch.complete)
plot(spe.ch.complete$height,nrow(spe):2,type='S',main='Fusion levels - Chord - Complete',ylab='k (number of clusters)',
     xlab='h (node height)',col='grey')
text(spe.ch.complete$height,nrow(spe):2,nrow(spe):2, col='red',cex=0.8)

#Plot the fusion level values of the UPGMA clustering
summary(spe.ch.UPGMA)
plot(spe.ch.UPGMA$height,nrow(spe):2,type='S',main='Fusion levels - Chord - UPGMA',ylab='k (number of clusters)',
     xlab='h (node height)',col='grey')
text(spe.ch.UPGMA$height,nrow(spe):2,nrow(spe):2, col='red',cex=0.8)

#Plot the fusion level values of the Ward clustering
summary(spe.ch.ward)
plot(spe.ch.ward$height,nrow(spe):2,type='S',main='Fusion levels - Chord - Ward',ylab='k (number of clusters)',
     xlab='h (node height)',col='grey')
text(spe.ch.ward$height,nrow(spe):2,nrow(spe):2, col='red',cex=0.8)

##Cut the trees to obtain k groups and compare the group contents using contingency tables

#choose a common number of groups
k<-4 #Number of groups where at least a small jump is present in all four graphs of fusion levels
#cut the dendrograms
spebc.single.g <- cutree(spe.ch.single,k)
spebc.complete.g <- cutree(spe.ch.complete,k)
spebc.UPGMA.g <- cutree(spe.ch.UPGMA,k)
spebc.ward.g <- cutree(spe.ch.ward,k)

#compare classifications by constructing contingency tables
#single vs. complete linkage
table(spebc.single.g,spebc.complete.g)
#single linkage vs. UPGMA
table(spebc.single.g,spebc.UPGMA.g)
#single linkage vs. Ward
table(spebc.single.g,spebc.ward.g)
#complete linkage vs. UPGMA
table(spebc.complete.g,spebc.UPGMA.g)
#complete linkage vs. Ward
table(spebc.complete.g,spebc.ward.g)
#UPGMA vs. Ward
table(spebc.UPGMA.g,spebc.ward.g)


##Two other methods for identifying an appropriate number of groups: silhouette widths and Mantel comparison
## silhouette widths - measure of the degree of membership of an object to its cluster, based on the average distance
## between this object and all objects of its cluster, compared to the same measure for the next closest cluster.
## Can be averaged over all objects of a partition

# Optimal number of clusters according to silhouette widths(Rousseeuw quality index)
# Plot average silhouette widths (using Ward clustering) for all partitions except for the trivial partition in a single group
# First, create empty vector in which the asw values will be written
asw <- numeric(nrow(spe))

#Retrieve and write the asw values into the vector
for (k in 2:(nrow(spe)-1)) {
  sil <- silhouette(cutree(spe.ch.ward,k=k),spe.ch)
  asw[k] <- summary(sil)$avg.width
}

# Best (largest) silhouette width
k.best <- which.max(asw)

#the plot is produced by plot.silhouette
plot(1:nrow(spe),asw,type='h',main="Silhouette - optimal number of clusters, Ward",xlab='k (number of groups)',
     ylab="Average silhouette width")
axis(1,k.best,paste('optimum',k.best,sep="\n"),col='red',font=2,col.axis='red')
points(k.best,max(asw),pch=16,col='red',cex=1.5)

cat("", "Silhouette - Optimal number of clusters k =", k.best,'\n', "with an average silhouette width of", max(asw),"\n")

## Optimal number of cluster according to Mantel statistic (Pearson)
# Function to compute a binary distance matrix from groups
grpdist <- function(X) {
  require(cluster)
  gr <- as.data.frame(as.factor(X))
  distgr <- daisy(gr,"gower")
  distgr
}

#Run based on Ward clustering
 kt <- data.frame(k=1:nrow(spe), r=0)
 
 for (i in 2:(nrow(spe)-1)) {
   gr<-cutree(spe.ch.ward,i)
   distgr <- grpdist(gr)
   mt <- cor(spe.ch,distgr,method="pearson")
   kt[i,2]<-mt
 }
 
kt
k.best <- which.max(kt$r)

plot(kt$k,kt$r,type='h',main='Mantel - optimal number of clusters - Ward',xlab='k (number of groups)',
    ylab="Pearson's correlation")
axis(1,k.best,paste('optimum',k.best,sep="\n"),col='red',font=2,col.axis='red')
points(k.best,max(kt$r),pch=16,col='red',cex=1.5)

## we choose k=4 as our final number of groups, and Ward clustering as our best-option dendrogram
# Silhouette plot of final partition
k<-4
# Silhouette plot
cutg <- cutree(spe.ch.ward,k=k)
sil <- silhouette(cutg,spe.ch)
silo<- sortSilhouette(sil)
rownames(silo) <- row.names(spe)[attr(silo,'iOrd')]
plot(silo,main='Silhouette plot - Chord - Ward',cex.names=0.8,col=silo+1,nmax.lab=100)

##Final dendrogram with the selected groups
# Reorder dendrogram from hclust(). reorder.hclust() reorders objects so that their order in the dissimilarity matrix
# is respected as much as possible. This does not affect the topology of the dendrogram
spe.chwo <- reorder.hclust(spe.ch.ward,spe.ch)

#plot reordered dendrogram with group labels
plot(spe.chwo,hang=-1,xlab="4 groups",sub='',ylab='Height',main='Chord - Ward (reordered)', labels=cutree(spe.chwo,k=k))
rect.hclust(spe.chwo,k=k)

#plot the final dendrogram with group colors (RGBCMY)
hcoplot(spe.ch.ward,spe.ch,k=4)


##Other representations of this final result
# Plot of the 4 Ward clusters on a map of the river
plot(spa,asp=1,type='n',main="Four Ward Groups",xlab='x coord(km)',ylab='y coord(km)')
lines(spa,col='lightblue')
text(50,10,"Upstream",cex=1.2,col='red')
text(25,115,"Downstream",cex=1.2,col='red')
#add the 4 groups
grw <- spebc.ward.g
k <- length(levels(factor(grw)))
for(i in 1:k) {
  points(spa[grw==i,1], spa[grw==i,2], pch=i+20,cex=3,col=i-1,bg=i+1)
}
text(spa, row.names(spa), cex=0.8, col='white',font=2)
legend('bottomright',paste('Group',1:k), pch=(1:k)+20,col=2:(k+1),pt.bg=2:(k+1),pt.cex=2,bty='n')

##Heat map of the distance matrix ordered with the dendrogram
dend <- as.dendrogram(spe.chwo)
heatmap(as.matrix(spe.ch), Rowv=dend,symm=TRUE,margin=c(3,3))
#Ordered community table
# species are ordered by the weighted averages on site scores
or <- vegemite(spe,spe.chwo)
#and a heatmap for the doubly ordered community table
heatmap(t(spe[rev(or$species)]), Rowv=NA,Colv=dend,col=c('white',brewer.pal(5,'Greens')), scale='none',
        margin=c(4,4),ylab='Species (weighted averages of sites)',xlab='Sites')


#### 4.8.1 K-means partitioning ####
## Note: for a description of the method and its reasoning, see the text
# ************************************

# k-means partitioning of the pre-transformed species data
# *******************************************************

spe.kmeans <- kmeans(spe.norm, centers=4, nstart=100)

# Comparison with the 4-group classification derived from Ward clustering:
table(spe.kmeans$cluster,spebc.ward.g)

# k-means partitioning, 2 to 10 groups
# **************************************

spe.KM.cascade <- cascadeKM(spe.norm, inf.gr=2, sup.gr=10, iter=100, criterion ="ssi")
plot(spe.KM.cascade,sortg=T)

summary(spe.KM.cascade)
spe.KM.cascade$results
spe.KM.cascade$partition

# Next, we should examine the contents of the defined k-means clusters
# Example
# Reorder the sites according to the k-means result
spe[order(spe.kmeans$cluster),]

# Reorder sites and species using fucntion vegemite()
ord.KM <- vegemite(spe, spe.kmeans$cluster)
spe[ord.KM$sites,ord.KM$species]


#### 4.8.2 Partitioning around medoids ####
# searches for k representative objects or medoids among the observations of the dataset. These observations should 
# represent the sturcture of the data. After finding a set of k medoids, k clusters are constructed by assigning
# each observation to the nearest medoids. The goals is to find k representative objects which minimize the sum
# of the dissimilarities of the observations to their closest representative object.

## Partitioning around medoids (PAM)
# Computed on the chord distance matrix
# *****************************************

# Choices of the number of clusters
asw <- numeric(nrow(spe))

# Loop to compute average silhouette width fo 2 to 28 groups.
# asw means "average silhouette width"

for (k in 2:(nrow(spe)-1)) asw[k] <- pam(spe.ch, k, diss=TRUE)$silinfo$avg.width
k.best<-which.max(asw)
cat("","Silhouette-optimal number of clusters k =",k.best,"\n", "with an average silhouette width of",max(asw),"\n")
plot(1:nrow(spe),asw,type="h",main="Choice of the number of clusters",xlab="k (number of groups)", 
     ylab="Average silhouette width")
axis(1, k.best, paste("optimum",k.best,sep="\n"),col="red",font=2, col.axis="red")
points(k.best,max(asw),pch=16,col="red",cex=1.5)

# PAM for k = 4 groups (even though it's not really the best choice)
spe.ch.pam <- pam(spe.ch, k=4, diss=TRUE)
summary(spe.ch.pam)
spe.ch.pam.g <- spe.ch.pam$clustering
spe.ch.pam$silinfo$widths
#Compare with classification from Ward clustering and from k-means
table(spe.ch.pam.g, spebc.ward.g)
table(spe.ch.pam.g,spe.kmeans$cluster)
#PAM result is quite different from ward clustering and k means
# Silhouette profile for k=4, k-means and PAM methods
par(mfrow=c(1,2))
plot(silhouette(spe.kmeans$cluster,spe.ch),main="Silhouette plot - k-means", cex.names=0.8, col=sort(spe.kmeans$cluster)+1)
plot(silhouette(spe.ch.pam),main="Silhouette plot - PAM", cex.names=0.8, col=sort(spe.ch.pam$silinfo$widths[,"cluster"])+1)

# Plot of the 4 k-means clusters on a map of the Doubs river
# **************************************************************
dev.off()
plot(spa, asp=1,type="n",main="Four k-means groups", xlab="x-coordinate (km)",ylab="y coordinate (km)")
lines(spa,col="light blue")
text(50,10,"Upstream",cex=1.2, col="red")
text(25,115,"Downsteam",cex=1.2,col="red")

grKM <- spe.kmeans$cluster
k <- length(levels(factor(grKM)))

for(i in 1:k) points(spa[grKM==i,1],spa[grKM==i,2], pch=i+20, cex=3,col=i+1, bg=i+1)

text(spa,row.names(spa), cex=0.8,col="white",font=2)
legend("bottomright",paste("Group", 1:k), pch=(1:k)+20, col=2:(k+1),pt.bg=2:(k+1),pt.cex=2,bty="n")

#### 4.9.1 Comparing a Typology with External Data (ANOVA approach) ####

# Relationships between fish clusters and 4 environmental variables based on the k-means clustering results
# ***********************************************************************

attach(env)
# Boxplots of quantitative environmental variables :
# Altitude, Slope, Oxygen, Ammonium

par(mfrow=c(2,2))

boxplot(sqrt(alt) ~ spe.kmeans$cluster, main="Altitude", las=1, ylab="sqrt (alt)",col=2:5,varwidth=TRUE)
boxplot(log(slo) ~ spe.kmeans$cluster, main="Slope", las=1, ylab="log(slo)",col=2:5,varwidth=TRUE)
boxplot(oxy ~ spe.kmeans$cluster, main="Oxygen", las=1, ylab="oxy",col=2:5,varwidth=TRUE)
boxplot(sqrt(amm) ~ spe.kmeans$cluster, main="Ammonium", las=1, ylab="sqrt (amm)",col=2:5,varwidth=TRUE)

# Test of ANOVA assumptions
# Normality of residuals
shapiro.test(resid(lm(sqrt(alt) ~ as.factor(spe.kmeans$cluster))))
shapiro.test(resid(lm(log(slo) ~ as.factor(spe.kmeans$cluster))))
shapiro.test(resid(lm(oxy ~ as.factor(spe.kmeans$cluster))))
shapiro.test(resid(lm(sqrt(amm) ~ as.factor(spe.kmeans$cluster))))

# Homogeneity of variances
bartlett.test(sqrt(alt),as.factor(spe.kmeans$cluster)) #This one has heterogeneous variances, and is not appr. for ANOVA
bartlett.test(log(slo),as.factor(spe.kmeans$cluster))
bartlett.test(oxy,as.factor(spe.kmeans$cluster))
bartlett.test(sqrt(amm),as.factor(spe.kmeans$cluster))

# ANOVA of the testable variables
summary(aov(log(slo) ~ as.factor(spe.kmeans$cluster)))
summary(aov(oxy ~ as.factor(spe.kmeans$cluster)))
summary(aov(sqrt(amm) ~ as.factor(spe.kmeans$cluster)))

# Kruskal-Wallis test of the other variable (altitude)
kruskal.test(alt ~ as.factor(spe.kmeans$cluster))

detach(env)

#### 4.9.2 Comparing Two Typologies (Contingency Table Approach) ####
# What if we produced two indepedent typologies, one from species data and one from environmental, and compare them?

# Contingency table of two typologies
# ***********************************

# Environment-based typology (similar to Ch. 2)
env2 <- env[,-1]
env.de <- vegdist(scale(env2),"euc") # Dissimilarity matrix
env.kmeans <- kmeans(env.de, centers=4, nstart=100)
env.KM.4 <- env.kmeans$cluster

# Table crossing the k-means and environment 4-group typologies
table(as.factor(spe.kmeans$cluster), as.factor(env.kmeans$cluster))

#### 4.10 Species ASsemblages ####
# Many approaches exist to the problem of identifying species associations in a data set. Here are some examples

#### 4.10.1 Simples statistics on group contents ####
# Using the previous sections, we can compute simple statistics from typologies (clustering) and look for 
# species that are more present, abundant, or specific in each cluster of sites

# Mean abundances on k-means site clusters
# ****************************************

groups <- as.factor(spe.kmeans$cluster)
spe.means <- matrix(0,ncol(spe), length(levels(groups)))
row.names(spe.means) <- colnames(spe)
for(i in 1:ncol(spe)) spe.means[i,] <- tapply(spe[,i], spe.kmeans$cluster,mean)

# Mean species abundances of the four groups
group1 <- round(sort(spe.means[,1], decreasing=TRUE),2)
group2 <- round(sort(spe.means[,2], decreasing=TRUE),2)
group3 <- round(sort(spe.means[,3], decreasing=TRUE),2)
group4 <- round(sort(spe.means[,4], decreasing=TRUE),2)

# Species with abundances greater than group mean species abundance
group1.domin <- which(group1>mean(group1))
group1
group1.domin
group2.domin <- which(group2>mean(group2))
group2
group2.domin
group3.domin <- which(group3>mean(group3))
group3
group3.domin
group4.domin <- which(group4>mean(group4))
group4
group4.domin

#### 4.10.2 Kendall's W Coefficient of Concordance ####
# "An overall test of independence of species is first carried out. If the null hypothesis is rejected, one looks
# for groups of correlated species and, within each group, tests the contribution of each species to the overall
# statistic, using a permutation test"
# Note: the search for assemblages is done without any reference to a typology of sites known a priori, and aims
# to find the most encompassing assemblages, i.e., the smallest number of groups containing the largest number of
# positively and significantly associated species

# Kendall's W
# ***************************************************

# Extraction of the most abundant species
sp.sum <- apply(spe, 2, sum)
spe.sorted <- spe[,order(sp.sum,decreasing=TRUE)]
spe.small <- spe.sorted[,1:20]

# Transformation of species data and transposition
spe.small.hel <- decostand(spe.small,"hellinger")
spe.small.std <- decostand(spe.small.hel, "standardize")
spe.small.t <- t(spe.small.std) # transpose (sites x species to species x sites)

# k-means partitioning of species
spe.t.kmeans.casc <- cascadeKM(spe.small.t, inf.gr=2,sup.gr=8,iter=100,criterion="calinski")
plot(spe.t.kmeans.casc,sortg=TRUE)
# result indicates that 2 groups may be a good choice

# The partition into 2 groups is found in column 1 of the object $partition
clusters <- spe.t.kmeans.casc$partition[,1]
clusters

# Now that we have two groups of species, let us run a global Kenall W test on these groups
spe.kendall.global <- kendall.global(spe.small.hel,clusters)
spe.kendall.global

# NOTE: If the corrected permutational p-values are equal to or small than 0.05, you can consider that all groups
# are globally significant, i.e. on the whole, they contain species that are concordant (not necessarily all
# species, but at least some). If the corrected p-values for some groups are not significant, it indicates that
# these groups have to subdivided into smaller groups.  To see WHICH species in each group are significantly
# concordant, we run a posteriori tests:

spe.kendall.post <- kendall.post(spe.small.hel,clusters, nperm=9999)
spe.kendall.post

#### 4.10.3 Species Assemblages in Presence-Absence Data ####
# (skipped for now, p. 95 in text)

#### 4.10.4 IndVal: Species Indicator Values ####
# From Dufrene and Legendre (1997), invdval index
# A high indicator value is obtained by a combination of large mean abundance within a group compared to
# the other groups (specificity), and presence in most sites of that group (fidelity)
# "The indval approach looks for species that are both necessary and sufficient, i.e. if you find that
# species you should be in that type, and if you in that type you should find that species."

# Various ways to group sites (e.g. many of the previous sections), but clustering by species is a bit
# circular.  Instead, maybe cluster sites on the basis of independent data (e.g. env. variables).
# Then indicator species can be considered indicator in a true sens of the word (indicative of the 
# ecological conditions of the group).

# Species indicator values example, using distance from source as the explanatory variable

# Divide sites into 4 groups depending on the distance to the source of the river
das.D1 <- dist(data.frame(das=env$dfs, row.names=rownames(env)))
dasD1.kmeans <- kmeans(das.D1, centers=4, nstart=100)
dasD1.kmeans$cluster

# Indicator species for this typology of the sites
(iva <- indval(spe, dasD1.kmeans$cluster))
# objects : relfrq is the relative freq of the species in each group
# relabu is the relative abundance of the species across groups
# indval is the indicator value of each species
# The two next items are extracted from the indval table
# the last contains the results of the permutation tests

# Table of the significant indicator species
gr <- iva$maxcls[iva$pval<=0.05]
iv <- iva$indcls[iva$pval<=0.05]
pv <- iva$pval[iva$pval<=0.05]
fr <- apply(spe>0, 2, sum)[iva$pval<=0.05]
fidg <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
