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

#Single Linkage Agglomerative Clustering
spe.ch.single <- hclust(spe.dc, method='single')
plot(spe.ch.single)

#Complete Linkage Agglomerative Clustering
spe.ch.complete <- hclust(spe.dc, method='complete')
plot(spe.ch.complete)

#UPGMA (average distance) agglomerative clustering
spe.ch.UPGMA <- hclust(spe.dc,method='average')
plot(spe.ch.UPGMA)

#centroid clustering
spe.ch.centroid <- hclust(spe.dc, method='centroid')
plot(spe.ch.centroid)

#Ward's minimum variance clustering
spe.ch.ward <- hclust(spe.dc, method='ward')
plot(spe.ch.ward)
spe.ch.ward$height <- sqrt(spe.ch.ward$height)
plot(spe.ch.ward)

## Cophenetic correlation for the above ##
# single linkage clustering
spe.ch.single.coph <- cophenetic(spe.ch.single)
cor(spe.dc,spe.ch.single.coph)
#complete linkage clustering
spe.ch.comp.coph <- cophenetic(spe.ch.complete)
cor(spe.dc,spe.ch.comp.coph)
#average clustering
spe.ch.UPGMA.coph <- cophenetic(spe.ch.UPGMA)
cor(spe.dc,spe.ch.UPGMA.coph)
#ward clustering
spe.ch.ward.coph <- cophenetic(spe.ch.ward)
cor(spe.dc,spe.ch.ward.coph)
#Shepard-like diagrams to illustrate the relationship between a distance matrix and a set of cophenetic matrices
par(mfrow=c(2,2))
plot(spe.dc,spe.ch.single.coph,xlab='Chord Distance',ylab='Cophenetic Distance',asp=1,xlim=c(0,sqrt(2)),ylim=c(0,sqrt(2)),
     main=c('Single linkage', paste('Cophenetic correlation',round(cor(spe.dc,spe.ch.single.coph),3))))
abline(0,1)
lines(lowess(spe.dc,spe.ch.single.coph),col='red')

plot(spe.dc,spe.ch.comp.coph,xlab='Chord Distance',ylab='Cophenetic Distance',asp=1,xlim=c(0,sqrt(2)),ylim=c(0,sqrt(2)),
     main=c('Complete linkage', paste('Cophenetic correlation',round(cor(spe.dc,spe.ch.comp.coph),3))))
abline(0,1)
lines(lowess(spe.dc,spe.ch.comp.coph),col='red')

plot(spe.dc,spe.ch.UPGMA.coph,xlab='Chord Distance',ylab='Cophenetic Distance',asp=1,xlim=c(0,sqrt(2)),ylim=c(0,sqrt(2)),
     main=c('UPGMA', paste('Cophenetic correlation',round(cor(spe.dc,spe.ch.UPGMA.coph),3))))
abline(0,1)
lines(lowess(spe.dc,spe.ch.UPGMA.coph),col='red')

plot(spe.dc,spe.ch.ward.coph,xlab='Chord Distance',ylab='Cophenetic Distance',asp=1,xlim=c(0,sqrt(2)),ylim=c(0,sqrt(2)),
     main=c('Ward clustering', paste('Cophenetic correlation',round(cor(spe.dc,spe.ch.ward.coph),3))))
abline(0,1)
lines(lowess(spe.dc,spe.ch.ward.coph),col='red')

#Gower distance (another choice for determining best clustering)
gow.dist.single <- sum((spe.dc-spe.ch.single.coph)^2)
gow.dist.comp <- sum((spe.dc-spe.ch.comp.coph)^2)
gow.dist.UPGMA <- sum((spe.dc-spe.ch.UPGMA.coph)^2)
gow.dist.ward <- sum((spe.dc-spe.ch.ward.coph)^2)
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
  sil <- silhouette(cutree(spe.ch.ward,k=k),spe.dc)
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
   mt <- cor(spe.dc,distgr,method="pearson")
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
sil <- silhouette(cutg,spe.dc)
silo<- sortSilhouette(sil)
rownames(silo) <- row.names(spe)[attr(silo,'iOrd')]
plot(silo,main='Silhouette plot - Chord - Ward',cex.names=0.8,col=cutg+1,nmax.lab=100)

##Final dendrogram with the selected groups
# Reorder dendrogram from hclust(). reorder.hclust() reorders objects so that their order in the dissimilarity matrix
# is respected as much as possible. This does not affect the topology of the dendrogram
spe.chwo <- reorder.hclust(spe.ch.ward,spe.dc)

#plot reordered dendrogram with group labels
plot(spe.chwo,hang=-1,xlab="4 groups",sub='',ylab='Height',main='Chord - Ward (reordered)', labels=cutree(spe.chwo,k=k))
rect.hclust(spe.chwo,k=k)

#plot the final dendrogram with group colors (RGBCMY)
hcoplot(spe.ch.ward,spe.dc,k=4)


##Other representations of this final result
# Plot of the 4 Ward clusters on a map of the river
spa <- doubs$xy
spa <- spa[-8,]
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
heatmap(as.matrix(spe.dc), Rowv=dend,symm=TRUE,margin=c(3,3))
#Ordered community table
# species are ordered by the weighted averages on site scores
or <- vegemite(spe,spe.chwo)
#and a heatmap for the doubly ordered community table
heatmap(t(spe[rev(or$species)]), Rowv=NA,Colv=dend,col=c('white',brewer.pal(5,'Greens')), scale='none',
        margin=c(4,4),ylab='Species (weighted averages of sites)',xlab='Sites')