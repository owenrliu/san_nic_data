### Numerical Ecology Chapter 6 ###


#### Packages and Helper Functions ####
library(ade4)
library(vegan)
library(packfor)
library(MASS)
library(ellipse)
library(FactoMineR)

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

#### 6.1 Objectives ####
# Simple ordination (chapter 5) is a pasdsive form of analysis, in that the user interprets the ordination results
# a posteriori.  Canonical ordination, on the contrary, associates two or more data sets in the ordination process
# itself.  Consequently, we can extract structures of a data set that are related to structure in other data sets,
# and/or formally test statistical hypotheses about the significance of these relatinoships.
# Symmetrical and assymetrical methods
# This chapter:

# Learn how to choose among various techniques:
# redundancy analysis (RDA)
# distance-based redundancy analysis (db-RDA)
# canonical correspondence analysis (CCA)
# linear discriminant analysis (LDA)
# canonical correlation analysis (CCorA)
# co-inerita analysis (CoIA)
# multiple factor analysis (MFA)

# Compute them and properly interpret them
# Explore particular applications of some methods


#### 6.3 Redundancy Analysis ####
# Method that comines regession and principal component analysis (PCA)
# Conceptually, it is a multivariate multiple linear regression followed by a PCA of the table of fitted values
# Steps (see text); on a matrix Y of centred responses data and a matrix X of centred (standardized) explantory vars:
  # Regress each y variable on explanatory table X and compute the fitted (yhat) and residual vectors.
  # Assemble all fitted vectors yhat in a matrix Yhat of fitted values
  # Compute a PCA of Yhat. produces a vector of canonical eigenvalues and matrix U of canonical eigenvectors
  # Use matrix U to compute two types of ordination site scores: use either the original centred matrix Y to
  # obtain an ordination in the space of the original variables Y (i.e., compute YU); or, use the matrix Yhat
  # of fitted values to obtain an ordination in the space of variables X (i.e. compute YhatU)
  # residual values can also be submitted to PCA to obtain an unconstrained ordination of the residuals

# In other words, this methods seeks, in successive order, a series of linear combinations of the explanatory
# variables that best explain the variation of the response matrix
# The difference from unconstrained ordination is important: the matrix of explanatory variables conditions
# the 'weights' (eigenvalues)., the orthogonality, and the direction of the ordination axes.
# therefore, in RDA, one can truly say that the axes explain or model the variation of the dependent matrix.
# Furthermore, a hypothesis of absence of linear relationship between Y and X can be tested in RDA
# Each of the canonical axes is a linear combination (multiple regression model) of ALL explanatory variables

### 6.3.2 RDA of the Doubs River Data ###
#data
data(doubs)
spe <- doubs$fish
env <- doubs$env
spa <- doubs$xy
spe <- spe[-8,]
env <- env[-8,]
spa <- spa[-8,]

# Set aside the variable 'das' for later use
das <- env$dfs

# Remove the 'dfs' variables from the nev dataset
env <- env[,-1]

# Recode the slope variable (slo) into a factor (qualitative) variable to show how these are handled in ordination
pen2 <- rep("very_steep", nrow(env))
pen2[env$slo <= quantile(env$slo)[4]] = "steep"
pen2[env$slo <= quantile(env$slo)[3]] = "moderate"
pen2[env$slo <= quantile(env$slo)[2]] = "low"
pen2 <- factor(pen2, levels=c('low','moderate','steep','very_steep'))
table(pen2)
# Create an env2 data frame with slope as a qualitative variables
env2 <- env
env2$slo <- pen2

# Create two subsets of explanatory variables
# Physiography (upstream-downstream gradient)
envtopo <- env[,c(1:3)]

# Water quality
envchem <- env[,c(4:10)]

# Hellinger-transform the species dataset
spe.hel <- decostand(spe,"hellinger")

# RDA of the Hellinger-transformed fish species data constrained by all the environmental variables contained in env2
# Observe the shortcut formula

spe.rda <- rda(spe.hel~., env2)
summary(spe.rda)

# see text for explanations and interpretations

# How to obtain canonical coefficients from an rda() object
coef(spe.rda)

# We may need an adjusted R2 (see text), which can be achieved in vegan
# Retrieval of the adjusted R^2
(R2 <- RsquareAdj(spe.rda)$r.squared)
(R2adj <- RsquareAdj(spe.rda)$adj.r.squared)

# "Triplot": sites, response variables, explanatory variables
# *******************************************************
# Scaling 1: distance triplot
plot(spe.rda, scaling=1, main="Triplot RDA spe.hel ~ env2 - scaling 1 - wa scores")
# arrows for species are missing. We can add them wihtout heads to make them different from the explanatory variables
spe.sc <- scores(spe.rda, choices=1:2, scaling=1, display="sp")
arrows(0,0,spe.sc[,1],spe.sc[,2],length=0,lty=1,col="red")

# Scaling 2 (default) : correlation triplot
plot(spe.rda, scaling=2, main="Triplot RDA spe.hel ~ env2 - scaling 2 - wa scores")
spe2.sc <- scores(spe.rda, choices=1:2, scaling=2, display="sp")
arrows(0,0,spe2.sc[,1],spe2.sc[,2],length=0,lty=1,col="red")

# weighted site scores vs. fitted site scores is a controversial issue: for fitted site scores, use display="lc"
# see text for explanation
# interpretation must be preceded by a test of statistical significance

# Site scores as linear combinations of the environmental variables (fitted site scores)
# Scaling 1
plot(spe.rda,scaling=1, display=c("sp","lc","cn"), main="Triplot RDA spe.hel ~ env2 - scaling 1 - lc scores")
arrows(0,0,spe.sc[,1],spe.sc[,2],length=0,lty=1,col="red")
# Scaling 2
plot(spe.rda,scaling=2, display=c("sp","lc","cn"), main="Triplot RDA spe.hel ~ env2 - scaling 1 - lc scores")
arrows(0,0,spe2.sc[,1],spe2.sc[,2],length=0,lty=1,col="red")

# Interpretation of the two scalings similar to PCA. However, the presence of vectors and centroids of
# explanatory variables call for additional interpretation rules (see the text for good explanation, p. 167)

# Global test of the RDA results
anova.cca(spe.rda, step=1000)
# Tests of all canonical axes
anova.cca(spe.rda, by="axis", step=1000)

# Apply Kaiser-Guttmna criterion to residual axes
spe.rda$CA$eig[spe.rda$CA$eig > mean(spe.rda$CA$eig)]

### 6.3.2.5 Partial RDA ###
# partial canonical ordination is the multivariate equivalent of partial linear regression
# for example, it is possible to run an RDA of species data, with explanatory variables, and also in the presence of
# covariables; allowing the user to display the patterns of species data uniquely explained by the explanatory
# variables, when the effect of the covariables is held constant.
# We'll do this with the physiographical variables from the Doubs set (as covariates), and the chem variables as explanatory

# Partial RDA: effect of water chemistry, holding physiography constat
# Simple interface; X and W may be separate tables of quantitative variables
spechem.physio <- rda(spe.hel, envchem, envtopo)
spechem.physio
summary(spechem.physio)
spechem.physio2 <- rda(spe.hel~pH+har+pho+nit+amm+oxy+bdo + Condition(alt + slo+ flo), data=env) # same result

# Test of the partial RDA (using the results with the formula interface to allow the tests of the axes to be run)
anova.cca(spechem.physio2,step=1000)
anova.cca(spechem.physio2,step=1000, by="axis")

# Partial RDA triplots (with fitted site scores)

# Scaling 1
plot(spechem.physio,scaling=1,display=c("sp","lc","cn"), main="Triplot RDA spe.hel ~ chem | Topo - scaling 1 - lc scores")
spe3.sc <- scores(spechem.physio, choices=1:2,scaling=1,display="sp")
arrows(0,0,spe3.sc[,1],spe3.sc[,2],length=0,lty=1,col="red")
# Scaling 2
plot(spechem.physio,scaling=2,display=c("sp","lc","cn"), main="Triplot RDA spe.hel ~ chem | Topo - scaling 2 - lc scores")
spe4.sc <- scores(spechem.physio, choices=1:2,scaling=2,display="sp")
arrows(0,0,spe4.sc[,1],spe4.sc[,2],length=0,lty=1,col="red")

# Variance inflation factors (VIF) in two RDAs, to reduce the number of explanatory variables
# if the VIFs are above 20, indicates strong collinearity. Ideally, VIFs above 10 should be at least examined, and
# avoided if possible
# *********************************************
# First RDA of the Chapter: all env vars
vif.cca(spe.rda) 
vif.cca(spechem.physio) # and the partial RDA

# How do we select the most important explanatory variables? Forward, backward, or stepwise selection of variables
# Forward selection of explanatory variables
# using a double stopping criterion (Blanchet et al. 2008a, see text p. 176-77)
# 1. RDA with all explanatory variables
spe.rda.all <- rda(spe.hel~.,data=env)
# 2. Global adjusted R^2
(R2a.all <- RsquareAdj(spe.rda.all)$adj.r.squared)
# 3. Forward selection using packfor's forward.sel()
forward.sel(spe.hel,env,adjR2thresh = R2a.all)
# the last variable exceeds the stopping criteria (see text), and should not be used
# another function, vegan's ordistep(), does a similar thing, but does not implement Blanchet's second stopping criterion
# Forward selection using ordistep()
# This function allows the use of factors. Options are also available for stepwise and backward selection of variables
step.forward <- ordistep(rda(spe.hel~1,data=env), scope=formula(spe.rda.all),direction="forward",pstep=1000)

# you can look manually for when the stopping criterion is exceeded (see text p. 178)

# But, now with just the parsimonious variables, we can re-do the analysis
# *******************************************
spe.rda.pars <- rda(spe.hel~alt+oxy+bdo, data=env)
spe.rda.pars
anova.cca(spe.rda.pars,step=1000)
anova.cca(spe.rda.pars,step=1000,by="axis")
vif.cca(spe.rda.pars)
(R2a.pars<-RsquareAdj(spe.rda.pars)$adj.r.squared)
# This process greatly improved the quality of the model (text p. 178)

# Triplots of the parsimonious RDA (with fitted site scores)
# ***********************************************************

#Scaling 1
plot(spe.rda.pars, scaling=1, display=c("sp","lc","cn"), main="Triplot RDA spe.hel ~ alt+oxy+bdo - scaling 1 - lc scores")
spe4.sc <- scores(spe.rda.pars, choices=1:2,scaling=1,display="sp")
arrows(0,0,spe4.sc[,1],spe4.sc[,2],length=0,lty=1,col="red")
#Scaling 2
plot(spe.rda.pars, scaling=2, display=c("sp","lc","cn"), main="Triplot RDA spe.hel ~ alt+oxy+bdo - scaling 2 - lc scores")
spe5.sc <- scores(spe.rda.pars, choices=1:2,scaling=2,display="sp")
arrows(0,0,spe5.sc[,1],spe5.sc[,2],length=0,lty=1,col="red")
# NOTE: since there is now a third significant canonical axis, you could plot other combinations: axes 1 and 3, axes 2 and 3

### 6.3.2.7 Variation Partitioning ###
# We may be interested in not only a partial analysis of our choice, but in quantifying the vation explained by
# all subset of the variables when controlling for the effect of the others (Fig. 6.5 in text)
# Conceptual steps:
  # if necessary, forward-select the explanatory variables separately in each subset
  # Run an RDA of the response data Y by X, yielding [a+b]
  # Run an RDA of the response data Y by W, yielding [b+c]
  # Run an RDA of the response data Y by X and W together, yielding [a+b+c]
  # Compute the adj R^2 of all of the above
  # Compute the fractions of adjusted variation by subtraction

# The whole procedure can be run in one R command with up to four explanatory matrices, with vegan's varpart()
# Variation partitioning with two set of explanatory variables
# *************************************************************

# Explanation of fraction labels
par(mfrow=c(1,3))
showvarparts(2)
showvarparts(3)
showvarparts(4) #... up to four explanatory matrices

# 1. variation partitioning with all explanatory variables
spe.part.all <- varpart(spe.hel,envchem, envtopo)
spe.part.all
plot(spe.part.all,digits=2) # venn circles not to scale

# 2. separate forward selection in each subset of environmental variables
spe.chem <- rda(spe.hel, envchem)
R2a.all.chem <- RsquareAdj(spe.chem)$adj.r.squared
forward.sel(spe.hel,envchem,adjR2thresh=R2a.all.chem,nperm=9999)
spe.topo <- rda(spe.hel,envtopo)
R2a.all.topo <- RsquareAdj(spe.topo)$adj.r.squared
forward.sel(spe.hel,envtopo,adjR2thresh=R2a.all.topo, nperm=9999)

# Parsimonious subsets of explanatory variables (based on forward selection)
envchem.pars <- envchem[,c(4,6,7)]
envtopo.pars <- envtopo[,c(1,2)]

# Variation partitioning
(spe.part <- varpart(spe.hel,envchem.pars,envtopo.pars))
plot(spe.part,digits=2)

# Tests of all testable fractions
# Test of fractions [a+b]
anova.cca(rda(spe.hel,envchem.pars),step=1000)
# Test of fractions [b+c]
anova.cca(rda(spe.hel,envtopo.pars),step=1000)
# Test of fractions [a+b+c]
env.pars<-cbind(envchem.pars,envtopo.pars)
anova.cca(rda(spe.hel,env.pars),step=1000)
# Test of fractions [a]
anova.cca(rda(spe.hel,envchem.pars,envtopo.pars),step=1000)
# Test of fraction [c]
anova.cca(rda(spe.hel,envtopo.pars,envchem.pars),step=1000)

# 3. Variation partitioning without the 'nit' variable
envchem.pars2 <- envchem[,c(6,7)]
(spe.part2 <- varpart(spe.hel,envchem.pars2, envtopo.pars))
plot(spe.part2,digits=2)
# Some interpretation and lessons--text p. 184-185 

## 6.3.2.8 RDA as a tool for multivariate ANOVA
# fictitious example using the fish data, coding factor variables as explanatory RDA variables
# first, a factor representing altitude, with three levels; and a fictitious, orthogonal pH variable

# Two-way MANOVA by RDA
# *******************************
# Creation of altitude factor (3 levels, 9 sites each)
alt.fac <- gl(3,9)
# Creation of factor mimicking pH
pH.fac <- as.factor(c(1,2,3,2,3,1,3,2,1,2,1,3,3,2,1,1,2,3,2,1,2,3,2,1,1,3,3))
# Are the factors balanced?
table(alt.fac,pH.fac)

# Creation of Helmert contrasts for the factors and their interaction
alt.pH.helm <- model.matrix(~alt.fac*pH.fac,contrasts=list(alt.fac="contr.helmert",pH.fac="contr.helmert"))
alt.pH.helm

# Check property 1 of Helmert contrasts: all variables sum to 0
apply(alt.pH.helm[,2:9],2,sum)
# Check property 2 of Helmert contrsats: variables are uncorreleated
cor(alt.pH.helm[,2:9])
# Verify multivariate homogeneity of within-group covariance matrices using the betadisper() fxn, implementing
# Marti Anderson's testing method
spe.hel.d1 <- dist(spe.hel[1:27,])
# Factor 'altitude'
(spe.hel.alt.MHV <- betadisper(spe.hel.d1,alt.fac))
anova(spe.hel.alt.MHV)
permutest(spe.hel.alt.MHV)
# The within-group coveriance matrices are homogeneous. We can proceed

# Test the interaction first. The factors alt and pH form the matrix of covariables
interaction.rda <- rda(spe.hel[1:27,],alt.pH.helm[,6:9],alt.pH.helm[,2:5])
anova(interaction.rda,step=1000,perm.max=1000)
# if the interaction is significant, it would preclude the analysis of the main factors since it would indicate 
# that the effect of a factor depends on the level of the other factor

# Test the main factor alt. The factor pH and the interaction form the matrix of covariables
factor.alt.rda <- rda(spe.hel[1:27,],alt.pH.helm[,2:3],alt.pH.helm[,4:9])
anova(factor.alt.rda,step=1000,perm.max=1000,strata=pH.fac)
# Test the factor pH.
factor.pH.rda <- rda(spe.hel[1:27,],alt.pH.helm[,4:5],alt.pH.helm[,c(2:3,6:9)])
anova(factor.pH.rda,step=1000,perm.max=1000,strata=alt.fac)

# RDA and triplot for the significant factor alt
alt.rda.out <- rda(spe.hel[1:27,]~., as.data.frame(alt.fac))
plot(alt.rda.out,scaling=1, display=c("sp","wa","cn"),main="MANOVA, factor altitude - scaling 1 - wa scores")
spe.manova.sc <- scores(alt.rda.out,choices=1:2,scaling=1,display="sp")
arrows(0,0,spe.manova.sc[,1],spe.manova.sc[,2],length=0,col="red")

### NOTE: Skipped some of the end of 6.3, including non-linear and db-RDA. Can return to these later if necessary

#### Canonical Correspondence Analysis ####
# The canonical counterpart of CA
# shares many characteristics with RDA
# weighted form of RDA applkied to the same matrix Q of contributions to the chi-sq statistic as used in CA
# shortcomings and warnings, text p. 198-99

### 6.4.2 CCA of the Doubs Data ###
## 6.4.2.1 CCA using vegan. Remember, here we use the raw, untransformed abundances

# Canonical Correspondence Analysis
# *********************************

# CCA of the raw fish species data, constrained by all the environmental variables in env2
spe.cca <- cca(spe~., env2)
spe.cca
summary(spe.cca) # scaling 2 (default)
# The differences from an RDA output:
  # variation now expressed as mean squared contingency coefficient
  # Species scores are represented by points in the triplot
  # site scores are weighted averages (instead of weighted sums) of species scores
# CCA Triplot using lc site scores
#****************************************

# Scaling 1: species cores scaled to relative eigenvalues, sites are weighted averages of the species
plot(spe.cca,scaling=1,display=c("sp","lc","cn"),main="Triplot CCA spe ~ env2 - scaling 1")

# Scaling 2: site scores scaled to relative eigenvalues, species are weighted averages of the sites
plot(spe.cca,display=c("sp","lc","cn"),main="Triplot CCA spe ~ env2 - scaling 2")

# CCA scaling 1 biplot without species (which makes the plot harder to see)
plot(spe.cca,scaling=1,display=c("lc","cn"),main="Biplot CCA spe ~ env2 - scaling 1")

# Similarly, CCA scaling 2 biplot without sites
plot(spe.cca,scaling=2,display=c("sp","cn"),main="Biplot CCA spe ~ env2 - scaling 2")

### 6.4.2.3 Permutation Tests in CCa, forward selection, and parsimonious CCA
# CCA results can be testesd for significance by permutations, in the same way as in RDA

# Permutation tests of CCA results
# *********************************

# of the global analysis
anova(spe.cca,step=1000)

# permutation test of each axis
anova(spe.cca,by="axis",step=1000)

# forward selection using vegan's ordistep(), since forward.sel only does RDA
# CCA-based forward selection using vegan's ordistep()
# ****************************************************
# This function allows the use of factors like 'slo' in env2
cca.step.forward <- ordistep(cca(spe~1, data=env2),scope=formula(spe.cca),direction="forward",pstep=1000)
# the result is the same as the most parsimonious one based on RDA. Therefore, we ca compute a parsimonious CCA
# on the basis of the same three varaibles: alt, oxy, and bdo

# Parsimonious CCA using alt, oxy and bdo
# ******************************************

(spe.cca.pars <- cca(spe~alt+oxy+bdo, data=env2))
anova.cca(spe.cca.pars,step=1000)
anova.cca(spe.cca.pars,step=1000,by="axis")

vif.cca(spe.cca)
vif.cca(spe.cca.pars)

### 6.4.2.4 Three-Dimensional Interactive Plots ###
# in the rgl package
#*************************************
library(rgl)
library(vegan3d)
# Plot of the sites only (wa scores)
# *************************************
ordirgl(spe.cca.pars,type="t",scaling=1)
# Connect weighted average scores to linear combination scores
orglspider(spe.cca.pars,scaling=1,col="purple")
# purple connections show how well the CCA model fits the data. The shorter the connections, the better the fit

# Plot the sites (wa scores) with a clustering result
# *******************************************************
# Colour the sites according to cluster membership
gr <- cutree(hclust(vegdist(spe.hel,"euc"),"ward"),4)
ordirgl(spe.cca.pars,type="t",scaling=1,ax.col="black",col=gr+1)
# Connect sites to cluster centoids
orglspider(spe.cca.pars,gr,scaling=1)

# Complete CCA 3D triplot
# ************************
ordirgl(spe.cca.pars,type="t",scaling=2)
orgltext(spe.cca.pars,display="species",type="t",scaling=2,col="cyan")

# Plot species groups (Jaccard similarity, useable in R mode)
# **********************************************************
gs <- cutree(hclust(vegdist(t(spe),method="jaccard"),"ward"),4)
ordirgl(spe.cca.pars,display="species",type="t",col=gs+1)


#### 6.5 Linear Discriminant Analysis ####

# LDA is different from RDA and CCA in that the response variable is agrouping of the sites; maybe from a clustering,
# maybe representing an ecological hypothesis.  LDA tries to determine to what extent an independent set of
# quantitative variables can explain this grouping.  Site typology MUST have been obtained independently from the
# explanatory variables used in teh LDA; otehrwise, the procedure is circular.

# LDA computes discriminant functions from standardized descriptors.  These coefficients quantify the relative
# contributions of the (standardized) explanatory variables to the discrimination of objects.  On the other hand,
# identification functions can be computed from the original descriptors and can be used to find the group to which
# a new object should be attributed.

# Linear Discriminant Analysis #
# ************************************

env.pars2 <- as.matrix(env[,c(1,9,10)])
# Verify multivariate homogeneity of within-group covariance matrices using the betadisper() fxn
env.pars2.d1 <- dist(env.pars2)
(env.MHV <- betadisper(env.pars2.d1,gr))
anova(env.MHV)
permutest(env.MHV)

# within-group covariance matrices are NOT homogeneous.  Let's try a log transformation of variables alt and bdo

env.pars3 <- cbind(log(env$alt),env$oxy, log(env$bdo))
colnames(env.pars3)<-c("alt.ln","oxy","bdo.ln")
row.names(env.pars3) <- row.names(env)
env.pars3.d1 <- dist(env.pars3)
(env.MHV2 <- betadisper(env.pars3.d1,gr))
permutest(env.MHV2)

# This time the within-group covariance matrices are homogeneous. We can proceed.

# Computation of LDA (discrimination)

env.pars3.df <- as.data.frame(env.pars3)
(spe.lda <- lda(gr ~ alt.ln + oxy + bdo.ln, data=env.pars3.df))

# The result object contains the informatino necessary to interpret the LDA
summary(spe.lda)
# Display the group means for the 3 variables
spe.lda$means
# Compute the normalized eigenvectors, which are the standardized discriminant function coefficients
(Cs <- spe.lda$scaling)
# Compute the canonical eigenvalues
spe.lda$svd^2
# Position the objects in the space of the canonical variates
(Fp <- predict(spe.lda)$x)
# Classification of the objects
(spe.class <- predict(spe.lda)$class)
# Posterior probabilities of the objects belonging to the groups
(spe.post <- predict(spe.lda)$posterior)
# Table of prior versus predicted classifications
(spe.table <- table(gr,spe.class))
# Proportion of correct classification
diag(prop.table(spe.table,1))
# Plot the objects in the spce of the canonical variates
# with colours according to their classifications
plot(Fp[,1],Fp[,2],type="n")
text(Fp[,2],Fp[,2],row.names(env),col=c(as.numeric(spe.class)+1))
abline(v=0,lty="dotted")
abline(h=0,lty="dotted")
# Draw 95% ellipses around the groups
for(i in 1:length(levels(as.factor(gr)))) {
  cov <- cov(Fp[gr==i,])
  centre <- apply(Fp[gr==i,],2,mean)
  lines(ellipse(cov,centre=centre,level=0.95))
}

# Classification of a new object (identification)
# a new object is created with ln(alt)=5.8,oxygen=90, and lb(bdo)=3.2
newo <- c(6.8,90,3.2)
newo <- as.data.frame(t(newo))
colnames(newo) <- colnames(env.pars3)
(predict.new <- predict(spe.lda,newdata=newo))

# LDA with jacknife-based classification (i.e., leave-one-out cross-validation)
spe.lda.jac <- lda(gr~alt.ln + oxy +bdo.ln,data=env.pars3.df,CV=TRUE)
summary(spe.lda.jac)

# Numbers and proportions of correct classification
spe.jac.class <- spe.lda.jac$class
spe.jac.table <- table(gr,spe.jac.class)
diag(prop.table(spe.jac.table,1))

#### Symmetrical Analysis of Two or More Data Sets ####
# - no "response" and "explanatory" variables, but rather all datasets are
# treated the same.  It is more exploratory or descriptive, and appropriate when no unidirectional causal hypothesis is
# embedded in the model. p. 211 in text.  Three methods here.

#### 6.8 Canonical Correlation Analysis ####
# aim of the method is to represent the observation along canonical axes that maximize the correlations between the two tables
# solutions is found by maximizing the between-set dispersion, expressed by the covariance matrix between the two sets of
# variables, with respect to the within-set dispersion.  Appropriate for exploratory purposes in cases where the two groups of
# varaibles are likely to influence each other.
# In this example, we can study how chemistry relates to physiography in the doubs data
# variables must be standardized since they have different physical dimensions, and some must be transformed for normality

# Canonical correlation analysis (CCorA)
#*******************************************

# Data prep
envchem2 <- envchem
envchem2$pho <- log(envchem$pho)
envchem2$nit <- sqrt(envchem$nit)
envchem2$amm <- log1p(envchem$amm)
envchem2$bdo <- log(envchem$bdo)

envtopo2 <- envtopo
envtopo2$alt <- log(envtopo$alt)
envtopo2$slo <- log(envtopo$slo)
envtopo2$flo <- sqrt(envtopo$flo)

#CCorA (on standardized variables)
chem.topo.ccora <- CCorA(envchem2,envtopo2,stand.Y=TRUE,stand.X = TRUE,permutations=999)
chem.topo.ccora
biplot(chem.topo.ccora,plot.type="biplot")

#### 6.9 Co-inertia Analysis ####
# Alternative to CCorA.  Compute the coariance matrix crossing the variables of the two tables. Sum of squared covariances is
# the total co-inertia.  Compute the eigenvalues and eigenvectors of that matrix.  The eigenvalues represent a partitioning of
# the total co-inertia.
# Project the points and variables of the two original data tables on the co-inertia axes.  By means of graphs, compare the 
# projections of the two data tables in the common co-inertia space.
# One advantage of CoIA is the possiblility to choose the type of ordination to apply to each data table prior to the analysis.
# We'll use PCA of the standardized data for both the chem and topo data.  Row weights must be equal in the two separate
# ordinations, and we check this.

# Co-inertia Analysis
# *****************************

# PCA on both matrices
dudi.chem <- dudi.pca(envchem2,scale=TRUE,scan=FALSE,nf=3)
dudi.topo <- dudi.pca(envtopo2,scale=TRUE,scan=FALSE,nf=2)
# Relative variation of eigenvalues
dudi.chem$eig/sum(dudi.chem$eig)
# Relative variation of eigenvalues
dudi.topo$eig/sum(dudi.topo$eig)
# Equal row weights in the two analyses?
all.equal(dudi.chem$lw,dudi.topo$lw)

# Co-inertia analysis
coia.chem.topo <- coinertia(dudi.chem,dudi.topo,scan=FALSE,nf=2)
coia.chem.topo
# Relative variation on the first eigenvalue
coia.chem.topo$eig[1]/sum(coia.chem.topo$eig)
summary(coia.chem.topo)
randtest(coia.chem.topo,nrepet=999)
plot(coia.chem.topo)

# See text p. 216 for interpretation

#### 6.10 Multiple Factor Analysis is skipped here, but may be useful (text p. 220)