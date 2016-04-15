## NUMERICAL ECOLOGY CHAPTER 5: UNCONSTRAINED ORDINATION ##

## While cluster analysis looks for discontinuities in a dataset, ordination extracts the main trends
## in the form of continuous axes.
## Objectives
## Learn how to choose among various ordination technique (PCA, CA, PCoA, NMDS), compute and interpret them
## Apply these techniques to the Doubs data
## Overlay the result of a cluster analysis on an ordination diagram to improve interpretation of clustering
## Interpret the sturcture in the species data using the env variables
## Write you own PCA function

#### 5.2 Ordination Overview ####
## 5.2.1 Multidimensional space
## variables = dimensions, which are often more than 2
# The aim of ordination methods, therefore, is to represent the data along a reduced number of
# orthogonal axes, contructed in such a way that they represent, in decreasing order, the main
# trends of the data.
# The methods in this chapter are descriptive: i.e. no statistical test is provided to assess 
# the significance of the structures detected. For that, see Ch. 6.

## 5.2.2 Ordination in Reduced Space
# Most methods are based on the extraction of the eigenvectos of an association matrix