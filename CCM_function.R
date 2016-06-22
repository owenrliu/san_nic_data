### CCM function, based on Sugihara et al. 2012 Science

predict_Y <- function(X,Y,E,tau=1) {
  require(dplyr)
  L <- length(X)
  tmin <- 1+(E-1)*tau
  tmax <- L
  tvec <- tmin:tmax
  
  
  Mx <- matrix(nrow=(tmax),ncol=E)
  for(i in 1:E) Mx[,i] <- lead(X,(tmin-((i-1)*tau)-1))
  Mx <- Mx[complete.cases(Mx),]
  
  # finding nearest neighbors
  n.neighbors <- E+1
  Mx.dist <- as.matrix(dist(Mx,upper=T))
  Ypred <- Y[tvec] # base vector of Y values from which predictions will be made
  Yhat <- Ypred # vector of predicted Y values

  for(i in 1:length(tvec)) {
    ordered.neighbors <- order(Mx.dist[i,]) # ordered indices of closest neighbors
    
    # If prediction set (Y) has missing values, remove these time indices from the library so it does not use
    # NA values for prediction
    ordered.neighbors <- ordered.neighbors[!(ordered.neighbors %in% which(is.na(Ypred)))]
    
    neighbors.t <- ordered.neighbors[2:(n.neighbors+1)] #removes closest neighbor (the point itself)
    
    # exponential weights for Yti
    ui <- exp(-Mx.dist[i,neighbors.t]/Mx.dist[i,neighbors.t[1]])
    wi <- ui/sum(ui)
    
    # predicted Y
    Y.ti <- Ypred[neighbors.t] # pull out analagous points on Y from which to estimate Yt
    Yhat[i] <- sum(wi*Y.ti) # the estimated Y value is the weighted sum of contemporaneous nearest neighbors
  }
  
  return(Yhat)
}

run_CCM <- function(X,Y,E,tau) {
  minL <- tau*(E-1)+3
  maxL <- length(X)-E+2
  tmin <- 1+(E-1)*tau
  Xlist <- lapply((minL:maxL),function(i) X[1:i])
  Ylist <- lapply((minL:maxL),function(i) X[1:i])
  corrs <- data_frame(L=minL:maxL,corr=NA)
  for(i in 1:(maxL-minL+1)) {
    lY <- length(Ylist[[i]])
    predicted.Y<- predict_Y(Xlist[[i]],Ylist[[i]],E,tau)
    corrs$corr[i] <- cor(Ylist[[i]][tmin:lY],predicted.Y)
  }
  return(corrs)
}