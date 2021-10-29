#### a function to calculate multidiversity
#### threshold can be a proportion of the maximum e.g. 0.5 or can be the median or mean
#### it can also be a vector of the same length as the number of diversities to
#### allow different thresholds and therefore different weightings for different diversities, 
#### setting the threshold to NA drops the group from the calculation of multidiversity
#### threshold = FALSE calculates the mean of the scaled diversities or functions
#### scaling can be by any function specified in "sc", i.e. max, mean, sd etc., max is the default
#### centering by the mean is possible with cent = TRUE, to take z-scores of the diversities/processes, 
#### use sc="sd" and cent = TRUE 
#### "by" specifies a vector of the same length as nrow(x) to be used to split the data and use different
#### thresholds for the groups in "by"
#### "weights" allows different weightings for the different groups and should be a vector of the same length 
#### as the number of groups


maxx <- function(x, ...){     
  return(mean(rev(sort(x))[1:5], ...))}
qtle <- function(x, ...){
  return(quantile(x, probs=.95, names = F, ...))} #0.95 or 0.975


multidiv <- function(x, threshold=FALSE, sc = "max", cent = FALSE, by =FALSE, weights = FALSE){
  
  result <- matrix(nrow=nrow(x),ncol=2) # 3 instead of 2
  
  if(any(by!=FALSE)){
    
    xs <- split(x, by)
    xst <- list()
    
    for(i in 1:length(unique(by))){
      xst[[i]] <- scale(xs[[i]], scale = apply(xs[[i]], 2, match.fun(sc), na.rm=T), center = cent) 
    }
    x.stand <- do.call("rbind", xst)
  }
  
  else{
    x.stand <- scale(x, scale = apply(x, 2, match.fun(sc), na.rm=T), center = cent)
  }
  
  if(length(threshold) == 1){
    gm <- apply(x, 1, function(x)(sum(complete.cases(x))))
  }
  
  else{
    x2 <- sweep(x, 2, threshold, "*")  ## remove diversities with NA threshold from calc. of how many measured
    gm <- apply(x2, 1, function(x)(sum(complete.cases(x))))
  }
  
  if(FALSE %in% threshold){   ### prevent error message if threshold is vector
    m <- apply(x.stand, 1, mean, na.rm=T)
  }
  
  else{
    
    if (any(c("median","mean") %in% threshold)){
      
      tf <- match.fun(threshold)
      
      x.thresh <- apply(x.stand, 2, function(x)(1*(x > tf(x, na.rm=T))))
      gg <- apply(x.thresh, 1, sum, na.rm=T)
      m <- gg/gm
      
    }
    
    else{ 
      
      x.thresh <- 1*sweep(x.stand,2,threshold, ">")  ### does each variable pass threshold?
      
      if(any(weights!=FALSE)){
        weights2 <- matrix(rep(weights, nrow(x.thresh)), nrow=nrow(x.thresh), byrow=TRUE)
        x.thresh <- x.thresh*weights2   ### multiply by weights
      }
      else{}
      gg <- apply(x.thresh, 1, sum,na.rm=T)
      m <- gg/gm
    }}
  
  result[, 1] <- m
  result[, 2] <- gm
  # result[, 3] <- gg # this is new
  
  colnames(result)<-c("m", "groups measured")#, "gg")
  
  return(result)
  # return(x.thresh)
}
