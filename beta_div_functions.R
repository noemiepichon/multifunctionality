###########################################################################################
########################## CqN ############################################################
###########################################################################################
# X is an S by N matrix which represent data for S species and N communities.
# Each element in X represents any measure of species importance.
# The parameters,"from","to", and "interval", are specified for the order q.
# Defaults are from=0,to=3,interval=0.25.
## Author: Anne Chao

CqN <- function(X,from=0,to=3,interval=0.25, plot = FALSE){
  # if(class(X)!="matrix") # changed 09/03/23
    {
    X=as.matrix(X,ncol=N)
  }
  
  if(any(is.na(X))){
    q=seq(from,to,interval)
    output <- matrix(rep(rep(NA, 7), length(q)), nrow=length(q))
    output <- cbind(q, output) }
  
  else{  
    
    I=which(colSums(X)==0)
    if(length(I)>0){
      output <- rep(NA,8)
    }
    
    N=ncol(X)
    q=seq(from,to,interval)
    if(length(colnames(X))==0){
      Name=paste0("Assemblage_",1:N)
    }else{
      Name=colnames(X)
    }
    
    ###### function ######
    qD_gamma=function(y,q){
      y=y/sum(y);y=y[y>0]
      if(q!=1){
        value=sum(y^q)^(1/(1-q))
      }else{
        value=exp(-sum(y*log(y)))
      }
      return(value)
    }
    qD_alpha=function(y,q){
      y=y/sum(y);y=y[y>0]
      if(q!=1){
        value=(sum(y^q)^(1/(1-q)))/N
      }else{
        value=exp(-sum(y*log(y))-log(N))
      }
      return(value)
    }
    
    ###### Calculation ######
    qD=sapply(1:N,function(j){
      sapply(1:length(q),function(k){
        qD_gamma(X[,j],q[k])
      })
    })
    colnames(qD)=paste("qD",Name,sep="_")
    X.pool=rowSums(X)
    qD.gamma=sapply(1:length(q),function(k){qD_gamma(X.pool,q[k])})
    qD.alpha=sapply(1:length(q),function(k){qD_alpha(X,q[k])})
    qD.beta=qD.gamma/qD.alpha
    CqN=ifelse(q!=1,(qD.beta^(1-q)-N^(1-q))/(1-N^(1-q)),
               1-log(qD.beta)/log(N))
    UqN=ifelse(q!=1,(qD.beta^(q-1)-N^(q-1))/(1-N^(q-1)),
               1-log(qD.beta)/log(N))
    
    ###### Output ######
    output=data.frame(q,qD,qD.gamma,qD.alpha,qD.beta,CqN,UqN)
    
    
    ###### Plot ######
    
    if(plot==TRUE){  
      par(mar=c(5,5,4,2)+0.1)
      par(mfrow=c(2,3))
      Max=max(qD)
      Min=min(qD)
      for(j in 1:N){
        if(j==1){
          plot(q,qD[,j],type="l",xlab="Order q",
               ylab=expression(paste(phantom()^q*D)),
               col=j,lty=j,ylim=c(Min,Max),
               main="Hill numbers")
        }else{
          points(q,qD[,j],type="l",col=j,lty=j)
        }
      }
      legend("topright",legend = Name,col=1:N,lty=1:N,bty="n")
      
      plot(q,qD.gamma,type="l",xlab="Order q",ylab=expression(paste(phantom()^q*D[gamma])),main=expression(paste(phantom()^q*D[gamma])))
      plot(q,qD.alpha,type="l",xlab="Order q",ylab=expression(paste(phantom()^q*D[alpha])),main=expression(paste(phantom()^q*D[alpha])))
      plot(q,qD.beta,type="l",xlab="Order q",ylab=expression(paste(phantom()^q*D[beta])),main=expression(paste(phantom()^q*D[beta])))
      plot(q,CqN,type="l",xlab="Order q",ylab=expression(paste("C"[qN])),main=expression(paste("C"[qN])))
      plot(q,UqN,type="l",xlab="Order q",ylab=expression(paste("U"[qN])),main=expression(paste("U"[qN])))}
    
    else{}  
  }  
  return(output)
}



###########################################################################################
########################## gdmplot.extract ################################################
###########################################################################################
### a function to extract fitted values for the relationships from a list of gdms 
### (much of the code is extracted from the function plot.gdm in the gdm package) 
### Author: Eric Allan

gdmplot.extract <- function(mod.list){
  
  ## number of predictors in model (all in list must have same number)
  preds <- length(mod.list[[1]]$predictors)
  
  fits <- rep(list(rep(list(1), preds)), length(mod.list))
  
  for(i in 1:length(mod.list)){
    
    PSAMPLE <- 200
    preddata <- rep(0,times=PSAMPLE)
    
    model <- mod.list[[i]]
    splineindex <- 1
    
    for (j in 1:preds) {
      
      ## only if the sum of the coefficients associated with this predictor is > 0.....
      numsplines <- model$splines[j]
      if ( sum(model$coefficients[splineindex:(splineindex+numsplines-1)]) > 0 ) {
        
        ## get predictor plot Y-data
        z <- .C( "GetPredictorPlotData",
                 pdata = as.double(preddata),
                 as.integer(PSAMPLE),
                 as.double(model$coefficients[splineindex:(splineindex+numsplines-1)]),
                 as.double(model$knots[splineindex:(splineindex+numsplines-1)]),
                 as.integer(numsplines), PACKAGE = "gdm"
        )
        
        fits[[i]][[j]] <- z$pdata
      }
      else{
        fits[[i]][[j]] <- NA}
      
      splineindex <- splineindex + numsplines
    }
  }
  
  ### create the x values, these are the same for each group
  xvals <- rep(list(1), preds)
  
  for(i in 1:preds){
    
    model <- mod.list[[1]]
    
    xvals[[i]] <- seq(from=model$knots[[(i*3)-2]],to=model$knots[[(i*3)]], length=200)
    
  }
  
  ### collect output and return a list of matrices, each element is for 1 predictor and the matrix has the x values and the predictions for each group
  
  xnames <- gsub("Plot1.","",mod.list[[1]]$predictors)
  
  output <- list()
  for(i in 1:preds){
    oo <- do.call(cbind,sapply(fits,"[",i))
    colnames(oo) <- names(mod.list)
    oo2 <- cbind("x"=xvals[[i]],oo)
    output[[i]] <- oo2
  }
  names(output) <- xnames
  
  return(output)
}