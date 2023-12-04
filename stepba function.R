#################################### stepba (stepWISE ELIMINATION OF NON-SIGNIFICANT TERMS BY baSCHI)
stepba <- function(mod=mod, nogo="", test="Chisq", print_results=FALSE, interactions="first", interaction.depth=1) {
  require(lme4)
  namex <- character(1)                       # name of the term with the highest P-value
  nogo.bak <- nogo
  if(interactions=="only"){
    nogo <- as.character(combn(unique(strsplit(paste(rownames(anova(mod)),collapse = ":"),split = ":")[[1]]), m = interaction.depth, FUN = function(x){paste(x, collapse = ":")}))
  }      # when main effects should be kept (interactions="only") replace "nogo" with the list of main effects or interactions to the degree of m!
  if(length(nogo)<1) nogo <- character(1)     # Check if nogo is defined
  z <- 1                                      # a running number
  # variables to store the results are defined
  listname <- character()
  listp <- character()
  listaic <- character()
  liststat <- character()
  for(i in 1:1000){                           # Start of the loop to remove non significant terms
    dr<-drop1(mod, test=test)                 # drop1: Which terms are (non-)significant
    if(print_results) print(dr)               # Print the drop1 table?
    
    namex <- rownames(dr)[order(dr[,ncol(dr)], decreasing = TRUE)[1]]   # Get the name of the term with the highest P-value
    
    
    if(dr[namex, ncol(dr) ]<=.05 ) break       # if there is nothing to be removed stop the loop
    
    
    nogox <- nogo                             # nogox is the nogo term(s), it might be expanded ... . It is newly defined in every iteration.
    
    # ... if interactions="first" we want to eliminate all non-significant interaction before touching the main-effects
    if( length(strsplit(namex, split = ":")[[1]]) == 1 & interactions=="first" ){  # is the namex a main-effect? and do we care?
      x <- row.names(dr)[-c(which(dr[,ncol(dr)]<.05), 1, which(row.names(dr)==namex)) ] # list of the terms which are not significant and are not the namex (which is a main-effect).
      if(length(x)>0){                        # is there any?
        xx <- numeric()                       # define a variable which will record the number of terms in interactions (or not)
        for(q in 1:length(x)){                # loop which will get xx
          xx[q] <- length ( strsplit(x[q], split = ":")[[1]])
        }
        if (sum(xx) > length(xx)) nogox <- c(nogo, namex, x[xx<2]) # when there is a non-significant interaction add the focal main-effect and the other main effects on the nogox list.
      }
    }
    
    if( namex %in% nogox ){  # is the term to be eliminated on the nogox list? If yes select another term
      x <- dr[!(rownames(dr) %in% nogox) & rownames(dr)!= "<none>",]  # exclude all terms on the nogox list
      x <- x[order(x[, ncol(dr)], decreasing=TRUE ),]             # order the list so that the non significant are on top.
      if(x[1,ncol(x)] > 0.05) {   # check if the one on top is non-significant
        namex <- rownames(x)[1]   # if yes -> make it the namex
      } else {break}              # if no -> break the loop
    }
    
    if(length(namex)==0){ break}      # nothing to eleminate?
    if(namex=="<none>"){ break }      # the same in other words
    z <- i                            # index z redefined
    #     for documentation: which term was eliminated because of which P-value, stat, AIC
    listname[i] <- namex
    listp[i] <- dr[row.names(dr)==namex, ncol(dr)]
    liststat[i] <- dr[row.names(dr)==namex, ncol(dr)-1]
    listaic[i] <- dr[row.names(dr)==namex,"AIC"]
    print(paste("deletation number:",z, ":",namex, ": p:",listp[i] ))   # give an impression what's going on
    #    require(beepr); beep()    # audio feedback (optional)
    r<- as.formula(paste(". ~. -", namex))  # The new formula is the old formula minus the term to be eliminated
    qq <- length( row.names(anova(mod)))    # length of the anova table.
    mod <- update(mod, r)                   # the model is updated with the new fromula
    modf <- mod                             # this new model might be the final one (modf)
    if(qq ==2 ){break}                      # if qq==2 then there is only 1 term left. We leave it at this...
  }
  
  if(length(listname)==0)listname <- listaic <- liststat <- listp <- 0; modf <- mod
  liste <- data.frame(number= (1:length(listname)), factor=listname, AIC=listaic, stat=liststat, P=listp)  # liste is a data frame with the sequence of deleted terms and their respective values
  
  # The result table is a mix between the summary of the final model and the iteration of the deleted terms...
  #define stat and p for Chisq and F tests
  if(test=="Chisq") {stat <- "LRT"; p <-"Pr(Chi)"}
  if(test=="F") {stat <- "F value"; p <- "Pr(>F)"}
  result.table <- data.frame(
    factor = c(rownames(summary(modf)$coef), rownames(dr), as.character(liste[rev(liste$number),"factor" ])),
    Estimate = c(as.data.frame(summary(modf)$coef)[,1], rep("NA", nrow(dr)+nrow(liste))),
    Std.Error = c(as.data.frame(summary(modf)$coef)[,2], rep("NA", nrow(dr)+nrow(liste))),
    AIC = c(rep("NA",nrow(summary(modf)$coef)),  round(dr[, "AIC"], 2), round(as.numeric(as.character(liste[rev(liste$number), "AIC" ])), 2)),
    stat = c(rep("NA",nrow(summary(modf)$coef)), round(dr[, stat], 3), round(as.numeric(as.character(liste[rev(liste$number), "stat" ])), 3)),
    p = c(rep("NA",nrow(summary(modf)$coef)), round(dr[,ncol(dr)], 3), round(as.numeric(as.character(liste[rev(liste$number), "P" ])), 3))
  )
  if(test=="Chisq") names(result.table)[5:6] <- c("LRT", "Pr(Chi)")
  if(test=="F")  names(result.table)[5:6] <- c("F value", "Pr(>F)")
  
  print(paste("number of deletations:", z)) # The final number of iterations
  print(getCall(modf)) # show the call of the final model
  
  # the final output
  return(list( modfinal=modf
               , list = liste
               , dr=dr
               , nogo=nogo
               , call=getCall(modf)
               , result.table = result.table))
}
##########################
