#4.5 mcRFop_cls:  backward variable elimination for RF classification using repeated importance averaging (nRep runs per iteration) to reduce stochastic bias in variable ranking caused by multicollinearity among climate predictors.
# x is predictor matrix; 
# y is dependent variable (classes)
# nTree= 100 as default, but can be changed
# nRep= 10 as default; number of repeated RF runs per iteration to stabilize importance ranking

mcRFop_cls <- function(x,y,nTree=100,nRep=10){
  
  library("foreach"); library("doSNOW");library(randomForest)
  nCore <- parallel::detectCores(); nCore = nCore-1
  cl <- makeCluster(nCore, type="SOCK")
  registerDoSNOW(cl)
  on.exit(stopCluster(cl))
  
  ## backward purging -----
  clmListAll <- mat.or.vec(ncol(x),2); colnames(clmListAll) <- c('Accy','variable');clmListAll
  varOut=""	  
  fewvarOut=""	   
  
  for(k in 1:ncol(x)-2){
    
    if(ncol(x) <= 2 && ncol(x) %% 2 == 0) {
      stop("Insufficient number of variables in the input data. Exiting.")
    } else if (ncol(x) <= 3 && ncol(x) %% 2 == 1) {
      stop("Insufficient number of variables in the input data. Exiting.")
    }
    
    x <- x[, !names(x) %in% varOut];names(x)
    
    print(paste("length of x: ", ncol(x)))
    print(paste("nTree: ", nTree))
    print(paste("nCore: ", nCore))
    print(paste("length of rep vector: ", length(rep(floor(nTree/nCore),nCore))))
    
    ## repeated RF to stabilize importance -----
    imp_accum <- NULL
    for(r in 1:nRep){
      rf <- foreach(ntree=rep(floor(nTree/nCore),nCore),.combine=combine,.packages="randomForest") %dopar% randomForest(x,y,ntree= ntree,importance=T)
      imp_r <- importance(rf)[,3:4]
      if(is.null(imp_accum)) imp_accum <- imp_r else imp_accum <- imp_accum + imp_r
    }
    imp <- imp_accum / nRep
    imp <- imp[order(imp[,1]),];head(imp)
    varOut <- row.names(imp)[1:2];varOut
    clmList <- rownames(imp);clmList
    clmListAll[nrow(imp),2] <- toString(clmList);clmListAll
    
    ##test-----
    td <- data.frame(y,p=predict(rf));head(td);str(td)
    td1 <- subset(td,y==p);dim(td1);dim(td);
    accy <- dim(td1)/dim(td);accy
    clmListAll[nrow(imp),1] <- accy[1];clmListAll
    
    if (ncol(x) <= 4 && ncol(x) %% 2 == 0) {break} 
    if (ncol(x) <= 3 && ncol(x) %% 2 == 1) {break}
  }
  return(clmListAll)
}