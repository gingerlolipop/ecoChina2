proSysSmpl <- function(x,byCol=3, minSz=2000){
  
  zfrq <- meanFrq(x,byCols=byCol,mean=F,frq=T)
  print(hd(zfrq)) #calculates the frequency of each class
  
  frqMin = min(zfrq[,2])
  print(hd(frqMin))  # find the class with the minimum frequency
  
  if(frqMin > minSz){  # if the minimum frequency is greater than the minimum sample size
    x <- rdSample(x,n=floor((frqMin-minSz)/minSz*nrow(x))) # resampling x
    zfrq <- meanFrq(x,byCols=byCol,mean=F,frq=T) #recalculate frequency after resampling
    print(hd(zfrq))
  }
  
  zfrq <- zfrq[,1:2] #keep only the first 2 columns
  colnames(zfrq) <- c('Zone','frq') 
  frqMax = max(zfrq[,2]) #find maximum freq
  print(hd(frqMax))
  
  frqMean = mean(zfrq[,2]) #mean freq
  print(hd(frqMean))
  
  dev <- log(1/frqMean*(frqMax-frqMin)^1.1) #calculate deviation
  print(hd(dev))
  
  zfrq$sr <- log(1/zfrq$frq*(frqMax-frqMin)^1.1)/dev #calculate sampling ratio for each class
  print(hd(zfrq)) #sr=sapling ratio
  
  vf2 = transform(zfrq,sr=sr/(log(1/minSz*(frqMax-minSz)^1.1)/dev)) #adjust sampling ratio based on minSz
  print(hd(vf2)) 
  
  vf2$sampled <- floor(vf2$frq*vf2$sr) #calc the sampled size for each class
  print(hd(vf2))
  
  vf3 <- vf2[order(vf2$frq),] #order df vf2 by frequency
  print(hd(vf3))
  
 # if(!any(is.nan(vf3$sr))) { 
  #  plot(frq~sr, data=vf3, col='black', xlab='Sampling rate', ylab='Sample size', xlim=c(min(vf3$sr),1)) 
 #   points(sampled~sr, data=vf2,pch=19)
 #   legend('topright', legend = c("Original", "Sampled"),pch=c(1,19))
 #   }
  
  
  sr <- vf3[,c('Zone','sr')]
  print(hd(sr));print(dim(sr)) # #sr=sampling ratio
  names(sr)[1] <- names(x)[byCol]
  
  sd <- dtJoin(x,sr)
  print(hd(sd));print(dim(sd))
  
  sd$rdn <- runif(nrow(sd))
  print(hd(sd))
  
  sp <- subset(sd,rdn<=sr)
  print(hd(sp)); print(dim(sp))
  
  sp2 <- sp[, 1:(ncol(sp)-2)]
  print(hd(sp2))
  
  return(sp2)
}

