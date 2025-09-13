rfModelMatrix=function(Z) {
  
  Z.class = class(Z)[1]
  
  if(Z.class=='factor') {
    Z.class='data.frame'
    Z=data.frame(Z=Z)
  }
  
  grp=NULL
  
  if(Z.class=='data.frame') {
    p=dim(Z)[2]
    znm = names(Z)
    is_factor = FALSE
    for(i in 1:p) {
      if(is.factor(Z[[i]])) {
        #delete last column
        Ztemp = fastDummies::dummy_cols(Z[[i]], remove_first_dummy = TRUE)
        Ztemp = as.data.frame(Ztemp[,2:ncol(Ztemp)])
        colnames(Ztemp) = paste(znm[i],1:ncol(Ztemp),sep='')
        Z[[i]]=Ztemp
        grp=c(grp, rep(i, ncol(Ztemp)))
        is_factor = TRUE
      } else {
        Z[[i]]=cbind(Z[[i]])
        colnames(Z[[i]])=znm[i]
        grp=c(grp, i)
      }
    }
    Ztemp=cbind(Z[[1]])
    if(p>1) for(i in 2:p) Ztemp=cbind(Ztemp, Z[[i]])
    #add a column of 1
    # intercept = rep(1, nrow(Ztemp))
    # Ztemp = cbind(Ztemp,intercept)
    # Ztemp = cbind(Ztemp,1)
    Z=Ztemp
  }
  else if(Z.class=='numeric' | Z.class=='integer') {
    Z=cbind(as.numeric(Z))
    grp=1
  }
  else if(Z.class=='NULL') return(Z)
  else if(Z.class!='matrix')
    stop('Expecting either a factor, a vector, a matrix or a data.frame')
  
  Z <- data.matrix(Z)
  
  return(list(Z=Z, grp=grp))
}
