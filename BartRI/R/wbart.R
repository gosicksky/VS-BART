wbart = function(x.train, 
                 y.train, 
                 z.train,
                 id,
                 x.test=matrix(0.0,0,0),
                 sparse=FALSE, 
                 theta=0, 
                 omega=1, 
                 a=0.5, 
                 b=1, 
                 augment=FALSE, 
                 rho=NULL,
                 xinfo=matrix(0.0,0,0), 
                 numcut=100L, 
                 usequants=FALSE, 
                 cont=FALSE, 
                 rm.const=TRUE, 
                 grp=NULL, 
                 xnames=NULL, 
                 categorical.idx=NULL,
                 power=2.0, 
                 base=-1.0, 
                 split.prob="polynomial",
                 k=2.0, 
                 sigmaf=NA, 
                 sigest=NA, 
                 sigdf=3, 
                 sigquant=.90, 
                 lambda=NA,
                 fmean=mean(y.train), 
                 w=rep(1,length(y.train)),
                 ntree=200L, 
                 ndpost=1000L, 
                 nskip=1000L, 
                 keepevery=1L,
                 nkeeptrain=ndpost, 
                 nkeeptest=ndpost, 
                 nkeeptestmean=ndpost, 
                 nkeeptreedraws=ndpost,
                 printevery=1000L, 
                 transposed=FALSE, 
                 z_c0, z_d0,
                 z_gamma_mean, z_gamma_cov,
                 z_lambda_mean, z_lambda_cov,
                 z_alpha,z_beta,
                 seed=1,
                 verbose=FALSE
                 ) {
  #--------------------------------------------------
  #data
  N = length(y.train)
  x_lm <- x.train
  
  #the start index for each group (id)
  start_index = as.numeric(c(1, cumsum(table(id))+1) - 1)
  #as.numeric for each element in start_index
  start_index = as.numeric(start_index)
  group_num = length(start_index) - 1
  
  print("Data loaded")
  
  if(!transposed) {
    if(is.null(dim(x.train))) {
      xnames = "X"
    } else {
      xnames = dimnames(x.train)[[2]]
    }
    
    print("start X transfer")
    
    temp = bartModelMatrix(x.train, numcut, usequants = usequants,
                           cont = cont, xinfo = xinfo, rm.const = rm.const)
    
    x.train = t(temp$X)
    numcut = temp$numcut
    xinfo = temp$xinfo
    if(length(x.test) > 0) {
      x.test = bartModelMatrix(x.test)
      x.test = t(x.test[ , temp$rm.const])
    }
    rm.const = temp$rm.const
    grp = temp$grp
    print(grp)
    rm(temp)
    
    if(length(grp) == 0){
      p0 = nrow(x.train)
      grp = 1:p0
    } else {
      p0 = length(unique(grp))  # number of predictors before dummification
    }
    categorical.idx = unique(grp)[which(sapply(unique(grp), function(s) sum(s==grp)) > 1)]
    
    print("X tranfer")
    
    z_temp = rfModelMatrix(z.train)
    z.train = t(z_temp$Z)
    z_grp = z_temp$grp
    
    if(length(z_grp) == 0){
      q0 = nrow(z.train)
      z_grp = 1:q0
    } else {
      q0 = length(unique(z_grp))
    }
    
    print("Z tranfer")
    
  } else {
    if(any(length(rm.const) == 0, length(grp) == 0, length(xnames) == 0))
      stop('Did not provide rm.const, grp and xnames for x.train after transpose!')
    if(is.logical(rm.const))
      stop('Did not provide rm.const for x.train after transpose!')
    if((length(grp) > length(unique(grp))) & (length(categorical.idx) <= 0))
      stop('Did not provide categorical.idx for x.train that contains categorical predictors!')
    
    p0 = length(unique(grp))
  }
  
  if(N != ncol(x.train))
    stop('The length of y.train and the number of rows in x.train must be identical')
  if(N != ncol(z.train))
    stop('The length of y.train and the number of rows in z.train must be identical')
  
  p = nrow(x.train)
  np = ncol(x.test)
  q = nrow(z.train)
  if(length(rho) == 0) rho=p
  
  if(!(split.prob %in% c("polynomial", "exponential"))) {
    stop("split.prob is either polynomial or exponential.")
  } else {
    if(split.prob == "polynomial") {
      if(base < 0)
        base = 0.95
    }
    if(split.prob == "exponential") {
      power = -1.0
      if(base < 0)
        base = 0.5
    }
  }
  
  #fmean should clutered by id
  
  #fmean = ave(y.train, id, FUN = mean)
  
  y.train = y.train - fmean
  
  #--------------------------------------------------
  #set nkeeps for thinning
  if((nkeeptrain != 0) & ((ndpost %% nkeeptrain) != 0)) {
    nkeeptrain = ndpost
    cat('*****nkeeptrain set to ndpost\n')
  }
  if((nkeeptest != 0) & ((ndpost %% nkeeptest) != 0)) {
    nkeeptest = ndpost
    cat('*****nkeeptest set to ndpost\n')
  }
  if((nkeeptestmean != 0) & ((ndpost %% nkeeptestmean) != 0)) {
    nkeeptestmean = ndpost
    cat('*****nkeeptestmean set to ndpost\n')
  }
  if((nkeeptreedraws != 0) & ((ndpost %% nkeeptreedraws) != 0)) {
    nkeeptreedraws = ndpost
    cat('*****nkeeptreedraws set to ndpost\n')
  }
  
  #--------------------------------------------------
  print("Calculate prior")
  prior_indicator = 0
  nu = sigdf
  if(is.na(lambda)) {
    if(is.na(sigest)) {
      if(p < N & q * length(unique(id)) < N) {
        #df = data.frame(t(x.train),t(z.train), y.train)
        #lmf = lm(y.train~., df)
        prior_indicator = 1
        colnames(x_lm) <- paste0("x", seq_len(ncol(x_lm)))
        z_lm <- t(z.train)
        colnames(z_lm) <- paste0("z", seq_len(ncol(z_lm)))
        data_lm <- as.data.frame(cbind(y.train, x_lm, z_lm, id = id))
        formula <- as.formula(paste("y.train ~", paste(colnames(x_lm), collapse = " + "), "+(", paste(colnames(z_lm), collapse = " + "), "-1 | id)"))
        lmf <- lmer(formula, data = data_lm)
        covariance_matrix <- as.matrix(bdiag(VarCorr(lmf)))
        print("covariance_matrix")
        print(covariance_matrix)
        # 
        # decomposed <- modified_cholesky_decomposition(covariance_matrix)
        # 
        # z.train <- z.train[order(decomposed$order), ]
        # 
        # Lambda <- decomposed$Lambda
        # z_lambda_mean <- diag(Lambda)
        # zero_index <- which(z_lambda_mean < 1e-16)
        # z_lambda_mean[zero_index] <- 0
        # 
        # Gamma <- decomposed$Gamma
        # #zero index column and row be zero for Gamma
        # for (i in zero_index) {
        #   Gamma[i, ] <- 0
        #   Gamma[, i] <- 0
        #   Gamma[i, i] <- 1
        # }
        # 
        # get_lower_triangular <- function(matrix) {
        #   element <- rep(0, nrow(matrix) * (nrow(matrix) - 1) / 2)
        #   k = 1
        #   for(i in 2:ncol(matrix)) {
        #     for(j in 1:(i - 1)) {
        #       element[k] <- matrix[i, j]
        #       k <- k + 1
        #     }
        #   }
        #   return(element)
        # }
        # z_gamma_mean <- get_lower_triangular(Gamma)
        
        # z_lambda_mean = merTools::REsdExtract(lmf)
        print("z_lambda_mean")
        print(z_lambda_mean)
        
        sigest = summary(lmf)$sigma

        # print("sigma")
        # print(sigest)
        # print("lambda")
        # print(z_lambda_mean)
        # print("gamma")
        # print(z_gamma_mean)
        
        ###
        lambda_indicator = which(diag(covariance_matrix) < 0.05)
      } else {
        sigest = sd(y.train)
      }
    }
    qchi = qchisq(1.0 - sigquant, nu)
    lambda = (sigest * sigest * qchi) / nu #lambda parameter for sigma prior
  } else {
    sigest = sqrt(lambda)
  }
  
  if(is.na(sigmaf)) {
    tau = (max(y.train) - min(y.train)) / (2 * k * sqrt(ntree))
  } else {
    tau = sigmaf / sqrt(ntree)
  }
  
  print("prior calculated")
  
  #--------------------------------------------------
  
  #random effect
  if (prior_indicator == 1) {
    z_b = matrix(rnorm(group_num*q,0,1),group_num,q,byrow=TRUE)
    z_gamma = rnorm(q*(q-1)/2,0,1)
    z_gamma = rep(0, q*(q-1)/2)
    #z_gamma = c(0,1,0,1,0,1,0,0,0,0,1,0,1,1,0)
   # z_gamma_mean = c(0,1,0,1,0,1,0,0,0,0,1,0,1,1,0)
    # z_gamma = c(0,1,0,0,0,0,1,0,1,0)
    # z_gamma_mean = c(0,1,0,0,0,0,1,0,1,0)
    z_lambda = merTools::REsdExtract(lmf)
    # z_lambda_mean = merTools::REsdExtract(lmf)
    # z_lambda = c(0,5,0,0,0,5)
    # z_lambda_mean =  c(0,5,0,0,0,5)
    # z_lambda = c(4,0,1,0.5,2)
    # z_lambda_mean = c(4,0,1,0.5,2)
    # z_lambda = c(1,0.1,2,0.1,1,0.5,2)
    # z_lambda_cov = rep(0.1,q)
    # z_lambda = rtruncnorm(q, a = 0, b = Inf, mean = 0, sd = 1)
    print("z_lambda")
    print(z_lambda)
    print("z_lambda_mean")
    print(z_lambda_mean)
    
    z_lambda_largest = max(z_lambda)
    z_alpha = 0.99 * z_lambda / (z_lambda_largest - 0.99 * z_lambda)
    
    print("z_alpha")
    print(z_alpha)
  }
  else {
    z_b = matrix(rnorm(group_num*q,0,1),group_num,q,byrow=TRUE)
    z_gamma = rnorm(q*(q-1)/2,0,1)
    z_lambda = rtruncnorm(q, a = 0, b = Inf, mean = 0, sd = 1)
  }
  
  #--------------------------------------------------
  
  print("z")
  print(z.train[1])
  
  #call c++ function
  print(getwd())
  ptm = proc.time()
  sourceCpp(paste0(getwd(),'/src/cwbart.cpp'))
  res = cwbart(
    N,  #number of observations in training data
    p,  #dimension of x
    np, #number of observations in test data
    x.train,   #pxn training data x
    y.train,   #pxn training data x
    x.test,   #p*np test data x
    ntree,
    numcut,
    ndpost*keepevery,
    nskip,
    power,
    base,
    tau,
    nu,
    lambda,
    sigest,
    w,
    sparse,
    theta,
    omega,
    grp,
    a,
    b,
    rho,
    augment,
    nkeeptrain,
    nkeeptest,
    nkeeptestmean,
    nkeeptreedraws,
    printevery,
    xinfo,
    z.train, #Nxq training data z, have been transposed
    q, #dimension of z
    group_num, #number of groups
    start_index, #start index for each group
    z_c0, #c0 for sigma
    z_d0, #d0 for sigma
    t(z_b), #random effect
    z_gamma, #start gamma for random effect
    z_gamma_mean, #mean for gamma
    t(z_gamma_cov), #cov for gamma
    z_lambda, #start lambda for random effect
    z_lambda_mean, #mean for lambda
    z_lambda_cov, #cov for lambda
    z_alpha, #probabiliy for lambda
    z_beta, #probabiliy for lambda
    seed
  )
  
  res$proc.time = proc.time() - ptm
  
  #--------------------------------------------------
  #fixed effect returns
  res$mu = fmean
  res$yhat.train.mean = res$yhat.train.mean + fmean
  res$yhat.train = res$yhat.train + fmean 
  res$yhat.test.mean = res$yhat.test.mean + fmean
  res$yhat.test = res$yhat.test + fmean
  res$fix.train.mean = res$fix.train.mean
  res$fix.train = res$fix.train
  res$random.train.mean = res$random.train.mean
  res$random.train = res$random.train
  
  if(nkeeptreedraws > 0)
    names(res$treedraws$cutpoints) = xnames
  
  #--------------------------------------------------
  #importance
  if(length(grp) == length(unique(grp))) {
    ## no dummy variables
    print("No dummy variables")
    print(dimnames(res$varcount)[[2]])
    print(xnames)
    dimnames(res$varcount)[[2]] = as.list(xnames)
    dimnames(res$varprob)[[2]] = as.list(xnames)
    
    ## vip: variable inclusion proportions
    res$vip = colMeans(t(apply(res$varcount, 1, function(s) s / sum(s))))
    
    ## (marginal) posterior variable inclusion probability
    res$pvip = colMeans(res$varcount > 0)
    
    ## posterior s_j's (only in DART)
    res$varprob.mean = colMeans(res$varprob)
    
    ## mi: Metropolis importance
    mr.vecs = lapply(res$mr_vecs, function(s) lapply(s, function(v) v[-1]))  # remove the meaningless first 0
    res$mr_vecs = NULL
    res$mr.vecs = mr.vecs
    mr.mean = matrix(unlist(lapply(mr.vecs, function(s) lapply(s, function(v) ifelse(length(v) > 0, mean(v), 0.0)))), 
                     ncol = p, byrow = TRUE)
    res$mr.mean = mr.mean
    res$mi = colMeans(t(apply(mr.mean, 1, function(s) s / sum(s))))
    names(res$mi) = as.list(xnames)
    dimnames(res$mr.mean)[[2]] = as.list(xnames)
  } else {
    print("Dummy variables detected")
    ## merge importance scores for dummy variables
    varcount = matrix(NA, nrow = nkeeptreedraws, ncol = p0)
    varprob = matrix(NA, nrow = nkeeptreedraws, ncol = p0)
    
    mr.vecs = lapply(res$mr_vecs, function(s) list(s[[1]][-1]))
    #mr0_vecs = lapply(res$mr0_vecs, function(s) list(s[[1]][-1]))
    varcount[, 1] = res$varcount[, 1]
    varprob[, 1] = res$varprob[, 1]
    
    j = 1
    for (l in 2:p) {
      if (grp[l] == grp[l-1]) {
        varcount[, j] = varcount[, j] + res$varcount[, l]
        varprob[, j] = varprob[, j] + res$varprob[, l]
        for (i in 1:nkeeptreedraws) {
          mr.vecs[[i]][[j]] = c(mr.vecs[[i]][[j]], res$mr_vecs[[i]][[l]][-1])
          #mr0_vecs[[i]][[j]] = c(mr0_vecs[[i]][[j]], res$mr0_vecs[[i]][[l]][-1])
        }
      } else {
        j = j + 1
        varcount[, j] = res$varcount[, l]
        varprob[, j] = res$varprob[, l]
        for (i in 1:nkeeptreedraws) {
          mr.vecs[[i]][[j]] = res$mr_vecs[[i]][[l]][-1]
          #mr0_vecs[[i]][[j]] = res$mr0_vecs[[i]][[l]][-1]
        }
      }
    }
    
    print(dimnames(varcount)[[2]])
    print(xnames)
    dimnames(varcount)[[2]] = as.list(xnames)
    dimnames(varprob)[[2]] = as.list(xnames)
    
    res$varcount = varcount
    res$varprob = varprob
    res$mr.vecs = mr.vecs
    res$mr_vecs = NULL
    #res$mr0_vecs = mr0_vecs
    
    ## vip
    res$vip <- colMeans(t(apply(varcount, 1, function(s) {
      total <- sum(s)
      if (total == 0) rep(0, length(s)) else s / total
    })))
    
    ## within-type vip
    within.type.vip = rep(0, p0)
    for (i in 1:nkeeptreedraws) {
      if (sum(varcount[i, categorical.idx]) != 0) {
        within.type.vip[categorical.idx] = within.type.vip[categorical.idx] + 
          varcount[i, categorical.idx] / sum(varcount[i, categorical.idx])
      }
      if (sum(varcount[i, -categorical.idx]) != 0) {
        within.type.vip[-categorical.idx] = within.type.vip[-categorical.idx] + 
          varcount[i, -categorical.idx] / sum(varcount[i, -categorical.idx])
      }
    }
    res$within.type.vip = within.type.vip / nkeeptreedraws
    names(res$within.type.vip) = xnames
    
    ## (marginal) posterior variable inclusion probability
    res$pvip = colMeans(varcount > 0)
    
    ## posterior s_j's (in DART)
    res$varprob.mean = colMeans(varprob)
    
    ## mi
    mr.mean = matrix(unlist(lapply(mr.vecs, function(s) lapply(s, function(v) ifelse(length(v) > 0, mean(v), 0.0)))), 
                     ncol = p0, byrow = TRUE)
    res$mr.mean = mr.mean
    res$mi = colMeans(t(apply(mr.mean, 1, function(s) {
      total <- sum(s)
      if (total == 0) rep(0, length(s)) else s / total
    })))
    dimnames(res$mr.mean)[[2]] = as.list(xnames)
    names(res$mi) = as.list(xnames)
  }
  
  #--------------------------------------------------
  #random effect
  #combine the lambda, becasue lambda >= 0
  if (length(z_grp) == length(unique(z_grp))){
    res$combined.lambda = res$rf_lambda
  }
  else{
    z_lambda = matrix(NA, nrow = nkeeptreedraws, ncol = q0)
    z_lambda[, 1] = res$rf_lambda[, 1]
    j = 1
    for (l in 2:(q-1)) {
      if (z_grp[l] == z_grp[l-1]) {
        z_lambda[, j] = z_lambda[, j] + res$rf_lambda[, l]
      } else {
        j = j + 1
        z_lambda[, j] = res$rf_lambda[, l]
      }
    }
    
    res$combined.lambda = z_lambda
  }
  
  res$original.lambda = res$rf_lambda
  res$rf_lambda = NULL
  res$original.gamma = res$rf_gamma
  res$rf_gamma = NULL
  # res$rf_p = res$rf_p
  
  #return list
  # res$D <- matrix(sapply(1:ndpost, function(i){
  #   Gamma = lower_matrix(res$original.gamma[i,])
  #   Lambda = diag(res$original.lambda[i,])
  #   return(diag(Lambda %*% Gamma %*% t(Gamma) %*% Lambda))
  # }), nrow = ndpost, ncol = q, byrow = TRUE)
  
  res$rm.const = rm.const
  
  attr(res, 'class') = 'wbart'
  return(res)
}
