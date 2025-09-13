permute.vs = function(x.train, 
                      y.train, 
                      z.train,
                      id,
                      probit=FALSE, 
                      npermute=100L,                      ## number of permutations 
                      nreps=10L,                          ## number of replicates 
                      alpha=0.05,                         ## local threshold
                      true.idx=NULL,
                      plot=TRUE, 
                      n.var.plot=Inf,
                      xinfo=matrix(0.0,0,0), 
                      numcut=100L,
                      usequants=FALSE, 
                      cont=FALSE, 
                      rm.const=TRUE, 
                      k=2.0,
                      power=2.0, 
                      base=0.95,
                      split.prob="polynomial",
                      ntree=20L, 
                      ndpost=1000, 
                      nskip=1000,
                      keepevery=1L, 
                      printevery=100L,
                      z_c0, z_d0,
                      z_gamma_mean, z_gamma_cov,
                      z_lambda_mean, z_lambda_cov,
                      z_alpha, z_beta,
                      verbose=FALSE) {
  
  x.train = as.data.frame(x.train)
  z.train = as.data.frame(z.train)
  
  res = list()
  varcounts = list() ## collect all the varcount matrices
  
  #------------------------------
  # timer starts
  start = Sys.time()
  
  #-----------------------------------------------------------
  # data
  categorical.idx = which(sapply(x.train, function(s) {is.factor(s)}))
  categorical.names = names(categorical.idx)
  
  #-----------------------------------------------------------
  # get avg/median variable importance from the original data
  if(verbose) cat("original data set...")
  
  avg.vip.mtx = matrix(NA, nrow = nreps, ncol = ncol(x.train))
  median.mi.mtx = matrix(NA, nrow = nreps, ncol = ncol(x.train))
  original.z_lambda = list()
  original.z_gamma = list()
  #calculate the average for each lambda to be zero, and if x1 is significant, and x2,x3,x4 not, then the 9th element will be the largest
  #if x1,x3 is significant, and x2,x4 not, then the 11th element will be the largest
  if(ncol(rfModelMatrix(z.train)$Z) > ncol(z.train))
    avg.lambda.mtx = matrix(NA, nrow = nreps, ncol = 2 ^ (ncol(z.train) + 1))
  else
    avg.lambda.mtx = matrix(NA, nrow = nreps, ncol = 2 ^ ncol(z.train))
  
  if (length(categorical.idx) > 0)
    avg.within.type.vip.mtx = matrix(NA, nrow = nreps, ncol = ncol(x.train))
  
  cnt = 0
  while (cnt < nreps) {
    print(cnt)
    bart = wbart(x.train = x.train, y.train = y.train, z.train = z.train, id = id,
                 sparse = FALSE,
                 xinfo = xinfo, numcut = numcut, usequants = usequants, cont = cont, rm.const = rm.const,
                 k = k, power = power, base = base, split.prob = split.prob,
                 ntree = ntree, ndpost = ndpost, nskip = nskip, keepevery = keepevery,
                 z_c0 = z_c0, z_d0 = z_d0,
                 z_gamma_mean = z_gamma_mean, z_gamma_cov = z_gamma_cov,
                 z_lambda_mean = z_lambda_mean, z_lambda_cov = z_lambda_cov, z_alpha = z_alpha, z_beta = z_beta,
                 verbose = verbose)
    cnt = cnt + 1
    varcounts[[cnt]] = bart$varcount
    avg.vip.mtx[cnt, ] = bart$vip
    median.mi.mtx[cnt, ] = bart$mi
    original.z_lambda[[cnt]] = bart$original.lambda
    original.z_gamma[[cnt]] = bart$original.gamma
    
    #random effect
    for (i in 1:(2 ^ ncol(bart$combined.lambda))) {
      index = as.binary(i - 1) #check the index[i] th variable probability
      if(length(index) != ncol(bart$combined.lambda)) {
        index = rev(c(rep(0, ncol(bart$combined.lambda) - length(index)), index))
      }
      else{
        index = rev(index)
      }
      temp = bart$combined.lambda
      for (j in 1:ncol(bart$combined.lambda)) {
        if(is.null(nrow(temp))) break
        if(as.numeric(index[j]) == 0) {
          temp = temp[which(temp[,j] == 0),]
        }
        else {
          temp = temp[which(temp[,j] != 0),]
        }
      }
      if(is.null(nrow(temp))) {
        avg.lambda.mtx[cnt, i] = 0
      }
      else {
        avg.lambda.mtx[cnt, i] = nrow(temp) / nrow(bart$combined.lambda)
      }
    }
    
    if (length(categorical.idx) > 0) {
      avg.within.type.vip.mtx[cnt, ] = bart$within.type.vip
    }
  }
    
  avg.vip = colMeans(avg.vip.mtx)
  names(avg.vip) = colnames(x.train)
  avg.vip = sort(avg.vip, decreasing = T)
  
  median.mi = apply(median.mi.mtx, 2, median)
  names(median.mi) = colnames(x.train)
  median.mi = sort(median.mi, decreasing = T)
  
  avg.lambda = colMeans(avg.lambda.mtx)
  
  if (length(categorical.idx) > 0) {
    avg.within.type.vip = colMeans(avg.within.type.vip.mtx)
    names(avg.within.type.vip) = colnames(x.train)
    avg.within.type.vip = sort(avg.within.type.vip, decreasing = T)
  }

  if(verbose) cat("complete! \n")
  
  
  #-----------------------------------------------------------
  # build null permutation
  if(verbose) cat("null data sets...")
  
  ## set up permute matrix
  permute.vips = matrix(NA, nrow = npermute, ncol = ncol(x.train))
  permute.mis = matrix(NA, nrow = npermute, ncol = ncol(x.train))
  permute.z_lambda = list()
  permute.z_gamma = list()
  #if there are categorical variables, then should add 1 because of intercept
  if(ncol(rfModelMatrix(z.train)$Z) > ncol(z.train))
    permute.lambda = matrix(NA, nrow = npermute, ncol = 2 ^ (ncol(z.train) + 1))
  else
    permute.lambda = matrix(NA, nrow = npermute, ncol = 2 ^ ncol(z.train))
  
  if (length(categorical.idx) > 0)
    permute.within.type.vips = matrix(NA, nrow = npermute, ncol = ncol(x.train))
  
  cnt = 0
  while (cnt < npermute) {
    print(cnt)
    y.permuted = ave(y.train, id, FUN = function(x) sample(x, replace = FALSE))
    #shuffle y
    # y.permuted = sample(y.permuted, replace = FALSE)
    #shuffle Z
    
    #shuffle X
    # x.permuted <- x.train
    # for (cluster in unique(id)) {
    #   cluster_idx <- which(id == cluster)
    #   perm <- sample(length(cluster_idx))  # generate permutation within cluster
    #   x.permuted[cluster_idx, ] <- x.train[cluster_idx, ][perm, , drop = FALSE]
    # }
    
    #check if z.train is constant within id or not
    z_constant <- all(tapply(seq_len(nrow(z.train)), id, function(rows) {
      all(apply(z.train[rows, , drop = FALSE], 2, function(col) length(unique(col)) == 1))
    }))
    
    if (z_constant) {
      z.permuted <- z.train 
      print("z.train is constant within id")
      z_unique <- z.train[!duplicated(id), , drop = FALSE]
      id_unique <- unique(id)
      
      # Shuffle the rows across IDs
      z_shuffled <- z_unique[sample(nrow(z_unique)), , drop = FALSE]
      
      # Create a permuted version of z.train
      z.permuted <- z.train  # same shape
      
      # Assign shuffled values back by ID
      for (i in seq_along(id_unique)) {
        rows <- which(id == id_unique[i])
        z.permuted[rows, ] <- z_shuffled[rep(i, length(rows)), , drop = FALSE]
      }
      
    } else {
      print("z.train is not constant within id")
      z.permuted <- z.train  # make a copy
      z.permuted[] <- NA     # initialize with NA
      
      for (cluster_id in unique(id)) {
        idx <- which(id == cluster_id)
        z.permuted[idx, ] <- z.train[sample(idx, length(idx), replace = FALSE), ]
      }
    }
    z.permuted <- as.data.frame(z.permuted)
    print("shuffle complete")
    
    bart = wbart(x.train = x.train, y.train = y.permuted, z.train = z.train, id = id,
                 sparse = FALSE,
                 xinfo = xinfo, numcut = numcut, usequants = usequants, cont = cont, rm.const = rm.const,
                 k = k, power = power, base = base, split.prob = split.prob,
                 ntree = ntree, ndpost = ndpost, nskip = nskip, keepevery = keepevery,
                 z_c0 = z_c0, z_d0 = z_d0,
                 z_gamma_mean = z_gamma_mean, z_gamma_cov = z_gamma_cov,
                 z_lambda_mean = z_lambda_mean, z_lambda_cov = z_lambda_cov, z_alpha = z_alpha, z_beta = z_beta,
                 verbose = verbose)
    
    cnt = cnt + 1
    
    varcounts[[nreps + cnt]] = bart$varcount
    permute.vips[cnt, ] = bart$vip
    permute.mis[cnt, ] = bart$mi
    permute.z_lambda[[cnt]] = bart$original.lambda
    permute.z_gamma[[cnt]] = bart$original.gamma
    
    for (i in 1:(2 ^ ncol(bart$combined.lambda))) {
      index = as.binary(i - 1)
      if(length(index) != ncol(bart$combined.lambda)) {
        index = rev(c(rep(0, ncol(bart$combined.lambda) - length(index)), index))
      }
      temp = bart$combined.lambda
      for (j in 1:ncol(bart$combined.lambda)) {
        if(is.null(nrow(temp))) break
        if(as.numeric(index[j]) == 0) {
          temp = temp[which(temp[,j] == 0),]
        }
        else {
          temp = temp[which(temp[,j] != 0),]
        }
      }
      if(is.null(nrow(temp))) {
        permute.lambda[cnt, i] = 0
      }
      else {
        permute.lambda[cnt, i] = nrow(temp) / nrow(bart$combined.lambda)
      }
    }
    
    if (length(categorical.idx) > 0) {
      permute.within.type.vips[cnt, ] = bart$within.type.vip
    }
  }
  
  if(verbose) cat("complete! \n")
  
  res$original.z_lambda = original.z_lambda
  res$original.z_gamma = original.z_gamma
  res$permute.z_lambda = permute.z_lambda
  res$permute.z_gamma = permute.z_gamma
  
  
  #-----------------------------------------------------------
  # sort permute mat and return results
  colnames(permute.vips) = colnames(x.train)
  permute.vips = permute.vips[, names(avg.vip)]
  
  colnames(permute.mis) = colnames(x.train)
  permute.mis = permute.mis[, names(median.mi)]
  
  if (length(categorical.idx) > 0) {
    colnames(permute.within.type.vips) = colnames(x.train)
    permute.within.type.vips = permute.within.type.vips[, names(avg.within.type.vip)]
  }
  
  #-----------------------------------------------------------
  # use local cutoff & returns
  vip.pointwise.cutoffs = apply(permute.vips, 2, quantile, probs = 1 - alpha)
  vip.imp.names = names(avg.vip[(avg.vip > vip.pointwise.cutoffs) & (avg.vip > 0)])
  vip.imp.cols = sapply(1:length(vip.imp.names), function(x) {which(vip.imp.names[x] == colnames(x.train))})
  res$vip.imp.cols = vip.imp.cols
  res$vip.imp.names = vip.imp.names
  res$avg.vip = avg.vip
  res$avg.vip.mtx = avg.vip.mtx
  res$permute.vips = permute.vips
  
  mi.pointwise.cutoffs = apply(permute.mis, 2, quantile, probs = 1 - alpha)
  mi.imp.names = names(median.mi[(median.mi > mi.pointwise.cutoffs) & (median.mi > 0)])
  mi.imp.cols = sapply(1:length(mi.imp.names), function(x) {which(mi.imp.names[x] == colnames(x.train))})
  res$mi.imp.cols = mi.imp.cols
  res$mi.imp.names = mi.imp.names
  res$median.mi = median.mi
  res$median.mi.mtx = median.mi.mtx
  res$permute.mis = permute.mis
  
  lambda.pointwise.cutoffs = apply(permute.lambda, 2, quantile, probs = 1 - alpha)
  res$avg.lambda = avg.lambda
  res$avg.lambda.mtx = avg.lambda.mtx
  res$permute.lambda = permute.lambda
  res$lambda.pointwise.cutoffs = lambda.pointwise.cutoffs
  
  if (length(categorical.idx) > 0) {
    within.type.vip.pointwise.cutoffs = apply(permute.within.type.vips, 2, quantile, probs = 1 - alpha)
    within.type.vip.imp.names = names(avg.within.type.vip[(avg.within.type.vip > within.type.vip.pointwise.cutoffs) & (avg.within.type.vip > 0)])
    within.type.vip.imp.cols = sapply(1:length(within.type.vip.imp.names), function(x) {which(within.type.vip.imp.names[x] == colnames(x.train))})
    
    res$within.type.vip.imp.cols = within.type.vip.imp.cols
    res$within.type.vip.imp.names = within.type.vip.imp.names
    res$avg.within.type.vip = avg.within.type.vip
    res$avg.within.type.vip.mtx = avg.within.type.vip.mtx
    res$permute.within.type.vips = permute.within.type.vips
  }
  
  res$varcounts = varcounts
  
  #-----------------------------------------------------------
  # score results
  if(length(true.idx) > 0) {
    
    true.len = length(true.idx)
    res$true.idx = true.idx
    
    ## vip
    tp = length(which(vip.imp.cols %in% true.idx))
    positive.len = length(vip.imp.cols)
    
    res$vip.precision = (tp * 1.0) / (positive.len * 1.0)
    res$vip.recall = (tp * 1.0) / (true.len * 1.0)
    res$vip.f1 = 2 * res$vip.precision * res$vip.recall / (res$vip.precision + res$vip.recall)
    
    ## mi
    tp = length(which(mi.imp.cols %in% true.idx))
    positive.len = length(mi.imp.cols)
    
    res$mi.precision = (tp * 1.0) / (positive.len * 1.0)
    res$mi.recall = (tp * 1.0) / (true.len * 1.0)
    res$mi.f1 = 2 * res$mi.precision * res$mi.recall / (res$mi.precision + res$mi.recall)
    
    if (length(categorical.idx) > 0) {
      ## within-type vip
      tp = length(which(within.type.vip.imp.cols %in% true.idx))
      positive.len = length(within.type.vip.imp.cols)
      
      res$wt.vip.precision = (tp * 1.0) / (positive.len * 1.0)
      res$wt.vip.recall = (tp * 1.0) / (true.len * 1.0)
      res$wt.vip.f1 = 2 * res$wt.vip.precision * res$wt.vip.recall / (res$wt.vip.precision + res$wt.vip.recall)
    }
  }
  
  #-----------------------------------------------------------
  if (plot) {
    
    if ((n.var.plot == Inf) | (n.var.plot > ncol(x.train))){
      n.var.plot = ncol(x.train)
    }
    
    ## vip
    non.zero.idx = which(avg.vip > 0)[1:min(n.var.plot, length(which(avg.vip > 0)))]
    plot.n = length(non.zero.idx)
    if(length(non.zero.idx) < length(avg.vip)) 
      warning(paste(length(which(avg.vip == 0)), "predictors with inclusion proportions of 0 omitted from plots."))
    max.cut = max(apply(permute.vips, 2, quantile, probs = 1 - alpha, na.rm = TRUE))
    
    plot(1:plot.n, avg.vip[non.zero.idx], 
         type = "n", xlab = "Predictors", xaxt = "n", ylim = c(0, max(max(avg.vip), max.cut * 1.1)),
         main = "Permutation-Based Variable Selection", ylab = "BART VIP")
    axis(1, at = 1:plot.n, labels = names(avg.vip[non.zero.idx]), las = 2)
    for (j in non.zero.idx){
      points(j, avg.vip[j], 
             pch = ifelse(avg.vip[j] <= quantile(permute.vips[, j], 1 - alpha), 1, 16),
             col = ifelse(names(avg.vip[j]) %in% categorical.names, 'green', 'red'))
    }
    sapply(non.zero.idx, function(s) {segments(s, 0, x1 = s, quantile(permute.vips[, s], 1 - alpha), col = "grey")})
    
    ## mi
    non.zero.idx = which(median.mi > 0)[1:min(n.var.plot, length(which(median.mi > 0)))]
    plot.n = length(non.zero.idx)
    if(length(non.zero.idx) < length(median.mi)) 
      warning(paste(length(which(median.mi == 0)), "predictors with Metropolis importance of 0 omitted from plots."))
    max.cut = max(apply(permute.mis, 2, quantile, probs = 1 - alpha, na.rm = TRUE))
    
    plot(1:plot.n, median.mi[non.zero.idx], 
         type = "n", xlab = "Predictors", xaxt = "n", ylim = c(0, max(max(median.mi), max.cut * 1.1)),
         main = "Permutation-Based Variable Selection", ylab = "BART MI")
    axis(1, at = 1:plot.n, labels = names(median.mi[non.zero.idx]), las = 2)
    for (j in non.zero.idx){
      points(j, median.mi[j], 
             pch = ifelse(median.mi[j] <= quantile(permute.mis[, j], 1 - alpha), 1, 16),
             col = ifelse(names(median.mi[j]) %in% categorical.names, 'green', 'red'))
    }
    sapply(non.zero.idx, function(s) {segments(s, 0, x1 = s, quantile(permute.mis[, s], 1 - alpha), col = "grey")})
    
    if (length(categorical.idx) > 0) {
      ## within-type vip
      non.zero.idx = which(avg.within.type.vip > 0)[1:min(n.var.plot, length(which(avg.within.type.vip > 0)))]
      plot.n = length(non.zero.idx)
      if(length(non.zero.idx) < length(avg.within.type.vip)) 
        warning(paste(length(which(avg.within.type.vip == 0)), 
                      "predictors with inclusion proportions (within-type) of 0 omitted from plots."))
      max.cut = max(apply(permute.within.type.vips, 2, quantile, probs = 1 - alpha, na.rm = TRUE))
      
      plot(1:plot.n, avg.within.type.vip[non.zero.idx], 
           type = "n", xlab = "Predictors", xaxt = "n", ylim = c(0, max(max(avg.within.type.vip), max.cut * 1.1)),
           main = "Permutation-Based Variable Selection", ylab = "BART Within-Type VIP")
      axis(1, at = 1:plot.n, labels = names(avg.within.type.vip[non.zero.idx]), las = 2)
      for (j in non.zero.idx){
        points(j, avg.within.type.vip[j], 
               pch = ifelse(avg.within.type.vip[j] <= quantile(permute.within.type.vips[, j], 1 - alpha), 1, 16),
               col = ifelse(names(avg.within.type.vip[j]) %in% categorical.names, 'green', 'red'))
      }
      sapply(non.zero.idx, function(s) {segments(s, 0, x1 = s, quantile(permute.within.type.vips[, s], 1 - alpha), col = "grey")})
    }
  }
  
  #------------------------------
  # timer ends
  end = Sys.time()
  if(verbose) cat("Elapsed", end-start, '\n')
  
  return(res)
}
               
