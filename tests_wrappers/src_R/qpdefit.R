#' file: qpdefit.R
#' author: Cristian Castiglione (with small changes applied by M. De Sanctis and I. Di Battista)
#' last change: 18/02/2023
#' content:
#' - function: qgam.fem-fit
#' - function: qgam.xfem-fit


### SQR-PDE (no covariates)
qgam.fem.fit_our <- function(
    y,
    Z,
    R0,
    R1,
    alpha = 0.5,
    lambda = 1.0,
    eps = 0.0,
    tune = 1.0,
    tol = 1e-6, 
    ftol = 1e-4, 
    xtol = 1e-4, 
    gcvout = TRUE,
    exact = TRUE,
    tol.weights = 1e-6, 
    mass_lumping = TRUE,   
    symmetry = TRUE, 
    sparsity = TRUE, 
    block = FALSE, 
    nsample = 1000,
    maxiter = 200,
    verbose = TRUE
) {
  
  "
  Nonparametric PDE-spatial quantile regression
  ---------------------------------------------
   
  y : n x 1 vector of response variables;
  Z : n x d matrix of FEM basis;
  R0 : d x d mass matrix;
  R1 : d x d stiffness matrix;
  alpha : quantile level lying in (0,1);
  lambda : positive smoothing parameter;
  eps : positive parameter controlling the smoothing factor in the GCV calculation;
  tune : positive tuning parameter rescaling the d.f. in the GCV calculation;
  tol: tolerance for the stopping criterion based on the functional J 
  ftol : tollerance parameter for the penalized likelihood;
  xtol : tollerance parameter for the regression parameters;
  gcvout : true for returning the GCV score;
  exact : true for an exact calculation of the GCV, false for an approximation;
  mass_lumping: TRUE to active mass lumping for the computation of the penalty matrix
  nsample : number of random vectors used for the stochastic approximation of the GCV
  maxiter : maximum number of iteration of the FPIRLS algorithm;
  verbose : true for printing some intermediate results during the execution;
  "
  
  # input dimensions
  n <- nrow(Z)
  d <- ncol(Z)
  
  # parameter check
  param.check(
    nobs = n, alpha = alpha,
    lambda = lambda, eps = eps, tune = tune,
    ftol = ftol, xtol = xtol, maxiter = maxiter,
    exact = exact, gcvout = gcvout, verbose = verbose
  )
  
  # transform to sparse matrix format
  Z  <- as(Z,  "sparseMatrix")
  R0 <- as(R0, "sparseMatrix")
  R1 <- as(R1, "sparseMatrix")
  
  # penalization matrix
  P = get.penalty.matrix(R0, R1, mass_lumping)
  R <- (2.0*n*lambda) * P     
  

  # initialization
  ZWZ <- crossprod(Z, Z)   
  ZWy <- crossprod(Z, y)
  
  if(!block) {   # direct implementation
    if(symmetry)
      A = as(ZWZ + R, "symmetricMatrix")  
    else 
      A <-  ZWZ + R 
    
    b <- as.vector(ZWy)
    f <- as.vector(solve(A, b, sparse = sparsity))  
    g <- as.vector(solve(R0, R1%*%f, sparse = sparsity))  
    
  }
  else {      # block implementation (here, symmetry not active)
    A = rbind(cbind(-ZWZ/(2.0*n) , lambda*t(R1)), cbind(lambda*R1, lambda*R0))
    # force the vectors to be columns to correctly implement the rhs 
    b1 = as.matrix(as.vector(-ZWy/(2*n)), nrow = d, ncol = 1)
    b2 = as.matrix(rep(0,d), nrow = d, ncol = 1) 
    b <- rbind(b1, b2)
    
    sol <- as.vector(solve(A, b, sparse = sparsity))  
    f <- sol[1:d]
    g <- sol[(d+1):(2*d)]
              
  }

  fn <- as.vector(Z %*% f)
  mu <- as.vector(fn)
  
  loss.J    <- as.numeric(sum( ( (1 / sqrt((2*n))) * (y-mu) )^2 ) ) 
  penalty <- lambda * as.numeric(crossprod(g, R0 %*% g))      
  
  sigma <- (as.numeric(sum(loss.fun(y, mu, alpha))) + 
              0.5 * lambda * as.numeric(crossprod(f, P %*% f))) / n
  
  niter <- maxiter
  check <- FALSE
  
  # memory allocation for storing the history of the estimation process
  est_mu    <- matrix(0, nrow = n, ncol = maxiter)
  est_fn    <- matrix(0, nrow = n, ncol = maxiter)
  est_obj   <- rep(0, maxiter)
  est_sigma <- rep(0, maxiter)
  
  # printing
  if (verbose) print.output(0, 1.0, 1.0, 1.0)
  
  # optimization loop
  for (iter in 1:maxiter) {
    
    # residuals, weights and pseudo-observations
    r <- get.abs.resid(y, mu)                   
    w <- get.weight.vector(r, tol.weights)      
    u <- get.pseudo.data(y, r, alpha)           

    # linear system construction
    ZWZ <- crossprod(Z, Z / w)             
    ZWy <- crossprod(Z, u / w)        
    
    # sparse linear system
    if(!block) {
      if(symmetry)
        A = forceSymmetric(ZWZ + R)  
      else 
        A <-  ZWZ + R 
      
      b <- as.vector(ZWy)
      f <- as.vector(solve(A, b, sparse = sparsity))  
      g <- as.vector(solve(R0, R1%*%f, sparse = sparsity))  
      
    }
    else {
      A = rbind(cbind(-ZWZ/(2.0*n) , lambda*t(R1)), cbind(lambda*R1, lambda*R0))
      
      # force the vectors to be columns to correctly implement the rhs 
      b1 = as.matrix(as.vector(-ZWy/(2*n)), nrow = d, ncol = 1)
      b2 = as.matrix(rep(0,d), nrow = d, ncol = 1) 
      b <- rbind(b1, b2)
      
      sol <- as.vector(solve(A, b, sparse = sparsity))  
      f <- sol[1:d]
      g <- sol[(d+1):(2*d)]
      
    }
    
    # fn and beta estimates
    fn <- as.vector(Z %*% f)
    mu <- as.vector(fn)
    
    loss.J    <- as.numeric(sum( ( (1 / sqrt((2*n)*w)) * (u-mu) )^2 ) )    
    penalty <- lambda * as.numeric(crossprod(g, R0 %*% g))       
    
    # storing
    est_mu[,iter] <- mu
    est_fn[,iter] <- fn
    est_obj[iter] <- loss.J + penalty   
    
    # evaluate if the convergence has been reached
    if (iter == 1) {
      delta  <- 1.0
    }
    if (iter > 1) {
      delta <- abs(est_obj[iter] - est_obj[iter - 1])
      if (delta < tol) check <- TRUE
    }
    
    # print the intermediate results
    if (verbose) print.output(iter, delta, est_obj[iter])
    
    # exit from the loop when convergence is reached
    if (check) {
      
      # total number of iterations
      niter <- iter 

      # exclude the unnecessary elements of the history matrices
      est_mu    <-  est_mu[, 1:(iter - 1)]
      est_fn    <-  est_fn[, 1:(iter - 1)]
      est_obj   <-  est_obj[1:(iter - 1)]

      break
    }
  }
  
  # print the end of loop
  if(verbose) print.output(-1, 1.0, 1.0, 1.0)
  
  # GCV calculation
  loss <- -1
  edf  <- -1
  gcv  <- -1
  reml <- -1
  
  res <- y - mu
  sigma <- mean(loss.fun(y, mu, alpha))
  
  if (gcvout) {
    tr   <- get.prod.trace(ZWZ + R, ZWZ, exact, nsample, verbose)  
    loss <- sum(loss.fun(y, mu, alpha, eps))                      
    edf  <- n - tune * tr                                        
    degree.gcv = tr 
    gcv  <- ifelse(edf > 1e-08, loss / edf, Inf)    
  }
  
  # output
  return(list(
    res = as.vector(res),
    mu = as.vector(mu),
    fn = as.vector(fn),
    f = as.vector(f),
    w = as.vector(w),
    P = P, 
    A = A, 
    b = b, 
    u = u, 
    sigma = sigma,
    alpha = alpha,
    lambda = lambda,
    eps = eps,
    penalty = penalty,
    trace.mu = est_mu,
    trace.fn = est_fn,
    est_obj = est_obj,
    edf = edf,
    tr = tr, 
    num_gcv = loss, 
    gcv = gcv,
    degree.gcv = degree.gcv, 
    niter = niter,
    converged = check
  ))
}



### SQR-PDE (with covariates)
qgam.xfem.fit_our <- function(
    y,
    X,
    Z,
    R0,
    R1,
    alpha = 0.5,
    lambda = 1.0,
    eps = 0.0,
    tune = 1.0,
    tol = 1e-6, 
    ftol = 1e-04,
    xtol = 1e-04,
    gcvout = TRUE,
    gcvmtd = "diag", # "diag", "dense"
    exact = TRUE,
    tol.weights = 1e-6, 
    mass_lumping = TRUE, 
    symmetry = TRUE, 
    sparsity = TRUE, 
    nsample = 1000,
    maxiter = 100,
    verbose = TRUE
) {
  
  "
  Semiparametric PDE-spatial quantile regression
  ----------------------------------------------
   
  y : n x 1 vector of response variables;
  X : n x q matrix of covariates;
  Z : n x d matrix of FEM basis;
  R0 : d x d mass matrix;
  R1 : d x d stiffness matrix;
  alpha : quantile level lying in (0,1);
  lambda : positive smoothing parameter;
  eps : positive parameter controlling the smoothing factor in the GCV calculation;
  tune : positive tuning parameter rescaling the d.f. in the GCV calculation;
  tol: FPIRLS tolerance for stopping criterion based on functional J 
  ftol : tollerance parameter for the penalized likelihood;
  xtol : tollerance parameter for the regression parameters;
  gcvout : true for returning the GCV score;
  gcvmtd : 'diag' for a sparse approximate calculation of the GCV, 'dense' for an exact calculation;
  exact : true for an exact calculation of the GCV, false for a stochastic approximation;
  mass_lumping: TRUE to active mass lumping for the computation of the penalty matrix 
  maxiter : maximum number of iteration of the FPIRLS algorithm;
  verbose : true for printing some intermediate results during the execution;
  "
  
  # input dimensions
  n <- nrow(Z) # sample size
  q <- ncol(X) # number of covariates
  d <- ncol(Z) # number of FEM knots
  
  # parameter check
  param.check(
    nobs = n, alpha = alpha,
    lambda = lambda, eps = eps, tune = tune,
    ftol = ftol, xtol = xtol, maxiter = maxiter,
    exact = exact, gcvout = gcvout, verbose = verbose
  )
  
  # transform to sparse matrix format
  Z  <- as(Z, "sparseMatrix")
  R0 <- as(R0, "sparseMatrix")
  R1 <- as(R1, "sparseMatrix")
  
  # penalization matrix (via mass lumping)
  P <- get.penalty.matrix(R0, R1, mass_lumping) 
  R <- (2.0*n*lambda) * P  
  
  # parameter indices
  idx_x <- 1:q
  idx_f <- (q + 1):(d + q)
  
  # completed design matrix
  ZWZ <- crossprod(Z, Z) + R
  XWX <- crossprod(X, X)
  ZWX <- crossprod(Z, X)
  XWZ <- t(ZWX)
  XWy <- crossprod(X, y)
  ZWy <- crossprod(Z, y)
  
  # initialization
  if(symmetry){
    A <- as(ZWZ, "symmetricMatrix")
  }
  else{
    A <- ZWZ
  }
  b <- as.vector(ZWy)
  f <- as.vector(solve(A, b, sparse = sparsity)) 
  g <- as.vector(solve(R0, R1%*%f, sparse = sparsity))
  

  beta <- rep(0, q)
  bn <- as.vector(X %*% beta) 
  fn <- as.vector(Z %*% f)
  mu <- as.vector(bn + fn) 
  
  sigma <- 1
  loss.J    <- as.numeric(sum( (y-mu)^2 )/(2*n))
  penalty.J <- lambda * as.numeric(crossprod(g, R0 %*% g))
  
  # initialization of the convergence parameters
  niter <- maxiter
  check <- FALSE
  
  # memory allocation for storing the history of the estimation process
  est_mu    <- matrix(0, nrow = n, ncol = maxiter)
  est_bn    <- matrix(0, nrow = n, ncol = maxiter)
  est_fn    <- matrix(0, nrow = n, ncol = maxiter)
  est_sigma <- rep(0, maxiter)
  est_obj   <- rep(0, maxiter)
  
  # printing
  if (verbose) print.output(0, 1.0, 1.0, 1.0)
  
  # optimization loop
  for(iter in 1:maxiter){
    
    cat("\n Iteration number:", iter)
    
    # residuals, weights and pseudo-observations
    r <- get.abs.resid(y, mu)
    w <- get.weight.vector(r, tol.weights)
    u <- get.pseudo.data(y, r, alpha)    
    
    # linear system construction
    ZWZ <- crossprod(Z, Z / w) + R
    XWX <- crossprod(X, X / w)
    ZWX <- crossprod(Z, X / w)
    XWZ <- t(ZWX)
    XWy <- crossprod(X, u / w)
    ZWy <- crossprod(Z, u / w)
    
    # regression parameter update
    v <- u - as.vector(Z %*% f)     
    if (q == 1) {
      # n x 1 covariate matrix
      beta <- sum(X * v / w) / as.numeric(XWX)
    } else {
      # n x q covariate matrix
      beta <- solve(XWX, crossprod(X, v / w)) 
    }
    
    # sparse linear system
    v <- u - as.vector(X %*% beta)            
    if(symmetry)
      A = forceSymmetric(ZWZ)   
    else 
      A <- ZWZ
    
    b <- crossprod(Z, v / w)
    f <- as.vector(solve(A, b, sparse = sparsity)) 
    g <- as.vector(solve(R0, R1%*%f, sparse = sparsity))  

    # fn and beta estimates
    bn <- as.vector(X %*% beta)                 
    fn <- as.vector(Z %*% f)
    mu <- as.vector(bn + fn)    
    
    loss.J    <- as.numeric(sum( (u-mu)^2/w )/(2*n))
    penalty.J <- lambda * as.numeric(crossprod(g, R0 %*% g))       

    # storing
    est_mu[, iter]  <- mu
    est_bn[, iter]  <- bn
    est_fn[, iter]  <- fn
    est_sigma[iter] <- sigma
    est_obj[iter] <- loss.J + penalty.J  

    # evaluate if the convergence has been reached
    if (iter == 1) {
      delta <- 1.0
    }
    if (iter > 1) {
      delta <- abs(est_obj[iter] - est_obj[iter - 1])
      if (delta < tol) check <- TRUE
    }
    
    # print the intermediate results
    if (verbose) print.output(iter, delta_par, delta_llik, pen_log_lik / n)
    
    # exit from the loop when convergence is reached
    if (check) {
      
      # total number of iterations
      niter <- iter
       
      # exclude the unnecessary elements of the history matrices
      est_mu    <- est_mu[, 1:(iter - 1)]
      est_bn    <- est_bn[, 1:(iter - 1)]
      est_fn    <- est_fn[, 1:(iter - 1)]
      est_sigma <- est_sigma[1:(iter - 1)]
      est_obj   <- est_obj[1:(iter - 1)]
      
      break
    }
  }
  
  # print the end of loop
  if (verbose) print.output(-1, 1.0, 1.0, 1.0)
  
  # GCV calculation
  loss <- -1
  edf  <- -1
  gcv  <- -1
  
  res <- y - mu 
  sigma <- mean(loss.fun(y, mu, alpha))  
  
  if (gcvout) {
    
    if (gcvmtd == "dense") {
      WX <- X / w
      Q <- diag(1 / w) - WX %*% solve(XWX, t(WX))
      ZQZ <- crossprod(Z, Q %*% Z)
      A <- forceSymmetric(ZQZ + 2.0*n*lambda*P)   
      B <- ZQZ
    } else if(gcvmtd == "diag"){
      ZWZ <- crossprod(Z, Z / w)
      A <- ZWZ + 2.0*n*lambda*P
      B <- ZWZ
    } else {
      stop("Only gcvmtd = 'diag' and gcvmtd = 'dense' methods are allowed for the GCV calculation.")
    }
    
    tr   <- get.prod.trace(A, B, exact, nsample, verbose)
    loss <- sum(loss.fun(y, mu, alpha, eps))
    edf  <- n - q - tune * tr
    degree.gcv = tr + q 
    gcv  <- ifelse(edf > 1e-08, loss / edf, Inf)
  }
  
  cat("\n J at last iteration: ", loss.J + penalty.J, "\n")
  
  # output
  return(list(
    res = as.vector(res),
    mu = as.vector(mu),
    bn = as.vector(bn),      # parametric part of the solution 
    fn = as.vector(fn),      # non-parametric part 
    beta = as.vector(beta),
    f = as.vector(f),
    w = as.vector(w),
    P = P, 
    A = A, 
    J_final = est_obj, 
    sigma = sigma,
    alpha = alpha,
    lambda = lambda,
    eps = eps,
    trace.mu = est_mu,
    trace.fn = est_fn,
    trace.sigma = est_sigma,
    est_obj = est_obj,
    edf = edf,
    tr = tr, 
    num_gcv = loss, 
    gcv = gcv,
    degree.gcv = degree.gcv, 
    niter = niter,
    converged = check
  ))
}










