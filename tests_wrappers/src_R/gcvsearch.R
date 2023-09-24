#' file: utilities.R
#' author: Cristian Castiglione (with small changes applied by M. De Sanctis and I. Di Battista)
#' last change: 18/02/2023
#' content: 
#' - function: qgam.fem
#' - function: qgam.xfem


### GCV optimization: non-parametric model
qgam.fem_our <- function(
    y,
    Z,
    R0,
    R1,
    alpha = 0.5,
    lambdas = 10^seq(-2, 3, by = 0.2),
    eps = 0.001,
    tune = 1.0,
    tol = 1e-6, 
    ftol = 1e-04,
    xtol = 1e-04,
    exact = FALSE,
    tol.weights = 1e-6, 
    mass_lumping = TRUE, 
    symmetry = TRUE, 
    sparsity = TRUE, 
    block = FALSE,
    nsample = 100,
    maxiter = 100,
    verbose = FALSE,
    show = TRUE,
    xlab = "log-lambda",
    ylab = "GCV score",
    main = "Smoothing parameter search",
    path = "",    
    measure.time = FALSE, 
    ...
) {
  
  "
  Nonparametric PDE-spatial quantile regression: Smoothing parameter search
  -------------------------------------------------------------------------
  
  y : n x 1 vector of response variables;
  Z : n x d matrix of FEM basis;
  R0 : d x d mass matrix;
  R1 : d x d stiffness matrix;
  alpha : quantile level lying in (0,1);
  lambdas : sequence of positive smoothing parameters;
  tune : positive tuning parameter rescaling the d.f. in the GCV calculation;
  ftol : tollerance parameter for the penalized likelihood;
  xtol : tollerance parameter for the regression parameters;
  sparse : true for using a sparse linear algebra routine;
  gcvout : true for returning the GCV score;
  maxiter : maximum number of iteration of the FPIRLS algorithm;
  verbose : true for printing some intermediate results during the execution;
  "
  
  cat("GCV optimization: start ")
  
  if(length(lambdas) == 1)
    save = FALSE
  else
    save = TRUE
  
  # estimation function
  fit <- function(x) {
    cat(".")
    qgam.fem.fit_our(
      y, Z, R0, R1, alpha = alpha, lambda = x,
      eps = eps, tune = tune, tol = tol, ftol = ftol, xtol = xtol,
      gcvout = TRUE, exact = exact, 
      tol.weights = tol.weights, 
      mass_lumping = mass_lumping, 
      symmetry = symmetry, 
      sparsity = sparsity, 
      block = block,
      nsample = nsample,
      maxiter = maxiter, verbose = verbose)
  }
  # GCV scores calculation
  if(measure.time){
    start_time <- Sys.time()
  }
  fits <- lapply(lambdas, fit)
  if(measure.time){
    end_time <- Sys.time()
    duration = difftime(end_time, start_time, units = "secs")
  }
  gcvs <- as.vector(unlist(lapply(fits, function(x) x$gcv)))
  edfs <- as.vector(unlist(lapply(fits, function(x) x$edf)))
  degrees.gcv <- as.vector(unlist(lapply(fits, function(x) x$degree.gcv)))
  num_gcvs <- as.vector(unlist(lapply(fits, function(x) x$num_gcv)))
  idx  <- which.min(gcvs)
  
  cat(" end \n")
  
  # plotting
  if (show) {
    plot(x = log(lambdas, base = 10), y = gcvs, type = "l",
         xlab = xlab, ylab = ylab, main = main, ...)
    points(x = log(lambdas, base = 10), y = gcvs, pch = 20)
    points(x = log(lambdas, base = 10)[idx], y = gcvs[idx], col = 2, pch = 20)
    abline(v = log(lambdas, base = 10)[idx], lty = 2, col = 2)
  }
  
  # warning
  if (idx == 1 || idx == length(lambdas)) {
    cat("Warning: minimum GCV at the limits. \n")
  }
  
  outfit <- fits[[idx]]
  outfit$lambdas <- lambdas
  outfit$gcvs <- gcvs
  outfit$edfs <- edfs
  outfit$degrees.gcv <- degrees.gcv 
  outfit$num_gcvs <- num_gcvs
  outfit$idx <- idx
  if(measure.time){
    outfit$time_fit <- duration
  }
   
  # Save results
  alpha_string = as.character(alpha*100)
  if(save){
    write.csv(format(lambdas, digits=16), 
                paste0(path, "/lambdasR.csv"))
    write.csv(format(gcvs, digits=16), 
                paste0(path, "/GCV_scoresR.csv"))
    write.csv(format(degrees.gcv, digits=16), 
                paste0(path, "/edfR.csv"))
    write.csv(format(num_gcvs, digits=16), 
                paste0(path, "/numR.csv"))

  }

  # output
  return(outfit)
  
}


### GCV optimization: semi-parametric model
qgam.xfem_our <- function(
    y,
    X,
    Z,
    R0,
    R1,
    alpha = 0.5,
    lambdas = 10^seq(-2, 3, by = 0.2),
    eps = 0.0,
    tune = 1.0,
    tol = 1e-6, 
    ftol = 1e-04,
    xtol = 1e-04,
    gcvmtd = "diag",
    exact = FALSE,
    tol.weights = 1e-6, 
    mass_lumping = TRUE,
    symmetry = TRUE, 
    sparsity = TRUE, 
    nsample = 100,
    maxiter = 100,
    verbose = FALSE,
    show = TRUE,
    xlab = "log-lambda",
    ylab = "GCV score",
    main = "Smoothing parameter search",
    path = "",  
    measure.time = FALSE, 
    ...
){
  "
  Semiparametric PDE-spatial quantile regression: Smoothing parameter search
  --------------------------------------------------------------------------
  
  y : n x 1 vector of response variables;
  Z : n x d matrix of FEM basis;
  R0 : d x d mass matrix;
  R1 : d x d stiffness matrix;
  alpha : quantile level lying in (0,1);
  lambdas : sequence of positive smoothing parameters;
  tune : positive tuning parameter rescaling the d.f. in the GCV calculation;
  tol: 
  ftol : tollerance parameter for the penalized likelihood;
  xtol : tollerance parameter for the regression parameters;
  sparse : true for using a sparse linear algebra routine;
  gcvout : true for returning the GCV score;
  maxiter : maximum number of iteration of the FPIRLS algorithm;
  verbose : true for printing some intermediate results during the execution;
  "
  
  cat("GCV optimization: start ")
  
  if((length(lambdas) == 1) || (path == ""))
    save = FALSE
  else
    save = TRUE
  
  
  # estimation function

  fit <- function(x) {
    cat(".")
    qgam.xfem.fit_our(
      y, X, Z, R0, R1, alpha = alpha, lambda = x,
      eps = eps, tune = tune, tol = tol, ftol = ftol, xtol = xtol,
      gcvout = TRUE, gcvmtd = gcvmtd, exact = exact, 
      tol.weights = tol.weights, 
      mass_lumping = mass_lumping, 
      symmetry = symmetry, 
      sparsity = sparsity, 
      nsample = nsample, maxiter = maxiter, verbose = verbose)
  }

  
  
  # GCV scores calculation
  if(measure.time){
    start_time <- Sys.time()
  }
  fits <- lapply(lambdas, fit)
  if(measure.time){
    end_time <- Sys.time()
    duration = difftime(end_time, start_time, units = "secs")
  } 
  gcvs <- as.vector(unlist(lapply(fits, function(x) x$gcv)))
  edfs <- as.vector(unlist(lapply(fits, function(x) x$edf)))
  degrees.gcv <- as.vector(unlist(lapply(fits, function(x) x$degree.gcv)))
  num_gcvs <- as.vector(unlist(lapply(fits, function(x) x$num_gcv)))
  idx  <- which.min(gcvs)
  
  cat(" end \n")
  
  # plotting
  if (show) {
    plot(x = log(lambdas, base = 10), y = gcvs, type = "l",
         xlab = xlab, ylab = ylab, main = main, ...)
    points(x = log(lambdas, base = 10), y = gcvs, pch = 20)
    points(x = log(lambdas, base = 10)[idx], y = gcvs[idx], col = 2, pch = 20)
    abline(v = log(lambdas, base = 10)[idx], lty = 2, col = 2)
  }
  
  # warning
  if (idx == 1 || idx == length(lambdas)) {
    cat("Warning: minimum GCV at the limits. \n")
  }
  
  outfit <- fits[[idx]]
  outfit$lambdas <- lambdas
  outfit$gcvs <- gcvs
  outfit$edfs <- edfs
  outfit$degrees.gcv <- degrees.gcv 
  outfit$num_gcvs <- num_gcvs
  outfit$idx <- idx
  if(measure.time){
    outfit$time_fit <- duration
  }
   
  
  # Save results
  alpha_string = as.character(alpha*100)
  if(eps > 1e-15){
    eps_string = "_eps"
  }else{
    eps_string = ""
  }
  
  if(save){
    write.csv(format(lambdas, digits=16), 
              paste0(path, "/lambdasR", eps_string, ".csv"))
    write.csv(format(gcvs, digits=16), 
              paste0(path, "/GCV_scoresR", eps_string, ".csv"))
    write.csv(format(degrees.gcv, digits=16), 
              paste0(path, "/edfR", eps_string, ".csv"))
    write.csv(format(num_gcvs, digits=16), 
              paste0(path, "/numR", eps_string, ".csv"))
  }
    

  
  # output
  return(outfit)
}

































