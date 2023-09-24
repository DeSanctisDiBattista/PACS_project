#' file: utilities.R
#' author: Cristian Castiglione
#' last change: 18/02/2023
#' content: 
#' - function: param.check
#' - function: get.penalty.matrix
#' - function: get.abs.resid
#' - function: get.weight.vector
#' - function: get.pseudo.data
#' - function: get.prod.trace
#' - function: loss.fun
#' - function: llik.fun
#' - function: absmax
#' - function: print.output


### parameter checks
param.check <- function(
  nobs, alpha, lambda, eps, tune,
  ftol, xtol, maxiter, exact, gcvout, verbose) {
  
  if ((alpha <= 1 / nobs) || (alpha >= 1 - 1 / nobs)) {
    stop("Invalid parameter: 'alpha'. \n")
  }
  if (lambda <= 0) {
    stop("Invalid parameter: 'lambda'. \n")
  }
  # if (eps <= 0) {
  #   stop("Invalid parameter: 'eps'. \n")
  # }
  if (tune <= 0) {
    stop("Invalid parameter: 'tune'. \n")
  }
  if (ftol <= 0) {
    stop("Invalid parameter: 'ftol'. \n")
  }
  if (xtol <= 0) {
    stop("Invalid parameter: 'xtol'. \n")
  }
  if (maxiter <= 0) {
    stop("Invalid parameter: 'maxiter'. \n")
  }
  if (!is.logical(gcvout)) {
    stop("Invalid parameter: 'gcvout'. \n")
  }
  if (!is.logical(verbose)) {
    stop("Invalid parameter: 'verbose'. \n")
  }
  if (!is.logical(exact)) {
    stop("Invalid parameter: 'exact'. \n")
  }
}

### penalization matrix
get.penalty.matrix <- function(R0, R1, lumping = TRUE) {
  if (!lumping) {
    P <- crossprod(R1, solve(R0, R1, sparse = TRUE))
    P <- forceSymmetric(P)
  } else {
    P <- crossprod(R1, R1 / colSums(R0))
    P <- forceSymmetric(P)
  }
  return(P)
}

### calculate the absolute residual vector
get.abs.resid <- function(y, mu) {
  return(as.vector(abs(y - mu)))
}

### calculate the weighting vector
get.weight.vector <- function(r, eps = 1e-06) {
  return(as.vector(ifelse(r > eps, r, r + eps)))  # I : se r troppo piccolo aumento un po'
                                                  #     per evitare wi = 0
}

### calculate the pseudo-data vector
get.pseudo.data <- function(y, w, alpha) {
  return(as.vector(y + (2 * alpha - 1) * w))      # I : restituisce la z del report (pag 8)
}

### compute tr(inv(A)*B) via either exact or stochastic approximations
get.prod.trace <- function(A, B, exact = FALSE, nsample = 100, verbose = FALSE){
  
  if (exact) {
    # exact calculation
    FF <- solve(A, B)
    tr <- sum(diag(FF))
  } else {
    # stochastic approximation
    nobs <- nrow(A)
    step <- as.integer(nsample / 20)
    
    # sparse Cholesky factorization
    L <- Cholesky(A)   # original
    #L <- chol(A)       # modifica
    if (verbose) cat(" GCV : Cholesky |")
    
    # forward substitution
    U <- sample(c(-1, +1), size = nobs * nsample, replace = TRUE)
    U <- matrix(U, nr = nobs, nc = nsample)
    V <- solve(L, U)            # original 
    #V <- solve(t(L)*L, U)      # In original, usa solve(L, U) (che va bene se uso "Cholesky", ma non se uso "chol")  
                                # NOTA: L ? l'upper se uso chol ! 
    if (verbose) cat(" Solve |")
    
    # incremental estimation
    tr <- 0
    for (k in 1:nsample) {
      uk <- as.vector(U[,k])
      vk <- as.vector(V[,k])
      tk <- as.numeric((uk %*% B %*% vk))
      tr <- tr + tk / nsample
      
      if (verbose && (k %% step == 0)) cat(".")
    }
    
    if (verbose) cat("| End \n")
    if (verbose) print.output(-1, 1.0, 1.0, 1.0)
  }
  
  # output
  return(tr)
}

### check loss function
loss.fun <- function(x, mu, alpha, eps = 0.0) {   # I : loss.fun ? rho_alpha (quantile check function)
  loss <- numeric(length(x))
  if (eps == 0.0) {
    # exact loss function
    loss <- 0.5 * abs(x - mu) + (alpha - 0.5) * (x - mu)
  } else {
    # smooth approximation
    loss <- (alpha - 1) * (x - mu) + eps * log1pexp((x - mu) / eps) #log(1 + exp(arg))
  }
  return(loss)
}

### log-likelihood function
llik.fun <- function(x, mu, sigma, alpha, eps = 0.0,
                     weights = rep(1, length(x))) {
  loss <- loss.fun(x, mu, alpha, eps)
  llik <- sum(log(weights / sigma) - loss / sigma)
  return(llik)
}

### absolute relative change
absmax <- function(x_new, x_old) {
  x_new <- as.vector(x_new)
  x_old <- as.vector(x_old)
  max_abs_diff <- as.numeric(max(abs(x_new - x_old)))
  max_abs_val  <- as.numeric(max(abs(c(x_new, x_old))))
  max_rel_diff <- as.numeric(max_abs_diff / max_abs_val)
  return(max_rel_diff)
}

### print output summary
print.output <- function(it, dx, df, ff) {
  dx_ <- round(dx, digits = 5)
  df_ <- round(df, digits = 5)
  ff_ <- round(ff, digits = 5)
  
  if (it == -1) {
    cat("------------------------------------------------------\n")
  } else if(it == 0) {
    cat(" iter: \t delta.par: \t delta.llik: \t pen.log.lik:   \n")
    cat("------------------------------------------------------\n")
  } else if(it == 1) {
    cat(paste("", it, "\t -.----- \t -.----- \t", ff_, "\n", sep=" ", collapse=""))
  } else if(it >= 2) {
    cat(gettextf(c(" %.0f \t", "%.5f \t", "%.5f \t", "%.5f \n"), c(it, dx_, df_, ff_)))
  }
}





