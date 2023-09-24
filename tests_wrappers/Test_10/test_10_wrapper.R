###------------------------TEST 10: R-C++ comparison-----------------------------------------------

# domain:       unit square [0,1] x [0,1]
# sampling:     locations != nodes
# penalization: laplacian 
# covariates:   yes
# BC:           no
# order FE:     1

# The two methods are compared with their default mass lumping options, namely 
# lumping = FALSE for C++ and lumping = TRUE for R 

## RMK: set the directory to the source file location 
##      ("Session" --> "Set Working Directory" --> "To source file location")

rm(list = ls())
graphics.off()

options(digits=16)
alpha = 0.5
alpha_string = as.character(alpha*100)

seed = 476813

source("../src_R/qpde.R")
source("../utilities.R")

library(Matrix)
library(plot3D)
library(mvnormtest)
library(geoR)
library(fdaPDE2)


## create mesh, unit square
nxx = 44
x <- seq(0, 1, length.out = nxx)
y <- x
mesh <- fdaPDE::create.mesh.2D(expand.grid(x, y))
nodes = mesh$nodes
nnodes = dim(nodes)[1]

mesh_data <- list(
  "nodes"    = mesh$nodes,
  "edges"    = mesh$edges,
  "elements" = mesh$triangles,
  "neigh"    = mesh$neighbors,
  "boundary" = mesh$nodesmarkers
)

## utility for data generation 
data_generation <- function(x, y, z = 1){
  return( 36*(1/3*x^3-1/2*x^2)*(1/3*y^3-1/2*y^2) ) 
}

## beta true 
beta_true = c(2, -1)

## define regularizing term
pde <- new(Laplacian_2D_Order1, mesh_data)

## set zero forcing term
quadrature_nodes <- pde$get_quadrature_nodes()
f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
pde$set_forcing_term(as.matrix(f))

## sequence of sample sizes
seq_nxx_loc =  c(16, 22, 32, 45, 63) 
seq_n = seq_nxx_loc^2 

M = 10  # number of simulations

# prepare room for results 
RMSE_Cpp = matrix(nrow = M, ncol = length(seq_n))
RMSE_R = matrix(nrow = M, ncol = length(seq_n))
beta_hat_1_Cpp = matrix(nrow = M, ncol = length(seq_n))
beta_hat_1_R = matrix(nrow = M, ncol = length(seq_n))
beta_hat_2_Cpp = matrix(nrow = M, ncol = length(seq_n))
beta_hat_2_R = matrix(nrow = M, ncol = length(seq_n))
time_Cpp = matrix(nrow = M, ncol = length(seq_n))
time_R = matrix(nrow = M, ncol = length(seq_n))

for(j in 1:length(seq_n)){
  
  n = seq_n[j]
  
  ## define SQRPDE model
  model <- new(SQRPDE_Laplacian_2D_GeoStatLocations, pde)
  model$set_alpha(alpha)
  
  # Locations 
  locations = as.matrix(expand.grid(seq(0,1,length.out = sqrt(n)), seq(0,1,length.out = sqrt(n))))
  model$set_locations(locations)

  for(m in 1:M){
    
    # Covariates 
    set.seed(seed*m)
    cov1 <- rnorm(n, mean = 0, sd = 1)
    set.seed(seed*m)
    cov2 <- rexp(n, rate = 3)
    X = cbind(cov1, cov2)
    
    ## data
    data_result = data_creation_skewed(data_generation = data_generation,
                                       alpha = alpha, locations = locations, nodes = nodes, 
                                       beta_true = beta_true, covariates = X, seed = seed*2*m)
    data = data_result$data
    fn_true = data_result$fn_true       # true field at locations 
    f_true = data_result$f_true         # true field at mesh nodes 
    

    ### C++ 
    
      ## set smoothing parameter (optimized through exact GCV strategy)

      suppressWarnings({
        lambda_Cpp <- read.csv(paste0(getwd(),"/n_", as.character(n), "/sim_", as.character(m), 
                                      "/LambdaCpp.csv"), header = FALSE)$V1
      })

        
      model$set_lambda_s(lambda_Cpp)
      
      ## set observations and covariates 
      model$set_observations(as.matrix(data))
      model$set_covariates(X)
      
      ## solve smoothing problem
      t0 = Sys.time()
      model$solve()
      t1 = Sys.time()
      
      ## extract the results 
      f_hat_Cpp = model$f()
      fn_hat_Cpp = model$fn()
      beta_hat_Cpp = model$beta()
      
      ## RMSE & SpRMSE
      RMSE_Cpp[m, j] = RMSE(f_true, f_hat_Cpp)
      beta_hat_1_Cpp[m, j] = beta_hat_Cpp[1]
      beta_hat_2_Cpp[m, j] = beta_hat_Cpp[2]
      ## time 
      time_Cpp[m, j] = difftime(t1, t0, units = "secs")   
    
    
    ### R 
      
      suppressWarnings({
        lambda_R <- read.csv(paste0(getwd(),"/n_", as.character(n), 
                                    "/sim_", as.character(m), "/LambdaR.csv"))$x
      })

      t2 = Sys.time()
      
      ## build the model 
      mesh <- fdaPDE::create.mesh.2D(nodes = nodes, order = 1)
      # PDE parameters: diffusion tensor, transport vector, reaction term
      K <- diag(2)
      b <- rep(0, 2)
      c <- 0
      # FEM, mass, and stiffness matrices
      coeff <- diag(dim(nodes)[1]) 
      basis <- fdaPDE::create.FEM.basis(mesh)         # FEM basis expansion 
      fun <- fdaPDE::FEM(coeff, basis)
      psi <- as(fdaPDE::eval.FEM(fun, locations), "sparseMatrix")
      

      R0 <- diag(colSums(get.FEM.Mass.Matrix(mesh)))  # lumping 
      R1 <- get.FEM.PDE.Matrix(mesh, K, b, c)         # stiffness matrix for a general elliptic PDE
      
      t3 = Sys.time()
      time_build_model = difftime(t3, t2, units = "secs")
      
      ## fit 
      fit_R <- qgam.xfem_our(y = data, X = X, Z = psi, R0 = R0, R1 = R1, alpha = alpha,
                              lambdas = lambda_R, eps = 0.0, tune = 1.0,
                              maxiter = 200, tol = 1e-6, 
                              exact = TRUE, tol.weights = 1e-6,
                              mass_lumping_sys = TRUE, mass_lumping_GCV = TRUE,
                              measure.time = TRUE)
      ## extract the results 
      f_hat_R = fit_R$f
      fn_hat_R = fit_R$fn
      beta_hat_R = fit_R$beta
      
      ## RMSE & SpRMSE
      RMSE_R[m, j] = RMSE(f_true, f_hat_R)
      beta_hat_1_R[m, j] = beta_hat_R[1]
      beta_hat_2_R[m, j] = beta_hat_R[2]
      ## time 
      time_R[m, j] = fit_R$time_fit + time_build_model  # in seconds
      
  }
}


## Plots

# RMSE_f
plot_boxplot_RvsCpp(matrix_R = RMSE_R, matrix_Cpp = RMSE_Cpp, 
                    seq_x = seq_n, names = as.character(seq_n), 
                    ylim = c(min(min(RMSE_R), min(RMSE_Cpp)), 
                             max(max(RMSE_R), max(RMSE_Cpp))), 
                    shift = TRUE,
                    title = "RMSE")

# Betas
plot_boxplot_RvsCpp(matrix_R = beta_hat_1_R, matrix_Cpp = beta_hat_1_Cpp, 
                    seq_x = seq_n, names = as.character(seq_n), 
                    true = beta_true[1],
                    ylim = c(min(min(beta_hat_1_R), min(beta_hat_1_Cpp)), 
                             max(max(beta_hat_1_R), max(beta_hat_1_Cpp))), 
                    shift = TRUE,
                    title = "Beta_1")

plot_boxplot_RvsCpp(matrix_R = beta_hat_2_R, matrix_Cpp = beta_hat_2_Cpp, 
                    seq_x = seq_n, names = as.character(seq_n), 
                    true = beta_true[2],
                    ylim = c(min(min(beta_hat_2_R), min(beta_hat_2_Cpp)), 
                             max(max(beta_hat_2_R), max(beta_hat_2_Cpp))), 
                    shift = TRUE,
                    title = "Beta_2")

# LogTime 
plot_boxplot_RvsCpp(matrix_R = log10(time_R), matrix_Cpp = log10(time_Cpp), 
                    title = "Log(Time)",
                    seq_x = seq_n, names = as.character(round(log10(seq_n), 1)), 
                    ylim = log10(c(min(min(time_R), min(time_Cpp)), 
                             max(max(time_R), max(time_Cpp)))), 
                    shift = TRUE)

