###------------------------TEST 9: R-C++ comparison-----------------------------------------------

# domain:       c-shaped
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
data(horseshoe2D)
options("scipen" = 999)
mesh_temp <- create.mesh.2D(
  nodes    = horseshoe2D$boundary_nodes,
  segments = horseshoe2D$boundary_segments
)

seq_max_area_mesh = c(0.008, 0.004, 0.002, 0.001)
seq_N = c()
for(area in seq_max_area_mesh) {
  nodes <- refine.mesh.2D(mesh_temp, maximum_area = area)$nodes  
  N <- nrow(nodes)
  seq_N = c(seq_N, N)
}

## generate locations
area_loc =  0.002
locations <- refine.mesh.2D(mesh_temp, maximum_area = area_loc )$nodes 
n <- nrow(locations)

## utility for data generation 

mean.function = function(locations){ 
  f = fs.test(locations[,1], locations[,2], b=1, exclude = FALSE)
  return(f)
} 

std.function = function(locations){ 
  f = fs.test(locations[,1], locations[,2], b=0.5, exclude = FALSE)
  min = min(f)
  f = 0.2*(f - 1.2*min(f))
  return(f)
}


## beta true 
beta_true = c(2, -1)


M = 10  # number of simulations

# prepare room for results 
RMSE_Cpp = matrix(nrow = M, ncol = length(seq_N))
RMSE_R = matrix(nrow = M, ncol = length(seq_N))
beta_hat_1_Cpp = matrix(nrow = M, ncol = length(seq_N))
beta_hat_1_R = matrix(nrow = M, ncol = length(seq_N))
beta_hat_2_Cpp = matrix(nrow = M, ncol = length(seq_N))
beta_hat_2_R = matrix(nrow = M, ncol = length(seq_N))
time_Cpp = matrix(nrow = M, ncol = length(seq_N))
time_R = matrix(nrow = M, ncol = length(seq_N))

for(j in 1:length(seq_N)){
  
  N = seq_N[j]
  
  mesh <- refine.mesh.2D(mesh_temp, maximum_area = seq_max_area_mesh[j])
  N_string = as.character(N)
  nodes = mesh$nodes
  
  mesh_data <- list(
    "nodes"    = mesh$nodes,
    "edges"    = mesh$edges,
    "elements" = mesh$triangles,
    "neigh"    = mesh$neighbors,
    "boundary" = mesh$nodesmarkers
  )
  
  ## define regularizing term
  pde <- new(Laplacian_2D_Order1, mesh_data)
  
  ## set zero forcing term
  quadrature_nodes <- pde$get_quadrature_nodes()
  f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
  pde$set_forcing_term(as.matrix(f))
  
  ## define SQRPDE model
  model <- new(SQRPDE_Laplacian_2D_GeoStatLocations, pde)
  model$set_alpha(alpha)
  
  # Locations 
  model$set_locations(locations)

  for(m in 1:M){
    
    # Covariates
    set.seed(21*m)
    cov1 <- (rnorm(n, mean = 0 , sd = 1))
    cov2 <- (rexp(n, rate = 1))
    X = cbind(cov1 ,cov2)
    
    # Data
    data = rnorm(n, mean = mean.function(locations), 
                 sd = std.function(locations)) + X%*%beta_true 
    
    # True solution
    fn_true = mean.function(locations) + qnorm(alpha) * std.function(locations)
    f_true = mean.function(nodes) + qnorm(alpha) * std.function(nodes)

    ### C++ 
      ## set smoothing parameter (optimized through exact GCV strategy)
      suppressWarnings({
        lambda_Cpp <- read.csv(paste0("N_", as.character(N), "/sim_", as.character(m), 
                                      "/LambdaCpp_massFALSE.csv"), header = FALSE)$V1
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
        lambda_R <- read.csv(paste0(getwd(),"/N_", as.character(N), "/sim_", as.character(m), 
                                    "/LambdaR_50.csv"))$x
      })

      t2 = Sys.time()
      
      ## build the model 
      
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
                    seq_x = seq_N, names = as.character(seq_N), 
                    ylim = c(min(min(RMSE_R), min(RMSE_Cpp)), 
                             max(max(RMSE_R), max(RMSE_Cpp))), 
                    xlab = "Sample size", shift = TRUE)

# Betas
plot_boxplot_RvsCpp(matrix_R = beta_hat_1_R, matrix_Cpp = beta_hat_1_Cpp, 
                    true = beta_true[1],
                    seq_x = seq_N, names = as.character(seq_N), 
                    ylim = c(min(min(beta_hat_1_R), min(beta_hat_1_Cpp)), 
                             max(max(beta_hat_1_R), max(beta_hat_1_Cpp))), 
                    xlab = "Sample size", shift = TRUE)

plot_boxplot_RvsCpp(matrix_R = beta_hat_2_R, matrix_Cpp = beta_hat_2_Cpp, 
                    true = beta_true[2],
                    seq_x = seq_N, names = as.character(seq_N), 
                    ylim = c(min(min(beta_hat_2_R), min(beta_hat_2_Cpp)), 
                             max(max(beta_hat_2_R), max(beta_hat_2_Cpp))), 
                    xlab = "Sample size", shift = TRUE)

# LogTime 
plot_boxplot_RvsCpp(matrix_R = log10(time_R), matrix_Cpp = log10(time_Cpp), 
                    seq_x = seq_N, names = as.character(round(log10(seq_N), 1)), 
                    ylim = log10(c(min(min(time_R), min(time_Cpp)), 
                             max(max(time_R), max(time_Cpp)))), 
                    xlab = paste0(expression(log10), "(sample size)"), shift = TRUE)

