###------------------------TEST 8: Mass lumping comparison----------------------

# domain:       unit square [0,1] x [0,1]
# sampling:     locations != nodes
# penalization: laplacian 
# covariates:   no
# BC:           no
# order FE:     1


## RMK: set the directory to the source file location 
##      ("Session" --> "Set Working Directory" --> "To source file location")

rm(list = ls())
graphics.off()

seed = 7893475
alpha = 0.5

source("../src_R/qpde.R")   # load the source files to launch the R version of the code 
source("../utilities.R")

library(Matrix)
library(plot3D)
library(geoR)
library(fdaPDE)
library(fdaPDE2)


## utility for data generation 
data_generation = function(x, y, z = 1){
  coe <- function(x,y) 1/2*sin(5*pi*x)*exp(-x^2)+1
  sin(2*pi*(coe(y,1)*x*cos(z-2)-y*sin(z-2)))*cos(2*pi*(coe(y,1)*x*cos(z-2+pi/2) + 
                                                         coe(x,1)*y*sin((z-2)*pi/2)))
}

## create locations
n = 784
locations = as.matrix(expand.grid(seq(0,1,length.out=sqrt(n)), seq(0,1,length.out=sqrt(n))))

M = 10  # number of simulations

## sequence of mesh nodes
seq_N =  c(256,484,1024,2025,3969)

mat_RMSE_R_lumpFALSE = matrix(ncol = M, nrow = length(seq_N))
mat_RMSE_Cpp_lumpFALSE = matrix(ncol = M, nrow = length(seq_N))
mat_RMSE_R_lumpTRUE = matrix(ncol = M, nrow = length(seq_N))
mat_RMSE_Cpp_lumpTRUE = matrix(ncol = M, nrow = length(seq_N))

mat_time_R_lumpFALSE = matrix(ncol = M, nrow = length(seq_N))
mat_time_Cpp_lumpFALSE = matrix(ncol = M, nrow = length(seq_N))
mat_time_R_lumpTRUE = matrix(ncol = M, nrow = length(seq_N))
mat_time_Cpp_lumpTRUE = matrix(ncol = M, nrow = length(seq_N))

## sequence of optimal lambdas (found with GCV exact strategy)
lambdas_lumpFALSE = 10^seq(-6.92, -6.8, by = 0.01)
lambdas_lumpTRUE = 10^seq(-6.68, -6.58, by = 0.01)

for(index.N in 1:length(seq_N)){
  ## create mesh, unit square
  nnodes = seq_N[index.N]
  nodes = expand.grid(seq(0,1,length.out=sqrt(nnodes)), seq(0,1,length.out=sqrt(nnodes)))
  mesh <- fdaPDE::create.mesh.2D(nodes = nodes, order = 1)
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
  
  RMSE_R_lumpFALSE = c()
  RMSE_Cpp_lumpFALSE = c()
  RMSE_R_lumpTRUE = c()
  RMSE_Cpp_lumpTRUE = c()
  
  time_R_lumpFALSE = c()
  time_Cpp_lumpFALSE = c()
  time_R_lumpTRUE = c()
  time_Cpp_lumpTRUE = c()
  
  for(m in 1:M){
    ## data
    data_result = data_creation_skewed(data_generation = data_generation, nodes = nodes, 
                                  locations = locations, alpha = alpha, seed = seed*m)
    data = data_result$data
    fn_true = data_result$fn_true       # true field at locations 
    f_true = data_result$f_true         # true field at mesh nodes 
    
    
    ### C++ lumping = FALSE
    pde$setMassLumpingSystem(FALSE)
    ## define SQRPDE model
    model_lumpFALSE <- new(SQRPDE_Laplacian_2D_GeoStatLocations, pde)
    model_lumpFALSE$set_alpha(alpha)
    
    model_lumpFALSE$set_observations(as.matrix(data))
    model_lumpFALSE$set_locations(locations)
    # Optimum lambda chosen as mean of the optimum sequence
    lambda_opt_Cpp_lumpFALSE = mean(lambdas_lumpFALSE)
    model_lumpFALSE$set_lambda_s(lambda_opt_Cpp_lumpFALSE)
    
    ## solve smoothing problem
    t0 = Sys.time()
    model_lumpFALSE$solve()
    t1 = Sys.time()
    
    ## RMSE & time
    RMSE_Cpp_lumpFALSE = c(RMSE_Cpp_lumpFALSE, RMSE(f_true, model_lumpFALSE$f()))
    time_Cpp_lumpFALSE = c(time_Cpp_lumpFALSE, difftime(t1, t0, units = "secs"))  
    
    
    ### R lumping = FALSE
    lambda_opt_R_lumpFALSE = lambda_opt_Cpp_lumpFALSE
    
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
    
    R0 <- get.FEM.Mass.Matrix(mesh)  # no lumping 
    
    R1 <- get.FEM.PDE.Matrix(mesh, K, b, c) # stiffness matrix for a general elliptic PDE
    
    t3 = Sys.time()
    time_build_model_lumpFALSE = difftime(t3, t2, units = "secs")
    
    ## fit 
    fit_R_lumpFALSE <- qgam.fem_our(y = data, Z = psi, R0 = R0, R1 = R1, alpha = alpha,
                                    lambdas = lambda_opt_R_lumpFALSE, eps = 0.0, tune = 1.0,
                                    maxiter = 200, tol = 1e-6, 
                                    exact = TRUE, tol.weights = 1e-6,
                                    mass_lumping = FALSE, 
                                    measure.time = TRUE)
    
    ## RMSE & time
    RMSE_R_lumpFALSE = c(RMSE_R_lumpFALSE, RMSE(f_true, fit_R_lumpFALSE$f))
    time_R_lumpFALSE= c(time_R_lumpFALSE, fit_R_lumpFALSE$time_fit + time_build_model_lumpFALSE)  # in seconds
    
    
    #---
    
    ### C++ lumping = TRUE
    pde$setMassLumpingSystem(TRUE)
    ## define SQRPDE model
    model_lumpTRUE <- new(SQRPDE_Laplacian_2D_GeoStatLocations, pde)
    model_lumpTRUE$set_alpha(alpha)
    model_lumpTRUE$set_observations(as.matrix(data))
    model_lumpTRUE$set_locations(locations)
    
    # Optimum lambda chosen as mean of the optimum sequence
    lambda_opt_Cpp_lumpTRUE = mean(lambdas_lumpTRUE)
    
    model_lumpTRUE$set_lambda_s(lambda_opt_Cpp_lumpTRUE)
    
    ## solve smoothing problem
    t0 = Sys.time()
    model_lumpTRUE$solve()
    t1 = Sys.time()
    
    ## RMSE & time
    RMSE_Cpp_lumpTRUE = c(RMSE_Cpp_lumpTRUE, RMSE(f_true, model_lumpTRUE$f()))
    time_Cpp_lumpTRUE = c(time_Cpp_lumpTRUE, difftime(t1, t0, units = "secs"))  
    
    
    ### R
    lambda_opt_R_lumpTRUE = lambda_opt_Cpp_lumpTRUE
    
    
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
    
    R1 <- get.FEM.PDE.Matrix(mesh, K, b, c) # stiffness matrix for a general elliptic PDE
    
    t3 = Sys.time()
    time_build_model_lumpTRUE = difftime(t3, t2, units = "secs")
    
    ## fit 
    fit_R_lumpTRUE <- qgam.fem_our(y = data,Z = psi, R0 = R0, R1 = R1, alpha = alpha,
                                   lambdas = lambda_opt_R_lumpTRUE, eps = 0.0, tune = 1.0,
                                   maxiter = 200, tol = 1e-6, 
                                   exact = TRUE, tol.weights = 1e-6,
                                   mass_lumping = TRUE, 
                                   measure.time = TRUE)
    
    ## RMSE & time
    RMSE_R_lumpTRUE = c(RMSE_R_lumpTRUE, RMSE(f_true, fit_R_lumpTRUE$f))
    time_R_lumpTRUE= c(time_R_lumpTRUE, fit_R_lumpTRUE$time_fit + time_build_model_lumpTRUE)  # in seconds
  }
  
  mat_RMSE_R_lumpFALSE[index.N, ] = RMSE_R_lumpFALSE
  mat_RMSE_Cpp_lumpFALSE[index.N, ] = RMSE_Cpp_lumpFALSE
  mat_RMSE_R_lumpTRUE[index.N, ] = RMSE_R_lumpTRUE
  mat_RMSE_Cpp_lumpTRUE[index.N, ] = RMSE_Cpp_lumpTRUE
  
  mat_time_R_lumpFALSE[index.N, ] = time_R_lumpFALSE
  mat_time_Cpp_lumpFALSE[index.N, ] = time_Cpp_lumpFALSE
  mat_time_R_lumpTRUE[index.N, ] = time_R_lumpTRUE
  mat_time_Cpp_lumpTRUE[index.N, ] = time_Cpp_lumpTRUE
  
}

# Plot results 

myred <- rgb(255, 0, 0, max = 255, alpha = 125)
myred_dark <- rgb(255, 0, 0, max = 255, alpha = 255)

myblue <- rgb(0, 0, 255, max = 255, alpha = 125)
myblue_dark <- rgb(0, 0, 255, max = 255, alpha = 255)
mygrey_dark = rgb(105, 105, 105, max = 255, alpha = 255)

mybarbiepink <- rgb(0.85*255, 0.09*255, 0.52*255, max = 255, alpha = 175)
mybarbiepink_dark <- rgb(85, 9, 52, max = 255, alpha = 255)

myturquoise <- rgb(0, 0.81*255, 0.82*255, max = 255, alpha = 175)
myturquoise_dark <- rgb(0, 0.81*255, 0.82*255, max = 255, alpha = 255)

ylim = c(min(mat_time_R_lumpFALSE, mat_time_R_lumpTRUE), 
         max(mat_time_R_lumpFALSE, mat_time_R_lumpTRUE))
boxplot(t(mat_time_R_lumpFALSE), 
        col = myturquoise, 
        medcol =  myturquoise_dark,
        boxcol =  myturquoise_dark,
        whiskcol = myturquoise_dark,
        staplecol = myturquoise_dark,
        outcol = myturquoise_dark,
        outpch = 19, ylim = ylim,
        main = "Computational times [secs] (R)")
boxplot(t(mat_time_R_lumpTRUE), 
        col = myblue, 
        medcol =  myblue_dark,
        boxcol =  myblue_dark,
        whiskcol = myblue_dark,
        staplecol = myblue_dark,
        outcol = myblue_dark,
        outpch = 19, add = TRUE)
legend(0.5, 1100, legend = c("lumpFALSE", "lumpTRUE"), col = c(myturquoise, myblue), 
       pch = 19)

ylim = c(min(mat_time_Cpp_lumpFALSE, mat_time_Cpp_lumpTRUE), 
         max(mat_time_Cpp_lumpFALSE, mat_time_Cpp_lumpTRUE))
boxplot(t(mat_time_Cpp_lumpFALSE), 
        col = myred, 
        medcol =  myred_dark,
        boxcol =  myred_dark,
        whiskcol = myred_dark,
        staplecol = myred_dark,
        outcol = myred_dark,
        outpch = 19, ylim = ylim, 
        main = "Computational times [secs] (C++)", xlab = "N", names = as.character(seq_N))
boxplot(t(mat_time_Cpp_lumpTRUE), add = TRUE,
        col = mybarbiepink, 
        medcol =  mybarbiepink_dark,
        boxcol =  mybarbiepink_dark,
        whiskcol = mybarbiepink_dark,
        staplecol = mybarbiepink_dark,
        outcol = mybarbiepink_dark,
        outpch = 19, ylim = ylim)
legend(0.5, 6, legend = c("lumpFALSE", "lumpTRUE"), col = c(myred, mybarbiepink), 
       pch = 19)














