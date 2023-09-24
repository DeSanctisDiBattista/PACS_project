###------------------------TEST 2-----------------------------------------------

# domain:       unit square [0,1]x[0,1] 
# sampling:     locations = nodes
# penalization: diffusion constant coefficient 
# covariates:   no
# BC:           no
# order FE:     1

## RMK: set the directory to the source file location 
##      ("Session" --> "Set Working Directory" --> "To source file location")

rm(list = ls())
graphics.off()

seed = 7893475
alphas = c(0.1, 0.5, 0.9)

source("../utilities.R")

library(Matrix)
library(plot3D)
library(geoR)
library(fdaPDE)
library(fdaPDE2)


## create mesh, unit square
nxx = 32
x <- seq(0, 1, length.out = nxx)
y <- x
locations <- expand.grid(x, y)
n = dim(locations)[1]
mesh <- fdaPDE::create.mesh.2D(locations)
nodes = mesh$nodes
nnodes = dim(nodes)[1]

mesh_data <- list(
  "nodes"    = mesh$nodes,
  "edges"    = mesh$edges,
  "elements" = mesh$triangles,
  "neigh"    = mesh$neighbors,
  "boundary" = mesh$nodesmarkers
)

pde_data <- list(
  "diffusion" = cbind(c(1,0), c(0,4)), 
  "transport" = rep(0,2), 
  "reaction" = 0 
)

## define regularizing term
pde <- new(ConstantCoefficients_2D_Order1, mesh_data)
pde$set_PDE_parameters(pde_data)

## set zero forcing term
quadrature_nodes <- pde$get_quadrature_nodes()
f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
pde$set_forcing_term(as.matrix(f))

## define SQRPDE model
model <- new(SQRPDE_ConstantCoefficients_2D_GeoStatNodes, pde)

M = 10  # number of simulations

## prepare room for results
RMSE_matrix = matrix(nrow = M, ncol=length(alphas))
SpRMSE_matrix = matrix(nrow = n, ncol=length(alphas))

alpha_index = 0
for(alpha in alphas){
  alpha_index = alpha_index + 1 
  RMSE_vector = vector(length = M)
  SpRMSE_vector = rep(0, n)
  
  model$set_alpha(alpha)
  for(m in 1:M){
    
    ## data
    data_generation <- function(x, y, z = 1){
      a1 <- 1
      a2 <- 4
      a1*sin(2*pi*(x+0.17))*cos(2*pi*y)+a2*sin(3*pi*(x+0.17))
    }   
    data_result = data_creation_skewed(data_generation = data_generation, 
                                           alpha = alpha, locations = locations,
                                           nodes = nodes, seed = seed*m)
    data = data_result$data
    fn_true = data_result$fn_true                # true field at locations 
    f_true = data_result$f_true                  # true field at mesh nodes 
    noise_true_sd = data_result$noise_true_sd    # true noise std deviation
    
    ## set smoothing parameter (optimized through exact GCV strategy)
    suppressWarnings({
      lambda <- read.csv(paste0(getwd(), "/alpha_", as.character(alpha*100), 
                                "/sim_", as.character(m), "/LambdaCpp.csv"), header = FALSE)$V1
    })
    model$set_lambda_s(lambda)
    
    ## set observations, tolerances and the quantile order 
    model$set_observations(as.matrix(data))

    ## solve smoothing problem
    model$solve()
    
    ## extract the results 
    f_hat = model$f()
    fn_hat = model$fn()
    
    ## RMSE & SpRMSE
    RMSE_vector[m] = RMSE(f_true, f_hat)
    SpRMSE_vector = SpRMSE_vector + (fn_hat - fn_true)^2
    
    
  }
  SpRMSE_vector = sqrt(SpRMSE_vector/M)
  
  SpRMSE_matrix[, alpha_index] = SpRMSE_vector
  RMSE_matrix[, alpha_index] = RMSE_vector
}


## Plot results 

ylim_rmse = c(min(min(RMSE_matrix), min(SpRMSE_matrix)), 
              max(max(RMSE_matrix), max(SpRMSE_matrix)))

boxplot(RMSE_matrix, ylim = ylim_rmse)
title(expression(RMSE), cex = 2)

boxplot(SpRMSE_matrix, ylim = ylim_rmse)
title(expression(SpRMSE), cex = 2)


































