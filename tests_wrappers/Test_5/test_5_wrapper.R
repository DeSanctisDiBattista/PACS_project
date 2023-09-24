###------------------------TEST 5-----------------------------------------------

# domain:       c-shaped network  
# sampling:     locations = nodes
# penalization: laplacian 
# covariates:   no
# BC:           no
# order FE:     1

## RMK: set the directory to the source file location 
##      ("Session" --> "Set Working Directory" --> "To source file location")

rm(list = ls())
graphics.off()

seed = 476813
set.seed(seed)

alphas = c(0.1, 0.5, 0.9)

source("../utilities.R")
library(Matrix)
library(plot3D)
library(purrr)
library(spatstat)
library(fdaPDE)
library(fdaPDE2)


## create mesh, C-shaped network

mesh_components = generate_mesh_Cnetwork(eps = 0.5)
nodes = 
mesh_data <- list(
  "nodes"    = mesh_components$nodes,
  "edges"    = mesh_components$edges,
  "elements" = mesh_components$edges,
  "neigh"    = mesh_components$neigh,
  "boundary" = matrix(mesh_components$boundary, nrow=length(mesh_components$boundary), ncol=1)
)

nodes = mesh_data$nodes

pde <- new(Laplacian_1_5D_Order1, mesh_data)

# smooth function generating the data 
data_generation = function(x,y,z=0){
  coe <- function(x,y) 1/2*sin(5*pi*x)*exp(-x^2)+1
  return(sin(2*pi*(coe(y,1)*x*cos(z-2)-y*sin(z-2)))*cos(2*pi*(coe(y,1)*x*cos(z-2+pi/2) + 
                                                                coe(x,1)*y*sin((z-2)*pi/2))))
}

## set zero forcing term
quadrature_nodes <- pde$get_quadrature_nodes()
f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
pde$set_forcing_term(as.matrix(f))

## define SQRPDE model
model <- new(SQRPDE_Laplacian_1_5D_GeoStatNodes, pde)

## simulations
M = 10
locations = nodes
n = dim(locations)[1]

## prepare room for results
RMSE_matrix = matrix(nrow = M, ncol=length(alphas))
SpRMSE_matrix = matrix(nrow = n, ncol=length(alphas))

alpha_index = 0
for(alpha in alphas){
  alpha_index = alpha_index + 1 
  model$set_alpha(alpha)
  
  RMSE_vector = rep(0, M)
  SpRMSE_vector = rep(0, n)
  for(m in 1:M){
    
    ## data
    data_unnnoised = data_generation(nodes[, 1], nodes[, 2])
    sd_noise = 0.1*(max(data_unnnoised)-min(data_unnnoised))
    set.seed(seed + m*100)
    noise = rnorm(n,0,sd_noise)
    data = data_unnnoised + noise
    
    set.seed(seed + m*100)
    f_true =  data_generation(nodes[, 1], nodes[, 2]) + qnorm(alpha,0,sd_noise)
    fn_true = f_true
    
    ## set smoothing parameter (optimized through exact GCV strategy)
    suppressWarnings({
      lambda <- read.csv(paste0(getwd(), "/alpha_", as.character(alpha*100), "/sim_", as.character(m), 
                                "/LambdaCpp.csv"), header = FALSE)$V1
    })
    model$set_lambda_s(lambda)
    
    ## setters
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

boxplot(RMSE_matrix)
title(expression(RMSE), cex = 2)

boxplot(SpRMSE_matrix)
title(expression(SpRMSE), cex = 2)





