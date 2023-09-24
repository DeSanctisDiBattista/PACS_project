###------------------------TEST 4-----------------------------------------------

# domain:       unit sphere
# sampling:     locations != nodes
# penalization: laplacian 
# covariates:   yes
# BC:           no
# order FE:     1

## RMK: set the directory to the source file location 
##      ("Session" --> "Set Working Directory" --> "To source file location")

rm(list = ls())
graphics.off()

seed = 21
set.seed(seed)

alphas = c(0.1, 0.5, 0.9)

source("../utilities.R")
library(rgl)
library(Matrix)
library(plot3D)
library(fdaPDE)
library(fdaPDE2)


## utility for random locations generation 
generate_random_locations = function(n, seed){
  
  set.seed(seed + 4) 
  phi = runif(n,0,2*pi)
  set.seed(seed + 100) 
  theta = runif(n,0,2*pi)
  set.seed(seed + 1234) 
  rho = runif(n,0,0.9)
  # (setting different seeds avoid undesired constraints like 'theta = phi')
  
  x = rho*sin(phi)*cos(theta)
  y = rho*sin(phi)*sin(theta)
  z = rho*cos(phi)
  
  return(cbind(x,y,z))
  
}


## create mesh, unit sphere
load("sphere3Ddata.RData")
mesh = create.mesh.3D(sphere3Ddata$nodes, sphere3Ddata$tetrahedrons)
nodes = mesh$nodes
mesh_data <- list(
  "nodes"    = mesh$nodes,
  "edges"    = mesh$faces,
  "elements" = mesh$tetrahedrons,
  "neigh"    = mesh$neighbors,
  "boundary" = mesh$nodesmarkers
)


n = 2000
locations = generate_random_locations(n = n, seed = seed)


## mean and standard deviation functions for the data generation 
mean.function = function(p){    
  
  x = p[,1]
  y = p[,2]
  z = p[,3]
  
  f = -sin(x)-2*sin(y)-3*sin(z)
  
  return(f)
}
std.function = function(p){    
  
  x = p[,1]
  y = p[,2]
  z = p[,3]
  
  center = rep(0,3)
  
  f =  1/10 * 1/( ((x-center[1])^2 + (y-center[2])^2 + (z-center[3])^2)^1/4 + 0.1 ) 
  
  return(f)
}

## true beta vector 
beta_true = c(2,1)
q = length(beta_true)

## covariates
set.seed(12345)
cov1 <- (rnorm(n, mean = 0 , sd = 1))
set.seed(12345)
cov2 <- (rexp(n, rate = 1))
X = cbind(cov1 ,cov2)

## define regularizing term
pde <- new(Laplacian_3D_Order1, mesh_data)

## set zero forcing term
quadrature_nodes <- pde$get_quadrature_nodes()
f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
pde$set_forcing_term(as.matrix(f))

## define SQRPDE model
model <- new(SQRPDE_Laplacian_3D_GeoStatLocations, pde)
model$set_locations(locations)

M = 10  # number of simulations

## prepare room for results
RMSE_matrix = matrix(nrow = M, ncol=length(alphas))
SpRMSE_matrix = matrix(nrow = n, ncol=length(alphas))
beta1_matrix = matrix(nrow = M, ncol=length(alphas))
beta2_matrix = matrix(nrow = M, ncol=length(alphas))

alpha_index = 0
for(alpha in alphas){
  alpha_index = alpha_index + 1 
  model$set_alpha(alpha)

  RMSE_vector = rep(0, M)
  SpRMSE_vector = rep(0, n)
  beta_hat_1_vector = rep(0, M)
  beta_hat_2_vector = rep(0, M)
  
  for(m in 1:M){
    
    ## data
    data_list = data_creation_hetero_fun(n, mean.function, std.function, 
                                     beta_true = beta_true, 
                                     locations = locations, nodes = nodes, 
                                     covariates = X, seed = seed*m)
    data = data_list$data
    f_true = data_list$f_true
    fn_true = data_list$fn_true

    ## set smoothing parameter (optimized through exact GCV strategy)
    suppressWarnings({
      lambda <- read.csv(paste0(getwd(), "/alpha_", 
                                as.character(alpha*100), "/sim_", 
                                m, "/LambdaCpp.csv"), header = FALSE)$V1
    })
    
    model$set_lambda_s(lambda)
    
    ## set observations
    model$set_observations(as.matrix(data))
    
    ## set covariates
    model$set_covariates(as.matrix(X))
    
    ## solve smoothing problem
    model$solve()
    
    ## extract the results 
    f_hat = model$f()
    fn_hat = model$fn()
    beta_hat = model$beta()
    
    ## RMSE & SpRMSE
    RMSE_vector[m] = RMSE(f_true, f_hat)
    SpRMSE_vector = SpRMSE_vector + (fn_hat - fn_true)^2
    
    beta_hat_1_vector[m] = beta_hat[1]
    beta_hat_2_vector[m] = beta_hat[2]
    
    
  }
  SpRMSE_vector = sqrt(SpRMSE_vector/M)
  
  SpRMSE_matrix[, alpha_index] = SpRMSE_vector
  RMSE_matrix[, alpha_index] = RMSE_vector
  beta1_matrix[, alpha_index] = beta_hat_1_vector
  beta2_matrix[, alpha_index] = beta_hat_2_vector
  
}

## Plot results 

ylim_rmse = c(min(min(RMSE_matrix), min(SpRMSE_matrix)), 
              max(max(RMSE_matrix), max(SpRMSE_matrix)))

boxplot(RMSE_matrix)
title(expression(RMSE), cex = 2)

boxplot(SpRMSE_matrix)
title(expression(SpRMSE), cex = 2)

boxplot(beta1_matrix)
title(expression(beta[1]), cex = 2)
abline(h = beta_true[1], col = "darkgreen", lwd = 2.5)

boxplot(beta2_matrix)
title(expression(beta[2]), cex = 2)
abline(h = beta_true[2], col = "darkgreen", lwd = 2.5)

