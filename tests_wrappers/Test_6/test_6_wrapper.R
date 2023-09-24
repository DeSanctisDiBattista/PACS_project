###------------------------TEST 6-----------------------------------------------

# domain:       Hub
# sampling:     locations = nodes
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


## create mesh, hub domain
data(hub2.5D)
mesh_base = create.mesh.2.5D(nodes = hub2.5D$hub2.5D.nodes, triangles = hub2.5D$hub2.5D.triangles)
mesh = refine.by.splitting.mesh.2.5D(mesh = mesh_base)
mesh_data <- list(
  "nodes"    = mesh$nodes,
  "edges"    = mesh$edges,
  "elements" = mesh$triangles,
  "neigh"    = mesh$neighbors,
  "boundary" = mesh$nodesmarkers
)

nodes = mesh$nodes
locations = nodes
n = dim(locations)[1]

## mean and standard deviation functions for the data generation 
mean.function = function(p){    
  
  x = p[,1]
  y = p[,2]
  z = p[,3]
  
  f = 1*sin(2*pi*x)+1*sin(2*pi*y)+3*sin(2*pi*z)
  
  return(f)
}
std.function = function(p){    
  
  x = p[,1]
  y = p[,2]
  z = p[,3]
  
  f = 2*abs(z-0.5)+0.1
  
  return(f)
}

## true beta vector 
beta_true = c(2,1)
q = length(beta_true)

## covariates
set.seed(seed)
cov1 <- (rnorm(n, mean = 0 , sd = 1))
set.seed(seed)
cov2 <- (rexp(n, rate = 1))
X = cbind(cov1 ,cov2)

## define regularizing term
pde <- new(Laplacian_2_5D_Order1, mesh_data)   

## set zero forcing term
quadrature_nodes <- pde$get_quadrature_nodes()
f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
pde$set_forcing_term(as.matrix(f))

## define SQRPDE model
model <- new(SQRPDE_Laplacian_2_5D_GeoStatNodes, pde)

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
    data_result = data_creation_hetero_fun(n, mean.function, std.function, beta_true = beta_true, 
                                           locations=nodes, nodes=nodes, covariates = X, seed = 37*m)
    data = data_result$data
    fn_true = data_result$fn_true       # true field at locations 
    f_true = data_result$f_true         # true field at mesh nodes 
    
    ## set smoothing parameter (optimized through exact GCV strategy)
    suppressWarnings({
      lambda <- read.csv(paste0(getwd(), "/alpha_", as.character(alpha*100), 
                                "/sim_", as.character(m), "/LambdaCpp.csv"), header = FALSE)$V1
    })
    model$set_lambda_s(lambda)
    
    ## setters
    model$set_observations(as.matrix(data))
    model$set_covariates(X)
    
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




