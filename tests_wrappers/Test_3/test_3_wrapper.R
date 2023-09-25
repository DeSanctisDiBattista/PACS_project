###------------------------TEST 3-----------------------------------------------

# domain:       c-shaped
# sampling:     areal
# penalization: laplacian 
# covariates:   yes
# BC:           no
# order FE:     1


## RMK: set the directory to the source file location 
##      ("Session" --> "Set Working Directory" --> "To source file location")

rm(list = ls())
graphics.off()

alphas = c(0.1,0.5,0.9)

source("../utilities.R")

library(Matrix)
library(plot3D)
library(geoR)
library(fdaPDE)
library(fdaPDE2)


## Load mesh
load("horseshoe2D_areal.RData")
nodes = mesh$nodes

mesh_data <- list(
  "nodes"    = mesh$nodes,
  "edges"    = mesh$edges,
  "elements" = mesh$triangles,
  "neigh"    = mesh$neighbors,
  "boundary" = mesh$nodesmarkers
)

incidence_matrix = as.matrix(read.csv("incidence_matrix.csv")[, -1])
n = dim(incidence_matrix)[1]    # number of subdomains
nnodes = dim(nodes)[1]

## define regularizing term
pde <- new(Laplacian_2D_Order1, mesh_data)

## set zero forcing term
quadrature_nodes <- pde$get_quadrature_nodes()
f <- as.matrix(rep(0., times = dim(quadrature_nodes)[1]))
pde$set_forcing_term(as.matrix(f))

## define SQRPDE model
model <- new(SQRPDE_Laplacian_2D_Areal, pde)

## Set subdomains
model$set_subdomains(incidence_matrix) # -> tolto t 

## beta true
beta_true = 3
## number of simulations 
M = 10
## prepare room for results
RMSE_matrix = matrix(nrow = M, ncol=length(alphas))
SpRMSE_matrix = matrix(nrow = n, ncol=length(alphas))
beta_hat_matrix = matrix(nrow = M, ncol=length(alphas))

alpha_index = 0
for(alpha in alphas){
  alpha_index = alpha_index + 1 
  RMSE_vec = c()
  beta_hat_vec = c()
  fn_hat_vec = data.frame(matrix(ncol = M, nrow = n))
  fn_true_vec = data.frame(matrix(ncol = M, nrow = n))
  
  model$set_alpha(alpha)
  
  for(m in 1:M){
    
    data_result = data_creation_areal(mesh = mesh, incidence_matrix = incidence_matrix,
                                      alpha = alpha, beta_true = beta_true, m = m)
    data = data_result$data
    X = data_result$covariates
    fn_true = data_result$fn_true
    f_true = data_result$f_true
    sol_integrand = data_result$sol_integrand
    sd_true = data_result$std.field
    
    ## set smoothing parameter
    suppressWarnings({
      lambda <- read.csv(paste0(getwd(), "/alpha_", 
                                as.character(alpha*100), "/sim_", 
                                m, "/LambdaCpp.csv"), header = TRUE)$x
    })
    
    model$set_lambda_s(lambda)
    
    ## set observations
    model$set_observations(as.matrix(data))
    
    ## set covariates
    model$set_covariates(as.matrix(X))
    
    ## solve smoothing problem
    model$solve()
    
    f_hat = model$f()
    fn_hat = model$fn()
    beta_hat = model$beta()
    
    RMSE_vec = c(RMSE_vec, RMSE(f_true, f_hat))
    beta_hat_vec = c(beta_hat_vec, beta_hat)
    fn_hat_vec[, m] = fn_hat
    fn_true_vec[, m] = fn_true
    
  }
  SpRMSE_vec = SpRMSE(fn_hat_vec, fn_true_vec)
  
  SpRMSE_matrix[, alpha_index] = SpRMSE_vec
  RMSE_matrix[, alpha_index] = RMSE_vec
  beta_hat_matrix[, alpha_index] = beta_hat_vec
  
}


# Plot boxplot 
ylim_rmse = c(min(min(RMSE_matrix), min(SpRMSE_matrix)), 
              max(max(RMSE_matrix), max(SpRMSE_matrix)))

boxplot(RMSE_matrix, ylim = ylim_rmse)
title(expression(RMSE), cex = 2)

boxplot(SpRMSE_matrix, ylim = ylim_rmse)
title(expression(SpRMSE), cex = 2)

boxplot(beta_hat_matrix)
abline(h=beta_true, col = "darkgreen", lwd=2.5, lty=1)
title(expression(beta), cex = 2)










