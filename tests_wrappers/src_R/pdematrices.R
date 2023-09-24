#' file: pdematrices.R
#' author: Cristian Castiglione
#' last change: 18/02/2023
#' content:
#' - function: get.FEM.Psi.Matrix
#' - function: get.FEM.Mass.Matrix
#' - function: get.FEM.Stiff.Matrix
#' - function: get.FEM.PDE.Matrix


### build the FEM design matrix
get.FEM.Psi.Matrix <- function(mesh, locs, sparse = TRUE) {
  coeff <- diag(dim(mesh$nodes)[1])
  basis <- fdaPDE::create.FEM.basis(mesh)
  fun <- fdaPDE::FEM(coeff, basis)
  psi <- fdaPDE::eval.FEM(fun, locs)
  if (sparse) psi <- as(psi, "sparseMatrix")
  return(psi)
}

### build the FEM mass matrix R0
get.FEM.Mass.Matrix <- function(mesh, sparse = TRUE) {
  basis <- fdaPDE::create.FEM.basis(mesh)
  mass <- fdaPDE:::CPP_get.FEM.Mass.Matrix(basis)
  if (!sparse) mass <- as.matrix(mass)
  return(mass)
}

### build the FEM stiffness matrix R1 for an isotropic diffusion PDE
get.FEM.Stiff.Matrix <- function(mesh, sparse = TRUE) {
  basis <- fdaPDE::create.FEM.basis(mesh)
  stiff <- fdaPDE:::CPP_get.FEM.Stiff.Matrix(basis)
  if (!sparse) stiff <- as.matrix(stiff)
  return(stiff)
}

### build the FEM stiffness matrix R1 for a general elliptic PDE
get.FEM.PDE.Matrix <- function(
  mesh, K = diag(2), b = rep(0, 2), c = 0, sparse = TRUE) {
  basis <- fdaPDE::create.FEM.basis(mesh)
  param <- list(K = K, b = b, c = c)
  obs <- rep(0, nrow(mesh$nodes))
  pdemat <- fdaPDE:::CPP_get.FEM.PDE.Matrix(obs, basis, param)
  if (!sparse) pdemat <- as.matrix(pdemat)
  return(pdemat)
}
