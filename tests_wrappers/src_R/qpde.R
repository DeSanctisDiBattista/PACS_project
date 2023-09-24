#' file: qpde.R
#' author: Cristian Castiglione
#' last change: 18/02/2023
#' content: wrapper

### Libraries
library(fdaPDE)
library(mgcv)
library(qgam)
library(Matrix)
library(plot3D)

### SQR-PDE functions
source("../src_R/pdematrices.R")
source("../src_R/utilities.R")
source("../src_R/qpdefit.R")
source("../src_R/gcvsearch.R")
