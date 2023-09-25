rm(list = ls())

## RMK: set the directory to the source file location 
#      ("Session" --> "Set Working Directory" --> "To source file location")

## install dependencies (execute this only once, if you not have Rcpp installed yet)
# install.packages("Rcpp")
# install.packages("RcppEigen")

# Check if the packages are installed 
system.file(package='Rcpp')
system.file(package='RcppEigen')

# load Rcpp library
library(Rcpp)

setwd(paste0(getwd(), "/../fdaPDE/wrappers/R/"))

# update RCppExports.cpp
compileAttributes(".")

# install fdaPDE
install.packages(".", type="source", repos=NULL)
