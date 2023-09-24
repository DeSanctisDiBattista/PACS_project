## Test characteristics
# domain:       unit square [0,1]x[0,1]  
# sampling:     locations != nodes  
# penalization: laplacian 
# covariates:   yes
# BC:           no
# order FE:     1

This test is meant to compare the accuracy and performances of R and C++ in the case where the two algorithms differ, namely the semiparametric scenario. 
The results are shown as function of an increasing sequence of sample size (256, 484, 1024, 2025, 3969)

Each "n_" folder denotes the sample size of that run. The "sim_" folders contain the results of each simulation.

The results of the two software have been obtained with the default mass lumping options, respectively. That is: 
mass lumping = TRUE, for R
mass lumping = FALSE, for C++ 