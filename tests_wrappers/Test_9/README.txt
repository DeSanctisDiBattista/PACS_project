## Test characteristics
# domain:       C-shaped  
# sampling:     locations != nodes  
# penalization: laplacian 
# covariates:   yes
# BC:           no
# order FE:     1

This test is meant to compare the accuracy and performances of R and C++ in the case where the two algorithms differ, namely the semiparametric scenario. 
The results are shown as function of an increasing sequence of mesh nodes (742, 1452, 2789, 5502). 

Each "N_" folder denotes number of nodes of that run. The "sim_" folders contain the results of each simulation.

The results of the two software have been obtained with the default mass lumping options, respectively. That is: 
mass lumping = TRUE, for R
mass lumping = FALSE, for C++ 