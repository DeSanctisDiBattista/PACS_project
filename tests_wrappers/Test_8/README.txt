## Test characteristics
# domain:       unit square [0,1]x[0,1] 
# sampling:     locations != nodes  
# penalization: laplacian 
# covariates:   no
# BC:           no
# order FE:     1

This test is meant to analyze the effect of the mass lumping option on both R and C++ software. The results are stored in "lumpTRUE" and "lumpFALSE" when mass lumping is activated and disactivated, respectively.  