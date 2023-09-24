### fdaPDE with SQR-PDE implementation 

This repository contains the development version of fdaPDE library. In particular we included the class for the Spatial Quantile Regression with PDE Regularization (SQR-PDE).

## Subfolder structure:
- /fdaPDE contains all C++ code (development version of the fdaPDE library forked from https://github.com/AlePalu/fdaPDE). The implemented class SQRPDE is contained in /fdaPDE/fdaPDE/models/regression. 

- /tests_wrappers contains: 
   -- the tests's folders /Test_1,...,/Test_6,/Test_8,...,/Test_10 to run the correspondent simulation tests from R, relying on the C++/R wrapper structure contained in ../fdaPDE/wrappers. 
     Each folder contains an .R file to run the test, a README file which briefly explains the structure of the test and a subfolders' structure to store the data for the tests. 
 - - /compile_fdaPDE2.R which contains the instructions to install the fdaPDE2 library required to run the test with R/C++ wrappers. 
 - - src_R which contains the source code made by C.Castiglione to run the R version algorithm of SQR-PDE. Some functions in src_R have been slightly modified to accomplish our purposes. 
 - - utilities.R which is an R script with some useful routines to run the tests


## Installation & test run:
Download the .zip file from the repository and unzip it to have at disposal the C++ source code.

To compile the C++ library, the following dependencies have to be installed:
 - a C++17 compliant compiler (gcc versions higher than 7 should be fine)
 - make
 - CMake
 - Eigen linear algebra library (version higher than 3.3)

 To generate the Makefile, compile and run the code, go to fdaPDE/test and execute from shell: 
   cmake CMakeLists.txt
   make
   ./fdaPDE_test

 This will run all the tests included in test/MainTest.cpp, in particular the GCV tests implemented in test/calibration/GCVTest.cpp   

To run the simulation tests relying on the R/C++ wrappers, the following R packages have to be installed: 
- fdaPDE2: to install it, follow the content in /compile_fdaPDE2.R. 
- fdaPDE: to install it run from R terminal the command install.packages('fdaPDE')
- minor packages (e.g. 'Matrix', 'plot3D', 'geoR'...): easly installed with the R command install.packages('PackageName')

