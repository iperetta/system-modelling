# system-modelling
Files from thesis


List of files in this project:

** Classes (C++):

cAAST.h
cexpAAST.h
cFitness.h
cMatrix.h
cPRNGXS.h
GenProg.h
myfun.h
prepareData.h

** GP Application (C++):

CASM_GP.cpp

** Support to plots (Matlab):

./graph/GP_graphs.m
./TGE/plotGMFitness_1D.m
./TGE/plotGMFitness_2D.m
./TGE/u_hat.m

** Testing the LPDE solver (Matlab):

solveLPDE.m

Hierarchical order is: 
*myfun* -> *prepareData* -> *cAAST/cexpAAST* -> *cMatrix* -> *cPRNGXS* -> *cFitness* -> *GenProg* -> Application (.cpp). 

To develop a Genetic Programming application to be used as a CASM tool, is enough to just include *GenProg.h*. Note that classes *myfun*, *cAAST*, *cMatrix* and *cPRNGXS* could be used independently for purposes other than the original ones.

To compile the C++ code, download all related files, go to the proper folder and use the command:

>> g++ -std=c++0x -o EXECNAME FILENAME.CPP -llapacke -lblas

Note that LAPACK and BLAS must be installed in the host system.

To run examples, use following datasets:

** To run GP, C++ examples:

data1D_poissonelectro30pts.txt
data1D_underdamposcil30pts.txt
data2D_concentration81pts.txt

** To plot information about GP run (GP_graphs):

./graph/GP_3759bestfit_vec.txt
./graph/GP_3759meanfit_vec.txt
./graph/GP_10439bestfit_vec.txt
./graph/GP_10439meanfit_vec.txt

** To plot TGE solutions with Matlab (plotGMFitness_1D):

./TGE/TGE_3476_data1D_poissonelectro30pts.txt
./TGE/TGE_3759_data1D_underdamposcil30pts.txt

** To plot TGE solutions with Matlab (plotGMFitness_2D):

./TGE/TGE_10439_data2D_concentration81pts.txt

Usage (after compile and build of CASM_GP.cpp): 

>> CASM_GP(.exe) {0:set_basic, 1:set_common} {namefile} {polynomial_degree} 
   {max_diff_order} {mcint powerof2} {population} {maxgen}
   
Example, linux: CASM_GP 0 data1D_underdamposcil30pts.txt 4 2 9 100 70

Example, windows: CASM_GP.exe 0 data1D_underdamposcil30pts.txt 4 2 9 100 70

Feel free to use any datasets available for you. Authors will appreciate
if you please could share results and point out possible problems.

(c)2015 - iperetta@ieee.org
