# IQSLEv1
IQSLEv1 released

Dependencies:
1. gcc && gfortran
2. OpenBlAS or intel mkl BLAS (Recommended)
3. SuiteSparse (need CHOLMOD modular)

Install:
1. If intel mkl is available
  make MKL_LIB=dir_path_contain_your_mkl_lib MKL_IOMP5=dir_path_contain_your_mkl_iomp5

2. If OpenBLAS is available
  make BLAS_LIB=dir_path_contain_your_openblas

3. Use OpenBLAS in the repo
  make
