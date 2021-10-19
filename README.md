# IQSLEv1
IQSLEv1 released

Current version is based on x86_64 linux system 

Dependencies:
1. gcc && gfortran
2. intel mkl BLAS (Recommended) or OpenBLAS or native BLAS (not Recommended)
3. SuiteSparse (need CHOLMOD module)

Install:
If intel mkl is available
```
make MKL_LIB=dir_path_contain_your_mkl_lib MKL_IOMP5=dir_path_contain_your_mkl_iomp5
```
Else If OpenBLAS/other BLAS is available
```
make BLAS_LIB=dir_path_contain_your_openblas
```
Else Use OpenBLAS in the repo
```
make
```
