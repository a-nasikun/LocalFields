#pragma once
#ifndef TEST_SOLVER_H
#define TEST_SOLVER_H


// For MKL
#include "mkl_lapacke.h"
#include "mkl_pardiso.h"
#include "mkl_types.h"

// For CUDA
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>

#include <Eigen/Sparse>

/* ================================== MKL ==============================================*/
/* Auxiliary routine: printing a matrix */
void print_matrix(char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda);
/* Auxiliary routine: printing a vector of integers */
void print_int_vector(char* desc, MKL_INT n, MKL_INT* a);
/* Testing LAPACKE LDLT Solver */
void testLAPACKEdsysv();

/* ================================== MKL PARDISO =======================================*/
void testMKL_Pardiso();
void solveLDLT_MKLPardiso(const Eigen::SparseMatrix<double> &A, const Eigen::MatrixXd& B, Eigen::MatrixXd& Xf);
void solveLDLT_EigenPardiso(const Eigen::SparseMatrix<double> &A, const Eigen::MatrixXd& B, Eigen::MatrixXd& Xf);

/* ================================== CUDA ==============================================*/
void testCUDA_LULinearSolver();
void printMatrix(int m, int n, const double*A, int lda, const char* name);
void testLDLTLinearSolver();
void solveLDLTinCUDA(const Eigen::SparseMatrix<double> &A, const Eigen::MatrixXd& B, Eigen::MatrixXd& Xf);

#endif // !TEST_SOLVER_H


