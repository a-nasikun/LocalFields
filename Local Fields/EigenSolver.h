#pragma once
#ifndef EIGEN_SOLVER_H
#define EIGEN_SOLVER_H

#include "Utility.h"

// For CUDA
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cusolverDn.h"
#include "cublas_v2.h"

// For MATLAB
#include <MatlabDataArray.hpp>
#include <MatlabEngine.hpp>
#include <tuple>
#include "engine.h"
#include "mex.h"

/* * How to compile (assume cuda is installed at /usr/local/cuda/) * nvcc -c -I/usr/local/cuda/include sygvd_example.cpp * g++ -o -fopenmp a.out sygvd_example.o -L/usr/local/cuda/lib64 -lcublas -lcusolver * */
#include <stdio.h> 
#include <stdlib.h> 
#include <assert.h> 
#include <time.h>
#include <iostream>


using namespace std;

/* Computing Eigenstructure in GPU */
void computeEigenGPU(cusolverDnHandle_t& cusolverH, Eigen::SparseMatrix<double> &S_, Eigen::SparseMatrix<double> &M_, Eigen::MatrixXd &LDEigVec, Eigen::VectorXd &LDEigVal);
void computeEigenGPU(Eigen::SparseMatrix<double> &S_, Eigen::SparseMatrix<double> &M_, Eigen::MatrixXd &LDEigVec, Eigen::VectorXd &LDEigVal);
void computeEigenGPU(Eigen::MatrixXd &S_, Eigen::MatrixXd &M_, Eigen::MatrixXd &LDEigVec, Eigen::VectorXd &LDEigVal);

/* Computing Eigenstructure in Matlab */
void computeEigenMatlab(Eigen::SparseMatrix<double> &S, Eigen::SparseMatrix<double> &M, Eigen::MatrixXd &EigVec, Eigen::VectorXd &EigVal);
void computeEigenMatlab(Eigen::SparseMatrix<double> &S, Eigen::SparseMatrix<double> &M, const int& numEigs, Eigen::MatrixXd &EigVec, Eigen::VectorXd &EigVal, const string& filename);
void computeEigenMatlab(Eigen::SparseMatrix<double> &S, const int& numEigs, Eigen::MatrixXd &EigVec, Eigen::VectorXd &EigVal, const string& filename);
void computeEigenMatlab(Engine*& ep, const int tid, Eigen::SparseMatrix<double> &S, Eigen::SparseMatrix<double> &M, const int& numEigs, Eigen::MatrixXd &EigVec, Eigen::VectorXd &EigVal, const string& filename);

/* Computing eigenstructure of 2x2 matrix explicitly/analytically */
void computeEigenExplicit(const Eigen::Matrix2d& M, Eigen::Vector2d& EigVal, Eigen::Matrix2d& EigVect);

///* Computing Eigenstructure in Spectra */
void computeEigenSpectra(Eigen::SparseMatrix<double> &S, Eigen::SparseMatrix<double> &M, const int& numEigs, Eigen::MatrixXd &EigVec, Eigen::VectorXd &EigVal, const string& filename);
void computeEigenSpectra_GenSym(Eigen::SparseMatrix<double> &S, Eigen::SparseMatrix<double> &M, const int& numEigs, Eigen::MatrixXd &EigVec, Eigen::VectorXd &EigVal, const string& filename);
void computeEigenSpectra_RegNSym(Eigen::SparseMatrix<double> &S, Eigen::SparseMatrix<double> &Minv, const int& numEigs, Eigen::MatrixXd &EigVec, Eigen::VectorXd &EigVal, const string& filename);
///
///void testViennaCL2();
///void testViennaCL2(const Eigen::SparseMatrix<double> &S, const Eigen::SparseMatrix<double> &Minv, Eigen::MatrixXd &EigVects, Eigen::VectorXd &EigVals);

#endif // !EIGEN_SOLVER_H

