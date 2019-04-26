#pragma once
#ifndef UTILITY_H
#define UTILITY_H



/* [STANDARD LIBRARIES] */
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

/* [ADDITIONAL LIBRARIES] */
#include <chrono>

/* [EIGEN] */
#include <Eigen/SparseCore>

/* [MATLAB] */
#include <MatlabDataArray.hpp>
#include <MatlabEngine.hpp>
#include <tuple>
#include "engine.h"
#include "mex.h"

/* [CUDA] */
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cusparse.h"
#include "cusolverDn.h"
#include "cusolverSp.h"

using namespace std;

/* Data structure for Priority Queue */
struct VertexPair
{
	int		vId;
	double	distance;
	bool	operator> (const VertexPair &ref) const { return this->distance > ref.distance; }
	bool	operator< (const VertexPair &ref) const { return this->distance < ref.distance; }
};

/* Data structure for Triangle Neighborhood (3 neighbors only) */
struct VtoFPair
{
	int		sId; 
	int		vId;
	int		fId; 
	bool	operator< (const VtoFPair &ref) const { return this->sId < ref.sId; }
	bool	operator> (const VtoFPair &ref) const { return this->sId > ref.sId; }
};

struct Edge_VPair
{
	int		id;
	int		v1, v2;
	bool	operator< (const Edge_VPair &ref) const { return this->id < ref.id; }
	bool	operator> (const Edge_VPair &ref) const { return this->id > ref.id; }
};

struct FacePair
{
	int		id;
	int		f1;
	bool	operator< (const FacePair &ref) const { return this->id < ref.id; }
	bool	operator> (const FacePair &ref) const { return this->id > ref.id; }
};

double SparseMatrixMaxValue(const Eigen::SparseMatrix<double> &M);
void WriteDenseMatrixToMatlab(const Eigen::MatrixXd& M, const string& filename);
void WriteSparseMatrixToMatlab(const Eigen::SparseMatrix<double>& M, const string& filename);
void ReadDenseMatrixFromMatlab(Eigen::MatrixXd& M, const string& filename);
void ReadSparseMatrixFromMatlab(Eigen::SparseMatrix<double>& M, const string& filename);
void ReadVectorFromMatlab(Eigen::VectorXd& v, const string& filename);
void visualizeSparseMatrixInMatlab(const Eigen::SparseMatrix<double> &M);
void ReadChristopherStiffnessMatrix(const string &filename, Eigen::SparseMatrix<double> &M);
double LoadSparseMatrixFromTxtFile(const string& filename, Eigen::SparseMatrix<double> &M);
double ConstructInverseMassMatrix(Eigen::SparseMatrix<double> &M, Eigen::SparseMatrix<double> &MInv);
void WriteSTDVectorToTxtFile(const vector<int>& vector, const string& filename);
void LoadSTDVectorFromTxtFile(const string& filename, vector<int>& vector);
void WriteEigenVectorToTxtFile(const Eigen::VectorXd& vector, const string& filename);
void LoadEigenVectorFromTxtFile(const string& filename, Eigen::VectorXd&);


void writeEigenSparseMatrixToBinary(Eigen::SparseMatrix<double> &m, const std::string &filename);
void readEigenSparseMatrixFromBinary(const std::string &filename, Eigen::SparseMatrix<double> &m);


//template<typename Scalar>
void manuallyDestroySparseMatrix(Eigen::SparseMatrix<double> &M);


//#include "MatlabDataArray.hpp"
//#include "MatlabEngine.hpp"
//#include <iostream>

void evalSurfaceGraph();

#endif // !UTILITY_H

