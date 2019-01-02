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
//#include <MatlabDataArray.hpp>
//#include <MatlabEngine.hpp>
//#include <tuple>
//#include "engine.h"
//#include "mex.h"

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

#endif // !UTILITY_H

