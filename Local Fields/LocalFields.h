#pragma once
#ifndef LOCAL_FIELDS_H
#define LOCAL_FIELDS_H

#include <set>
#include <vector>
#include <queue>

#include <Eigen/Sparse>
#include <Eigen/PardisoSupport>

#include "Utility.h"

using namespace std; 

class LocalFields
{
public:
	LocalFields(const int &sampleID);
	void constructSubdomain(const int &sampleID, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const double &avgEdgeLength, const Eigen::MatrixXi &AdjMF3N);
	void constructSubdomain(const int &sampleID, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const double &avgEdgeLength, const Eigen::MatrixXi &AdjMF3N, const double& distRatio);
	void constructBoundary(const Eigen::MatrixXi& F, const Eigen::MatrixXi &AdjMF3N, const vector<set<int>> &AdjMF2Ring);
	void constructLocalElements(const Eigen::MatrixXi &F);
	void constructMatrixBLocal(const Eigen::SparseMatrix<double>& B2D);
	void constructMatrixBLocal(const Eigen::SparseMatrix<double>& B2D, const vector<set<int>>& AdjMF2Ring);
	void constructMatrixBLocalDirectInsert(const Eigen::SparseMatrix<double>& B2D, const vector<set<int>>& AdjMF2Ring);
	void constructLocalConstraints();
	void setupRHSLocalProblemMapped();
	void setupLHSLocalProblemMapped();
	void solveLocalSystemMappedLDLT(vector<Eigen::Triplet<double>> &BTriplet);
	void measureXF(const Eigen::VectorXd& doubleArea, const Eigen::SparseMatrix<double>& J);

private:
	int id;// , sampleID;
	Eigen::MatrixXd					XfLoc, cLoc, bLoc, gLoc, hLoc;
	Eigen::SparseMatrix<double>		BLoc, ALoc, CLoc, SF2DLoc, SF_Curl, SF_Div;	
	vector<int>						LocalElements, GlobToLocMap;
	Eigen::VectorXd					vEstimateLoc;

public:
	int								sampleID;
	set<int>						SubDomain, Boundary, BeyondBoundary;

};


#endif // !LOCAL_FIELDS_H

