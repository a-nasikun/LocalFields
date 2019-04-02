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
	void computeDijkstraFaceDistance(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &FC, const Eigen::MatrixXi &AdjMF3N);
	void constructMatrixBLocal(const Eigen::SparseMatrix<double>& B2D);
	//void constructMatrixBLocal(const Eigen::SparseMatrix<double>& B2D, const vector<set<int>>& AdjMF2Ring);
	void constructMatrixBLocal(const Eigen::SparseMatrix<double>& B2D, const vector<set<int>>& AdjMF2Ring, vector<Eigen::Triplet<double>>& BTriplet);
	void constructMatrixBLocalDirectInsert(const Eigen::SparseMatrix<double>& B2D, const vector<set<int>>& AdjMF2Ring);
	void constructLocalConstraints(vector<Eigen::Triplet<double>>& C1Triplet, vector<Eigen::Triplet<double>>& C2Triplet);	
	void constructLocalConstraintsWithLaplacian(const Eigen::VectorXd& doubleArea, const vector<set<int>>& AdjMF2Ring, const Eigen::SparseMatrix<double>& SF2D, vector<Eigen::Triplet<double>>& C1Triplet, vector<Eigen::Triplet<double>>& C2Triplet);
	void constructLocalConstraintsWithLaplacian(const Eigen::VectorXd& doubleArea, const Eigen::SparseMatrix<double>& SF2D, vector<Eigen::Triplet<double>>& C1Triplet, vector<Eigen::Triplet<double>>& C2Triplet);
	//void constructLocalConstraintsWithLaplacian(const Eigen::VectorXd& doubleArea, const Eigen::SparseMatrix<double>& SF2D, vector<Eigen::Triplet<double>>& C1Triplet, vector<Eigen::Triplet<double>>& C2Triplet);
	void constructLocalEigenProblem(const Eigen::SparseMatrix<double>& SF2D, const vector<set<int>>& AdjMF2Ring, Eigen::VectorXd& doubleArea, Eigen::MatrixXd &EigLocal);
	void constructLocalEigenProblemWithSelector(const Eigen::SparseMatrix<double>& SF2D, const vector<set<int>>& AdjMF2Ring, Eigen::VectorXd& doubleArea, Eigen::MatrixXd &EigLocal);
	/* For basis construction */
	void constructLocalEigenProblem(const Eigen::SparseMatrix<double>& SF2D, const vector<set<int>>& AdjMF2Ring, Eigen::VectorXd& doubleArea, vector<Eigen::Triplet<double>>& BTriplet);
	void constructLocalEigenProblem(const Eigen::SparseMatrix<double>& SF2D, const Eigen::MatrixXd &AdjMF3N, Eigen::VectorXd& doubleArea, vector<Eigen::Triplet<double>>& BTriplet);
	void constructLocalEigenProblemWithSelector(const Eigen::SparseMatrix<double>& SF2D, const vector<set<int>>& AdjMF2Ring, Eigen::VectorXd& doubleArea, vector<Eigen::Triplet<double>>& BTriplet);
	//void constructLocalConstraints();
	void setupRHSLocalProblemMapped();
	void setupLHSLocalProblemMapped(const vector<Eigen::Triplet<double>>& BTriplet, const vector<Eigen::Triplet<double>>& C1Triplet, const vector<Eigen::Triplet<double>>& C2Triplet);
	void solveLocalSystemMappedLDLT(vector<Eigen::Triplet<double>> &BTriplet);
	void measureXF(const Eigen::VectorXd& doubleArea, const Eigen::SparseMatrix<double>& J);

private:
	int id;// , sampleID;
	Eigen::MatrixXd					XfLoc, cLoc, bLoc, gLoc, hLoc;
	Eigen::SparseMatrix<double>		BLoc, ALoc, CLoc, SF2DLoc, SF_Curl, SF_Div;	
	Eigen::SparseMatrix<double>		SelectorA; 
	vector<int>						LocalElements, GlobToLocMap, GlobToLocInnerMap;// , InnerElements;
	Eigen::VectorXd					vEstimateLoc, dijksFaceDist;

public:
	int								sampleID;
	set<int>						SubDomain, Boundary, BeyondBoundary;
	vector<int>						InnerElements;
	Eigen::VectorXd					dijksFaceDistMapped, scalingFactor;

};


#endif // !LOCAL_FIELDS_H

