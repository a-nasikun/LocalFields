#pragma once
#ifndef VECTOR_FIELDS_H
#define VECTOR_FIELDS_H

#include "Utility.h"
#include "EigenSolver.h"
#include "LocalFields.h"

#include <igl/edges.h>
#include <igl/grad.h>
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/adjacency_matrix.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/opengl/glfw/Viewer.h>

#include <igl/jet.h>
#include <igl/parula.h>

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/PardisoSupport>
//#include <Eigen/PaStiXSupport>

// For Matlab
//#include "engine.h"
//#include "mex.h"

#include "TestSolver.h"

using namespace std;

enum MassMatrixToShow	{MASS_MV, MASS_MVinv, MASS_MF2D, MASS_MF2Dinv, MASS_MF3D, MASS_MF3Dinv};
enum GradientToShow		{GRAD_3D, GRAD_2D};

class VectorFields
{
public:
	// MESH-related Functions
	void readMesh(const string &meshFile);
	void getVF(Eigen::MatrixXd &V, Eigen::MatrixXi &F);
	void computeFaceCenter();

	// SETTING UP UTILITY MATRICES
	void computeAverageEdgeLength();
	void computeFaceNormal();
	void constructFaceAdjacency2RingMatrix();
	void constructFaceAdjacency3NMatrix();
	void constructFaceAdjacencyMatrix_IGL();
	void constructVertexAdjacencyMatrix();
	void constructNeigborRings(const int &idx);
	void computeDijkstraDistanceVertex(const int &source, Eigen::VectorXd &D);
	void computeDijkstraDistanceFace(const int &source, Eigen::VectorXd &D);
	void computeDijkstraDistanceFaceForSampling(const int &source, Eigen::VectorXd &D);
	void computeDijkstraDistanceFaceMultSource(const Eigen::VectorXi &source, Eigen::VectorXd &D);
	void computeEigenLaplace2D();
	void computeEigenLaplace3D();
	void computeEigenstructureGradient3D();
	//void computeEigenstructureGradient2D();
	void constructLaplace2D();
	void constructLaplace3D();
	void computeEdges();
	void constructVFNeighbors();
	void constructVFNeighborsFull();

	// SETTING UP MATRICES
	void constructGlobalMatrices();
	void constructMassMatrices();
	void constructMassMatrixMV();
	void constructMassMatrixMVinv();
	void constructMassMatrixMF2D();
	void constructMassMatrixMF2Dinv();
	void constructMassMatrixMF3D();
	void constructMassMatrixMF3Dinv();
	void constructStiffnessMatrices();
	void constructStiffnessMatrixSV();
	void constructStiffnessMatrixSF2D();
	void constructStiffnessMatrixSF3D();
	void constructStiffnessMatrixCurlPart3D();
	void constructStiffnessMatrixCurlPart3DandCurl4F();
	void constructStiffnessMatrixCurlPart2D();
	void constructStiffnessMatrixCurlPart2D_Direct();
	void constructStiffnessMatrixDivPart3D();
	void constructStiffnessMatrixDivPart3D_Implicit();
	void constructStiffnessMatrixDivPart3D_Explicit();
	void constructStiffnessMatrixDivPart3DandDiv4F_Explicit();
	void constructStiffnessMatrixDivPart2D();
	void constructStiffnessMatrixDivPart2D_Direct();
	void constructSF2DPacked();
	void constructGradient3D();
	void rearrangeGradient3D();
	void constructGradient2D();
	void computeDivergent3D();
	void computeDivergent2D();
	void computeCurl2D();
	void computeCurl3D();

	// Deal with GLOBAL Problem
	void constructRotationMatrix();
	void constructMappingMatrix();
	void constructMatrixB();
	void constructConstraints();
	void construct1CentralConstraint();
	void constructRingConstraints();
	void constructSpecifiedConstraints();
	void constructSingularities();
	void constructSpecifiedConstraintsWithSingularities();
	void setupGlobalProblem();
	void setupRHSGlobalProblem();
	void setupRHSGlobalProblemMapped();
	void setupLHSGlobalProblem();
	void setupLHSGlobalProblemMapped();
	void solveGlobalSystem();
	void solveGlobalSystemMappedLDLT();

	// LOCAL SYSTEM
	void constructSamples(const int &n);
	void farthestPointSampling();

	void constructBasis();	
	void gatherBasisElements(const vector<vector<Eigen::Triplet<double>>> &UiTriplet);
	void normalizeBasis();
	void normalizeBasisAbs();

	// GLOBAL SYSTEM BASED ON BASIS
	void setAndSolveUserSystem();
	void setupUserBasis();
	void getUserConstraints();
	void getUserConstraintsRandom();
	void getUserConstraintsGUI(const vector<vector<int>> &selectedFaces);
	void getUserConstraintsSpecified();
	void setupRHSUserProblemMapped();
	void setupLHSUserProblemMapped();
	void solveUserSystemMappedLDLT();
	void mapSolutionToFullRes();
	void obtainUserVectorFields();

	// COMPARING RESULTS
	void measureApproxAccuracyL2Norm();

	void constructSubdomainSingle(const int &source);
	void constructBoundary();
	void constructLocalElements();
	void constructMatrixBLocal();
	void constructLocalConstraints();
	void setupRHSLocalProblem();
	void setupLHSLocalProblem();
	void solveLocalSystem();

	// ITEMS FOR TESTING ONLY
	void constructArbitraryField();
	void computeGradArbField3D();
	void computeGradArbField2D();
	void computeCoGradArbField2D();
	void computeCoGradArbField3D();
	void computeCurlGradArbField3D();
	void computeCurlGradArbField2D();
	void computeCurlCoGradArbField3D();
	void computeCurlCoGradArbField2D();
	void computeDivGradArbField3D();
	void computeDivGradArbField2D();
	void computeDivCoGradArbField3D();
	void computeDivCoGradArbField2D();
	void testMappingMatrix();
	void testAdjMV();
	void testAdjacencyAndEdges();
	void testDijkstraFace();
	void testCurlEnergy();
	int selectRandomFace();
	void checkB2DStructure();
	
	// VISUALIZATION of TESTING
	//void visualizeSparseMatrixInMatlab(const Eigen::SparseMatrix<double> &M);
	void visualizeGradient3DArbField(igl::opengl::glfw::Viewer &viewer);
	void visualizeGradient2DArbField(igl::opengl::glfw::Viewer &viewer);
	void visualizeCoGrad3DArbField(igl::opengl::glfw::Viewer &viewer);
	void visualizeCoGrad2DArbField(igl::opengl::glfw::Viewer &viewer);
	void visualizeDivGrad3DArbField(igl::opengl::glfw::Viewer &viewer);
	void visualizeDivCoGrad2DArbField(igl::opengl::glfw::Viewer &viewer);
	void visualizeDivCoGrad3DArbField(igl::opengl::glfw::Viewer &viewer);
	void visualizeDivGrad2DArbField(igl::opengl::glfw::Viewer &viewer);
	void visualizeCurlGrad3DArbField(igl::opengl::glfw::Viewer &viewer);
	void visualizeCurlGrad2DArbField(igl::opengl::glfw::Viewer &viewer);
	void visualizeCurlCoGrad3DArbField(igl::opengl::glfw::Viewer &viewer);
	void visualizeCurlCoGrad2DArbField(igl::opengl::glfw::Viewer &viewer);
	void visualizeLaplaceGrad2DArbField(igl::opengl::glfw::Viewer &viewer);
	void visualizeLaplaceGrad3DArbField(igl::opengl::glfw::Viewer &viewer);
	void visualizeFaceNeighbors(igl::opengl::glfw::Viewer &viewer, const int &idx);
	void visualizeVertexFacesNeighbors(igl::opengl::glfw::Viewer &viewer, const int &idx);
	void visualizeNeighboringRings(igl::opengl::glfw::Viewer &viewer);
	void visualizeDijkstra(igl::opengl::glfw::Viewer &viewer);
	void visualizeEigenfields(igl::opengl::glfw::Viewer &viewer);
	void visualizeArbField(igl::opengl::glfw::Viewer &viewer);
	void visualizeEigFieldsDiv(igl::opengl::glfw::Viewer &viewer, const int &eigID);
	void visualizeRandomFace(igl::opengl::glfw::Viewer &viewer, const int &faceID);
	void visualizeDijkstraFace(igl::opengl::glfw::Viewer &viewer);
	void visualizeSubdomain(igl::opengl::glfw::Viewer &viewer);
	void visualizeSamples(igl::opengl::glfw::Viewer &viewer);

	// VISUALIZATION of IMPORTANT ELEMENTS
	void visualizeMassMatrix(igl::opengl::glfw::Viewer &viewer, const MassMatrixToShow &type);
	void visualizeGradient(igl::opengl::glfw::Viewer &viewer, const GradientToShow &type);
	void visualizeLocalFrames(igl::opengl::glfw::Viewer &viewer);
	void visualizeApproximatedFields(igl::opengl::glfw::Viewer &viewer);
	void visualize2Dfields(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color);
	void visualize3Dfields(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field3D, const Eigen::RowVector3d &color);
	void visualize2DfieldsNormalized(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color);
	void visualize2DfieldsScaled(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color);
	void visualize2DfieldsRegular(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color);
	void visualizeBasis(igl::opengl::glfw::Viewer &viewer, const int &id);
	void visualizeBasisNormalized(igl::opengl::glfw::Viewer &viewer, const int &id);
	void visualizeBasisSum(igl::opengl::glfw::Viewer &viewer, const int &id);
	void visualizeApproxResult(igl::opengl::glfw::Viewer &viewer, const int &id);
	void visualizeUserConstraints(igl::opengl::glfw::Viewer &viewer);
	void visualizeGlobalConstraints(igl::opengl::glfw::Viewer &viewer);
	void visualizeSingularitiesConstraints(igl::opengl::glfw::Viewer &viewer);

protected:
	Eigen::MatrixXd					V, FC, NF, EigVects, EigVectsDiv, EigFieldsDiv3D, EigFieldsDiv2D, Xf, c, cBar, b, bBar, g, h;
	Eigen::MatrixXd					XfLoc, XfLocSum, XfNorm, BasisSum, BasisSumN, cLoc, bLoc, XLowDim, XFullDim, Curl3DPacked, Div3DPacked, SF2DPacked;
	Eigen::MatrixXi					F, E, AdjMF3N, EdgePairMatrix;
	Eigen::SparseMatrix<double>		MV, MVinv, MF2D, MF2Dinv, MF3D, MF3Dinv, BasisTemp, Basis, SV, SF2D, SF3D, L2D, L3D, B2D, B2Dbar, B3D, LapCurl3D, LapCurl2D, LapDiv3D, LapDiv2D;
	Eigen::SparseMatrix<double>		BLoc, ALoc, CLoc; 
	Eigen::SparseMatrix<double>		GF3D, GF2D, Div3D, Div2D, Curl3D, Curl2D, A, J, C, Cbar, A_LHSbar, A_LHS, AjdMF_igl;
	Eigen::VectorXd					doubleArea, arbField, gradArbField3D, gradArbField2D, eigVals, eigValsDiv, lambda, vFields;
	Eigen::VectorXd					curlGradArbField3D, curlGradArbField2D, curlCoGradArbField3D, curlCoGradArbField2D;
	Eigen::VectorXd					divGradArbField3D, divGradArbField2D, coGradArbField3D, coGradArbField2D;
	Eigen::VectorXd					divCoGradArbField3D, divCoGradArbField2D, vEstUser, vEst, gbar, hbar,pbar;
	vector<set<int>>				AdjMV, AdjMF2Ring, NeighRing/*, AdjMF3N_temp*/;
	vector<set<VtoFPair>>			VFNeighbors, VFNeighFull;
	vector<set<Edge_VPair>>			EdgePairsList; 
	vector<set<FacePair>>			AdjMF3N_temp;
	vector<int>						LocalElements, LocToGlobMap, GlobToLocMap, Sample, userConstraints, globalConstraints;
	vector<int>						neighCCW, singularities;
	vector<chrono::duration<double>>durations;
	vector<vector<int>>				SingNeighCC;
	set<int>						SubDomain, Boundary; 
	int								sample, numSample; 
	

	// FOR TESTING ONLY
	Eigen::VectorXd					laplaceGradArbField2D, laplaceGradArbField3D, dijkstraFace;

private:
	double avgEdgeLength; 
};


#endif // !VECTOR_FIELDS_H

