#pragma once
#ifndef VECTOR_FIELDS_H
#define VECTOR_FIELDS_H

#define EIGEN_USE_MKL_ALL

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
	void constructVertexAdjacencyMatrix();
	void computeDijkstraDistanceVertex(const int &source, Eigen::VectorXd &D);
	void computeDijkstraDistanceFace(const int &source, Eigen::VectorXd &D);
	void computeDijkstraDistanceFaceForSampling(const int &source, Eigen::VectorXd &D);
	void computeDijkstraDistanceFaceMultSource(const Eigen::VectorXi &source, Eigen::VectorXd &D);
	void computeEdges();
	void constructVFNeighbors();
	void constructVFNeighborsFull();
	void constructVFAdjacency();
	void testAdjacency();

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
	void constructStiffnessMatrixSF2D();
	void constructStiffnessMatrixSF3D();
	void constructStiffnessMatrixCurlPart3D();
	void constructStiffnessMatrixCurlPart3DandCurl4F();
	void constructStiffnessMatrixCurlPart2D();
	void constructStiffnessMatrixDivPart3D();
	void constructStiffnessMatrixDivPart3D_Implicit();
	void constructStiffnessMatrixDivPart3D_Explicit();
	void constructStiffnessMatrixDivPart3DandDiv4F_Explicit();
	void constructStiffnessMatrixDivPart2D();
	void constructGradient3D();
	void rearrangeGradient3D();
	void constructGradient2D();
	void computeDivergent3D();
	void computeDivergent2D();
	void computeCurl3D();
	void computeCurl2D();

	// Deal with GLOBAL Problem
	void constructRotationMatrix();
	void constructMappingMatrix();
	void constructMatrixB();

	void setupGlobalProblem();
	void constructConstraints();
	void construct1CentralConstraint();
	void constructRingConstraints();
	void constructSpecifiedConstraints();
	void constructSingularities();
	void constructSpecifiedConstraintsWithSingularities();
	void setupRHSGlobalProblemMapped();
	void setupLHSGlobalProblemMapped();
	void solveGlobalSystemMappedLDLT();
	void solveGlobalSystemMappedLU_GPU();

	// LOCAL SYSTEM
	void constructSamples(const int &n);
	void farthestPointSampling();
	void constructBasis();	
	void gatherBasisElements(const vector<vector<Eigen::Triplet<double>>> &UiTriplet);
	void normalizeBasis();
	void normalizeBasisAbs();

	// REDUCED-GLOBAL SYSTEM BASED ON BASIS
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
	void measureU1andJU0();

	// ITEMS FOR TESTING ONLY
	void constructArbitraryField();
	void testMappingMatrix();
	void testAdjMV();
	void testAdjacencyAndEdges();
	void testDijkstraFace();
	void testCurlEnergy();
	int selectRandomFace();
	void checkB2DStructure();
	
	// VISUALIZATION of TESTING
	//void visualizeSparseMatrixInMatlab(const Eigen::SparseMatrix<double> &M);	
	void visualizeFaceNeighbors(igl::opengl::glfw::Viewer &viewer, const int &idx);
	void visualizeVertexFacesNeighbors(igl::opengl::glfw::Viewer &viewer, const int &idx);
	void visualizeNeighboringRings(igl::opengl::glfw::Viewer &viewer);
	void visualizeDijkstra(igl::opengl::glfw::Viewer &viewer);
	void visualizeEigenfields(igl::opengl::glfw::Viewer &viewer);
	void visualizeArbField(igl::opengl::glfw::Viewer &viewer);
	void visualizeRandomFace(igl::opengl::glfw::Viewer &viewer, const int &faceID);
	void visualizeDijkstraFace(igl::opengl::glfw::Viewer &viewer);
	void visualizeSubdomain(igl::opengl::glfw::Viewer &viewer);
	void visualizeSamples(igl::opengl::glfw::Viewer &viewer);
	void visualizeSharedEdges(igl::opengl::glfw::Viewer &viewer);

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
	// Variable (Matrix, Vector or regular variables) for Matrix construction
	Eigen::MatrixXd					V, FC, NF;
	Eigen::MatrixXd					Curl3DPacked, Div3DPacked;
	Eigen::MatrixXi					F, E, AdjMF3N, EdgePairMatrix;
	Eigen::SparseMatrix<double>		MV, MVinv, MF2D, MF2Dinv, MF3D, MF3Dinv, SF2D, SF3D, B2D, B2Dbar, LapCurl3D, LapCurl2D, LapDiv3D, LapDiv2D;
	Eigen::SparseMatrix<double>		GF3D, GF2D, Div3D, Div2D, Curl3D, Curl2D, A, J;
	Eigen::SparseMatrix<bool>		VFAdjacency;
	Eigen::VectorXd					doubleArea;
	vector<set<int>>				AdjMV, AdjMF2Ring, NeighRing;
	vector<set<VtoFPair>>			VFNeighbors, VFNeighFull;
	vector<set<Edge_VPair>>			EdgePairsList;
	vector<set<FacePair>>			AdjMF3N_temp;
	double							avgEdgeLength;

	// Variable related to global problem
	Eigen::MatrixXd					b, g, h;
	Eigen::SparseMatrix<double>		C, A_LHS;
	Eigen::VectorXd					vEst, Xf, c;
	vector<int>						LocalElements, userConstraints, globalConstraints;
	vector<int>						singularities;
	vector<vector<int>>				SingNeighCC;
	set<int>						SubDomain, Boundary;

	// Variable related to subspace construction
	Eigen::SparseMatrix<double>		BasisTemp, Basis;
	Eigen::MatrixXd					BasisSum, BasisSumN;
	vector<int>						Sample;
	vector<chrono::duration<double>>durations;
	int								numSample;

	// Variable related to manipulation within the subspace
	Eigen::MatrixXd					cBar, bBar;
	Eigen::MatrixXd					XLowDim, XFullDim;
	Eigen::SparseMatrix<double>		Cbar, A_LHSbar;
	Eigen::VectorXd					vEstUser, gbar, hbar, pbar;	

	// FOR TESTING ONLY
	Eigen::VectorXd					dijkstraFace, arbField;
	vector<vector<int>>				sharedEdgesVect; 

private:
	
};


#endif // !VECTOR_FIELDS_H

