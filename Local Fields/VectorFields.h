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
#include <igl/principal_curvature.h>
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
	void constructStiffnessMatrixSF2D(Eigen::SparseMatrix<double>& LapCurl3D, Eigen::SparseMatrix<double>& LapCurl2D, Eigen::SparseMatrix<double>& LapDiv3D, Eigen::SparseMatrix<double>& LapDiv2D);
	void constructStiffnessMatrixSF3D(Eigen::SparseMatrix<double>& LapCurl3D, Eigen::SparseMatrix<double>& LapDiv3D);
	void constructStiffnessMatrixCurlPart3D(Eigen::SparseMatrix<double>& LapCurl3D);
	void constructStiffnessMatrixCurlPart2D(Eigen::SparseMatrix<double>& LapCurl3D, Eigen::SparseMatrix<double>& LapCurl2D);
	void constructStiffnessMatrixDivPart3D(Eigen::SparseMatrix<double>& LapDiv3D);
	void constructStiffnessMatrixDivPart3D_Implicit(Eigen::SparseMatrix<double>& LapDiv3D);
	void constructStiffnessMatrixDivPart3D_Explicit(Eigen::SparseMatrix<double>& LapDiv3D);
	void constructStiffnessMatrixDivPart2D(Eigen::SparseMatrix<double>& LapDiv3D, Eigen::SparseMatrix<double>& LapDiv2D);
	void constructGradient3D();
	void rearrangeGradient3D();
	void rearrangeGradient3D(Eigen::SparseMatrix<double>& Grad3D);
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
	void pushNewUserConstraints(const int& fInit, const int& fEnd);
	void constructSpecifiedHardConstraints();
	void constructInteractiveConstraints();
	void constructSingularities();
	void constructHardConstraintsWithSingularities();
	void constructHardConstraintsWithSingularities_Cheat();
	void constructHardConstraintsWithSingularitiesWithGauss();
	void constructSoftConstraints();
	void constructCurvesAsConstraints(const int& init, const int& end, vector<int>& curve);
	void projectCurvesToFrame();
	void setupRHSGlobalProblemMapped(Eigen::VectorXd& g, Eigen::VectorXd& h, Eigen::VectorXd& vEst, Eigen::VectorXd& b);
	void setupLHSGlobalProblemMapped(Eigen::SparseMatrix<double>& A_LHS);
	void setupRHSGlobalProblemSoftConstraints(const double& lambda, Eigen::VectorXd& b);
	void setupLHSGlobalProblemSoftConstraints(const double& lambda, Eigen::SparseMatrix<double>& A_LHS);
	void solveGlobalSystemMappedLDLT(Eigen::VectorXd& vEst, Eigen::SparseMatrix<double>& A_LHS, Eigen::VectorXd& b);
	void solveGlobalSystemMappedLU_GPU(Eigen::VectorXd& vEst, Eigen::SparseMatrix<double>& A_LHS, Eigen::VectorXd& b);
	void solveGlobalSystemMappedLDLTSoftConstraints(Eigen::SparseMatrix<double>& A_LHS, Eigen::VectorXd& b);

	// APPLICATIONS ON GLOBAL SYSTEM
	void computeSmoothing(const double& mu, const Eigen::VectorXd& v_in, Eigen::VectorXd& v_out);

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
	void setupRHSUserProblemMapped(Eigen::VectorXd& gBar, Eigen::VectorXd& hBar, Eigen::VectorXd& vEstBar, Eigen::VectorXd& bBar);
	void setupLHSUserProblemMapped(Eigen::SparseMatrix<double>& A_LHSBar);
	void setupRHSUserProblemMappedSoftConstraints(const double& lambda, Eigen::VectorXd& bBar);
	void setupLHSUserProblemMappedSoftConstraints(const double& lambda, Eigen::SparseMatrix<double>& A_LHSBar);
	void solveUserSystemMappedLDLT(Eigen::VectorXd& vEstBar, Eigen::SparseMatrix<double>& A_LHSBar, Eigen::VectorXd& bBar);
	void solveUserSystemMappedLDLTSoftConstraints(Eigen::SparseMatrix<double>& A_LHSBar, Eigen::VectorXd& bBar);
	void mapSolutionToFullRes();
	void obtainUserVectorFields();

	// COMPARING RESULTS
	void measureApproxAccuracyL2Norm();
	void measureDirichletEnergy();
	void testBasis();
	void measureU1andJU0();

	// APPLICATIONS ON REDUCED BASIS
	void computeSmoothingApprox(const double& mu, const Eigen::VectorXd& v_in, Eigen::VectorXd& v_out);

	void ConstructCurvatureTensor(igl::opengl::glfw::Viewer &viewer);
	void ComputeCurvatureFields();

	// ITEMS FOR TESTING ONLY
	void constructArbitraryField();
	void constructArbitraryField2D();
	void testMappingMatrix();
	void testAdjMV();
	void testAdjacencyAndEdges();
	void testDijkstraFace();
	void testCurlEnergy();
	int selectRandomFace();
	void checkB2DStructure();
	void testEdgesAddition(igl::opengl::glfw::Viewer &viewer);
	
	// VISUALIZATION of TESTING
	//void visualizeSparseMatrixInMatlab(const Eigen::SparseMatrix<double> &M);	
	void visualizeFaceNeighbors(igl::opengl::glfw::Viewer &viewer, const int &idx);
	void visualizeVertexFacesNeighbors(igl::opengl::glfw::Viewer &viewer, const int &idx);
	//void visualizeNeighboringRings(igl::opengl::glfw::Viewer &viewer);
	void visualizeDijkstra(igl::opengl::glfw::Viewer &viewer);
	//void visualizeEigenfields(igl::opengl::glfw::Viewer &viewer);
	void visualizeArbField(igl::opengl::glfw::Viewer &viewer);
	void visualizeRandomFace(igl::opengl::glfw::Viewer &viewer, const int &faceID);
	void visualizeDijkstraFace(igl::opengl::glfw::Viewer &viewer);
	void visualizeSubdomain(igl::opengl::glfw::Viewer &viewer);
	void visualizeSamples(igl::opengl::glfw::Viewer &viewer);
	void visualizeSharedEdges(igl::opengl::glfw::Viewer &viewer);
	void visualizeCurveConstraints(igl::opengl::glfw::Viewer &viewer);
	void visualizeSoftConstraints(igl::opengl::glfw::Viewer &viewer);
	void visualize1FieldOnCenter(igl::opengl::glfw::Viewer &viewer, const bool& even);	

	// VISUALIZATION of IMPORTANT ELEMENTS
	void selectFaceToDraw(const int& numFaces);
	void visualizeMassMatrix(igl::opengl::glfw::Viewer &viewer, const MassMatrixToShow &type);
	void visualizeGradient(igl::opengl::glfw::Viewer &viewer, const GradientToShow &type);
	void visualizeLocalFrames(igl::opengl::glfw::Viewer &viewer);
	void visualizeApproximatedFields(igl::opengl::glfw::Viewer &viewer);
	void visualize2Dfields(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color, const double& scale, const bool& normalized = false);
	void visualize3Dfields(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field3D, const Eigen::RowVector3d &color);
	void visualize2DfieldsNormalized(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color, const int &numFaces);
	void visualize2DfieldsScaled(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color, const double &scale);
	void visualize2DfieldsRegular(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color);
	void visualize2DfieldsScaled(igl::opengl::glfw::Viewer &viewer, const Eigen::SparseMatrix<double> &Field2D, const int &idx, const Eigen::RowVector3d &color);
	void visualizeBasis(igl::opengl::glfw::Viewer &viewer, const int &id);
	void visualizeBasisNormalized(igl::opengl::glfw::Viewer &viewer, const int &id);
	void visualizeBasisSum(igl::opengl::glfw::Viewer &viewer, const int &id);
	void visualizeApproxResult(igl::opengl::glfw::Viewer &viewer);
	void visualizeUserConstraints(igl::opengl::glfw::Viewer &viewer);
	void visualizeGlobalConstraints(igl::opengl::glfw::Viewer &viewer);
	void visualizeSingularitiesConstraints(igl::opengl::glfw::Viewer &viewer);

	// VISUALIZATION of APPLICATIONs
	void visualizeSmoothing(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd& v);
	void visualizeCurvatureTensor(igl::opengl::glfw::Viewer &viewer);

protected:
	// Variable (Matrix, Vector or regular variables) for Matrix construction
	//Eigen::MatrixXd					V, FC, NF;
	Eigen::MatrixXd					V, NF;
	Eigen::MatrixXi					F, E, AdjMF3N, EdgePairMatrix;
	Eigen::SparseMatrix<double>		MV, MVinv, MF2D, MF2Dinv, MF3D, MF3Dinv, SF2D, SF3D, B2D;
	Eigen::SparseMatrix<double>		GF3D, GF2D, Div3D, Div2D, Curl3D, Curl2D, A, J;
	Eigen::SparseMatrix<bool>		VFAdjacency;
	Eigen::VectorXd					doubleArea;
	vector<set<int>>				AdjMV, AdjMF2Ring, NeighRing;
	vector<set<VtoFPair>>			VFNeighbors, VFNeighFull;
	vector<set<Edge_VPair>>			EdgePairsList;
	vector<set<FacePair>>			AdjMF3N_temp;
	double							avgEdgeLength;
	vector<int>						FaceToDraw;

	// Variable related to global problem
	Eigen::SparseMatrix<double>		C; 
	Eigen::VectorXd					c, Xf;
	vector<int>						LocalElements, userConstraints, globalConstraints, userVisualConstraints;
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
	Eigen::MatrixXd					cBar;
	Eigen::VectorXd					XLowDim, XFullDim;
	Eigen::SparseMatrix<double>		CBar, B2DBar;

	// Variables related to Applications
	Eigen::MatrixXd					CurvatureTensorField2D;
	Eigen::SparseMatrix<double>		CurvatureTensor2D; 

	// FOR TESTING ONLY
public: 
	Eigen::VectorXd					dijkstraFace, arbField, arbField2D, wb;
	Eigen::VectorXd					sampleDistance;
	vector<vector<int>>				sharedEdgesVect, curvesConstraints; 
	vector<vector<Eigen::Vector2d>>	constraintVect2D;
	Eigen::MatrixXd					FC;
private:
	
};


#endif // !VECTOR_FIELDS_H

