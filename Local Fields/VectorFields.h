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
#include "engine.h"
#include "mex.h"

#include "TestSolver.h"

using namespace std;

enum MassMatrixToShow	{MASS_MV, MASS_MVinv, MASS_MF2D, MASS_MF2Dinv, MASS_MF3D, MASS_MF3Dinv};
enum GradientToShow		{GRAD_3D, GRAD_2D};

class VectorFields
{
public:
	// MESH-related Functions
	void readMesh(const string &meshFile);
	void scaleMesh();
	void readArrowMesh(const string &meshFile);
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
	void computeDijkstraForParallelTransport(const int &source, const int& target);
	void computeEdges();
	void constructVFNeighbors();
	void constructVFNeighborsFull();
	void constructVFAdjacency();
	void testAdjacency();
	void constructEVList();
	void constructEFList();

	// SETTING UP MATRICES
	void constructGlobalMatrices();
	void constructMassMatrices();
	void constructMassMatrixMV();
	void constructMassMatrixMVinv();
	void constructMassMatrixMF2D();
	void constructMassMatrixMF2Dinv();
	void constructMassMatrixMStarAndInv();
	void constructMassMatrixMF3D();
	void constructMassMatrixMF3Dinv();
	void constructStiffnessMatrices();
	void constructStiffnessMatrices_Implicit();
	void loadStiffnessMatrices();
	void constructStiffnessMatrixSF2D(Eigen::SparseMatrix<double>& Matrix3D, Eigen::SparseMatrix<double>& Matrix2D);
	void constructStiffnessMatrixSF2D(Eigen::SparseMatrix<double>& LapCurl3D, Eigen::SparseMatrix<double>& LapCurl2D, Eigen::SparseMatrix<double>& LapDiv3D, Eigen::SparseMatrix<double>& LapDiv2D);
	//void constructStiffnessMatrixSF3D(Eigen::SparseMatrix<double>& LapCurl3D, Eigen::SparseMatrix<double>& LapDiv3D);
	//void constructStiffnessMatrixCurlPart3D(Eigen::SparseMatrix<double>& LapCurl3D);
	//void constructStiffnessMatrixCurlPart2D(Eigen::SparseMatrix<double>& LapCurl3D, Eigen::SparseMatrix<double>& LapCurl2D);
	//void constructStiffnessMatrixDivPart3D(Eigen::SparseMatrix<double>& LapDiv3D);
	//void constructStiffnessMatrixDivPart3DFromCurl3D(Eigen::SparseMatrix<double>& LapCurl3D, Eigen::SparseMatrix<double>& LapDiv3D);
	//void constructStiffnessMatrixDivPart3D_Implicit(Eigen::SparseMatrix<double>& LapDiv3DAsym);
	//void constructStiffnessMatrixDivPart3D_Implicit(Eigen::SparseMatrix<double>& LapDiv3D);
	//void constructStiffnessMatrixDivPart3D_Explicit(Eigen::SparseMatrix<double>& LapDiv3D);
	void constructStiffnessMatrixSF2DAsym(Eigen::SparseMatrix<double>& LapCurl3D, Eigen::SparseMatrix<double>& LapDiv3D);
	void constructStiffnessMatrixSF3D(Eigen::SparseMatrix<double>& LapCurl3D_Conform, Eigen::SparseMatrix<double>& LapDiv3D_Conform, Eigen::SparseMatrix<double>& LapCurl3D_NonConform, Eigen::SparseMatrix<double>& LapDiv3D_NonConform);
	void constructStiffnessMatrixCurlPart3D_Conform(Eigen::SparseMatrix<double>& LapCurl3D_Conform);
	void constructStiffnessMatrixCurlPart3D_NonConform(Eigen::SparseMatrix<double>& LapCurl3D_NonConform);
	void constructStiffnessMatrixCurlPart2D(Eigen::SparseMatrix<double>& LapCurl3D, Eigen::SparseMatrix<double>& LapCurl2D);
	void constructStiffnessMatrixDivPart3D(Eigen::SparseMatrix<double>& LapDiv3D);
	void constructStiffnessMatrixDivPart3DFromCurl3D(Eigen::SparseMatrix<double>& LapCurl3D, Eigen::SparseMatrix<double>& LapDiv3D);
	void constructStiffnessMatrixDivPart3D_Conform(Eigen::SparseMatrix<double>& LapDiv3D);
	void constructStiffnessMatrixDivPart3D_NonConform(Eigen::SparseMatrix<double>& LapDiv3D);
	void constructStiffnessMatrixDivPart2D(Eigen::SparseMatrix<double>& LapDiv3D, Eigen::SparseMatrix<double>& LapDiv2D);
	void constructGradient3D();
	void rearrangeGradient3D();
	void rearrangeGradient3D(Eigen::SparseMatrix<double>& Grad3D);
	void constructGradient2D();
	void constructGradientStar3D();
	void constructGradientStar2D();
	void computeDivergent3D();
	void computeDivergent2D();
	void computeCurl3D();
	void computeCurl2D();

	// Deal with GLOBAL Problem
	void constructRotationMatrix();
	void constructMappingMatrix();
	void constructMatrixB();

	void setupGlobalProblem(const Eigen::Vector3d& lambda);
	void setupGlobalProblem(const Eigen::Vector3d& lambda, Eigen::MatrixXd& M);
	void constructConstraints();
	void construct1CentralConstraint();
	void constructRingConstraints();
	void pushNewUserConstraints(const int& fInit, const int& fEnd);
	void constructSpecifiedHardConstraints();
	void constructRandomHardConstraints();
	void constructInteractiveConstraints();
	void constructInteractiveConstraintsWithLaplacian();
	void resetInteractiveConstraints();
	void constructSingularities();
	void constructHardConstraintsWithSingularities();
	void constructHardConstraintsWithSingularities_Cheat();
	void constructHardConstraintsWithSingularitiesWithGauss();
	void constructSoftConstraints();
	void constructCurvesAsConstraints(const int& init, const int& end, vector<int>& curve);
	void measureSoftConstraintError(const Eigen::Vector3d& lambda);
	void projectCurvesToFrame();
	void setupRHSGlobalProblemMapped(Eigen::VectorXd& g, Eigen::VectorXd& h, Eigen::VectorXd& vEst, Eigen::VectorXd& b);
	void setupLHSGlobalProblemMapped(Eigen::SparseMatrix<double>& A_LHS);
	void setupRHSGlobalProblemSoftConstraints(const Eigen::Vector3d& lambda, Eigen::VectorXd& b);
	void setupLHSGlobalProblemSoftConstraints(const Eigen::Vector3d& lambda, Eigen::SparseMatrix<double>& A_LHS);
	void solveGlobalSystemMappedLDLT(Eigen::VectorXd& vEst, Eigen::SparseMatrix<double>& A_LHS, Eigen::VectorXd& b);
	void solveGlobalSystemMappedLU_GPU(Eigen::VectorXd& vEst, Eigen::SparseMatrix<double>& A_LHS, Eigen::VectorXd& b);
	void solveGlobalSystemMappedLDLTSoftConstraints(Eigen::SparseMatrix<double>& A_LHS, Eigen::VectorXd& b);

	
	// Alignment fields (Maximal curvature direction) 
	void computeMaximalPrincipalCurvature(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::VectorXd &PD, Eigen::VectorXd& PV);

	// APPLICATIONS ON GLOBAL SYSTEM
	void computeSmoothing(const double& mu, const Eigen::VectorXd& v_in, Eigen::VectorXd& v_out);
	void computeEigenFields(const int &numEigs, const string& filename );
	void retrieveEigenFields(const string& filename);

	// LOCAL SYSTEM
	void constructSamples(const int &n);
	void farthestPointSampling();
	void constructMultiBasis();
	void constructBasis();
	void constructBasis_LocalEigenProblem();
	void constructBasis_OptProblem();
	void constructBasis_LocalEigenProblem10();
	void constructBasis_GradOfLocalFunction(Eigen::SparseMatrix<double>& BasisFunctions);
	void constructBasis_EigenPatch(Eigen::SparseMatrix<double>& BasisFunctions);
	void constructBasisEigenVects();
	void gatherBasisElements(const vector<vector<Eigen::Triplet<double>>> &UiTriplet, const int& NUM_EIGEN);
	void loadAndConstructBasis();
	void writeBasisElementsToFile(const vector<vector<Eigen::Triplet<double>>> &UiTriplet, const int& NUM_EIGEN);
	void normalizeBasis();
	void normalizeBasisAbs(const int& stride);
	void storeBasis(const string& filename);
	void retrieveBasis(const string& filename);

	// REDUCED-GLOBAL SYSTEM BASED ON BASIS
	void setAndSolveUserSystem(const Eigen::Vector3d& lambda);
	void setupReducedBiLaplacian();
	void getUserConstraints();
	void setupRHSUserProblemMapped(Eigen::VectorXd& gBar, Eigen::VectorXd& hBar, Eigen::VectorXd& vEstBar, Eigen::VectorXd& bBar);
	void setupLHSUserProblemMapped(Eigen::SparseMatrix<double>& A_LHSBar);
	void setupRHSUserProblemMappedSoftConstraints(const Eigen::Vector3d& lambda, Eigen::VectorXd& bBar);
	void setupLHSUserProblemMappedSoftConstraints(const Eigen::Vector3d& lambda, Eigen::SparseMatrix<double>& A_LHSBar);
	void solveUserSystemMappedLDLT(Eigen::VectorXd& vEstBar, Eigen::SparseMatrix<double>& A_LHSBar, Eigen::VectorXd& bBar);
	void solveUserSystemMappedLDLTSoftConstraints(Eigen::SparseMatrix<double>& A_LHSBar, Eigen::VectorXd& bBar);
	void mapSolutionToFullRes();
	void obtainUserVectorFields();

	// COMPARING RESULTS
	void measureApproxAccuracyL2Norm();
	void measureDirichletEnergy();
	void testProjection_MyBasis_NoRegularizer(const Eigen::SparseMatrix<double>& Basis, const Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> &sparseSolver, const Eigen::SparseMatrix<double>& B, const Eigen::VectorXd& a, const Eigen::VectorXd& inputFields, double &error);
	void testProjection_EigenBasis_NoRegularizer(const Eigen::MatrixXd& Basis, const Eigen::LDLT<Eigen::MatrixXd>& denseSolver, const Eigen::VectorXd& a, const Eigen::VectorXd& inputFields, double &error);
	void testProjection_MyBasis_WithRegularizer(const Eigen::SparseMatrix<double>& Basis, const Eigen::VectorXd& inputFields, const Eigen::SparseMatrix<double>& MReg, double &error);
	void testProjection_EigenBasis_WithRegularizer(const Eigen::MatrixXd& Basis, const Eigen::VectorXd& inputFields, const Eigen::SparseMatrix<double>& MReg, double &error);
	void testProjection_MyBasis_WithRegularizer(const Eigen::SparseMatrix<double>& Basis, const Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> &sparseSolver_Ref, const Eigen::SparseMatrix<double>& B_Ref, const Eigen::VectorXd& a_Ref, const Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> &sparseSolver_Red, const Eigen::SparseMatrix<double>& B_Red, const Eigen::VectorXd& a_Red, const Eigen::VectorXd& inputFields, const Eigen::SparseMatrix<double>& MReg, double &error);
	void testProjection_EigenBasis_WithRegularizer(const Eigen::MatrixXd& Basis, const Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> &sparseSolver_Ref, const Eigen::SparseMatrix<double>& B_Ref, const Eigen::VectorXd& a_Ref, const Eigen::LDLT<Eigen::MatrixXd>& denseSolver_Red, const Eigen::MatrixXd& B_Red, const Eigen::VectorXd& a_Red, const Eigen::VectorXd& inputFields, const Eigen::SparseMatrix<double>& MReg, double &error);
	//void testProjection_MyBasis_WithRegularizer(const Eigen::SparseMatrix<double>& Basis, const Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> &sparseSolver_Ref, const Eigen::SparseMatrix<double>& B_Ref, const Eigen::VectorXd& a_Ref, const Eigen::LDLT<Eigen::MatrixXd>& denseSolver_Red, const Eigen::MatrixXd& B_Red, const Eigen::VectorXd& a_Red, const Eigen::VectorXd& inputFields, const Eigen::SparseMatrix<double>& MReg, double &error);
	void projectionTest();
	void projectionTest(bool &readDesFieldsFromFile, bool &readPertFieldsFromFile, bool &useEigenBasis, int start, int nTests);
	void compareModalBasis_SamePerformance();
	void compareModalBasis_SameStorage();
	void convergenceTest();
	void measureU1andJU0();
	void measureL2NormEigVectors();
	void vectorFieldsDesignTest();
	void vectorFieldsDesignTest_Normalized();

	// APPLICATIONS ON REDUCED BASIS
	void computeSmoothingApprox(const double& mu, const Eigen::VectorXd& v_in, Eigen::VectorXd& v_out);

	void computeApproxEigenFields(const int &numEigs, const string& filename);
	void retrieveApproxEigenFields();
	void ConstructCurvatureTensor(igl::opengl::glfw::Viewer &viewer);
	void ComputeCurvatureFields();

	// ITEMS FOR TESTING ONLY
	void TEST_VECTOR(igl::opengl::glfw::Viewer &viewer, const string& meshFile);
	void constructParallelTransport();
	void writeBasisToFile();
	void writeField3DToFile();
	void printDataForVTK();
	void writeEigenFieldsForVTK();
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
	void testEnergyOfLocalPatch(igl::opengl::glfw::Viewer &viewer);
	void testGradients();
	void testRotation();
	void testMassMatrix();
	void projectionMatrixTest();
	void testCurlAndDiv();
	void perturbVectorFields(Eigen::VectorXd& inputFields);
	void perturbVectorFieldsRegular(Eigen::VectorXd& inputFields);
	void testSpectra();
	void testSparseMatrix();
	
	
	// VISUALIZATION of TESTING
	//void visualizeSparseMatrixInMatlab(const Eigen::SparseMatrix<double> &M);
	void visualizeFaceNeighbors(igl::opengl::glfw::Viewer &viewer, const int &idx);
	void visualizeVertexFacesNeighbors(igl::opengl::glfw::Viewer &viewer, const int &idx);
	//void visualizeNeighboringRings(igl::opengl::glfw::Viewer &viewer);
	void visualizeDijkstra(igl::opengl::glfw::Viewer &viewer);
	void visualizeEigenfields(igl::opengl::glfw::Viewer &viewer, int i);
	void visualizeApproxEigenfields(igl::opengl::glfw::Viewer &viewer, int i, int iRef);
	void visualizeArbField(igl::opengl::glfw::Viewer &viewer);
	void visualizeRandomFace(igl::opengl::glfw::Viewer &viewer, const int &faceID);
	void visualizeDijkstraFace(igl::opengl::glfw::Viewer &viewer);
	void visualizeSubdomain(igl::opengl::glfw::Viewer &viewer);
	void visualizeSamples(igl::opengl::glfw::Viewer &viewer);
	void visualizeSharedEdges(igl::opengl::glfw::Viewer &viewer);
	void visualizeLocalSubdomain(igl::opengl::glfw::Viewer &viewer);
	void visualizeParallelTransportPath(igl::opengl::glfw::Viewer &viewer);
	void visualizeParallelTransport(igl::opengl::glfw::Viewer &viewer);
	void visualizeCurveConstraints(igl::opengl::glfw::Viewer &viewer);
	void visualizeSoftConstraints(igl::opengl::glfw::Viewer &viewer);
	void visualize1FieldOnCenter(igl::opengl::glfw::Viewer &viewer, const bool& even);	
	void visualizePatchDijkstra(igl::opengl::glfw::Viewer &viewer);
	void visualizeAreaOfLaplaceConstraint(igl::opengl::glfw::Viewer &viewer);
	void visualizeGradientFields(igl::opengl::glfw::Viewer &viewer);
	void visualizeRotationTests(igl::opengl::glfw::Viewer &viewer);

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
	void write2DFieldsForVTK(const Eigen::VectorXd &field2D, const string& dataName, const string& filename);

	// VISUALIZATION of APPLICATIONs
	void visualizeSmoothing(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd& v);
	void visualizeCurvatureTensor(igl::opengl::glfw::Viewer &viewer);

	void writeVectorFieldsToFile(const Eigen::VectorXd &vfields, const string& filename);
	
	// GETTER AND SETTER of IMPORTANT ELEMENTS
	Eigen::VectorXd getRefFields() const; 

protected:
	// Variable (Matrix, Vector or regular variables) for Matrix construction
	//Eigen::MatrixXd					V, FC, NF;
	Eigen::MatrixXd					V, NF, VArrow;
	Eigen::MatrixXi					F, FArrow, E, AdjMF3N, EdgePairMatrix;
	Eigen::SparseMatrix<double>		MV, MVinv, MF2D, MF2Dinv, MStar, MStarInv, MF3D, MF3Dinv, MF2DhNeg, MF2DhPos, SF2D, SF3D, B2D, B2DAsym;
	Eigen::SparseMatrix<double>		SF2DAsym;
	Eigen::SparseMatrix<double>		GF3D, GF2D, GFStar3D, GFStar2D, Div3D, Div2D, Curl3D, Curl2D, A, AT2R, J, J3D;
	Eigen::SparseMatrix<bool>		VFAdjacency;
	Eigen::VectorXd					doubleArea;
	vector<set<int>>				AdjMV, AdjMF2Ring, NeighRing;
	vector<set<VtoFPair>>			VFNeighbors, VFNeighFull;
	vector<set<Edge_VPair>>			EdgePairsList;
	vector<set<FacePair>>			AdjMF3N_temp;
	vector<set<int>>				VENeighbors;
	Eigen::MatrixXi					FE, EF;
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
	double							numSupport; 

	// Variable related to manipulation within the subspace
	Eigen::MatrixXd					cBar;
	//Eigen::VectorXd					XLowDim, XFullDim;
	Eigen::SparseMatrix<double>		CBar, B2DBar;

	// Variables related to Applications
	Eigen::MatrixXd					CurvatureTensorField2D;
	Eigen::SparseMatrix<double>		CurvatureTensor2D; 

	// FOR TESTING ONLY
public: 
	//Eigen::VectorXd					dijkstraFace, arbField, arbField2D, wb;
	Eigen::VectorXd					localSystem;
	Eigen::VectorXd					eigValuesFull, eigValuesReduced;
	Eigen::MatrixXd					eigFieldReduced2D, eigFieldFull2D;
	Eigen::MatrixXd					eigFieldsLocal;
	vector<vector<Eigen::Vector2d>> mappedBasis; 
	vector<vector<Eigen::Vector2d>> mappedBasis2;
	vector<int>						PTpath, PTsharedEdges, localPatchElements;
	vector<Eigen::Vector2d>			parallelTransport;
	vector<vector<int>>				sharedEdgesVect, curvesConstraints;
	Eigen::VectorXd					XLowDim, XFullDim;
	//Eigen::VectorXd					dijkstraFace, arbField, arbField2D, wb;
	Eigen::VectorXd					dijkstraFace, arbField, arbField2D, arbFieldE3D, arbFieldE2D, wb, wbEigen, projRef, projApprox, pertFields;
	Eigen::VectorXd					sampleDistance, patchDijkstraDist;
	vector<vector<Eigen::Vector2d>>	constraintVect2D;
	Eigen::MatrixXd					FC;
	int								testID;
private:
	
};


#endif // !VECTOR_FIELDS_H

