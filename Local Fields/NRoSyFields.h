#ifndef NROSY_FIELDS_H
#define NROSY_FIELDS_H

#include "VectorFields.h"

using namespace std;

struct NRoSy
{
	Eigen::VectorXd					magnitude;
	Eigen::VectorXd					theta;
};

class NRoSyFields
{
public:
	/* Reading data*/
	void readMesh(const string &filename);
	void readMesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);
	void scaleMesh();
	void getVF(Eigen::MatrixXd &V, Eigen::MatrixXi &F);
	void computeFaceCenter();
	void constructFaceAdjacency3NMatrix();
	void constructFaceAdjacency2RingMatrix();
	void constructEVList();
	void constructEFList();
	void computeAverageEdgeLength();


	/* Utilities */
	void constructFrameBasis();
	void constructMappingMatrix();
	void selectFaceToDraw(const int& numFaces);
	void computeDijkstraDistanceFaceForSampling(const int &source, Eigen::VectorXd &D);
	void constructMassMatrixMF3D();
	void buildStiffnessMatrix_Combinatorial();
	void buildStiffnessMatrix_Geometric();
	void computeFrameRotation(igl::opengl::glfw::Viewer &viewer);
	void computeEigenFields_generalized(const int &numEigs, const string& filename);
	void computeEigenFields_regular(const int &numEigs, const string& filename);
	void storeBasis(const string& filename);
	void retrieveBasis(const string& filename);

	/* Creating NRoSyFields */
	void representingNRoSyFields(const Eigen::MatrixXd& NFields);
	void constructNRoSyFields(const int& nRoSy, const Eigen::MatrixXd& NFields);
	void constructNRoSyFields(const int& nRoSy, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);

	/* Creating 2 fields from 2 (maximal) curvature vector fields */
	void computeMaximalPrincipalCurvature(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::VectorXd &PD, Eigen::VectorXd& PV);
	void smoothNRoSyFields(double lambda, Eigen::PardisoLDLT<Eigen::SparseMatrix<double>>& sparseSolver, const Eigen::VectorXd& inputFields, Eigen::VectorXd& outputFields);

	/* Rep. Vectors and N-RoSy Fields interface */
	void createNRoSyFromVectors(const Eigen::VectorXd& vectorFields);
	void createNRoSyFromVectors(const Eigen::VectorXd& vectorFields, NRoSy& nRoSyFields);
	void convertNRoSyToRepVectors(const NRoSy& nRoSyFields, Eigen::VectorXd& repVect);
	void convertRepVectorsToNRoSy(const Eigen::VectorXd& repVect, NRoSy& nRoSyFields);

	/* Visualizing the NRoSyFields */
	void visualizeNRoSyFields(igl::opengl::glfw::Viewer &viewer, const NRoSy& nRoSyFields, const Eigen::RowVector3d& color);
	void visualizeRepVectorFields(igl::opengl::glfw::Viewer &viewer, const NRoSy& nRoSyFields, const Eigen::RowVector3d& color);
	void visualizeRepVectorFields(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd& repVector, const Eigen::RowVector3d& color);
	void visualize2Dfields(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color, const double& scale, const bool& normalized = false);
	void visualizeEigenFields(igl::opengl::glfw::Viewer &viewer, const int id);
	void visualizeBasis(igl::opengl::glfw::Viewer &viewer, const int &id);
	void visualizeConstrainedFields(igl::opengl::glfw::Viewer &viewer);
	void visualizeConstrainedFields_Reduced(igl::opengl::glfw::Viewer &viewer);
	void visualizeConstraints(igl::opengl::glfw::Viewer &viewer);
	void visualizeSoftConstraints(igl::opengl::glfw::Viewer &viewer);

	/* N-FIELDS DESIGN */
	void nRoSyFieldsDesignRef();
	void nRoSyFieldsDesignRef_HardConstraints();
	void createAlignmentField(Eigen::VectorXd& v);
	void constructRandomHardConstraints(Eigen::SparseMatrix<double>& C, Eigen::VectorXd& c);
	void setupRHSBiharmSystemRef(const Eigen::SparseMatrix<double>& B2F, const Eigen::SparseMatrix<double>& C, const Eigen::VectorXd& c, Eigen::VectorXd& g, Eigen::VectorXd& h, Eigen::VectorXd& vEst,  const vector<double>& lambda, Eigen::VectorXd& b);
	void setupLHSBiharmSystemRef(const Eigen::SparseMatrix<double>& B2F, const Eigen::SparseMatrix<double>& C, const Eigen::VectorXd& c, const vector<double>& lambda, Eigen::SparseMatrix<double>& A_LHS);
	void solveBiharmSystemRef(const Eigen::VectorXd& vEst, const Eigen::SparseMatrix<double>& A_LHS, const Eigen::VectorXd& b, Eigen::VectorXd& Xf);

	// Interactive constraints for spline editing
	void nRoSyFieldsDesignRef_Splines();
	void pushNewUserConstraints(const int& fInit, const int& fEnd);
	void constructInteractiveConstraints();
	void setupWeight(double mu, vector<double>& lambda);

	void resetInteractiveConstraints();
	void setupRHSBiharmSystemRef_Chris(const Eigen::SparseMatrix<double>& B2F, const vector<double>& lambda, const Eigen::VectorXd& c, Eigen::VectorXd& g, Eigen::VectorXd& h, Eigen::VectorXd& b);
	void setupLHSBiharmSystemRef_Chris(const Eigen::SparseMatrix<double>& B2F, const vector<double>& lambda, const Eigen::SparseMatrix<double>& C, Eigen::SparseMatrix<double>& A_LHS);
	void solveBiharmSystemRef_Chris(const Eigen::SparseMatrix<double>& A_LHS, const Eigen::VectorXd& b, Eigen::VectorXd& Xf);

		// Soft constraints
	void nRoSyFieldsDesignRef_SoftConstraints();
	void constructSoftConstraints();
	void constructCurvesAsConstraints(const int& init, const int& end, vector<int>& curve);
	void measureSoftConstraintError(const Eigen::Vector3d& lambda);
	void projectCurvesToFrame();
	void setupRHSGlobalProblemSoftConstraints(const Eigen::Vector3d& lambda, Eigen::VectorXd& b);
	void setupLHSGlobalProblemSoftConstraints(const Eigen::Vector3d& lambda, Eigen::SparseMatrix<double>& A_LHS);
	void solveGlobalSystemMappedLDLTSoftConstraints(const Eigen::VectorXd& vEst, Eigen::SparseMatrix<double>& A_LHS, Eigen::VectorXd& b);


	/* SUBSPACE CONSTRUCTION */
	void constructBasis();
	void constructSamples(const int &n);
	void farthestPointSampling();
	void constructBasis_LocalEigenProblem();
	void gatherBasisElements(const vector<vector<Eigen::Triplet<double>>> &UiTriplet, const int& NUM_EIGEN);

	/* REDUCED N-FIELDS DESIGN */
	void nRoSyFieldsDesign_Reduced();
			// Hard constraints
	void nRoSyFieldsDesign_Reduced_HardConstraints();
	void constructRandomHardConstraints_Reduced();
	void setupRHSBiharmSystem_Reduced(const Eigen::SparseMatrix<double>& B2FBar, const Eigen::SparseMatrix<double>& MFBar, Eigen::VectorXd& gBar, Eigen::VectorXd& hBar, Eigen::VectorXd& vEstBar, const vector<double>& lambda, Eigen::VectorXd& bBar);
	void setupLHSBiharmSystem_Reduced(const Eigen::SparseMatrix<double>& B2FBar, const Eigen::SparseMatrix<double>& MFBar, const vector<double>& lambda, Eigen::SparseMatrix<double>& A_LHSBar);
	void solveBiharmSystem_Reduced(const Eigen::VectorXd& vEstBar, const Eigen::SparseMatrix<double>& A_LHSBar, const Eigen::VectorXd& bBar);
	
	void nRoSyFieldsDesign_Reduced_Splines();
	void getReducedConstraints();
	void setupRHSSplines_Reduced(const Eigen::SparseMatrix<double>& B2FBar, const vector<double>& lambda, Eigen::VectorXd& gBar, Eigen::VectorXd& hBar, Eigen::VectorXd& bBar);
	void setupLHSSplines_Reduced(const Eigen::SparseMatrix<double>& B2FBar, const vector<double>& lambda, Eigen::SparseMatrix<double>& A_LHSBar);
	void solveSplines_Reduced(const Eigen::SparseMatrix<double>& A_LHSBar, const Eigen::VectorXd& bBar);
	
	// SOFT Constraints
	void nRoSyFieldsDesign_Reduced_SoftConstraints();
	void constructSoftConstraints_Reduced();
	void setupRHSSoftConstraints_Reduced(const Eigen::Vector3d& lambda, Eigen::VectorXd& bBar);
	void setupLHSSoftConstraints_Reduced(const Eigen::Vector3d& lambda, Eigen::SparseMatrix<double>& A_LHSBar);
	void solveUserSystemMappedLDLTSoftConstraints(Eigen::SparseMatrix<double>& A_LHSBar, Eigen::VectorXd& bBar);

	void measureAccuracy();

	/* Testing stuff */
	void TEST_NROSY(igl::opengl::glfw::Viewer &viewer, const string& meshFile);
	void writeNRoSyFieldsToFile(const NRoSy& nRoSy, const string& filename);
	void writeConstraintsToFile(const string& filename);

	/* PROJECTION ON REDUCED FIELDS */
	void testProjection_MyBasis_NoRegularizer(const Eigen::SparseMatrix<double>& Basis, const Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> &sparseSolver, const Eigen::SparseMatrix<double>& B, const Eigen::VectorXd& a, const Eigen::VectorXd& inputFields, double &error);

public:
	NRoSy							nRoSy;
	Eigen::MatrixXd					V, FC;
	Eigen::MatrixXi					F, E;
	int								nRot;
	Eigen::VectorXd					frameBasis, repVector; 
	Eigen::SparseMatrix<double>		A; 
	vector<int>						FaceToDraw;
	double							avgEdgeLength;
	Eigen::SparseMatrix<double>		MF, MFinv;				// Triangle/face-based mass matrices (3 values per face)
	Eigen::SparseMatrix<double>		MF2DhNeg, MF2DhPos;		// M^(-1/2)
	Eigen::SparseMatrix<double>		SF;						//  Harmonic/Dirichlet Energy
	Eigen::VectorXd					doubleArea;				// (double) Area of each triangle
	Eigen::MatrixXi					FE, EF;					// Face-Edge and Edge-Face neighboring information matrix
	vector<set<int>>				VENeighbors;			// Vertex-Edge neighboring information
	vector<set<int>>				AdjMF2Ring;				// 2-ring neighborhood of triangles
	Eigen::MatrixXi					AdjMF3N;				// List of 3 neighbors of a triangle
	Eigen::MatrixXd					FrameRot;				// Rotation angle on each frame to the shared edge of two neighboring triangles
	Eigen::MatrixXd					eigFieldsNRoSyRef;
	Eigen::VectorXd					eigValuesNRoSyRef;

	// Variable related to subspace construction
	Eigen::SparseMatrix<double>		Basis;
	vector<int>						Sample;
	int								numSample;
	double							numSupport;
	Eigen::VectorXd					sampleDistance;
	Eigen::VectorXd					localSystem;
	set<int>						SubDomain, Boundary;

	/* Variable related to n-RoSy fields design */
	Eigen::VectorXd					Xf;
	Eigen::VectorXd					c;										// representation vector of the constraints
	Eigen::SparseMatrix<double>		C;										// selector matrix
	vector<int>						reducedConstraints, globalConstraints, userVisualConstraints;
	Eigen::VectorXd					alignFields;
	Eigen::SparseMatrix<double>		BF, BM;

	/* Variable related to REDUCED n-RoSy fields design */
	Eigen::VectorXd					XfBar;
	Eigen::VectorXd					cBar;										// representation vector of the constraints
	Eigen::SparseMatrix<double>		CBar;										// selector matrix
	Eigen::SparseMatrix<double>		BFBar, BMBar, MFBar;
	
	/* Variable on projection */
	Eigen::VectorXd					wb;											// projected representation fields
	vector<vector<int>>				curvesConstraints;
	vector<vector<Eigen::Vector2d>>	constraintVect2D;

	/* Testing variables */
	int								testID = 5;

};

#endif