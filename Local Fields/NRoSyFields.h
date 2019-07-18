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

	/* Creating NRoSyFields */
	void representingNRoSyFields(const Eigen::MatrixXd& NFields);
	void constructNRoSyFields(const int& nRoSy, const Eigen::MatrixXd& NFields);
	void constructNRoSyFields(const int& nRoSy, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);

	/* Rep. Vectors and N-RoSy Fields interface */
	void createNRoSyFromVectors(const Eigen::VectorXd& vectorFields);
	void createNRoSyFromVectors(const Eigen::VectorXd& vectorFields, NRoSy& nRoSyFields);
	void convertNRoSyToRepVectors(const NRoSy& nRoSyFields, Eigen::VectorXd& repVect);
	void convertRepVectorsToNRoSy(const Eigen::VectorXd& repVect, NRoSy& nRoSyFields);

	/* Visualizing the NRoSyFields */
	void visualizeNRoSyFields(igl::opengl::glfw::Viewer &viewer, const NRoSy& nRoSyFields);
	void visualizeRepVectorFields(igl::opengl::glfw::Viewer &viewer, const NRoSy& nRoSyFields);
	void visualizeRepVectorFields(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd& repVector);
	void visualize2Dfields(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color, const double& scale, const bool& normalized = false);
	void visualizeEigenFields(igl::opengl::glfw::Viewer &viewer, const int id);
	void visualizeBasis(igl::opengl::glfw::Viewer &viewer, const int &id);
	void visualizeConstrainedFields(igl::opengl::glfw::Viewer &viewer);
	void visualizeConstraints(igl::opengl::glfw::Viewer &viewer);

	/* N-FIELDS DESIGN */
	void nRoSyFieldsDesignRef();
	void nRoSyFieldsDesignRef_HardConstraints();
	void constructRandomHardConstraints(Eigen::SparseMatrix<double>& C, Eigen::VectorXd& c);
	void setupRHSBiharmSystemRef(const Eigen::SparseMatrix<double>& B2F, const Eigen::SparseMatrix<double>& C, const Eigen::VectorXd& c, Eigen::VectorXd& g, Eigen::VectorXd& h, Eigen::VectorXd& vEst, Eigen::VectorXd& b);
	void setupLHSBiharmSystemRef(const Eigen::SparseMatrix<double>& B2F, const Eigen::SparseMatrix<double>& C, const Eigen::VectorXd& c, Eigen::SparseMatrix<double>& A_LHS);
	void solveBiharmSystemRef(const Eigen::VectorXd& vEst, const Eigen::SparseMatrix<double>& A_LHS, const Eigen::VectorXd& b, Eigen::VectorXd& Xf);

	/* SUBSPACE CONSTRUCTION */
	void constructBasis();
	void constructSamples(const int &n);
	void farthestPointSampling();
	void constructBasis_LocalEigenProblem();
	void gatherBasisElements(const vector<vector<Eigen::Triplet<double>>> &UiTriplet, const int& NUM_EIGEN);

	/* Testing stuff */
	void TEST_NROSY(igl::opengl::glfw::Viewer &viewer, const string& meshFile);

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
	vector<int>						reducedConstraints, globalConstraints;

	/* Testing variables */
	int								testID;

};

#endif