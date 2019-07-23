#ifndef TENSOR_FIELDS_H
#define TENSOR_FIELDS_H

#include "Utility.h"

#include <igl/opengl/glfw/Viewer.h>

using namespace std;

class TensorFields
{
public:
	TensorFields();
	~TensorFields();

	/* Mesh-related items*/
	void readMesh(const string &meshFile);
	void getVF(Eigen::MatrixXd &V, Eigen::MatrixXi &F);
	void scaleMesh();
	void computeFaceCenter();
	void computeEdges();
	void computeAverageEdgeLength();
	void computeFaceNormal();
	void constructEVList();
	void constructEFList();

	/* SETTING UP UTILITY MATRICES */
	void constructMassMatrixMF3D();
	void constructMappingMatrix();
	void constructFaceAdjacency3NMatrix();
	void constructFaceAdjacency2RingMatrix();
	void selectFaceToDraw(const int& numFaces);
	void computeDijkstraDistanceFaceForSampling(const int &source, Eigen::VectorXd &D);
	void computeFrameRotation(igl::opengl::glfw::Viewer &viewer);
	void obtainTransformationForLaplacian(double tetha, Eigen::Matrix3d& G);
	void obtainTransformationForLaplacian(double cT, double sT, double cF, double sF, Eigen::Matrix3d& G);
	void buildStiffnessMatrix_Combinatorial();
	void buildStiffnessMatrix_Geometric();
	void convertTensorToVoigt(const Eigen::MatrixXd& tensor, Eigen::VectorXd& voigt);
	void convertTensorToVoigt_Elementary(const Eigen::Matrix2d& tensor, Eigen::Vector3d& voigt);
	void convertVoigtToTensor(const Eigen::VectorXd& voigt, Eigen::MatrixXd& tensor);
	void convertVoigtToTensor_Elementary(const Eigen::Vector3d& voigt, Eigen::Matrix2d& tensor);
	void constructTensorRepFields(const Eigen::MatrixXd& tensor, Eigen::MatrixXd& matrixRep);
	void storeBasis(const string& filename);
	void retrieveBasis(const string& filename);

	/* ADDITIONAL STUFF */
	void computeEigenFields_regular(const int &numEigs, const string& filename);
	void computeEigenFields_generalized(const int &numEigs, const string& filename);
	void loadEigenFields(const string& filename);

	/* SUBSPACE CONSTRUCTION */
	void constructBasis();
	void constructSamples(const int &n);
	void farthestPointSampling();
	void constructBasis_LocalEigenProblem();
	void gatherBasisElements(const vector<vector<Eigen::Triplet<double>>> &UiTriplet, const int& NUM_EIGEN);

	/* TENSOR STUFF */
	void computeTensorFields();
	void constructCurvatureTensor(igl::opengl::glfw::Viewer &viewer);
	void constructVoigtVector();

	/* VISUALIZATION */
	void visualize2Dfields(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color, const double& scale, const bool& normalized = false);
	void visualizeTensorFields(igl::opengl::glfw::Viewer &viewer, const Eigen::MatrixXd& tensorFields_);
	void visualizeSmoothedTensorFields(igl::opengl::glfw::Viewer &viewer);
	void visualizeSmoothedAppTensorFields(igl::opengl::glfw::Viewer &viewer);
	void visualizeSamples(igl::opengl::glfw::Viewer &viewer);
	void visualizeEigenTensorFields(igl::opengl::glfw::Viewer &viewer, int id);
	void visualizeBasis(igl::opengl::glfw::Viewer &viewer, const int &id);
	void visualizeReducedTensorFields(igl::opengl::glfw::Viewer &viewer);

	/* TESTING STUFF*/
	void TEST_TENSOR(igl::opengl::glfw::Viewer &viewer, const string& meshFile);
	void testDirichletAndLaplace();
	void testSmoothing(igl::opengl::glfw::Viewer &viewer, const Eigen::MatrixXd& inputTensor, Eigen::MatrixXd& outputTensor);
	void testTransformation(igl::opengl::glfw::Viewer &viewer);
	void tensorConvertNConvert(igl::opengl::glfw::Viewer &viewer);

	/* APPLICATION :: SMOOTHING */
	void smoothingRef(igl::opengl::glfw::Viewer &viewer, const Eigen::MatrixXd& inputTensor, Eigen::MatrixXd& outputTensor);
	void smoothing_Explicit_Combinatorial(igl::opengl::glfw::Viewer &viewer, const Eigen::MatrixXd& inputTensor, Eigen::MatrixXd& outputTensor);
	void smoothing_Explicit_Geometric(igl::opengl::glfw::Viewer &viewer, const Eigen::MatrixXd& inputTensor, Eigen::MatrixXd& outputTensor);
	void smoothing_Implicit_Combinatorial(igl::opengl::glfw::Viewer &viewer, const Eigen::MatrixXd& inputTensor, Eigen::MatrixXd& outputTensor);
	void smoothing_Implicit_Geometric(igl::opengl::glfw::Viewer &viewer, const Eigen::MatrixXd& inputTensor, Eigen::MatrixXd& outputTensor);

	void smoothingRed(igl::opengl::glfw::Viewer &viewer, const Eigen::MatrixXd& inputTensor, Eigen::MatrixXd& outputTensor);
	void smoothingRed_Explicit_Geometric(igl::opengl::glfw::Viewer &viewer, const Eigen::MatrixXd& inputTensor, Eigen::MatrixXd& outputTensor);
	void smoothingRed_Implicit_Geometric(igl::opengl::glfw::Viewer &viewer, const Eigen::MatrixXd& inputTensor, Eigen::MatrixXd& outputTensor);

	/* APPLICATION :: Sub-space Projection */
	void subspaceProjection(const Eigen::VectorXd& refField);


/* For convenience, all variables that should be private will be declared protected in this prototyping stage */
public:
	Eigen::MatrixXd					Tensor, tensorFields;	// 2-by-2 tensor and the representing vectors (using principal curvatures)	
	Eigen::VectorXd					voigtReps;				// a 3-by-1 representation of 2-by-2 tensor
	Eigen::MatrixXd					V, FC, NF;				// Vertex, Face-center, and Face-normals
	Eigen::MatrixXd					FrameRot;				// Rotation angle on each frame to the shared edge of two neighboring triangles
	Eigen::MatrixXi					F, E, AdjMF3N;			// Faces, Edges, and Face-Face adjacency matrices
	Eigen::VectorXd					doubleArea;				// (double) Area of each triangle
	Eigen::SparseMatrix<double>		A;						// A map from local to to world coordinate
	Eigen::SparseMatrix<double>		MF, MFinv;				// Triangle/face-based mass matrices (3 values per face)
	Eigen::SparseMatrix<double>		SF;						// Laplace matrix for the tensor fields (finite difference approach)
	double							avgEdgeLength;			// avg edge length -> to scale the fields
	vector<int>						FaceToDraw;				// Indices to faces that we'll draw the fields upon
	vector<set<int>>				VENeighbors;			// Vertex-Edge neighboring information
	vector<set<int>>				AdjMF2Ring;				// 2-ring neighborhood of triangles
	Eigen::MatrixXi					FE, EF;					// Face-Edge and Edge-Face neighboring information matrix
	//double							scale = 10000000;
	//double							scale = 100;		// regular eigenfields => arma 10k
	//double							scale = 1000.0;		// regular eigenfields => arma 43k
	//double								scale = 10;
	//double								scale = 0.01;		// smoothing
	double								scale = 0.1;		// smoothing torus
	//double							scale = 2.0; 

	//
	Eigen::MatrixXd eigFieldsTensorRef;
	Eigen::VectorXd eigValuesTensorRef;

	// Variable related to subspace construction
	Eigen::SparseMatrix<double>		Basis;
	vector<int>						Sample;
	int								numSample;
	double							numSupport;
	Eigen::VectorXd					sampleDistance;
	Eigen::VectorXd					localSystem;
	set<int>						SubDomain, Boundary;

	/* Application */
	Eigen::MatrixXd					smoothedTensorRef;			// Smoothed tensor, application;
	Eigen::MatrixXd					smoothedTensorRed;			// Smoothed tensor, application;

	/* Projection */
	Eigen::MatrixXd					TensorRed;
private:

};
#endif // !TENSOR_FIELDS_H

