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
	void scaleMesh();
	void computeFaceCenter();
	void computeEdges();
	void computeAverageEdgeLength();
	void computeFaceNormal();

	/* SETTING UP UTILITY MATRICES */
	void constructMappingMatrix();
	void constructFaceAdjacency3NMatrix();
	void selectFaceToDraw(const int& numFaces);
	void computeDijkstraDistanceFaceForSampling(const int &source, Eigen::VectorXd &D);
	void computeTensorFields();
	

	/* TENSOR STUFF */
	void constructCurvatureTensor(igl::opengl::glfw::Viewer &viewer);

	/* VISUALIZATION */
	void visualize2Dfields(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color, const double& scale, const bool& normalized = false);
	void visualizeTensorFields(igl::opengl::glfw::Viewer &viewer);

/* For convenience, all variables that should be private will be declared protected in this prototyping stage */
protected:
	Eigen::MatrixXd					Tensor, tensorFields;
	Eigen::VectorXd					voigtReps;
	Eigen::MatrixXd					V, FC, NF;
	Eigen::MatrixXi					F, E, AdjMF3N;
	Eigen::VectorXd					doubleArea;
	Eigen::SparseMatrix<double>		A;
	double							avgEdgeLength;
	vector<int>						FaceToDraw;

private:

};
#endif // !TENSOR_FIELDS_H

