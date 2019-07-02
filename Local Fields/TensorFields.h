#ifndef TENSOR_FIELDS_H
#define TENSOR_FIELDS_H

#include "Utility.h"

using namespace std;

class TensorFields
{
public:
	TensorFields();
	~TensorFields();
	void constructCurvatureTensor();

	/* Mesh-related items*/
	void readMesh(const string &meshFile);
	void scaleMesh();
	void computeFaceCenter();
	void computeEdges();

	/* SETTING UP UTILITY MATRICES */
	void computeAverageEdgeLength();
	void constructMappingMatrix();
	void constructFaceAdjacency3NMatrix();

/* For convenience, all variables that should be private will be declared protected in this prototyping stage */
protected:
	Eigen::MatrixXd					Tensor;
	Eigen::VectorXd					vectorReps;
	Eigen::MatrixXd					V, FC, NF;
	Eigen::MatrixXi					F, E, AdjMF3N;
	Eigen::VectorXd					doubleArea;
	Eigen::SparseMatrix<double>		A;
	double							avgEdgeLength;
private:

};
#endif // !TENSOR_FIELDS_H

