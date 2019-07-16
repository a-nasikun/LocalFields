#ifndef NROSY_FIELDS_H
#define NROSY_FIELDS_H

#include "VectorFields.h"

using namespace std;

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
	void computeAverageEdgeLength();


	/* Utilities */
	void constructFrameBasis();
	void constructMappingMatrix();
	void selectFaceToDraw(const int& numFaces);
	void computeDijkstraDistanceFaceForSampling(const int &source, Eigen::VectorXd &D);

	/* Creating NRoSyFields */
	void representingNRoSyFields(const Eigen::MatrixXd& NFields);
	void constructNRoSyFields(const int& nRoSy, const Eigen::MatrixXd& NFields);
	void constructNRoSyFields(const int& nRoSy, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);

	/* Rep. Vectors and N-RoSy Fields interface */
	void convertNRoSyToRepVectors(Eigen::VectorXd& vectorFields);
	void convertRepVectorsToNRoSy(const Eigen::VectorXd& vectorFields);
	void createNRoSyFromVectors(const Eigen::VectorXd& vectorFields);

	/* Visualizing the NRoSyFields */
	void visualizeNRoSyFields(igl::opengl::glfw::Viewer &viewer);
	void visualizeRepVectorFields(igl::opengl::glfw::Viewer &viewer);
	void visualize2Dfields(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color, const double& scale, const bool& normalized = false);

	/* Testing stuff */
	void TEST_NROSY(igl::opengl::glfw::Viewer &viewer, const string& meshFile);

private:
	Eigen::MatrixXd					V, FC;
	Eigen::MatrixXi					F, E, AdjMF3N;
	int								NRoSy;
	Eigen::VectorXd					frameBasis, magnitude, theta, repVector; 
	Eigen::SparseMatrix<double>		A; 
	vector<int>						FaceToDraw;
	double							avgEdgeLength;
	
};

#endif