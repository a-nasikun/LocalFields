#include "VectorFields.h"

/* ====================== ITEMS FOR TESTING ONLY ============================*/
void VectorFields::constructArbitraryField()
{
	// Construct the Field
	arbField.resize(V.rows());

	// Random Arbitrary field
	//for (int i = 0; i < V.rows(); i++) {
	//	arbField(i) = 0.0;
	//}
	//
	//srand(time(NULL));
	//int pID = rand() % V.rows();
	//int pID = 0;
	//arbField(pID) = 1.0;

	// Dijstra-based Arbitrary Field
	int pID = *(NeighRing[0].begin());
	//computeDijkstraDistanceVertex(pID, arbField);

	// Position based arbitrary scalar field
	for (int i = 0; i < V.rows(); i++) {
		arbField(i) = V(i, 0) *  V(i, 1) *  V(i, 2);
	}

}

void VectorFields::testMappingMatrix()
{
	// Should be identity
	Eigen::SparseMatrix<double> ATA = A.transpose()*A;
	//visualizeSparseMatrixInMatlab(ATA);
	cout << "Block of ATA " << endl << ATA.block(5, 5, 5, 5) << endl;

	// Orthogonality of the matrix
	Eigen::Vector3d e, f;
	srand(time(NULL));
	const int idx = rand() % F.rows();
	e << A.coeffRef(3 * idx, 2 * idx), A.coeffRef(3 * idx + 1, 2 * idx), A.coeffRef(3 * idx + 2, 2 * idx);
	f << A.coeffRef(3 * idx, 2 * idx + 1), A.coeffRef(3 * idx + 1, 2 * idx + 1), A.coeffRef(3 * idx + 2, 2 * idx + 1);
	cout << idx << " => e*f = " << e.dot(f) << endl;
}

void VectorFields::testAdjMV()
{
	for (int i = 0; i < V.rows(); i++) {
		cout << i << ":";
		for (std::set<int, double>::iterator it = AdjMV[i].begin(); it != AdjMV[i].end(); ++it) {
			cout << *it << ", ";
		}
		cout << endl;
	}
}

void VectorFields::testAdjacencyAndEdges()
{
	for (int i = 0; i < F.rows(); i++) {
		//for (int i = 0; i < min(100, (int)F.rows()); i++) {
		printf("F(%d) [%d, %d, %d] :=>"
			"(0) F(%d)[%d, %d, %d] on edges(%d, %d),"
			"(1) F(%d)[%d, %d, %d] on edges(%d, %d),"
			"(2) F(%d)[%d, %d, %d] on edges(%d, %d)\n",
			i, F(i, 0), F(i, 1), F(i, 2),
			AdjMF3N(i, 0), F(AdjMF3N(i, 0), 0), F(AdjMF3N(i, 0), 1), F(AdjMF3N(i, 0), 2), EdgePairMatrix(i, 0), EdgePairMatrix(i, 1),
			AdjMF3N(i, 1), F(AdjMF3N(i, 1), 0), F(AdjMF3N(i, 1), 1), F(AdjMF3N(i, 1), 2), EdgePairMatrix(i, 2), EdgePairMatrix(i, 3),
			AdjMF3N(i, 2), F(AdjMF3N(i, 2), 0), F(AdjMF3N(i, 2), 1), F(AdjMF3N(i, 2), 2), EdgePairMatrix(i, 4), EdgePairMatrix(i, 5));
	}
}

void VectorFields::testDijkstraFace()
{
	dijkstraFace.resize(F.rows());

	// For single-sourced Dijkstra
	//const int source = *(NeighRing[0].begin());
	//computeDijkstraDistanceFace(source, dijkstraFace);

	// For multiple-sourced Dijkstra
	const int numSource = 5;
	Eigen::VectorXi source(numSource);
	srand(time(NULL));
	for (int i = 0; i < numSource; i++) {
		source(i) = rand() % F.rows();
	}
	computeDijkstraDistanceFaceMultSource(source, dijkstraFace);
}

void VectorFields::testCurlEnergy()
{
	//double gradEnergy = gradArbField3D.transpose() * LapCurl3D * gradArbField3D;
	//cout << "The energy is " << gradEnergy << endl;
}

int VectorFields::selectRandomFace()
{
	srand(time(NULL));
	int randID = rand() % F.rows();
	cout << "Selected face: " << randID << endl;
	return randID;
}

void VectorFields::checkB2DStructure()
{
	for (int i = 0; i < B2D.outerSize(); i++) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(B2D, i); it; ++it) {
			if (it.row() == 200 && it.col() % 2 == 0) {
				cout << "VALUES IN B2D :: N(0)=";
				cout << it.col() / 2 << ", " << endl;
			}
		}
	}

	cout << "NEIGHBORS :: N(0)=";
	for (int i = 0; i < AdjMF3N.cols(); i++) {
		cout << (AdjMF3N(100, i)) << ", ";
	}
	cout << endl;

	cout << "2RING:: N(0)=";
	for (int i : AdjMF2Ring[100]) {
		cout << i << ", ";
	}
	cout << endl;

}
