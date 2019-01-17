#include "VectorFields.h"

/* ====================== ITEMS FOR TESTING ONLY ============================*/
void VectorFields::constructArbitraryField()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Constructing arbitrary field...";

	// Construct the Field
	arbField.resize(V.rows());
		
	// Dijstra-based Arbitrary Field
	//int pID = *(NeighRing[0].begin());
	//computeDijkstraDistanceVertex(pID, arbField);

	// Position based arbitrary scalar field
	for (int i = 0; i < V.rows(); i++) {
		//arbField(i) = V(i, 0) *  V(i, 1) *  V(i, 2);
		arbField(i) = V(i, 1);// *V(i, 1) *  V(i, 2);
	}

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;
}

void VectorFields::constructArbitraryField2D()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Construting gradient 2D of arbitrary field...";


	Eigen::SparseMatrix<double> Grad3D, Grad2D;
	igl::grad(V, F, Grad3D);
	rearrangeGradient3D(Grad3D);
	Grad2D = A.transpose()*Grad3D; 
	arbField2D = Grad2D * arbField;

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;
}

void VectorFields::testBasis()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Testing basis 2D of arbitrary field...\n";

	// Construct matrices for Test
	cout << "____Assigning variables\n";
	Eigen::SparseMatrix<double> BB = Basis;// BasisTemp;
	//Eigen::VectorXd				v = Xf; 
	Eigen::VectorXd				v =  arbField2D;
	Eigen::VectorXd				a = (v.transpose()*MF2D*BB).transpose();
	Eigen::SparseMatrix<double> B = BB.transpose() * MF2D * BB;

	cout << "____Solving linear system variables\n";
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> sparseSolver(B);
	Eigen::VectorXd w = sparseSolver.solve(a);

	if (sparseSolver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system. " << endl;
		if (sparseSolver.info() == Eigen::NumericalIssue)
			cout << "NUMERICAL ISSUE. " << endl;
		cout << sparseSolver.info() << endl;
		return;
	}

	//Eigen::VectorXd wb;
	wb.resize(BB.rows());

	for (int i = 0; i < BB.rows(); i++) {
		wb(i) = 0.0;
	}

	cout << "____Getting total SUM(wi*bi) \n";
	for (int i = 0; i < w.rows(); i++) {
		wb += w(i)*BB.col(i);
	}

	// Compare their L2-Norm
	cout << "____Computing L2-norm \n";
	Eigen::VectorXd diff = v - wb; 
	double norm1 = diff.transpose()*MF2D*diff;
	double norm2 = v.transpose()*MF2D*v;
	double normL2 = sqrt(norm1 / norm2); 

	cout << "The L-2 Norm is << " << normL2 << ".\n" << endl; 

	//cout << "Weight: \n" << w << endl; 

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;
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

void VectorFields::constructParallelTransport()
{
	Eigen::Vector2d PTvalue;
	PTvalue << sqrt(2.0)/2.0, sqrt(2.0)/2.0; 
	//PTvalue << 1.0, 0.0;

	parallelTransport.resize(PTpath.size());
	parallelTransport[0] = PTvalue;
	PTsharedEdges.resize(2 * PTpath.size());

	for (int i = 0; i < PTpath.size() - 1; i++) {
		int face1 = PTpath[i];
		int face2 = PTpath[i + 1];
		
		// 1. Find shared edge (naively)
		enum class SharedEdgeCase { Case1, Case2, Case3 };
		SharedEdgeCase edgeCase1, edgeCase2;
		Eigen::RowVector3d es;
		for (int f1 = 0; f1 < F.cols(); f1++) {
			for (int f2 = 0; f2 < F.cols(); f2++) {
				bool b1 = F(face1, (f1 + 1) % F.cols()) == F(face2, f2);
				bool b2 = F(face1, f1) == F(face2, (f2 + 1) % F.cols());
				if (b1 && b2) {						
					es = V.row(F(face1, (f1 + 1) % F.cols())) - V.row(F(face1, f1));
					PTsharedEdges[2 * i + 0] = F(face1, f1);
					PTsharedEdges[2 * i + 1] = F(face1, (f1 + 1) % F.cols());
					
					if (f1 == 0)		edgeCase1 = SharedEdgeCase::Case1;	// => edge V0->V1 is the shared edge => it takes 0 step to reach v0
					else if (f1 == 1)	edgeCase1 = SharedEdgeCase::Case3;	// => edge V1->V2 is the shared edge => it takes 2 step to reach v0
					else if (f1 == 2)	edgeCase1 = SharedEdgeCase::Case2;	// => edge V2->V0 is the shared edge => it takes 1 step to reach v0

					if (f2 == 0)		edgeCase2 = SharedEdgeCase::Case1;
					else if (f2 == 1)	edgeCase2 = SharedEdgeCase::Case3;
					else if (f2 == 2)	edgeCase2 = SharedEdgeCase::Case2;
				}
			}
		}

		// 2. Find angles between basis1 and es
		Eigen::VectorXd eVect;
		Eigen::RowVector3d b11, b12;
		eVect = A.block(3 * face1, 2 * face1 + 0, 3, 1);
		b11 << eVect(0), eVect(1), eVect(2);

		double cosR12 = (b11.dot(es)) / (b11.norm()*es.norm());
		if (cosR12 > 1.0) cosR12 = 1.0;
		if (cosR12 <-1.0) cosR12 = -1.0;
		const double angleR12_1 = (edgeCase1 == SharedEdgeCase::Case2 ? 2 * M_PI - acos(cosR12) : acos(cosR12));
		const double cosR12_1 = cos(angleR12_1);
		const double sinR12_1 = sin(angleR12_1);
		//printf("______[%.2f] Rotation matrix R12_1 = [%.3f,%.3f; %.3f, %.3f]\n", angleR12_1*180.0 / M_PI, cosR12_1, -sinR12_1, sinR12_1, cosR12_1);

		Eigen::RowVector3d b21, b22;
		eVect = A.block(3 * face2, 2 * face2, 3, 1);
		b21 << eVect(0), eVect(1), eVect(2);

		double cosR21 = (b21.dot(es)) / (b21.norm()*es.norm());
		if (cosR21 > 1.0) cosR21 = 1.0;
		if (cosR21 < -1.0) cosR21 = -1.0;
		double angleR21_1 = (edgeCase2 == SharedEdgeCase::Case3 ? 2 * M_PI - acos(cosR21) : acos(cosR21));
		const double cosR21_1 = cos(angleR21_1);
		const double sinR21_1 = sin(angleR21_1);
		//printf("______[%.2f] Rotation matrix R22_1 = [%.2f,%.2f; %.2f, %.2f]\n", angleR21_1*180.0 / M_PI, cosR21_1, -sinR21_1, sinR21_1, cosR21_1);
		
		double anglePhi1;
		Eigen::Vector2d bas(1.0, 0.0);
		if (parallelTransport[i](1) < 0)
			anglePhi1 = 2*M_PI - acos(bas.dot(parallelTransport[i])/parallelTransport[i].norm());
		else
			anglePhi1 = acos(bas.dot(parallelTransport[i]) / parallelTransport[i].norm());

		//printf("______basis to field : [%.2f]\n", anglePhi1*180.0/M_PI);
		const double anglePhi2 = anglePhi1 - angleR12_1 + angleR21_1;
		//printf("______NEW angle : [%.2f]\n", anglePhi2*180.0 / M_PI);
		const double cosBasis = cos(anglePhi2);
		const double sinBasis = sin(anglePhi2);

		// Rotate basis
		// Obtain the mapped basis
		Eigen::Matrix2d RotMat2D;
		RotMat2D << cosBasis, -sinBasis, sinBasis, cosBasis;
		Eigen::Vector2d transported = RotMat2D * bas;// parallelTransport[i];
		parallelTransport[i + 1] = transported;		
	}
}
