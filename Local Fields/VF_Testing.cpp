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
	cout << "> Constructing gradient 2D of arbitrary field...";


	Eigen::SparseMatrix<double> Grad3D, Grad2D;
	igl::grad(V, F, Grad3D);
	printf("Size of grad3d: %dx%d\n", Grad3D.rows(), Grad3D.cols());
	rearrangeGradient3D(Grad3D);
	printf("Size of grad3d [2]: %dx%d\n", Grad3D.rows(), Grad3D.cols());
	Grad2D = A.transpose()*Grad3D; 
	printf("Grad3d [2]: %dx%d || arbFields: %d. \n", Grad2D.rows(), Grad2D.cols(), arbField.rows());
	arbField2D = Grad2D * arbField;

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;
}

void VectorFields::testBasis_NoRegularizer(double &error)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Testing basis 2D of arbitrary field without regularizer...\n";

	// Construct matrices for Test
	cout << "____Assigning variables\n";
	Eigen::SparseMatrix<double> U = Basis;// BasisTemp;
	Eigen::VectorXd				v = Xf; 
	//Eigen::VectorXd				v =  arbField2D;
	Eigen::VectorXd				a = (v.transpose()*MF2D*U).transpose();
	Eigen::SparseMatrix<double> B = U.transpose() * MF2D * U;

	cout << "____Solving linear system variables\n";
	//Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> sparseSolver(B);
	Eigen::PardisoLLT<Eigen::SparseMatrix<double>> sparseSolver(B);

	Eigen::VectorXd w = sparseSolver.solve(a);

	if (sparseSolver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system. " << endl;
		if (sparseSolver.info() == Eigen::NumericalIssue)
			cout << "NUMERICAL ISSUE. " << endl;
		cout << sparseSolver.info() << endl;
		return;
	}

	//Eigen::VectorXd wb;
	wb.resize(U.rows());

	for (int i = 0; i < U.rows(); i++) {
		wb(i) = 0.0;
	}

	cout << "____Getting total SUM(wi*bi) \n";
	//for (int i = 0; i < w.rows(); i++) {
	//	wb += w(i)*U.col(i);
	//}
	wb = U*w; 

	// Compare their L2-Norm
	cout << "____Computing L2-norm \n";
	Eigen::VectorXd diff = v - wb; 
	double norm1 = diff.transpose()*MF2D*diff;
	double norm2 = v.transpose()*MF2D*v;
	double normL2 = sqrt(norm1 / norm2); 
	error = normL2; 

	cout << "The L-2 Norm is << " << normL2 << ".\n" << endl; 

	//cout << "Weight: \n" << w << endl; 

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;
}

void VectorFields::testBasis_WithRegularizer()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Testing basis 2D of arbitrary field (with regularizer)...\n";

	// Construct matrices for Test
	cout << "____Assigning variables\n";
	Eigen::SparseMatrix<double> U = Basis;// BasisTemp;
	Eigen::VectorXd				v = Xf; 
	//Eigen::VectorXd				v = arbField2D;
	Eigen::VectorXd				a = (U.transpose()*MF2D*v).transpose();
	//Eigen::SparseMatrix<double> B1 = U.transpose() * MF2D * U;
	//Eigen::SparseMatrix<double> B2 = U.transpose() * SF2D * MF2D * U;
	const double lambda = 0.005; 
	//Eigen::SparseMatrix<double> B = B1 + lambda*B2;
	Eigen::SparseMatrix<double> B = U.transpose()*(MF2D+lambda*B2D)*U;

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
	wb.resize(U.rows());

	for (int i = 0; i < U.rows(); i++) {
		wb(i) = 0.0;
	}

	cout << "____Getting total SUM(wi*bi) \n";
	//for (int i = 0; i < w.rows(); i++) {
	//	wb += w(i)*U.col(i);
	//}
	wb = U*w; 

	// Compare their L2-Norm
	cout << "____Computing L2-norm \n";
	Eigen::VectorXd diff = v - wb;
	double norm1 = diff.transpose()*MF2D*diff;
	double norm2 = v.transpose()*MF2D*v;
	double normL2 = sqrt(norm1 / norm2);

	cout << "The L-2 Norm is << " << normL2 << ".\n" << endl;

	double energyInput = v.transpose()*B2D*v; 
	double energyOutput = wb.transpose()*B2D*v;
	cout << "Energy: Input= " << energyInput << ", \t\t Output= " << energyOutput << endl; 

	//cout << "Weight: \n" << w << endl; 

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;
}

void VectorFields::projectionTest()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;	


	const int NUM_TEST = 1;
	Eigen::VectorXd errors(NUM_TEST);

	for (int i = 0; i < NUM_TEST; i++)
	{
		t1 = chrono::high_resolution_clock::now();

		/* Reference results */
		setupGlobalProblem();

		/* Projection to the subspace */
		//testBasis_NoRegularizer(errors(i));
		testBasis_WithRegularizer();

		t2 = chrono::high_resolution_clock::now();
		duration = t2 - t1;
		printf("[%d] run => Error=%.10f (in %.3f seconds) \n", i, errors(i), duration.count());		
	}

	cout << "ERRORS: \n" <<  errors << endl; 
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
	//std::rand() / ((RAND_MAX) / 6);
	//int randID = rand() % F.rows();
	int randID = rand() / (RAND_MAX/F.rows());
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

void VectorFields::testEdgesAddition(igl::opengl::glfw::Viewer &viewer)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Adding edges... ";

	/* Adding edges */
	Eigen::MatrixXd transformedFields(F.rows(), F.cols());
	Eigen::VectorXd fields3D = A*arbField2D; 
	for (int i = 0; i < F.rows(); i++)
	{
		transformedFields.row(i) = (fields3D.block(3 * i, 0, 3, 1)).transpose();
	}

	viewer.data().add_edges(FC, FC + transformedFields*avgEdgeLength, Eigen::RowVector3d(0.0, 0.8, 0.1));


	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;

}
