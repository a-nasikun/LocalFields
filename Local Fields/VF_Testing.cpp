#include "VectorFields.h"
#include <fstream>

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
		arbField(i) = V(i, 0) *  V(i, 1) *  V(i, 2);
		//arbField(i) = V(i, 2);
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

	cout << "____Getting total SUM(wi*bi) \n";
	wb = U*w; 

	// Compare their L2-Norm
	cout << "____Computing L2-norm \n";
	Eigen::VectorXd diff = v - wb; 
	double norm1 = diff.transpose()*MF2D*diff;
	double norm2 = v.transpose()*MF2D*v;
	double normL2 = sqrt(norm1 / norm2); 
	error = normL2; 

	cout << "The L-2 Norm is << " << normL2 << ".\n" << endl; 

	/* Measuring the energy */
	double energy1 = v.transpose()*B2D*v;
	double energy2 = wb.transpose()*B2D*wb;
	cout << "Bi-Harmonic ENERGY => Ref=" << energy1 << ", Approx:" << energy2 << endl;
	cout << "Relative energy: " << abs(energy1 - energy2) / energy1 << endl; 

	/* Measuring the 'length' of each vector */
	energy1 = v.transpose()*SF2D*v;
	energy2 = wb.transpose()*SF2D*wb;
	cout << "Harmonic Energy => Ref=" << energy1 << ", Approx:" << energy2 << endl;
	cout << "Relative energy: " << abs(energy1 - energy2) / energy1 << endl;

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;
}

void VectorFields::testBasis_WithRegularizer(const Eigen::SparseMatrix<double>& MReg, double &error)
{

	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Testing basis 2D of arbitrary field WITH regularizer...\n";

	// Construct matrices for Test	
	cout << "__[APPROXIMATION]....\n";
	//const double				lambda = 10000 * MF2D.coeff(0,0) / MReg.coeff(0,0);

	Eigen::SparseMatrix<double> U = Basis;// BasisTemp;
	Eigen::VectorXd				v = XFullDim;
	const double				inputEnergy = v.transpose()*MReg*v; 
	const double				lambda = 50000/inputEnergy;
	Eigen::VectorXd				a = U.transpose()*MF2D*v;
	Eigen::SparseMatrix<double> B = U.transpose() * (MF2D + lambda*MReg) * U;

	cout << "____Solving linear system variables (with lambda=" << lambda << ")\n";
	Eigen::PardisoLLT<Eigen::SparseMatrix<double>> sparseSolver(B);
	

	Eigen::VectorXd w = sparseSolver.solve(a);

	if (sparseSolver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system. " << endl;
		if (sparseSolver.info() == Eigen::NumericalIssue)
			cout << "NUMERICAL ISSUE. " << endl;
		cout << sparseSolver.info() << endl;
		return;
	}
	wb = U*w;
	
	cout << "__[REFERENCE]....\n";
	a = MF2D*v;
	B = MF2D + lambda*MReg; 
	Eigen::VectorXd  wRef;
	Eigen::PardisoLLT<Eigen::SparseMatrix<double>> sparseSolver2(B);
	wRef = sparseSolver2.solve(a);
	projRef = wRef;


	// Compare their L2-Norm
	cout << "____Computing L2-norm \n";
	Eigen::VectorXd diff = wRef - wb;
	double norm1 = diff.transpose()*MF2D*diff;
	double norm2 = wRef.transpose()*MF2D*wRef;
	double normL2 = sqrt(norm1 / norm2);
	error = normL2;
	cout << "The L-2 Norm is << " << normL2  << endl;

	/* Measuring the energy */
	double energy1 = wRef.transpose()*B2D*wRef;
	double energy2 = wb.transpose()*B2D*wb;
	cout << "BIHARMONIC ENERGY => Ref=" << energy1 << ", Approx:" << energy2 << ",initial: " << v.transpose()*B2D*v << endl;
	cout << "Relative energy: " << abs(energy1 - energy2) / energy1 << endl;

	/* Measuring the energy */
	energy1 = wRef.transpose()*SF2D*wRef;
	energy2 = wb.transpose()*SF2D*wb;
	cout << "HARMONIC ENERGY => Ref=" << energy1 << ", Approx:" << energy2 << ",initial: " << v.transpose()*SF2D*v << endl;
	cout << "Relative energy: " << abs(energy1 - energy2) / energy1 << endl;


	/* Measuring the 'length' of each vector */
	double length1 = wRef.transpose()*MF2D*wRef;
	double length2 = wb.transpose()*MF2D*wb;
	cout << "Length => Ref=" << length1 << ", Approx:" << length2 << ", initial: " << v.transpose()*MF2D*v << endl;
	cout << "Relative length: " << length2 / length1 * 100.0 << "%" <<  endl;

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;	
}

void VectorFields::projectionTest()
{
	cout << "PROJECTION TEST! \n";
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;	


	const int NUM_TEST = 1;
	Eigen::VectorXd errors(NUM_TEST);

	//Xf = XFullDim; 

	for (int i = 0; i < NUM_TEST; i++)
	{
		t1 = chrono::high_resolution_clock::now();

		

		/* Projection to the subspace */
		/* Reference results */
		//setupGlobalProblem(Eigen::Vector3d(1,1,1));
		//testBasis_NoRegularizer(errors(i));

		testBasis_WithRegularizer(B2D, errors(i));
		//testBasis_WithRegularizer(SF2D, errors(i));

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

void VectorFields::writeBasisToFile()
{
	printf("Basis=%dx%d || #F=%d\n", Basis.rows(), Basis.cols(), F.rows());
	Eigen::SparseMatrix<double> Basis3D = A*Basis;
	printf("Basis=%dx%d || Basis3D=%dx%d || #F=%d\n", Basis.rows(), Basis.cols(), Basis3D.rows(), Basis3D.cols(), F.rows());

	WriteSparseMatrixToMatlab(Basis3D, "hello");
}

void VectorFields::writeField3DToFile()
{
	Eigen::VectorXd Field3D = A*XFullDim;
	WriteDenseMatrixToMatlab(Field3D, "hello");
}

void VectorFields::printDataForVTK()
{
	/* PRint vertices*/
	cout << "POINTS " << V.rows() << " float\n";
	cout << V << endl;

	/* Print faces */
	cout << "POLYGONS " << F.rows() << " " << F.rows() * (F.cols()+1) << endl; 
	for (int i = 0; i < F.rows(); i++)
	{
		cout << F.cols() << " " << F(i, 0) << " " << F(i, 1) << " " << F(i, 2) << endl;
	}

	/* Print eigenfields*/
	Eigen::VectorXd eig3D;
	eig3D = A * eigFieldFull2D.col(0);
	cout << "CELL_DATA " << F.rows() << endl; 
	cout << "VECTORS EigenField float\n";
	for (int i = 0; i < F.rows(); i++)
	{
		cout << F.cols() << " " << eig3D(3 * i + 0) << " " << eig3D(3 * i + 1) << " " << eig3D(3 * i + 2) << endl;
	}
}

void VectorFields::writeEigenFieldsForVTK()
{
	for (int id = 0; id < 1; id++)
	{
		string filename = "D:/4_SCHOOL/TU Delft/Research/Projects/LocalFields/VTK and ParaView/Test Data/Torus_4k_EigFields_face_ref_"+ to_string(id+1) +".vtk";
		
		ofstream file(filename);
		if (file.is_open())
		{
			file << "# vtk DataFile Version 2.0\n";
			file << "Torus Eigenfields\n";
			file << "ASCII\n";
			file << "DATASET POLYDATA\n";

			/* PRint vertices*/
			file << "POINTS " << V.rows() << " double\n";
			for (int i = 0; i < V.rows(); i++)
			{
				file << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << "\n";
			}

			/* Print faces */
			file << "POLYGONS " << F.rows() << " " << F.rows() * (F.cols() + 1) << "\n";
			for (int i = 0; i < F.rows(); i++)
			{
				file << F.cols() << " " << F(i, 0) << " " << F(i, 1) << " " << F(i, 2) << "\n";
			}

			/* Print eigenfields*/
			Eigen::VectorXd eig3D;
			Eigen::MatrixXd EigFields;

			//EigFields = Basis.block(0, 0, Basis.rows(), eigFieldReduced2D.cols())*eigFieldReduced2D;
			EigFields = eigFieldFull2D;
			double const scale = 1.0;
			eig3D = A * EigFields.col(id);

			/* FACE-base DATA */
			file << "CELL_DATA " << F.rows() << "\n";
			file << "VECTORS EigenField_" << id+1 << " double\n";
			for (int i = 0; i < F.rows(); i++)
			{
				file << scale*eig3D(3 * i + 0) << " " << scale*eig3D(3 * i + 1) << " " << scale*eig3D(3 * i + 2) << "\n";
				cout << "face=" << i << ": " << scale*eig3D(3 * i + 0) << " " << scale*eig3D(3 * i + 1) << " " << scale*eig3D(3 * i + 2) << "\n";
			}

			/* POINT-base DATA */
			//Eigen::VectorXd VEigFields;
			//VEigFields.setZero(3 * V.rows());
			//Eigen::VectorXi VNumNeighFaces;
			//VNumNeighFaces.setZero(V.rows());
			//
			//for (int i = 0; i < F.rows(); i++)
			//{
			//	VEigFields.block(3 * F(i, 0), 0, 3, 1) += eig3D.block(3 + i, 0, 3, 1);
			//	VEigFields.block(3 * F(i, 1), 0, 3, 1) += eig3D.block(3 + i, 0, 3, 1);
			//	VEigFields.block(3 * F(i, 2), 0, 3, 1) += eig3D.block(3 + i, 0, 3, 1);
			//	//VEigFields.row(F(i, 0)) += EigFields.row(i);
			//	//VEigFields.row(F(i, 1)) += EigFields.row(i);
			//	//VEigFields.row(F(i, 2)) += EigFields.row(i);
			//	VNumNeighFaces(F(i, 0)) += 1;
			//	VNumNeighFaces(F(i, 1)) += 1;
			//	VNumNeighFaces(F(i, 2)) += 1;
			//}
			//
			//for (int i = 0; i < V.rows(); i++)
			//{
			//	VEigFields.block(3 * i, 0, 3, 1) /= (double)VNumNeighFaces(i);
			//	VEigFields.block(3 * i, 0, 3, 1) = VEigFields.block(3 * i, 0, 3, 1).normalized();
			//}
			//
			//file << "POINT_DATA " << V.rows() << "\n";
			//file << "VECTORS EigenField_" << id+1 << " double\n";
			//for (int i = 0; i < V.rows(); i++)
			//{
			//	file << scale*VEigFields(3 * i + 0) << " " << scale*VEigFields(3 * i + 1) << " " << scale*VEigFields(3 * i + 2) << "\n";
			//}

			file.close();
		}
	}
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

void VectorFields::testEnergyOfLocalPatch(igl::opengl::glfw::Viewer &viewer)
{
	// Identifying elements of local patch
	for (int i = 0; i < localPatchElements.size(); i++)
	{
		viewer.data().add_points(FC.row(localPatchElements[i]), Eigen::RowVector3d(0.0, 0.2, 0.9));
	}

	// Map from Global to Local
	vector<int> GlobToLocMap; 
	GlobToLocMap.resize(F.rows());
	for (int i = 0; i < F.rows(); i++)
	{
		GlobToLocMap[i] = -1;
	}

	int counter = 0;
	for (int face : localPatchElements) {
		GlobToLocMap[face] = counter;
		counter++;
	}

	// Construct Mass Matrix
	vector<Eigen::Triplet<double>> MTriplet;
	MTriplet.reserve(2 * localPatchElements.size());
	Eigen::SparseMatrix<double> MLocal(2 * localPatchElements.size(), 2 * localPatchElements.size());
	for (int i = 0; i < localPatchElements.size(); i++)
	{
		MTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, doubleArea(localPatchElements[i]) / 2.0));
		MTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, doubleArea(localPatchElements[i]) / 2.0));
	}
	MLocal.setFromTriplets(MTriplet.begin(), MTriplet.end());

	// Obtain the Stiffness matrix on the patch
	vector<Eigen::Triplet<double>> SFTriplet;
	SFTriplet.reserve(2 * localPatchElements.size());
	Eigen::SparseMatrix<double> SFLocal(2 * localPatchElements.size(), 2 * localPatchElements.size());
	for (int i = 0; i < localPatchElements.size(); i++)
	{
		// Diagonal elements
		int li = localPatchElements[i];
		SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, SF2D.coeff(2 * li + 0, 2 * li + 0)));
		SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 1, SF2D.coeff(2 * li + 0, 2 * li + 1)));
		SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 0, SF2D.coeff(2 * li + 1, 2 * li + 0)));
		SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, SF2D.coeff(2 * li + 1, 2 * li + 1)));

		// Non-diagonal elements
		for (int j : AdjMF2Ring[li]) {
			const int neigh = j;
			if (GlobToLocMap[neigh] >= 0) {
				int neighLoc = GlobToLocMap[neigh];
				/* Non diagonal elements*/
				SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * neighLoc + 0, SF2D.coeff(2 * li + 0, 2 * neigh + 0)));
				SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * neighLoc + 1, SF2D.coeff(2 * li + 0, 2 * neigh + 1)));
				SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * neighLoc + 0, SF2D.coeff(2 * li + 1, 2 * neigh + 0)));
				SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * neighLoc + 1, SF2D.coeff(2 * li + 1, 2 * neigh + 1)));

				/* Diagonal elements */
				//SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, -1.0 * SF2D.coeff(2 * li + 0, 2 * neigh + 0)));
				//SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 1, -1.0 * SF2D.coeff(2 * li + 0, 2 * neigh + 1)));
				//SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 0, -1.0 * SF2D.coeff(2 * li + 1, 2 * neigh + 0)));
				//SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, -1.0 * SF2D.coeff(2 * li + 1, 2 * neigh + 1)));
			}
		}		
	}

	SFLocal.setFromTriplets(SFTriplet.begin(), SFTriplet.end());
	visualizeSparseMatrixInMatlab(SFLocal);

	// Get the components of the harmonic energy on the patch
	Eigen::VectorXd localHarmEnergy(2 * localPatchElements.size());
	for (int i = 0; i < localPatchElements.size(); i++)
	{
		Eigen::Vector2d vectLoc;
		vectLoc(0) = eigFieldFull2D(2 * localPatchElements[i] + 0, 0);
		vectLoc(1) = eigFieldFull2D(2 * localPatchElements[i] + 1, 0);
		//vectLoc.normalize();
		localHarmEnergy(2 * i + 0) = vectLoc(0);
		localHarmEnergy(2 * i + 1) = vectLoc(1);
	}

	// Compute the dirichlet energy
	double energy = localHarmEnergy.transpose()*SFLocal*localHarmEnergy;
	printf("Local Energy = %.5f\n", energy);

	double totalEnergy = eigFieldFull2D.col(0).transpose() * SF2D * eigFieldFull2D.col(0);
	printf("Total Energy = %.5f\n", totalEnergy);
}