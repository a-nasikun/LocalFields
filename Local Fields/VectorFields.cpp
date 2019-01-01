#include "VectorFields.h"

/* ====================== SETTING UP MATRICES ============================*/
void VectorFields::constructGlobalMatrices()
{
	// Mass Matrices
	constructMassMatrices();
}

void VectorFields::constructMassMatrices()
{
	// Vertex-based Mass Matrices
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Constructing Mass matrices... \n";
		constructMassMatrixMV();
		constructMassMatrixMVinv();

	// Face-based Mass Matrices
		constructMassMatrixMF2D();
		constructMassMatrixMF2Dinv();
		constructMassMatrixMF3D();
		constructMassMatrixMF3Dinv();

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "..in Total of " << duration.count() << " seconds" << endl;
}

void VectorFields::constructMassMatrixMV() 
{
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "....Constructing vertex-based Mass matrix... ";

	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, MV);

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;

	//cout << MV << endl << endl;
}

void VectorFields::constructMassMatrixMVinv()
{
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "....Constructing Inverse of vertex-based Mass matrix... ";

	MVinv.resize(MV.rows(), MV.cols());
	vector<Eigen::Triplet<double>> MVTriplet;
	MVTriplet.reserve(MV.rows());

	/* Getting the sum of every non-zero elements in a row */
	for (int k = 0; k < MV.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(MV, k); it; ++it) {
			MVTriplet.push_back(Eigen::Triplet<double>(it.row(), it.col(), 1.0 / it.value()));
		}
	}
	MVinv.setFromTriplets(MVTriplet.begin(), MVTriplet.end());

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}

void VectorFields::constructMassMatrixMF2D()
{
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "....Constructing Face-based Mass matrix (2D)... ";

	igl::doublearea(V, F, doubleArea);

	MF2D.resize(2 * F.rows(), 2 * F.rows());
	vector<Eigen::Triplet<double>> MFTriplet;
	MFTriplet.reserve(2 * F.rows());

	for (int i = 0; i < F.rows(); i++) {
		double area = doubleArea(i) / 2;
		MFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, area));
		MFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, area));
	}
	MF2D.setFromTriplets(MFTriplet.begin(), MFTriplet.end());

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}

void VectorFields::constructMassMatrixMF2Dinv()
{
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "....Constructing inverse of Face-based Mass matrix (2D)... ";

	MF2Dinv.resize(MF2D.rows(), MF2D.cols());
	vector<Eigen::Triplet<double>> MFTriplet;
	MFTriplet.reserve(MF2D.rows());

	/* Getting the sum of every non-zero elements in a row */
	for (int k = 0; k < MF2D.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(MF2D, k); it; ++it) {
			MFTriplet.push_back(Eigen::Triplet<double>(it.row(), it.col(), 1.0 / it.value()));
		}
	}

	MF2Dinv.setFromTriplets(MFTriplet.begin(), MFTriplet.end());

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}

void VectorFields::constructMassMatrixMF3D()
{
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "....Constructing Face-based Mass matrix (3D)... ";


	MF3D.resize(3 * F.rows(), 3 * F.rows());
	vector<Eigen::Triplet<double>> MFTriplet;
	MFTriplet.reserve(2 * F.rows());

	for (int i = 0; i < F.rows(); i++) {
		double area = doubleArea(i) / 2;
		MFTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * i + 0, area));
		MFTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * i + 1, area));
		MFTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * i + 2, area));
	}
	MF3D.setFromTriplets(MFTriplet.begin(), MFTriplet.end());
	
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
	//cout << "The maximum element of MF3D is " << SparseMatrixMaxValue(MF3D) << endl; 
	
}

void VectorFields::constructMassMatrixMF3Dinv()
{
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "....Constructing inverse of Face-based Mass matrix (3D)... ";


	MF3Dinv.resize(MF3D.rows(), MF3D.cols());
	vector<Eigen::Triplet<double>> MFTriplet;
	MFTriplet.reserve(3 * F.rows());

	/* Getting the sum of every non-zero elements in a row */
	for (int k = 0; k < MF3D.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(MF3D, k); it; ++it) {
			MFTriplet.push_back(Eigen::Triplet<double>(it.row(), it.col(), 1.0 / it.value()));
		}
	}

	MF3Dinv.setFromTriplets(MFTriplet.begin(), MFTriplet.end());

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
	//cout << "The maximum element of MF3Dinv is " << SparseMatrixMaxValue(MF3Dinv) << endl;
}

void VectorFields::constructStiffnessMatrices()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Constructing stiffness matrices...\n";

	//constructStiffnessMatrixSV();
	constructStiffnessMatrixSF3D();
	constructStiffnessMatrixSF2D();

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "..in total of" << duration.count() << " seconds" << endl;

}

void VectorFields::constructStiffnessMatrixSV()
{
	igl::cotmatrix(V, F, SV);
	printf("Size of SV=%dx%d\n", SV.rows(), SV.cols());
	//cout << SV.block(0, 0, 5, 5) << endl;
	//cout << SV << endl;
}

void VectorFields::constructStiffnessMatrixSF2D()
{
	// Explicit construction
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;

	t1 = chrono::high_resolution_clock::now();
	cout << "....Constructing Stiffness Matrix (2D) Curl part ";
		constructStiffnessMatrixCurlPart2D();	
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;


	t1 = chrono::high_resolution_clock::now();
	cout << "....Constructing Stiffness Matrix (2D) Divergent part ";
		constructStiffnessMatrixDivPart2D();
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
		
	SF2D = LapDiv2D - LapCurl2D;

	// Implicit Construction
	//Eigen::SparseMatrix<double> GMG, JGMGJ;
	//printf("Dim check: G=%dx%d, M=%dx%d\n", GF2D.rows(), GF2D.cols(), MVinv.rows(), MVinv.cols());
	//GMG		= GF2D*MVinv*GF2D.transpose();
	//JGMGJ	= J*GMG*J;
	//SF2D	= MF2D*(GMG - JGMGJ)*MF2D;

}

void VectorFields::constructStiffnessMatrixSF3D() 
{	
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;

	t1 = chrono::high_resolution_clock::now();
	cout << "....Constructing Stiffness Matrix (3D) Curl part ";
		//constructStiffnessMatrixCurlPart3D();
		constructStiffnessMatrixCurlPart3DandCurl4F();
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1; 
	cout << "in " << duration.count() << " seconds" << endl;


	t1 = chrono::high_resolution_clock::now();
	cout << "....Constructing Stiffness Matrix (3D) Divergent part ";
		constructStiffnessMatrixDivPart3D();
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
		
	SF3D = LapDiv3D - LapCurl3D; 	
}

void VectorFields::constructStiffnessMatrixCurlPart3D()
{
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;

	Eigen::SparseMatrix<double> LapCurl3D_temp(3 * F.rows(), 3 * F.rows());
	LapCurl3D.resize(3 * F.rows(), 3 * F.rows());
	LapCurl3D.reserve(3 * F.rows() * 4);
	vector<Eigen::Triplet<double>> LTriplet;
	LTriplet.reserve(12 * LapCurl3D.rows());

	t1 = chrono::high_resolution_clock::now();
	for (int i = 0; i < F.rows(); i++) {
		double				area1 = doubleArea(i);
		Eigen::RowVector3d	n1 = NF.row(i);
		for (int j = 0; j < F.cols(); j++) {
			int				neigh = AdjMF3N(i, j);

			if (neigh > i) {
				double				area2 = doubleArea(neigh);
				double				area = area1 + area2;
				Eigen::RowVector3d	n2 = NF.row(neigh);
				Eigen::RowVector3d	n = (n1 + n2) / 2.0;
				Eigen::Vector3d		edge = V.row(EdgePairMatrix(i, 2 * j + 1)) - V.row(EdgePairMatrix(i, 2 * j));
				//edge = n.cross(edge);
				Eigen::Matrix3d		block = (-3.0 / area) * edge * edge.transpose();


				// ITS BLOCK
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * neigh + 0, block(0, 0)));	// row 1
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * neigh + 1, block(0, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * neigh + 2, block(0, 2)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * neigh + 0, block(1, 0)));	// row 2
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * neigh + 1, block(1, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * neigh + 2, block(1, 2)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * neigh + 0, block(2, 0)));	// row 3
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * neigh + 1, block(2, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * neigh + 2, block(2, 2)));

				// THE BLOCK that's the Transpose of this BLOCK
				block.transposeInPlace();
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 0, 3 * i + 0, block(0, 0)));	// row 1
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 0, 3 * i + 1, block(0, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 0, 3 * i + 2, block(0, 2)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 1, 3 * i + 0, block(1, 0)));	// row 2
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 1, 3 * i + 1, block(1, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 1, 3 * i + 2, block(1, 2)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 2, 3 * i + 0, block(2, 0)));	// row 3
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 2, 3 * i + 1, block(2, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 2, 3 * i + 2, block(2, 2)));
			}
		}
	}
	LapCurl3D_temp.setFromTriplets(LTriplet.begin(), LTriplet.end());

	/*
	Eigen::SparseMatrix<double> LapCurl3D_temp(3 * F.rows(), 3 * F.rows());
	LapCurl3D.resize(3 * F.rows(), 3 * F.rows());
	LapCurl3D.reserve(3 * F.rows() * 4);
	vector<Eigen::Triplet<double>> LTriplet;

	t1 = chrono::high_resolution_clock::now();
	for (int i = 0; i < F.rows(); i++){
		double area1 = doubleArea(i);
		for (int j = 0; j < F.cols(); j++) {
			int				neigh = AdjMF3N(i, j);

			if (neigh > i) {
				double				area2	= doubleArea(neigh);
				double				area	= area1 + area2;
				Eigen::VectorXd		edge	= V.row(EdgePairMatrix(i, 2 * j + 1)) - V.row(EdgePairMatrix(i, 2 * j));
				Eigen::MatrixXd		block	= (-3.0 / area) * edge * edge.transpose();

				// ITS BLOCK
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * neigh + 0, block(0, 0)));	// row 1
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * neigh + 1, block(0, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * neigh + 2, block(0, 2)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * neigh + 0, block(1, 0)));	// row 2
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * neigh + 1, block(1, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * neigh + 2, block(1, 2)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * neigh + 0, block(2, 0)));	// row 3
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * neigh + 1, block(2, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * neigh + 2, block(2, 2)));

				// THE BLOCK that's the Transpose of this BLOCK
				block.transposeInPlace();
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 0, 3 * i + 0, block(0, 0)));	// row 1
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 0, 3 * i + 1, block(0, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 0, 3 * i + 2, block(0, 2)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 1, 3 * i + 0, block(1, 0)));	// row 2
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 1, 3 * i + 1, block(1, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 1, 3 * i + 2, block(1, 2)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 2, 3 * i + 0, block(2, 0)));	// row 3
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 2, 3 * i + 1, block(2, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 2, 3 * i + 2, block(2, 2)));
			}
		}
	}		
	LapCurl3D_temp.setFromTriplets(LTriplet.begin(), LTriplet.end());	*/
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	printf("To set-up non-diagonal elements =%.4f seconds\n", duration.count());

	t1 = chrono::high_resolution_clock::now();
	// INSERTING the DIAGONAL Elements
	for (int i = 0; i < F.rows(); i++) {
		for (int j = 0; j < F.cols(); j++) {		// Working in 1 block (3x3 size) 
			for (int k = 0; k < F.cols(); k++) {	// Inserting in column order
				// Obtaining an element from the other 3 neighbor (same location of 3 different blocks)
				double v0, v1, v2;				
				v0 = LapCurl3D_temp.coeff(3 * i + k, 3 * AdjMF3N(i, 0) + j);
				v1 = LapCurl3D_temp.coeff(3 * i + k, 3 * AdjMF3N(i, 1) + j);
				v2 = LapCurl3D_temp.coeff(3 * i + k, 3 * AdjMF3N(i, 2) + j);

				double value = -(v0+v1+v2);
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + k, 3 * i + j, value));	// row 1
			}
		}
	}
	LapCurl3D.setFromTriplets(LTriplet.begin(), LTriplet.end());	
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	printf("To set-up diagonal block elements =%.4f seconds\n", duration.count());
}

void VectorFields::constructStiffnessMatrixCurlPart3DandCurl4F()
{
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;

	Eigen::SparseMatrix<double> LapCurl3D_temp(3 * F.rows(), 3 * F.rows());
	LapCurl3D.resize(3 * F.rows(), 3 * F.rows());
	LapCurl3D.reserve(3 * F.rows() * 4);
	vector<Eigen::Triplet<double>> LTriplet;
	LTriplet.reserve(12 * LapCurl3D.rows());

	// For Curl3DPacked
	Curl3DPacked.resize(3 * F.rows(), 4 * F.cols());

	t1 = chrono::high_resolution_clock::now();
	for (int i = 0; i < F.rows(); i++) {
		double				area1 = doubleArea(i);
		Eigen::RowVector3d	n1 = NF.row(i);
		for (int j = 0; j < F.cols(); j++) {
			int				neigh = AdjMF3N(i, j);

			double				area2 = doubleArea(neigh);
			double				area = area1 + area2;
			Eigen::RowVector3d	n2 = NF.row(neigh);
			Eigen::RowVector3d	n = (n1 + n2) / 2.0;
			Eigen::Vector3d		edge = V.row(EdgePairMatrix(i, 2 * j + 1)) - V.row(EdgePairMatrix(i, 2 * j));
			//edge = n.cross(edge);
			Eigen::Matrix3d		block = (-3.0 / area) * edge * edge.transpose();

			// ITS BLOCK
			LTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * neigh + 0, block(0, 0)));	// row 1
			LTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * neigh + 1, block(0, 1)));
			LTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * neigh + 2, block(0, 2)));
			LTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * neigh + 0, block(1, 0)));	// row 2
			LTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * neigh + 1, block(1, 1)));
			LTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * neigh + 2, block(1, 2)));
			LTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * neigh + 0, block(2, 0)));	// row 3
			LTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * neigh + 1, block(2, 1)));
			LTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * neigh + 2, block(2, 2)));

			// Structure of the packed data
			Curl3DPacked(3 * i + 0, 3 * j + 0) = block(0, 0);
			Curl3DPacked(3 * i + 0, 3 * j + 1) = block(0, 1);
			Curl3DPacked(3 * i + 0, 3 * j + 2) = block(0, 2);
			Curl3DPacked(3 * i + 1, 3 * j + 0) = block(1, 0);
			Curl3DPacked(3 * i + 1, 3 * j + 1) = block(1, 1);
			Curl3DPacked(3 * i + 1, 3 * j + 2) = block(1, 2);
			Curl3DPacked(3 * i + 2, 3 * j + 0) = block(2, 0);
			Curl3DPacked(3 * i + 2, 3 * j + 1) = block(2, 1);
			Curl3DPacked(3 * i + 2, 3 * j + 2) = block(2, 2);
		}
	}
	LapCurl3D_temp.setFromTriplets(LTriplet.begin(), LTriplet.end());

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	//printf("To set-up non-diagonal elements =%.4f seconds\n", duration.count());

	// Getting the diagonal elements, from and to packed ata
	for (int i = 0; i < F.rows(); i++) {
		for (int j = 0; j < F.cols(); j++) {
			for (int k = 0; k < F.cols(); k++) {
				Curl3DPacked(3 * i + j, 3 * F.cols() + k) = -(Curl3DPacked(3 * i + j, 0 * F.cols() + k) + Curl3DPacked(3 * i + j, 1 * F.cols() + k) + Curl3DPacked(3 * i + j, 2 * F.cols() + k));
			}
		}
	}

	t1 = chrono::high_resolution_clock::now();
	// INSERTING the DIAGONAL Elements
	for (int i = 0; i < F.rows(); i++) {
		for (int j = 0; j < F.cols(); j++) {		// Working in 1 block (3x3 size) 
			for (int k = 0; k < F.cols(); k++) {	// Inserting in column order => same column, increasing row => k=row
													// Obtaining an element from the other 3 neighbor (same location of 3 different blocks)
				double v0, v1, v2;
				//v0 = LapCurl3D_temp.coeff(3 * i + k, 3 * AdjMF3N(i, 0) + j);
				//v1 = LapCurl3D_temp.coeff(3 * i + k, 3 * AdjMF3N(i, 1) + j);
				//v2 = LapCurl3D_temp.coeff(3 * i + k, 3 * AdjMF3N(i, 2) + j);

				//double value = -(v0 + v1 + v2);
				double value = Curl3DPacked(3 * i + k, 3 * F.cols() + j);
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + k, 3 * i + j, value));	// row 1
			}
		}
	}
	LapCurl3D.setFromTriplets(LTriplet.begin(), LTriplet.end());
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	//printf("To set-up diagonal block elements =%.4f seconds\n", duration.count());
}

void VectorFields::constructStiffnessMatrixCurlPart2D()
{
	LapCurl2D = A.transpose() * LapCurl3D * A;
}

void VectorFields::constructStiffnessMatrixCurlPart2D_Direct()
{
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;

	Eigen::SparseMatrix<double> LapCurl2D_temp(2 * F.rows(), 2 * F.rows());
	LapCurl2D.resize(2 * F.rows(), 2 * F.rows());
	LapCurl2D.reserve(2 * F.rows() * 4);
	vector<Eigen::Triplet<double>> LTriplet;

	t1 = chrono::high_resolution_clock::now();
	for (int i = 0; i < F.rows(); i++) {
		double				area1 = doubleArea(i);
		for (int j = 0; j < F.cols(); j++) {
			int				neigh = AdjMF3N(i, j);

			if (neigh > i) {
				double				area2 = doubleArea(neigh);
				double				area = area1 + area2;
				Eigen::Vector3d		edge = V.row(EdgePairMatrix(i, 2 * j + 1)) - V.row(EdgePairMatrix(i, 2 * j));
				//edge = n.cross(edge);
				Eigen::Matrix3d		block = (-3.0 / area) * edge * edge.transpose();
				Eigen::MatrixXd		Aa(3, 2);
				Aa = A.block(3 * i, 2 * i, 3, 2);
				Eigen::Matrix2d     block2D = Aa.transpose() * block * Aa; 

				// ITS BLOCK
				LTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * neigh + 0, block2D(0, 0)));	// row 1
				LTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * neigh + 1, block2D(0, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * neigh + 0, block2D(1, 0)));	// row 2
				LTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * neigh + 1, block2D(1, 1)));
				// THE BLOCK that's the Transpose of this BLOCK
				block2D.transposeInPlace();
				LTriplet.push_back(Eigen::Triplet<double>(2 * neigh + 0, 2 * i + 0, block2D(0, 0)));	// row 1
				LTriplet.push_back(Eigen::Triplet<double>(2 * neigh + 0, 2 * i + 1, block2D(0, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(2 * neigh + 1, 2 * i + 0, block2D(1, 0)));	// row 2
				LTriplet.push_back(Eigen::Triplet<double>(2 * neigh + 1, 2 * i + 1, block2D(1, 1)));
			}
		}
	}
	LapCurl2D_temp.setFromTriplets(LTriplet.begin(), LTriplet.end());

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	//printf("To set-up non-diagonal elements =%.4f seconds\n", duration.count());

	t1 = chrono::high_resolution_clock::now();
	// INSERTING the DIAGONAL Elements
	for (int i = 0; i < F.rows(); i++) {
		for (int j = 0; j < 2; j++) {		// Working in 1 block (2x2 size of local frame) 
			for (int k = 0; k < 2; k++) {	// Inserting in column order
													// Obtaining an element from the other 3 neighbor (same location of 3 different blocks)
				double v0, v1, v2;
				v0 = LapCurl2D_temp.coeff(2 * i + k, 2 * AdjMF3N(i, 0) + j);
				v1 = LapCurl2D_temp.coeff(2 * i + k, 2 * AdjMF3N(i, 1) + j);
				v2 = LapCurl2D_temp.coeff(2 * i + k, 2 * AdjMF3N(i, 2) + j);

				double value = -(v0 + v1 + v2);
				LTriplet.push_back(Eigen::Triplet<double>(2 * i + k, 2 * i + j, value));	// row 1
			}
		}
	}
	LapCurl2D.setFromTriplets(LTriplet.begin(), LTriplet.end());
	
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	//printf("To set-up diagonal block elements =%.4f seconds\n", duration.count());
}

void VectorFields::constructStiffnessMatrixDivPart3D()
{
	// IMPLICIT Construction for Divergent Part
	//constructStiffnessMatrixDivPart3D_Implicit();

	// EXPLICIT Construction
	//constructStiffnessMatrixDivPart3D_Explicit();
	constructStiffnessMatrixDivPart3DandDiv4F_Explicit();
}

void VectorFields::constructStiffnessMatrixDivPart3D_Implicit()
{
	LapDiv3D = MF3D*GF3D*MVinv*GF3D.transpose()*MF3D;
}

void VectorFields::constructStiffnessMatrixDivPart3D_Explicit()
{
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;

	Eigen::SparseMatrix<double> LapDiv3D_temp(3 * F.rows(), 3 * F.rows());
	LapDiv3D.resize(3 * F.rows(), 3 * F.rows());
	LapDiv3D.reserve(3 * F.rows() * 4);
	vector<Eigen::Triplet<double>> LTriplet;
	LTriplet.reserve(12 * 3 * F.rows());

	t1 = chrono::high_resolution_clock::now();
	for (int i = 0; i < F.rows(); i++) {
		double				area1 = doubleArea(i);
		Eigen::RowVector3d	n1 = NF.row(i);
		for (int j = 0; j < F.cols(); j++) {
			int				neigh = AdjMF3N(i, j);

			if (neigh > i) {
				double				area2 = doubleArea(neigh);
				double				area = area1 + area2;
				Eigen::RowVector3d	n2 = NF.row(neigh);
				Eigen::RowVector3d	n = (n1 + n2) / 2.0;
				Eigen::Vector3d		edge = V.row(EdgePairMatrix(i, 2 * j + 1)) - V.row(EdgePairMatrix(i, 2 * j));
				
				edge = n.cross(edge);
				Eigen::Matrix3d		block = (3.0 / area) * edge * edge.transpose();
							

				// ITS BLOCK
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * neigh + 0, block(0, 0)));	// row 1
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * neigh + 1, block(0, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * neigh + 2, block(0, 2)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * neigh + 0, block(1, 0)));	// row 2
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * neigh + 1, block(1, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * neigh + 2, block(1, 2)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * neigh + 0, block(2, 0)));	// row 3
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * neigh + 1, block(2, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * neigh + 2, block(2, 2)));

				// THE BLOCK that's the Transpose of this BLOCK
				block.transposeInPlace();
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 0, 3 * i + 0, block(0, 0)));	// row 1
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 0, 3 * i + 1, block(0, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 0, 3 * i + 2, block(0, 2)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 1, 3 * i + 0, block(1, 0)));	// row 2
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 1, 3 * i + 1, block(1, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 1, 3 * i + 2, block(1, 2)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 2, 3 * i + 0, block(2, 0)));	// row 3
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 2, 3 * i + 1, block(2, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 2, 3 * i + 2, block(2, 2)));
			}
		}
	}
	LapDiv3D_temp.setFromTriplets(LTriplet.begin(), LTriplet.end());
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	printf("To set-up non-diagonal elements =%.4f seconds\n", duration.count());

	t1 = chrono::high_resolution_clock::now();
	// INSERTING the DIAGONAL Elements
	for (int i = 0; i < F.rows(); i++) {
		//LapDiv3D.reserve(F.cols()*F.cols());
		for (int j = 0; j < F.cols(); j++) {		// Working in 1 block (3x3 size) 
			for (int k = 0; k < F.cols(); k++) {	// Inserting in column order
													// Obtaining an element from the other 3 neighbor (same location of 3 different blocks)
				double v0, v1, v2;
				v0 = LapDiv3D_temp.coeff(3 * i + k, 3 * AdjMF3N(i, 0) + j);
				v1 = LapDiv3D_temp.coeff(3 * i + k, 3 * AdjMF3N(i, 1) + j);
				v2 = LapDiv3D_temp.coeff(3 * i + k, 3 * AdjMF3N(i, 2) + j);

				double value = -(v0 + v1 + v2);
				//LapDiv3D.insert(3 * i + k, 3 * i + j) = value;
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + k, 3 * i + j, value));
			}
		}
	}

	LapDiv3D.setFromTriplets(LTriplet.begin(), LTriplet.end());
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	printf("To set-up diagonal block elements =%.4f seconds\n", duration.count());
}

void VectorFields::constructStiffnessMatrixDivPart3DandDiv4F_Explicit()
{
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;

	Eigen::SparseMatrix<double> LapDiv3D_temp(3 * F.rows(), 3 * F.rows());
	LapDiv3D.resize(3 * F.rows(), 3 * F.rows());
	LapDiv3D.reserve(3 * F.rows() * 4);
	vector<Eigen::Triplet<double>> LTriplet;
	LTriplet.reserve(12 * 3 * F.rows());

	// For Curl3DPacked
	Div3DPacked.resize(3 * F.rows(), 4 * F.cols());

	t1 = chrono::high_resolution_clock::now();
	for (int i = 0; i < F.rows(); i++) {
		double				area1 = doubleArea(i);
		Eigen::RowVector3d	n1 = NF.row(i);
		for (int j = 0; j < F.cols(); j++) {
			int				neigh = AdjMF3N(i, j);

			double				area2 = doubleArea(neigh);
			double				area = area1 + area2;
			Eigen::RowVector3d	n2 = NF.row(neigh);
			Eigen::RowVector3d	n = (n1 + n2) / 2.0;
			Eigen::Vector3d		edge = V.row(EdgePairMatrix(i, 2 * j + 1)) - V.row(EdgePairMatrix(i, 2 * j));

			edge = n.cross(edge);
			Eigen::Matrix3d		block = (3.0 / area) * edge * edge.transpose();


			// ITS BLOCK
			LTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * neigh + 0, block(0, 0)));	// row 1
			LTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * neigh + 1, block(0, 1)));
			LTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * neigh + 2, block(0, 2)));
			LTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * neigh + 0, block(1, 0)));	// row 2
			LTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * neigh + 1, block(1, 1)));
			LTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * neigh + 2, block(1, 2)));
			LTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * neigh + 0, block(2, 0)));	// row 3
			LTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * neigh + 1, block(2, 1)));
			LTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * neigh + 2, block(2, 2)));		

			// Structure of the packed data
			Div3DPacked(3 * i + 0, 3 * j + 0) = block(0, 0);
			Div3DPacked(3 * i + 0, 3 * j + 1) = block(0, 1);
			Div3DPacked(3 * i + 0, 3 * j + 2) = block(0, 2);
			Div3DPacked(3 * i + 1, 3 * j + 0) = block(1, 0);
			Div3DPacked(3 * i + 1, 3 * j + 1) = block(1, 1);
			Div3DPacked(3 * i + 1, 3 * j + 2) = block(1, 2);
			Div3DPacked(3 * i + 2, 3 * j + 0) = block(2, 0);
			Div3DPacked(3 * i + 2, 3 * j + 1) = block(2, 1);
			Div3DPacked(3 * i + 2, 3 * j + 2) = block(2, 2);
		}
	}
	LapDiv3D_temp.setFromTriplets(LTriplet.begin(), LTriplet.end());
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	//printf("To set-up non-diagonal elements =%.4f seconds\n", duration.count());

	// Getting the diagonal elements, from and to packed ata
	for (int i = 0; i < F.rows(); i++) {
		for (int j = 0; j < F.cols(); j++) {
			for (int k = 0; k < F.cols(); k++) {
				Div3DPacked(3 * i + j, 3 * F.cols() + k) = -(Div3DPacked(3 * i + j, 0 * F.cols() + k) + Div3DPacked(3 * i + j, 1 * F.cols() + k) + Div3DPacked(3 * i + j, 2 * F.cols() + k));
			}
		}
	}

	t1 = chrono::high_resolution_clock::now();
	// INSERTING the DIAGONAL Elements
	for (int i = 0; i < F.rows(); i++) {
		//LapDiv3D.reserve(F.cols()*F.cols());
		for (int j = 0; j < F.cols(); j++) {		// Working in 1 block (3x3 size) 
			for (int k = 0; k < F.cols(); k++) {	// Inserting in column order
													// Obtaining an element from the other 3 neighbor (same location of 3 different blocks)
				double v0, v1, v2;
				//v0 = LapDiv3D_temp.coeff(3 * i + k, 3 * AdjMF3N(i, 0) + j);
				//v1 = LapDiv3D_temp.coeff(3 * i + k, 3 * AdjMF3N(i, 1) + j);
				//v2 = LapDiv3D_temp.coeff(3 * i + k, 3 * AdjMF3N(i, 2) + j);

				//double value = -(v0 + v1 + v2);
				double value = Div3DPacked(3 * i + k, 3 * F.cols() + j);
				//LapDiv3D.insert(3 * i + k, 3 * i + j) = value;
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + k, 3 * i + j, value));
			}
		}
	}

	LapDiv3D.setFromTriplets(LTriplet.begin(), LTriplet.end());
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	//printf("To set-up diagonal block elements =%.4f seconds\n", duration.count());
}

void VectorFields::constructStiffnessMatrixDivPart2D()
{
	// IMPLICT
	//LapDiv2D = MF2D*GF2D*MVinv*GF2D.transpose()*MF2D;

	// EXPLICIT
	LapDiv2D = A.transpose() * LapDiv3D * A;
}

void VectorFields::constructStiffnessMatrixDivPart2D_Direct()
{
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;

	Eigen::SparseMatrix<double> LapDiv2D_temp(2 * F.rows(), 2 * F.rows());
	LapDiv2D.resize(2 * F.rows(), 2 * F.rows());
	LapDiv2D.reserve(2 * F.rows() * 4);
	vector<Eigen::Triplet<double>> LTriplet;

	t1 = chrono::high_resolution_clock::now();
	for (int i = 0; i < F.rows(); i++) {
		double				area1 = doubleArea(i);
		Eigen::RowVector3d	n1 = NF.row(i);
		for (int j = 0; j < F.cols(); j++) {
			int				neigh = AdjMF3N(i, j);

			if (neigh > i) {
				double				area2 = doubleArea(neigh);
				double				area = area1 + area2;
				Eigen::RowVector3d	n2 = NF.row(neigh);
				Eigen::RowVector3d	n = (n1 + n2) / 2.0;
				Eigen::Vector3d		edge = V.row(EdgePairMatrix(i, 2 * j + 1)) - V.row(EdgePairMatrix(i, 2 * j));

				edge = n.cross(edge);
				Eigen::Matrix3d		block = (3.0 / area) * edge * edge.transpose();
				Eigen::MatrixXd		Aa = A.block(3 * i, 2 * i, 3, 2);
				Eigen::Matrix2d     block2D = Aa.transpose() * block * Aa;

				// ITS BLOCK
				LTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * neigh + 0, block2D(0, 0)));	// row 1
				LTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * neigh + 1, block2D(0, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * neigh + 0, block2D(1, 0)));	// row 2
				LTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * neigh + 1, block2D(1, 1)));

				// THE BLOCK that's the Transpose of this BLOCK
				block2D.transposeInPlace();
				LTriplet.push_back(Eigen::Triplet<double>(2 * neigh + 0, 2 * i + 0, block2D(0, 0)));	// row 1
				LTriplet.push_back(Eigen::Triplet<double>(2 * neigh + 0, 2 * i + 1, block2D(0, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(2 * neigh + 1, 2 * i + 0, block2D(1, 0)));	// row 2
				LTriplet.push_back(Eigen::Triplet<double>(2 * neigh + 1, 2 * i + 1, block2D(1, 1)));
			}
		}
	}
	LapDiv2D_temp.setFromTriplets(LTriplet.begin(), LTriplet.end());
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	printf("To set-up non-diagonal elements =%.4f seconds\n", duration.count());

	t1 = chrono::high_resolution_clock::now();
	// INSERTING the DIAGONAL Elements
	for (int i = 0; i < F.rows(); i++) {
		//LapDiv3D.reserve(F.cols()*F.cols());
		for (int j = 0; j < 2; j++) {		// Working in 1 block (2x2 size of local frame) 
			for (int k = 0; k < 2; k++) {	// Inserting in column order
													// Obtaining an element from the other 3 neighbor (same location of 3 different blocks)
				double v0, v1, v2;
				v0 = LapDiv2D_temp.coeff(2 * i + k, 2 * AdjMF3N(i, 0) + j);
				v1 = LapDiv2D_temp.coeff(2 * i + k, 2 * AdjMF3N(i, 1) + j);
				v2 = LapDiv2D_temp.coeff(2 * i + k, 2 * AdjMF3N(i, 2) + j);

				double value = -(v0 + v1 + v2);
				//LapDiv3D.insert(3 * i + k, 3 * i + j) = value;
				LTriplet.push_back(Eigen::Triplet<double>(2 * i + k, 2 * i + j, value));
			}
		}
	}

	LapDiv2D.setFromTriplets(LTriplet.begin(), LTriplet.end());
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	printf("To set-up diagonal block elements =%.4f seconds\n", duration.count());
}

void VectorFields::constructSF2DPacked() 
{
	SF2DPacked.resize(2 * F.rows(), 4 * 2); // 2*|F| by 3 neighbors (+1, the total)

	for (int i = 0; i < F.rows(); i++) {
		for (int j = 0; j < 4; j++) {
			Eigen::Matrix3d S3 = Div3DPacked.block(3 * i, j*F.cols(), 3, 3) - Curl3DPacked.block(3 * i, j*F.cols(), 3, 3);
			Eigen::Matrix3d Aa = A.block(3 * i, 2 * i, 3, 2);
			SF2DPacked.block(2 * i, j * 2, 2, 2) = Aa.transpose() * S3 * Aa; 
		}
	}
}
void VectorFields::constructGradient3D()
{
	// Construct Gradient of World-/Global-Coordidate
	igl::grad(V, F, GF3D);
	rearrangeGradient3D();
	//visualizeSparseMatrixInMatlab(GF3D);
	//cout << GF3D.block(0, 0, 7, GF3D.cols()) << endl;
}

void VectorFields::constructGradient2D()
{
	GF2D = A.transpose()*GF3D;
	//visualizeSparseMatrixInMatlab(GF2D);
}

void VectorFields::computeDivergent3D()
{
	Div3D = -MVinv*(GF3D.transpose()*MF3D);
}

void VectorFields::computeDivergent2D()
{
	Div2D = -MVinv*(GF2D.transpose()*MF2D);
}

void VectorFields::computeCurl2D()
{
	Curl2D = MVinv*GF2D.transpose()*J*MF2D;
	printf("Curl(2D)=%dx%d with %d non-zero entries.\n", Curl2D.rows(), Curl2D.cols(), Curl2D.nonZeros());
}

void VectorFields::computeCurl3D()
{
	if (Curl2D.nonZeros() == 0) {
		computeCurl2D();
		printf("Curl(2D)=%dx%d with %d non-zero entries.\n", Curl2D.rows(), Curl2D.cols(), Curl2D.nonZeros());
		Curl3D = Curl2D * A.transpose();
	}
	else {
		Curl3D = Curl2D * A.transpose();
	}
	
}

void VectorFields::rearrangeGradient3D()
{
	//MFinv.resize(MF.rows(), MF.cols());
	vector<Eigen::Triplet<double>> GTriplet;
	GTriplet.reserve(GF3D.nonZeros());
	int nRows = F.rows(), nCols = GF3D.cols();

	/* Getting the sum of every non-zero elements in a row */
	for (int k = 0; k < GF3D.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(GF3D, k); it; ++it) {
			if (it.row() < nRows) {					// values for x-coordinates
				GTriplet.push_back(Eigen::Triplet<double>(3 * it.row(), it.col(), it.value()));
			}
			else if (it.row() < 2 * nRows) {		// values for y-coordinates
				GTriplet.push_back(Eigen::Triplet<double>(3 * (it.row() - nRows) + 1, it.col(), it.value()));
			}
			else {									// values for z-coordinates
				GTriplet.push_back(Eigen::Triplet<double>(3 * (it.row() - 2 * nRows) + 2, it.col(), it.value()));
			}
		}
	}
	GF3D.resize(0, 0);
	GF3D.resize(3 * nRows, nCols);
	GF3D.setFromTriplets(GTriplet.begin(), GTriplet.end());
}

void VectorFields::constructRotationMatrix()
{
	J.resize(2 * F.rows(), 2 * F.rows());
	vector<Eigen::Triplet<double>> JTriplet;
	JTriplet.reserve(2 * 2*F.rows());
	const double cosT = 0.0, sinT = 1.0;

	for (int i = 0; i < F.rows(); i++) {
		// Constructing the triplet
		JTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, cosT));		// column 1
		JTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 0, sinT));
		JTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 1, -sinT));	// column 2
		JTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, cosT));
	}

	J.setFromTriplets(JTriplet.begin(), JTriplet.end());
	//cout << "Matrix J" << endl << J.block(0, 0, 10, 10) << endl << endl;
	//visualizeSparseMatrixInMatlab(J);
}

void VectorFields::constructMappingMatrix()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Constructing Mapping matrices (Global/World-Coord to Local Frame)... ";


	A.resize(3 * F.rows(), 2 * F.rows());
	vector<Eigen::Triplet<double>> ATriplet;
	ATriplet.reserve(3 * 2 * F.rows());
	Eigen::Vector3d e, f, n;

	for (int i = 0; i < F.rows(); i++) {
		e = V.row(F(i, 1)) - V.row(F(i, 0));
		e.normalize();

		n = NF.row(i);
		n.normalize();

		f = n.cross(e);
		f.normalize();

		ATriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 2 * i + 0, e(0)));
		ATriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 2 * i + 0, e(1)));
		ATriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 2 * i + 0, e(2)));
		ATriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 2 * i + 1, f(0)));
		ATriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 2 * i + 1, f(1)));
		ATriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 2 * i + 1, f(2)));
	}

	A.setFromTriplets(ATriplet.begin(), ATriplet.end());

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}

void VectorFields::constructMatrixB()
{
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;

	cout << "> Constructing Bi-Laplacian Matrix B... ";
	t1 = chrono::high_resolution_clock::now();
	B2D = SF2D * MF2Dinv * SF2D;
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;		
	cout << "in " << duration.count() << " seconds" << endl; 
}

void VectorFields::constructConstraints()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Constructing constraints... ";

	//construct1CentralConstraint();
	//constructRingConstraints();
	//constructSpecifiedConstraints();
	
	constructSingularities();
	constructSpecifiedConstraintsWithSingularities();

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;

	// Information about constraints
	printf("....Num of Constraints = %d\n", globalConstraints.size());
	printf("....Matrix C size = %dx%d \n", C.rows(), C.cols());
	printf("....Matrix c size = %dx%d \n", c.rows(), c.cols());
}

void VectorFields::construct1CentralConstraint()
{
	vector<Eigen::Triplet<double>>	CTriplet;
	CTriplet.reserve(2);
	const int constNum = 1;
	//srand(time(NULL));
	//const int centralID = rand()%F.rows(); 
	const int centralID = *(NeighRing[0].begin());
	cout << "Face ID is : " << centralID << endl;

	// Setting up matrix C
	C.resize(2 * constNum, B2D.cols());

	CTriplet.push_back(Eigen::Triplet<double>(0, 2 * centralID + 0, 1.0));
	CTriplet.push_back(Eigen::Triplet<double>(1, 2 * centralID + 1, 1.0));

	C.setFromTriplets(CTriplet.begin(), CTriplet.end());

	// Setting up vector c (There are 2 vector c)
	c.resize(2 * constNum, 2);
	c.col(0) << 1.0, 0.0;
	c.col(1) << 0.0, 1.0;
}

void VectorFields::constructRingConstraints()
{
	// Define necessary data/variables
	vector<Eigen::Triplet<double>>	CTriplet;
	
	const int outerRingID = 9;
	const int outerBoundaryID = min(outerRingID, (int) NeighRing.size()-3);
	const int constNum = 1 + (int) NeighRing[outerBoundaryID+1].size() + (int) NeighRing[outerBoundaryID + 2].size();
	const int centralID = *(NeighRing[0].begin());
	int counter = 0;

	// Setting up matrix C
	C.resize(2 * constNum, B2D.cols());
	CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * centralID + 0, 1.0));
	CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * centralID + 1, 1.0));
	for (int i = outerBoundaryID + 1; i <= outerBoundaryID + 2; i++) {
		for (std::set<int, double>::iterator it = NeighRing[i].begin(); it != NeighRing[i].end(); ++it) {
			CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * (*it) + 0, 1.0));
			CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * (*it) + 1, 1.0));
		}
	}
	C.setFromTriplets(CTriplet.begin(), CTriplet.end());

	// Setting up vector c (There are 2 vector c)
	c.resize(2 * constNum, 2);
	Eigen::VectorXd zeroElements(2 * constNum - 2);
	for (int i = 0; i < zeroElements.size(); i++) zeroElements(i) = 0.0;
	c.col(0) << 1.0, 0.0, zeroElements;
	c.col(1) << 0.0, 1.0, zeroElements;
}

void VectorFields::constructSpecifiedConstraints()
{
	// Define the constraints
	const int numConstraints = 2;
	set<int> constraints;
	//vector<int> globalConstraints(numConstraints);
	globalConstraints.resize(numConstraints);
	Eigen::VectorXd D;
	D.resize(F.rows());

	// Initialize the value of D
	for (int i = 0; i < F.rows(); i++) {
		D(i) = numeric_limits<double>::infinity();
	}

	constraints.insert(0);
	int curPoint = 0;

	do {
		Eigen::VectorXi::Index maxIndex;
		computeDijkstraDistanceFaceForSampling(curPoint, D);
		D.maxCoeff(&maxIndex);
		constraints.insert(maxIndex);
		curPoint = maxIndex;
	} while (constraints.size() <= numConstraints);

	int counter1 = 0;
	for (int i : constraints) {
		globalConstraints[counter1++] = i;
	}
	printf("Constraints = %d\n", globalConstraints.size());

	// Setting up matrix C
	Eigen::SparseMatrix<double> CTemp;
	vector<Eigen::Triplet<double>> CTriplet;
	CTriplet.reserve(2 * globalConstraints.size());

	int counter = 0;
	for (int i = 0; i < globalConstraints.size(); i++) {
		CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * globalConstraints[i] + 0, 1.0));
		CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * globalConstraints[i] + 1, 1.0));
	}
	C.resize(2 * globalConstraints.size(), B2D.rows());
	C.setFromTriplets(CTriplet.begin(), CTriplet.end());
	printf("Cp=%dx%d\n", C.rows(), C.cols());


	// Setting up vector c (There are 2 vector c)
	srand(time(NULL));
	c.resize(2 * globalConstraints.size(), 2);
	for (int i = 0; i < globalConstraints.size(); i++) {
		c(2 * i + 0, 0) = sqrt(2.0);
		c(2 * i + 1, 0) = sqrt(2.0);
		c(2 * i + 0, 1) = sqrt(2.0);
		c(2 * i + 1, 1) = sqrt(2.0);
	}
	printf("cBar=%dx%d\n", c.rows(), c.cols());
}

void VectorFields::constructSingularities()
{
	const int NUM_SINGS = 2;
	singularities.resize(NUM_SINGS);
	SingNeighCC.resize(NUM_SINGS);

	time_t t;
	srand((unsigned)time(&t));
	//srand(time(NULL));
	for (int id = 0; id < NUM_SINGS; id++) {
		const int SingLocation = rand() % V.rows();
		//const int SingLocation = id * (int)(V.rows() / NUM_SINGS) + 5;
		singularities[id] = SingLocation;
		const int SingNeighNum = VFNeighFull[SingLocation].size();
		const int firstNeigh = VFNeighFull[SingLocation].begin()->fId;

		SingNeighCC[id].resize(SingNeighNum);
		SingNeighCC[id][0] = firstNeigh;
		int curNeigh = firstNeigh;
		int vertex1 = SingLocation;

		//printf("Vertex %d has %d neighbors.\n", SingLocation, SingNeighNum);
		
		for (int i2 = 1; i2<SingNeighNum; i2++) {
			//printf("number of neighbors of vertex %d = %d (%d)\n", SingLocation, i2 /*SingNeighCC[id].size()*/, VFNeighFull[SingLocation].size());
			int vertex2;
			for (int i = 0; i < F.cols(); i++) {
				if (F(curNeigh, i%F.cols()) == vertex1) {
					vertex2 = F(curNeigh, (i + F.cols() - 1) % F.cols());
				}
			}
			for (std::set<VtoFPair>::iterator it1 = next(VFNeighFull[SingLocation].begin(), 1); it1 != VFNeighFull[SingLocation].end(); ++it1) {
				for (int i = 0; i < F.cols(); i++) {
					if (F(it1->fId, i) == vertex1 && F(it1->fId, (i + 1) % F.cols()) == vertex2) {
						SingNeighCC[id][i2] = it1->fId;
						//printf("     Inserting %d as neighbor\n", it1->fId);
						curNeigh = it1->fId;
					}
				}
			}
		}
	}
}

void VectorFields::constructSpecifiedConstraintsWithSingularities()
{
	// Define the constraints
	const int numConstraints = 30;
	set<int> constraints;

	globalConstraints.resize(numConstraints);
	Eigen::VectorXd D;
	D.resize(F.rows());

	// Initialize the value of D
	for (int i = 0; i < F.rows(); i++) {
		D(i) = numeric_limits<double>::infinity();
	}

	constraints.insert(0);
	int curPoint = 0;

	// Creating random constraints
	do {
		Eigen::VectorXi::Index maxIndex;
		computeDijkstraDistanceFaceForSampling(curPoint, D);
		D.maxCoeff(&maxIndex);
		constraints.insert(maxIndex);
		curPoint = maxIndex;
	} while (constraints.size() <= numConstraints);

	int counter1 = 0;
	for (int i : constraints) {
		globalConstraints[counter1++] = i;
	}
	//printf("Constraints = %d\n", globalConstraints.size());

	//constructSingularities();

	int numSingConstraints = 0;
	for (int i = 0; i < SingNeighCC.size(); i++) {
		// Use only n-1 neighboring faces as constraints
		for (int j = 0; j < (SingNeighCC[i].size()-1); j++) {
			numSingConstraints++;
		}
	}

	// Setting up matrix C and vector c
	c.resize(2 * (globalConstraints.size()+numSingConstraints), 2);
	//printf("cBar=%dx%d\n", c.rows(), c.cols());


	// HARD CONSTRAINTS
	Eigen::SparseMatrix<double> CTemp;
	vector<Eigen::Triplet<double>> CTriplet;
	CTriplet.reserve(2 * globalConstraints.size() + 2 * 3 * 7 * SingNeighCC.size());
	int counter = 0;
	for (int i = 0; i < globalConstraints.size(); i++) {
		// Matrix C
		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * globalConstraints[i] + 0, 1.0));
		c(counter, 0) = sqrt(2.0);
		c(counter++, 1) = sqrt(2.0);
		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * globalConstraints[i] + 1, 1.0));
		c(counter, 0) = sqrt(2.0);
		c(counter++, 1) = sqrt(2.0);
	}

	
	// SINGULARITIES CONSTRAINTS
	for (int id = 0; id < SingNeighCC.size(); id++) {
		const double rotAngle = 2 * M_PI / (double)SingNeighCC[id].size();
		const double cosA = cos(rotAngle);
		const double sinA = sin(rotAngle);
		for (int i = 0; i < (SingNeighCC[id].size() - 1); i++) {
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, cosA));
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, -sinA));
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, -1.0));
			c(counter, 0) = 0.0;
			c(counter, 1) = 0.0;
			counter++;

			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, sinA));
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, cosA));
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, -1.0));
			c(counter, 0) = 0.0;
			c(counter, 1) = 0.0;
			counter++;
		}
	}

	C.resize(2 * (globalConstraints.size()+numSingConstraints), B2D.rows());
	C.setFromTriplets(CTriplet.begin(), CTriplet.end());
	//printf("Cp=%dx%d\n", C.rows(), C.cols());

	// Setting up vector c (There are 2 vector c)	
	
}

void VectorFields::setupGlobalProblem()
{	
	constructConstraints();	
	setupRHSGlobalProblemMapped();
	setupLHSGlobalProblemMapped();
	//setupRHSGlobalProblem();
	//setupLHSGlobalProblem();	

	//solveGlobalSystem();
	solveGlobalSystemMappedLDLT();
}

void VectorFields::setupRHSGlobalProblem()
{
	Eigen::VectorXd		zeroElements(B2D.rows());
	for (int i = 0; i < B2D.rows(); i++) zeroElements(i) = 0.0;

	b.resize(B2D.rows() + c.rows(),2);
	b.col(0) << zeroElements, c.col(0); 
	b.col(1) << zeroElements, c.col(1);

	//cout << "b: " << endl << b << endl << endl; 
}

void VectorFields::setupRHSGlobalProblemMapped()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Constructing RHS... ";

	vEst.resize(B2D.cols());
	for (int i = 0; i < vEst.rows(); i++) {
		vEst(i) = 0.5;
	}

	g = B2D * vEst;
	b.resize(B2D.rows() + c.rows(), 2);

	// First column of b
	h = C * vEst - c.col(0);
	b.col(0) << g, h;

	// Second Column of b
	h = C * vEst - c.col(1);
	b.col(1) << g, h;

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}

void VectorFields::setupLHSGlobalProblem()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Constructing LHS... ";

	A_LHS.resize(B2D.rows() + C.rows(), B2D.cols() + C.rows());

	vector<Eigen::Triplet<double>>	ATriplet;
	ATriplet.reserve(10 * B2D.rows());

	for (int k = 0; k < B2D.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(B2D, k); it; ++it) {
			ATriplet.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
		}
	}

	for (int k = 0; k < C.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(C, k); it; ++it) {
			ATriplet.push_back(Eigen::Triplet<double>(B2D.rows() + it.row(), it.col(), it.value()));
			ATriplet.push_back(Eigen::Triplet<double>(it.col(), B2D.cols() + it.row(), -it.value()));
		}
	}
	A_LHS.setFromTriplets(ATriplet.begin(), ATriplet.end());

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}

void VectorFields::setupLHSGlobalProblemMapped()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Constructing LHS... ";

	A_LHS.resize(B2D.rows() + C.rows(), B2D.cols() + C.rows());

	vector<Eigen::Triplet<double>>	ATriplet;
	ATriplet.reserve(10 * B2D.rows());		// It should be #rows x 4 blocks @ 2 elements (8) + #constraints,
											// but made it 10 for safety + simplicity

	for (int k = 0; k < B2D.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(B2D, k); it; ++it) {
			ATriplet.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
		}
	}

	for (int k = 0; k < C.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(C, k); it; ++it) {
			ATriplet.push_back(Eigen::Triplet<double>(B2D.rows() + it.row(), it.col(), it.value()));
			ATriplet.push_back(Eigen::Triplet<double>(it.col(), B2D.cols() + it.row(), it.value()));
		}
	}
	A_LHS.setFromTriplets(ATriplet.begin(), ATriplet.end());

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}

void VectorFields::solveGlobalSystem()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Solving the global system... \n";

	Xf.resize(B2D.rows(), 2);

	// Setting up the solver
	Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
	Eigen::VectorXd x; 
	//solver.compute(A);
	
	solver.analyzePattern(A_LHS);

	cout << "....Factorizing LHS..." << endl;
	solver.factorize(A_LHS);
	if (solver.info() != Eigen::Success) {
		cout << "Factorization Failed." << endl;
		return;
	}

	//printf("Size of B=%d\n", b.rows());
	// FIRST BASIS
	cout << "....Solving first problem (first frame)." << endl;
	x = solver.solve(b.col(0));
	if (solver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system." << endl;
		return;
	}

	Xf.col(0) = x.block(0, 0, B2D.rows(), 1);	

	// SECOND BASIS
	cout << "....Solving second problem. (second frame)" << endl;
	x = solver.solve(b.col(1));
	if (solver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system." << endl;
		return;
	}
	Xf.col(1) = x.block(0, 0, B2D.rows(), 1);
	
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "..in " << duration.count() << " seconds" << endl;
}

void VectorFields::solveGlobalSystemMappedLDLT()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Solving the global system (Pardiso LDLT)... \n";

	//cout << "Starting to solve problem." << endl;
	Xf.resize(B2D.rows(), 2);
	
	// Setting up the solver
	//Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> sparseSolver(A_LHS);
	Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> sparseSolver(A_LHS);	
	//Eigen::PastixLDLT<Eigen::SparseMatrix<double>,1> sparseSolver(A_LHS);

	// FIRST BASIS
	cout << "....Solving first problem (first frame)..." << endl;
	Eigen::VectorXd x = sparseSolver.solve(b.col(0));
	
	if (sparseSolver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system. " << endl;
		if (sparseSolver.info() == Eigen::NumericalIssue)
			cout << "NUMERICAL ISSUE. " << endl;
		cout << sparseSolver.info() << endl;
		return;
	}
	
	Xf.col(0) = -x.block(0, 0, B2D.rows(), 1) + vEst;
	
	// SECOND BASIS
	cout << "....Solving second problem (second frame)." << endl;
	x = sparseSolver.solve(b.col(1));
	if (sparseSolver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system." << endl;
		return;
	}
	
	Xf.col(1) = -x.block(0, 0, B2D.rows(), 1) + vEst;
	
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "..in " << duration.count() << " seconds" << endl;
}

void VectorFields::constructSamples(const int &n)
{
	numSample = n; 

	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;

	t1 = chrono::high_resolution_clock::now();
	farthestPointSampling();
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;

	cout << "> Constructing " << n << " samples in " << duration.count() << "seconds" << endl;
}

void VectorFields::farthestPointSampling()
{
	Sample.resize(numSample);
	Eigen::VectorXd D;
	D.resize(F.rows());

	// Initialize the value of D
	for (int i = 0; i < F.rows(); i++) {
		D(i) = numeric_limits<double>::infinity();
	}

	srand(time(NULL));
	Sample[0] = rand() % F.rows();
	//Sample[0] = 0;

	//computeDijkstraDistanceFaceForSampling(Sample[0], D);
	//Eigen::VectorXi::Index maxIndex1;
	//D.maxCoeff(&maxIndex1);
	//Sample[1] = maxIndex1;

	for (int i = 1; i < numSample; i++) {
		Eigen::VectorXi::Index maxIndex;
		computeDijkstraDistanceFaceForSampling(Sample[i-1], D);
		D.maxCoeff(&maxIndex);
		Sample[i] = maxIndex;
	}
}

void VectorFields::constructBasis()
{	
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Construting Basis...\n";

	double	coef = sqrt(pow(1.1, 2) + pow(0.7, 2));
	double distRatio = coef * sqrt((double)V.rows() / (double) Sample.size());

	// Setup sizes of each element to construct basis
	try {
		BasisTemp.resize(2 * F.rows(), 2 * Sample.size());
	}
	catch(string &msg) {
		cout << "Cannot allocate memory for basis.." << endl;
	}
	
	Basis.resize(BasisTemp.rows(), BasisTemp.cols());
	vector<vector<Eigen::Triplet<double>>> UiTriplet(Sample.size());

	
	cout << "....Constructing and solving local systems...";
	const int NUM_PROCESS = 8;
	durations.resize(NUM_PROCESS);

	for (int i = 0; i < NUM_PROCESS; i++) {
		durations[i] = t1 - t1;
		//cout << "Init dur " << i<< " = " << durations[i].count() << " seconds" << endl;
	}
	
	int id, tid, ntids, ipts, istart, iproc;
	

#pragma omp parallel private(tid,ntids,ipts,istart,id)	
	{		
		iproc = omp_get_num_procs();
		//iproc = 1; 
		tid		= omp_get_thread_num();
		ntids	= omp_get_num_threads();
		ipts	= (int)ceil(1.00*(double)Sample.size() / (double)ntids);
		istart	= tid * ipts;
		if (tid == ntids - 1) ipts = Sample.size() - istart;
		if (ipts <= 0) ipts = 0;

		Eigen::VectorXd				D(F.rows());
		for (int i = 0; i < F.rows(); i++) {
			D(i) = numeric_limits<double>::infinity();
		}
		
		//cout << "[" << tid << "] Number of processors " << iproc << ", with " << ntids << " threads." << endl;

		UiTriplet[tid].reserve(2.0 * ((double)ipts / (double)Sample.size()) * 2 * 10.0 * F.rows());

		// Computing the values of each element
		for (id = istart; id < (istart + ipts); id++) {
			if (id >= Sample.size()) break;

			LocalFields localField(id);
				t1 = chrono::high_resolution_clock::now();
			localField.constructSubdomain(Sample[id], V, F, avgEdgeLength, AdjMF3N, distRatio);
				t2 = chrono::high_resolution_clock::now();
				durations[0] += t2 - t1;

				t1 = chrono::high_resolution_clock::now();
			localField.constructBoundary(F, AdjMF3N, AdjMF2Ring);
				t2 = chrono::high_resolution_clock::now();
				durations[1] += t2 - t1;

				t1 = chrono::high_resolution_clock::now();
			localField.constructLocalElements(F);
				t2 = chrono::high_resolution_clock::now();
				durations[2] += t2 - t1;

				t1 = chrono::high_resolution_clock::now();
			//localField.constructMatrixBLocal(B2D);
			localField.constructMatrixBLocal(B2D, AdjMF2Ring);
				t2 = chrono::high_resolution_clock::now();
				durations[3] += t2 - t1;

				t1 = chrono::high_resolution_clock::now();
			localField.constructLocalConstraints();
				t2 = chrono::high_resolution_clock::now();
				durations[4] += t2 - t1;

				t1 = chrono::high_resolution_clock::now();
			localField.setupRHSLocalProblemMapped();
				t2 = chrono::high_resolution_clock::now();
				durations[5] += t2 - t1;

				t1 = chrono::high_resolution_clock::now();
			localField.setupLHSLocalProblemMapped();
				t2 = chrono::high_resolution_clock::now();
				durations[6] += t2 - t1;

				t1 = chrono::high_resolution_clock::now();
			localField.solveLocalSystemMappedLDLT(UiTriplet[id]);
				t2 = chrono::high_resolution_clock::now();
				durations[7] += t2 - t1;
			//cout << "System " << id << " ( " << XfLoc.rows() << ") is solved." << endl; 
			//printf("System %d (%d) is solved.\n", id, XfLoc.rows());
		}

	}
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl; 

	cout << "....Gathering local elements as basis matrix... ";
	t1 = chrono::high_resolution_clock::now();
	gatherBasisElements(UiTriplet);
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
	
	cout << "....Partition of unity of the basis matrix... ";
	t1 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	//normalizeBasis();
	normalizeBasisAbs();
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "..in Total of " << duration.count() << " seconds" << endl;

	//for (int i = 0; i < NUM_PROCESS; i++) {
		//printf("Process %d takes %.8f seconds.\n", i, durations[i].count());
		//cout << "BasisTemp matrix is normalized in " << duration.count() << " seconds" << endl;
	//}

	// Information about the process timing
	printf("> Basis Timing information \n");
	printf("....[0] Constructing internal elements: %.8f seconds.\n", durations[0].count());
	printf("....[1] Constructing boundary: %.8f seconds.\n", durations[1].count());
	printf("....[2] Constructing local subdomains: %.8f seconds.\n", durations[2].count());
	printf("....[3] Constructing matrix B local: %.8f seconds.\n", durations[3].count());
	printf("....[4] Constructing local constraints: %.8f seconds.\n", durations[4].count());
	printf("....[5] Constructing RHS (mapped): %.8f seconds.\n", durations[5].count());
	printf("....[6] Constructing LHS (mapped): %.8f seconds.\n", durations[6].count());
	printf("....[7] Solving local systems (mapped, Pardiso LDLT): %.8f seconds.\n", durations[7].count());

	// Information about Basis
	printf("> Basis Structure information \n");
	printf("....Size = %dx%d\n", Basis.rows(), Basis.cols());
	printf("....NNZ per row = %.2f\n", (double) Basis.nonZeros() / (double) Basis.rows());
}

void VectorFields::gatherBasisElements(const vector<vector<Eigen::Triplet<double>>> &UiTriplet)
{
	vector<Eigen::Triplet<double>> BTriplet;
	BasisSum.resize(2 * F.rows(), 2);
	for (int i = 0; i < BasisSum.rows(); i++) {
		BasisSum(i, 0) = 0.0;
		BasisSum(i, 1) = 0.0;
	}

	int totalElements = 0;
	for (int j = 0; j < Sample.size(); j++) {
		totalElements += UiTriplet[j].size();
	}

	BTriplet.resize(totalElements);
	for (int j = 0; j < Sample.size(); j++) {
		int tripSize = 0;
		for (int k = 0; k < j; k++) {
			tripSize += UiTriplet[k].size();
		}
		std::copy(UiTriplet[j].begin(), UiTriplet[j].end(), BTriplet.begin() + tripSize);
	}
	BasisTemp.setFromTriplets(BTriplet.begin(), BTriplet.end());

	//printf("A basis matrix (%dx%d) is constructed.\n", BasisTemp.rows(), BasisTemp.cols());

	for (int k = 0; k < BasisTemp.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(BasisTemp, k); it; ++it) {
			BasisSum(it.row(), it.col() % 2) += it.value();
		}
	}

}

void VectorFields::normalizeBasis()
{
	Eigen::MatrixXd BasisSum(BasisTemp.rows(), 2);
	BasisSumN.resize(BasisTemp.rows(), 2);
	Eigen::MatrixXd normSum(BasisTemp.rows(), 2);
	Eigen::MatrixXd BasisNorm(F.rows(), 2);
	Eigen::MatrixXi nonZeros(BasisTemp.rows(), 2);
		for (int i = 0; i < nonZeros.rows(); i++) {
			for (int j = 0; j < nonZeros.cols(); j++) {
				nonZeros(i, j) = 0;
				BasisSum(i, j) = 0.0;
				BasisSumN(i, j) = 0.0;
			}
		}
	vector<Eigen::Triplet<double>> BNTriplet;
	BNTriplet.reserve(BasisTemp.nonZeros());

	// Getting the sum of each pair on each frame AND
	// Counting the non-zeros per rows
	for (int k = 0; k < BasisTemp.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(BasisTemp, k); it; ++it) {
			BasisSum(it.row(), it.col() % 2) += it.value();
			nonZeros(it.row(), it.col() % 2) += 1;
		}
	}
	
	// Computing the norm
	for (int i = 0; i < F.rows(); i++) {
		double frame1Norm, frame2Norm, a,b;
		a = BasisSum(2 * i + 0, 0);
		b = BasisSum(2 * i + 1, 0);
		frame1Norm = sqrt(a*a + b*b);

		a = BasisSum(2 * i + 0, 1);
		b = BasisSum(2 * i + 1, 1);
		frame2Norm = sqrt(a*a + b*b);

		BasisNorm(i, 0) = frame1Norm;
		BasisNorm(i, 1) = frame2Norm;		
	}

	// Constructing normalized basis each element of basis
	for (int k = 0; k < BasisTemp.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(BasisTemp, k); it; ++it) {
			//double newValue = it.value() / BasisNorm(it.row() / 2, it.col() % 2);
			double newValue = it.value();
			BNTriplet.push_back(Eigen::Triplet<double>(it.row(), it.col(), newValue));
		}
	}
	Basis.setFromTriplets(BNTriplet.begin(), BNTriplet.end());	

	// Check Normalization
	// Getting the sum of each pair on each frame
	for (int k = 0; k < Basis.outerSize(); ++k) {		
		for (Eigen::SparseMatrix<double>::InnerIterator it(Basis, k); it; ++it) {
			BasisSumN(it.row(), it.col() % 2) += it.value();
		}
	}

	// Computing the norm
	for (int i = 0; i < F.rows(); i++) {
		double frame1Norm, frame2Norm, a, b;
		a = BasisSumN(2 * i + 0, 0);
		b = BasisSumN(2 * i + 1, 0);
		frame1Norm = sqrt(a*a + b*b);

		a = BasisSumN(2 * i + 0, 1);
		b = BasisSumN(2 * i + 1, 1);
		frame2Norm = sqrt(a*a + b*b);

		BasisNorm(i, 0) = frame1Norm;
		BasisNorm(i, 1) = frame2Norm;
	}
	
	// Show result (The sum=> should all be equal to 1
	//for (int i = 0; i < F.rows(); i++) {
	//	printf("[%.6f][%.6f]\n", BasisSumN.block(2*i,0,2,1).norm(), BasisSumN.block(2 * i, 1, 2, 1).norm());
	//}

	int numNonZeroes = Basis.nonZeros();
	int numElements = Basis.rows();
	cout << "Average non-zeros is " << (double)numNonZeroes / (double)numElements << endl; 
}

//void VectorFields::normalizeBasisAbs()
//{
//	Eigen::MatrixXd BasisSum(BasisTemp.rows(), 2);
//	BasisSumN.resize(BasisTemp.rows(), 2);
//	Eigen::MatrixXd normSum(BasisTemp.rows(), 2);
//	Eigen::MatrixXd BasisNorm(F.rows(), 2);
//	Eigen::MatrixXi nonZeros(BasisTemp.rows(), 2);
//	for (int i = 0; i < nonZeros.rows(); i++) {
//		for (int j = 0; j < nonZeros.cols(); j++) {
//			nonZeros(i, j) = 0;
//			BasisSum(i, j) = 0.0;
//			BasisSumN(i, j) = 0.0;
//		}
//	}
//	vector<Eigen::Triplet<double>> BNTriplet;
//
//
//	// Getting the sum of each pair on each frame AND
//	// Counting the non-zeros per rows
//	for (int k = 0; k < BasisTemp.outerSize(); ++k) {
//		for (Eigen::SparseMatrix<double>::InnerIterator it(BasisTemp, k); it; ++it) {
//			BasisSum(it.row(), it.col() % 2) += abs(it.value());
//			nonZeros(it.row(), it.col() % 2) += 1;
//		}
//	}
//
//	// Computing the norm
//	for (int i = 0; i < F.rows(); i++) {
//		double frame1Norm, frame2Norm, a, b;
//		a = BasisSum(2 * i + 0, 0);
//		b = BasisSum(2 * i + 1, 0);
//		frame1Norm = sqrt(a*a + b*b);
//
//		a = BasisSum(2 * i + 0, 1);
//		b = BasisSum(2 * i + 1, 1);
//		frame2Norm = sqrt(a*a + b*b);
//
//		BasisNorm(i, 0) = frame1Norm;
//		BasisNorm(i, 1) = frame2Norm;
//	}
//
//	// Constructing normalized basis each element of basis
//	for (int k = 0; k < BasisTemp.outerSize(); ++k) {
//		for (Eigen::SparseMatrix<double>::InnerIterator it(BasisTemp, k); it; ++it) {
//			double newValue = it.value() / BasisNorm(it.row() / 2, it.col() % 2);
//			//double newValue = it.value();
//			BNTriplet.push_back(Eigen::Triplet<double>(it.row(), it.col(), newValue));
//		}
//	}
//	Basis.setFromTriplets(BNTriplet.begin(), BNTriplet.end());
//
//	// Check Normalization
//	// Getting the sum of each pair on each frame
//	for (int k = 0; k < Basis.outerSize(); ++k) {
//		for (Eigen::SparseMatrix<double>::InnerIterator it(Basis, k); it; ++it) {
//			BasisSumN(it.row(), it.col() % 2) += it.value();
//		}
//	}
//
//	// Computing the norm
//	for (int i = 0; i < F.rows(); i++) {
//		double frame1Norm, frame2Norm, a, b;
//		a = BasisSumN(2 * i + 0, 0);
//		b = BasisSumN(2 * i + 1, 0);
//		frame1Norm = sqrt(a*a + b*b);
//
//		a = BasisSumN(2 * i + 0, 1);
//		b = BasisSumN(2 * i + 1, 1);
//		frame2Norm = sqrt(a*a + b*b);
//
//		BasisNorm(i, 0) = frame1Norm;
//		BasisNorm(i, 1) = frame2Norm;
//	}
//
//	// Show result (The sum=> should all be equal to 1
//	//for (int i = 0; i < F.rows(); i++) {
//	//	printf("[%.6f][%.6f]\n", BasisSumN.block(2*i,0,2,1).norm(), BasisSumN.block(2 * i, 1, 2, 1).norm());
//	//}
//
//	int numNonZeroes = Basis.nonZeros();
//	int numElements = Basis.rows();
//	cout << "Average non-zeros is " << (double)numNonZeroes / (double)numElements << endl;
//}

void VectorFields::normalizeBasisAbs()
{
	Eigen::MatrixXd normSum(F.rows(), 2);
	BasisSumN.resize(BasisTemp.rows(), 2);
	vector<Eigen::Triplet<double>> BNTriplet;
	BNTriplet.reserve(BasisTemp.nonZeros());

	Eigen::MatrixXd BasisNorm(F.rows(), 2);

	for (int i = 0; i < normSum.rows(); i++) {
		for (int j = 0; j < normSum.cols(); j++) {
			normSum(i, j) = 0.0;
			BasisSumN(2 * i + 0, j) = 0.0;
			BasisSumN(2 * i + 1, j) = 0.0;
		}
	}

	// Getting the sum of norm on each pair on each frame 
	for (int k = 0; k < BasisTemp.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(BasisTemp, k); it; ++it) {
			if (it.row() % 2 == 1) continue;

			double a = it.value();
			double b = BasisTemp.coeff(it.row() + 1, it.col());
			double norm = sqrt(a*a + b*b);
			normSum(it.row() / 2, it.col() % 2) += norm;
		}
	}

	// Normalize the system
	for (int k = 0; k < BasisTemp.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(BasisTemp, k); it; ++it) {
			double newValue = it.value() / normSum(it.row() / 2, it.col() % 2);
			// To have the basis with norm 2.0
			BNTriplet.push_back(Eigen::Triplet<double>(it.row(), it.col(), newValue/2.0));
			BasisSumN(it.row(), it.col() % 2) += newValue;
		}
	}

	Basis.setFromTriplets(BNTriplet.begin(), BNTriplet.end());	
}
void VectorFields::setAndSolveUserSystem()
{
	setupUserBasis();
	getUserConstraints();
	setupRHSUserProblemMapped();
	setupLHSUserProblemMapped();
	solveUserSystemMappedLDLT();
	mapSolutionToFullRes();
	//obtainUserVectorFields();
}

void VectorFields::setupUserBasis()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Constructing Local Basis...";

	B2Dbar = Basis.transpose() * B2D * Basis; 

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;
	printf(".... Local Basis = %dx%d\n", B2Dbar.rows(), B2Dbar.cols());
}
void VectorFields::getUserConstraints()
{	
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Obtaining user constraints ";

	//getUserConstraintsRandom();

	//const vector<vector<int>> selectedFaces;
	//getUserConstraintsGUI(selectedFaces);
	//getUserConstraintsSpecified();
	userConstraints = globalConstraints; 
	Cbar = C * Basis;
	cBar = c;

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;

	printf(".... C_Lobal = %dx%d\n", Cbar.rows(), Cbar.cols());
	printf(".... c_Lobal = %dx%d\n", cBar.rows(), cBar.cols());
}

void VectorFields::getUserConstraintsRandom()
{
	// Define the constraints
	const int numConstraints = 20;

	srand(time(NULL));
	set<int> constraints;
	userConstraints.resize(numConstraints);
	do {
		int c = rand() % F.rows();
		constraints.insert(c);
	} while (constraints.size() <= numConstraints);

	int counter1 = 0;
	for (int i : constraints) {
		userConstraints[counter1++] = i;
	}
	printf("Constraints = %d\n", userConstraints.size());

	// Setting up matrix C
	Eigen::SparseMatrix<double> CTemp;
	vector<Eigen::Triplet<double>> CTriplet;
	CTriplet.reserve(2 * userConstraints.size());

	int counter = 0;
	for (int i = 0; i < userConstraints.size(); i++) {
		CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * userConstraints[i] + 0, 1.0));
		CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * userConstraints[i] + 1, 1.0));
	}
	CTemp.resize(2 * userConstraints.size(), Basis.rows());
	CTemp.setFromTriplets(CTriplet.begin(), CTriplet.end());
	printf("CTemp=%dx%d\n", CTemp.rows(), CTemp.cols());

	Cbar = CTemp * Basis;
	printf("Cbar=%dx%d\n", Cbar.rows(), Cbar.cols());
	// Setting up vector c (There are 2 vector c)
	srand(time(NULL));
	cBar.resize(2 * userConstraints.size(), 2);
	for (int i = 0; i < userConstraints.size(); i++) {
		cBar(2 * i + 0, 0) = ((double)(rand() % 1000) - 500) / 1000.0;
		cBar(2 * i + 1, 0) = ((double)(rand() % 1000) - 500) / 1000.0;

		cBar(2 * i + 0, 1) = ((double)(rand() % 1000) - 500) / 1000.0;
		cBar(2 * i + 1, 1) = ((double)(rand() % 1000) - 500) / 1000.0;
	}
	printf("cBar=%dx%d\n", cBar.rows(), cBar.cols());
}

void VectorFields::getUserConstraintsGUI(const vector<vector<int>> &selectedFaces)
{

}

void VectorFields::getUserConstraintsSpecified()
{
	// Define the constraints
	const int numConstraints = 20;
	set<int> constraints;
	userConstraints.resize(numConstraints);
	Eigen::VectorXd D;
	D.resize(F.rows());

	// Initialize the value of D
	for (int i = 0; i < F.rows(); i++) {
		D(i) = numeric_limits<double>::infinity();
	}

	constraints.insert(0);	
	int curPoint = 0;
		
	do {
		Eigen::VectorXi::Index maxIndex;
		computeDijkstraDistanceFaceForSampling(curPoint, D);
		D.maxCoeff(&maxIndex);
		constraints.insert(maxIndex);
		curPoint = maxIndex;
	} while (constraints.size() <= numConstraints);

	int counter1 = 0;
	for (int i : constraints) {
		userConstraints[counter1++] = i;
	}
	printf("Constraints = %d\n", userConstraints.size());

	// Setting up matrix C
	Eigen::SparseMatrix<double> CTemp;
	vector<Eigen::Triplet<double>> CTriplet;
	CTriplet.reserve(2 * userConstraints.size());
	int counter = 0;
	for (int i = 0; i < userConstraints.size(); i++) {
		CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * userConstraints[i] + 0, 1.0));
		CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * userConstraints[i] + 1, 1.0));
	}
	CTemp.resize(2 * userConstraints.size(), Basis.rows());
	CTemp.setFromTriplets(CTriplet.begin(), CTriplet.end());
	printf("CTemp=%dx%d\n", CTemp.rows(), CTemp.cols());

	Cbar = CTemp * Basis;
	printf("Cbar=%dx%d\n", Cbar.rows(), Cbar.cols());

	// Setting up vector c (There are 2 vector c)
	srand(time(NULL));
	cBar.resize(2 * userConstraints.size(), 2);
	for (int i = 0; i < userConstraints.size(); i++) {
		cBar(2 * i + 0, 0) = sqrt(2.0);
		cBar(2 * i + 1, 0) = sqrt(2.0);
		cBar(2 * i + 0, 1) = sqrt(2.0);
		cBar(2 * i + 1, 1) = sqrt(2.0);
	}
	printf("cBar=%dx%d\n", cBar.rows(), cBar.cols());
}

void VectorFields::setupRHSUserProblemMapped()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Constructing RHS (mapped)...";
	
	vEstUser.resize(B2Dbar.rows());
	for (int i = 0; i < vEstUser.rows(); i++) {
		vEstUser(i) = 0.5;
	}

	gbar = B2Dbar * vEstUser;
	//printf("gbar=%dx%d\n", gbar.rows(), gbar.cols());
	bBar.resize(B2Dbar.rows() + cBar.rows(), 2);
	//printf("bbar=%dx%d\n",bBar.rows(), bBar.cols());

	// First column of b
	hbar = Cbar * vEstUser - cBar.col(0);
	//printf("hbar=%dx%d\n", hbar.rows(), hbar.cols());

	bBar.col(0) << gbar, hbar;
	//printf("bbar=%dx%d\n", bBar.rows(), bBar.cols());
	// Second Column of b
	hbar = Cbar * vEstUser - cBar.col(1);
	//printf("hbar=%dx%d\n", hbar.rows(), hbar.cols());

	bBar.col(1) << gbar, hbar;


	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;
}

void VectorFields::setupLHSUserProblemMapped()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Constructing LHS (mapped)...";


	A_LHSbar.resize(B2Dbar.rows() + Cbar.rows(), B2Dbar.cols() + Cbar.rows());
	vector<Eigen::Triplet<double>>	ATriplet;
	ATriplet.reserve(B2Dbar.nonZeros() + 2 * Cbar.nonZeros());

	for (int k = 0; k < B2Dbar.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(B2Dbar, k); it; ++it) {
			ATriplet.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
		}
	}

	for (int k = 0; k < Cbar.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(Cbar, k); it; ++it) {
			ATriplet.push_back(Eigen::Triplet<double>(B2Dbar.rows() + it.row(), it.col(), it.value()));
			ATriplet.push_back(Eigen::Triplet<double>(it.col(), B2Dbar.cols() + it.row(), it.value()));
		}
	}
	A_LHSbar.setFromTriplets(ATriplet.begin(), ATriplet.end());

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;
	printf("....Local LHS = %dx%d\n", A_LHSbar.rows(), A_LHSbar.cols());
}

void VectorFields::solveUserSystemMappedLDLT()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Solving reduced system...\n";

	XLowDim.resize(B2Dbar.rows(), 2);
	XFullDim.resize(Basis.rows(), 2);
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> sparseSolver(A_LHSbar);

	cout << "....Solving for the first frame.\n";
	Eigen::VectorXd x = sparseSolver.solve(bBar.col(0));

	if (sparseSolver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system. " << endl;
		if (sparseSolver.info() == Eigen::NumericalIssue)
			cout << "NUMERICAL ISSUE. " << endl;
		cout << sparseSolver.info() << endl;
		return;
	}

	XLowDim.col(0) = -x.block(0, 0, B2Dbar.rows(), 1) + vEstUser;

	// SECOND BASIS
	cout << "....Solving for the second frame.\n";
	x = sparseSolver.solve(bBar.col(1));
	if (sparseSolver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system." << endl;
		return;
	}
	XLowDim.col(1) = -x.block(0, 0, B2Dbar.rows(), 1) + vEstUser;
	
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "..in Total of " << duration.count() << " seconds." << endl;

	//cout << "Solution (LowDim) \n" << XLowDim << endl; 	
}

void VectorFields::mapSolutionToFullRes()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Mapping to full-resolution...";

	XFullDim = Basis * XLowDim;


	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << " in " << duration.count() << " seconds." << endl;

	printf("....XFull (%dx%d) =  Basis (%dx%d) * XLowDim (%dx%d) \n", XFullDim.rows(), XFullDim.cols(), Basis.rows(), Basis.cols(), XLowDim.rows(), XLowDim.cols());
}

void VectorFields::obtainUserVectorFields()
{
	cout << "Hello there \n" << endl; 
}

void VectorFields::measureApproxAccuracyL2Norm()
{
	Eigen::VectorXd diff = Xf.col(0) - XFullDim.col(0);
	double xf = Xf.col(0).transpose() * MF2D * Xf.col(0);
	double L2norm = diff.transpose() * MF2D * diff; 
	printf("Diff 0 = %.10f\n", sqrt(L2norm / xf)); 

	diff = Xf.col(1) - XFullDim.col(1);
	xf = Xf.col(1).transpose() * MF2D * Xf.col(1);
	L2norm = diff.transpose() * MF2D * diff;
	printf("Diff 1 = %.10f\n", sqrt(L2norm / xf));
}


void VectorFields::constructLocalElements()
{
	LocalElements.resize(SubDomain.size() + Boundary.size());
	GlobToLocMap.resize(F.rows());
	for (int i = 0; i < F.rows(); i++) GlobToLocMap[i] = -1; 

	int counter = 0;
	for (int face : SubDomain) {
		LocalElements[counter] = face; 
		GlobToLocMap[face] = counter;
		counter++;
	}

	for (int face : Boundary) {
		LocalElements[counter] = face; 
		GlobToLocMap[face] = counter;
		counter++;
	}
}

void VectorFields::constructMatrixBLocal()
{
	BLoc.resize(2 * LocalElements.size(), 2 * LocalElements.size());
	vector<Eigen::Triplet<double>> BTriplet;
	BTriplet.reserve(10 * BLoc.rows());

	for (int i = 0; i < LocalElements.size(); i++) {
		int li = LocalElements[i];
		// Get the DIAGONAL Elements from B2D Matrix
		BTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, B2D.coeff(2 * li + 0, 2 * li + 0)));
		BTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 1, B2D.coeff(2 * li + 0, 2 * li + 1)));
		BTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 0, B2D.coeff(2 * li + 1, 2 * li + 0)));
		BTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, B2D.coeff(2 * li + 1, 2 * li + 1)));

		for (int j = i + 1; j < LocalElements.size(); j++) {
			// Get the HORIZONTAL Elements
			int lj = LocalElements[j];
			BTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * j + 0, B2D.coeff(2 * li + 0, 2 * lj + 0)));
			BTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * j + 1, B2D.coeff(2 * li + 0, 2 * lj + 1)));
			BTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * j + 0, B2D.coeff(2 * li + 1, 2 * lj + 0)));
			BTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * j + 1, B2D.coeff(2 * li + 1, 2 * lj + 1)));

			// Get the VERTICAL Elements
			BTriplet.push_back(Eigen::Triplet<double>(2 * j + 0, 2 * i + 0, B2D.coeff(2 * lj + 0, 2 * li + 0)));
			BTriplet.push_back(Eigen::Triplet<double>(2 * j + 0, 2 * i + 1, B2D.coeff(2 * lj + 0, 2 * li + 1)));
			BTriplet.push_back(Eigen::Triplet<double>(2 * j + 1, 2 * i + 0, B2D.coeff(2 * lj + 1, 2 * li + 0)));
			BTriplet.push_back(Eigen::Triplet<double>(2 * j + 1, 2 * i + 1, B2D.coeff(2 * lj + 1, 2 * li + 1)));
		}
	}

	BLoc.setFromTriplets(BTriplet.begin(), BTriplet.end());
}

void VectorFields::constructLocalConstraints()
{
	// Setting up matrix C
	CLoc.resize(2 * (1+LocalElements.size()), BLoc.cols());
	vector<Eigen::Triplet<double>> CTriplet;
	CTriplet.reserve(2 * (1 + LocalElements.size()));

	int counter = 0; 
	CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * GlobToLocMap[sample] + 0, 1.0));
	CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * GlobToLocMap[sample] + 1, 1.0));
	for(int bound : Boundary){
		CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * GlobToLocMap[bound] + 0, 1.0));
		CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * GlobToLocMap[bound] + 1, 1.0));
	}
	CLoc.setFromTriplets(CTriplet.begin(), CTriplet.end());
	cout << "Selection of faces is done." << endl;

	// Setting up vector c (There are 2 vector c)
	cLoc.resize(2 * (1 + LocalElements.size()), 2);
	Eigen::VectorXd zeroElements(2 * LocalElements.size());
	for (int i = 0; i < zeroElements.size(); i++) zeroElements(i) = 0.0;
	cLoc.col(0) << 1.0, 0.0, zeroElements;
	cLoc.col(1) << 0.0, 1.0, zeroElements;
	cout << "Definition of constraints is done." << endl;
}

void VectorFields::setupRHSLocalProblem()
{
	Eigen::VectorXd		zeroElements(BLoc.rows());
	for (int i = 0; i < zeroElements.size(); i++) zeroElements(i) = 0.0;

	bLoc.resize(BLoc.rows() + cLoc.rows(), 2);
	bLoc.col(0) << zeroElements, cLoc.col(0);
	bLoc.col(1) << zeroElements, cLoc.col(1);
	cout << "RHS is constructed." << endl;
}

void VectorFields::setupLHSLocalProblem()
{
	ALoc.resize(BLoc.rows() + CLoc.rows(), BLoc.cols() + CLoc.rows());	
	vector<Eigen::Triplet<double>>	ATriplet;
	ATriplet.reserve(BLoc.nonZeros() + 2 * CLoc.nonZeros());

	for (int k = 0; k < BLoc.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(BLoc, k); it; ++it) {
			ATriplet.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
		}
	}

	for (int k = 0; k < CLoc.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(CLoc, k); it; ++it) {
			ATriplet.push_back(Eigen::Triplet<double>(BLoc.rows() + it.row(), it.col(), it.value()));
			ATriplet.push_back(Eigen::Triplet<double>(it.col(), BLoc.cols() + it.row(), -it.value()));
		}
	}
	ALoc.setFromTriplets(ATriplet.begin(), ATriplet.end());
	cout << "LHS is constructed." << endl;
}

void VectorFields::solveLocalSystem()
{
	cout << "Starting to solve problem." << endl;
	XfLoc.resize(BLoc.rows(), 2);

	// Setting up the solver
	Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
	Eigen::VectorXd x;
	//solver.compute(A);

	solver.analyzePattern(ALoc);

	cout << "Factorizing LHS." << endl;
	solver.factorize(ALoc);
	if (solver.info() != Eigen::Success) {
		cout << "Factorization Failed." << endl;
		return;
	}

	printf("Size of B=%d\n", bLoc.rows());
	// FIRST BASIS
	cout << "Solving first problem." << endl;
	x = solver.solve(bLoc.col(0));
	if (solver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system." << endl;
		return;
	}
	

	XfLoc.col(0) = x.block(0, 0, BLoc.rows(), 1);

	// SECOND BASIS
	cout << "Solving second problem." << endl;
	x = solver.solve(bLoc.col(1));
	if (solver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system." << endl;
		return;
	}

	XfLoc.col(1) = x.block(0, 0, BLoc.rows(), 1);
	//for (int i = 0; i < B.rows(); i++) {
	//	//Xf.block(columnOrderMap(i), 0,1,1) = xx.block(i, 0, 1, 1);
	//}
	cout << "Problems solved." << endl;
	//cout << XfLoc << endl;
	printf("XfLoc=%dx%d\n", XfLoc.rows(), XfLoc.cols());
}

/* ====================== SETTING UP UTILITY MATRICES ============================*/ 
void VectorFields::computeAverageEdgeLength()
{
	Eigen::Vector3d e;
	double totalLength = 0.0;
	for (int i = 0; i < F.rows(); i++) {
		for (int j = 0; j < F.cols(); j++){
			e = V.row(F(i, (j+1)%F.cols())) - V.row(F(i, j));
			totalLength += e.norm();
		}		
	}
	avgEdgeLength = totalLength / (double)(F.rows()*F.cols());
}

void VectorFields::computeFaceNormal() 
{
	igl::per_face_normals(V, F, NF);
}

void VectorFields::constructFaceAdjacency2RingMatrix()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Building \"Face-to-Face\" Adjacency (All connected neigbors)... ";


	AdjMF2Ring.clear();
	AdjMF2Ring.resize(F.rows());

	vector<set<int>> faceNeighOnVert;
	faceNeighOnVert.resize(V.rows());
	
	//cout << "Construct face adjacency." << endl;

	//Getting the faces on which each vertex resides
	for (int i = 0; i < F.rows(); i++)
	{
		for (int j = 0; j < F.cols(); j++) {
			int a1 = F(i, j);
			faceNeighOnVert[a1].insert(i);
		}
	}

	// Getting the neighborhood structure
	for (int i = 0; i < V.rows(); i++) {
		for (std::set<int>::iterator it1 = faceNeighOnVert[i].begin(); it1 != faceNeighOnVert[i].end(); ++it1) {
			for (std::set<int>::iterator it2 = next(it1, 1); it2 != faceNeighOnVert[i].end(); ++it2) {
				AdjMF2Ring[*it1].insert(*it2);
				AdjMF2Ring[*it2].insert(*it1);
			}
		}
	}

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}

void VectorFields::constructFaceAdjacency3NMatrix()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Building \"Face-to-Face\" Adjacency (3 Neigbors)... ";

	// TWO RING
	AdjMF2Ring.clear();
	AdjMF2Ring.resize(F.rows());

	// ONE RING
	AdjMF3N_temp.resize(F.rows());
	EdgePairsList.resize(F.rows());
	int counter = 0, counter1=0;

	for (int i = 0; i < V.rows(); i++) {
		//printf("Vertex %d has %d neighbors, which arew: \n", i, VFNeighbors[i].size());
		for (std::set<VtoFPair>::iterator it1 = VFNeighbors[i].begin(); it1 != VFNeighbors[i].end(); ++it1) {
			//printf("%d (in F=%d) : ", it1->vId, it1->fId);
			for (std::set<VtoFPair>::iterator it2 = next(it1,1); it2 != VFNeighbors[i].end(); ++it2) {
				// ONE RING Neighbor
				if (it1->vId == it2->vId) {
					FacePair fp1{ counter1++,it2->fId };
					FacePair fp2{ counter1++,it1->fId };
					AdjMF3N_temp[it1->fId].insert(fp1);
					AdjMF3N_temp[it2->fId].insert(fp2);
					//AdjMF3N_temp[it1->fId].insert(it2->fId);
					//AdjMF3N_temp[it2->fId].insert(it1->fId);

					Edge_VPair ep1{ counter++, i, it1->vId };
					Edge_VPair ep2{ counter++, i, it2->vId };
					EdgePairsList[it1->fId].insert(ep1);
					EdgePairsList[it2->fId].insert(ep2);
				}
			}
		}
		//cout << endl; 
	}	

	// MOVING THE ADJACENCY it to matrix format
	AdjMF3N.resize(F.rows(), F.cols());
	for (int i = 0; i < F.rows(); i++) {
		int counter = 0;
		//for (std::set<int>::iterator it = AdjMF3N_temp[i].begin(); it != AdjMF3N_temp[i].end(); ++it) {
		for (std::set<FacePair>::iterator it = AdjMF3N_temp[i].begin(); it != AdjMF3N_temp[i].end(); ++it) {
			AdjMF3N(i, counter++) = it->f1;
		}
	}

	// Checking the result of EdgeList, print it
	//for (int i = 0; i < F.rows(); i++) {
	//	printf("F=%d has %d neighbors, which are", i, AdjMF3N_temp[i].size());
	//	int eCounter = 0; 
	//	for (std::set<Edge_VPair>::iterator it = EdgePairsList[i].begin(); it != EdgePairsList[i].end(); ++it) {
	//		cout << AdjMF3N(i,eCounter++) << " in edge (" << it->v1 << ", " << it->v2 << ") ";
	//	}
	//	cout << endl;
	//}

	// MOVING EDGE Adjacency to MATRIX format
	EdgePairMatrix.resize(F.rows(), 3 * F.cols());
	for (int i = 0; i < F.rows(); i++) {
		int eCounter = 0;
		for (std::set<Edge_VPair>::iterator it = EdgePairsList[i].begin(); it != EdgePairsList[i].end(); ++it) {
			EdgePairMatrix(i, eCounter++) = it->v1;
			EdgePairMatrix(i, eCounter++) = it->v2;
		} 
	}

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;

}

void VectorFields::constructFaceAdjacencyMatrix_IGL()
{
	igl::adjacency_matrix(F, AjdMF_igl);
	//visualizeSparseMatrixInMatlab(AjdMF_igl);
}

void VectorFields::constructVertexAdjacencyMatrix()
{	
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Building \"Vertex-to-Vertex\" Adjacency... ";
	
	AdjMV.clear();
	AdjMV.resize(V.rows());

	for (int i = 0; i < F.rows(); i++)
	{
		int a1 = F(i, 0);
		int a2 = F(i, 1);
		int a3 = F(i, 2);

		AdjMV[a1].insert(a2);
		AdjMV[a1].insert(a3);
		AdjMV[a2].insert(a1);
		AdjMV[a2].insert(a3);
		AdjMV[a3].insert(a1);
		AdjMV[a3].insert(a2);
	}

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}

void VectorFields::constructNeigborRings(const int &idx) 
{
	// Define unvisited triangles
	set<int> VisitedTrianglesSet;

	// Initialize the neighboring ring
	NeighRing.clear();
	set<int> oneRing;
	oneRing.insert(idx);
	VisitedTrianglesSet.insert(idx);
	NeighRing.push_back(oneRing);

	oneRing.clear();
	for (std::set<int, double>::iterator it = AdjMF2Ring[idx].begin(); it != AdjMF2Ring[idx].end(); ++it) {
		oneRing.insert(*it);
		VisitedTrianglesSet.insert(*it);
	}
	NeighRing.push_back(oneRing);

	// Constructing the rings via "set" data structure
	while (VisitedTrianglesSet.size()<F.rows())
	{
		int k = NeighRing.size() - 1;
		oneRing.clear();
		for (std::set<int, double>::iterator jt = NeighRing[k].begin(); jt != NeighRing[k].end(); ++jt) {
			for (std::set<int, double>::iterator kt = AdjMF2Ring[*jt].begin(); kt != AdjMF2Ring[*jt].end(); ++kt) {
				if (VisitedTrianglesSet.find(*kt) == VisitedTrianglesSet.end()) {
					oneRing.insert(*kt);
					VisitedTrianglesSet.insert(*kt);
				}
			}
		}
		NeighRing.push_back(oneRing);
	}

	// Displaying the list of the rings
	/*for (int i = 0; i < NeighRing.size(); i++) {
		printf("Ring %d: ", i);
		for (std::set<int, double>::iterator it = NeighRing[i].begin(); it != NeighRing[i].end(); ++it) {
			printf("%d =>", *it);
		}
		printf("\n");
	}*/
}

void VectorFields::computeDijkstraDistanceVertex(const int &source, Eigen::VectorXd &D)
{
	priority_queue<VertexPair, std::vector<VertexPair>, std::greater<VertexPair>> DistPQueue;
	
	// Computing distance for initial sample points S
	for (int i = 0; i < V.rows(); i++) {
		D(i) = numeric_limits<double>::infinity();
	}

	D(source) = 0.0f;
	VertexPair vp{ source,D(source) };
	DistPQueue.push(vp);

	double distFromCenter = numeric_limits<double>::infinity();

	// For other vertices in mesh
	do {
		if (DistPQueue.size() == 0) break;
		VertexPair vp1 = DistPQueue.top();
		distFromCenter = vp1.distance;
		DistPQueue.pop();

		// Updating the distance for neighbors of vertex of lowest distance in priority queue
		auto const& elem = AdjMV[vp1.vId];
		for (auto it = elem.begin(); it != elem.end(); ++it) {
			
			/* Regular Dikjstra */
			double dist = (V.row(vp1.vId) - V.row(*it)).norm();
			//VtoVDist(V.row(vp1.vId), V.row(*it), dist);
			double tempDist = distFromCenter + dist;

			/* Correct distance using Euclidean Distance */
			//double tempDist = (V.row(source) - V.row(*it)).norm();

			/* updating the distance */
			if (tempDist < D(*it)) {
				D(*it) = tempDist;
				VertexPair vp2{ *it,tempDist };
				DistPQueue.push(vp2);
			}
		}
	} while (!DistPQueue.empty());
	//} while (distFromCenter <= nDist);
}

void VectorFields::computeDijkstraDistanceFace(const int &source, Eigen::VectorXd &D)
{
	priority_queue<VertexPair, std::vector<VertexPair>, std::greater<VertexPair>> DistPQueue;

	// Computing distance for initial sample points S
	for (int i = 0; i < F.rows(); i++) {
		D(i) = numeric_limits<double>::infinity();
	}

	D(source) = 0.0f;
	VertexPair vp{ source,D(source) };
	DistPQueue.push(vp);

	double distFromCenter = numeric_limits<double>::infinity();

	// For other vertices in mesh
	do {
		if (DistPQueue.size() == 0) break;
		VertexPair vp1 = DistPQueue.top();
		distFromCenter = vp1.distance;
		DistPQueue.pop();

		// Updating the distance for neighbors of vertex of lowest distance in priority queue
		//auto const& elem = AdjMF3N_temp[vp1.vId];
		int const elem = vp1.vId;
		Eigen::Vector3d const c1 = (V.row(F(elem, 0)) + V.row(F(elem, 1)) + V.row(F(elem, 2))) / 3.0; 
		for (auto it = 0; it != F.cols(); ++it) {
			/* Regular Dikjstra */
			const int neigh = AdjMF3N(elem, it); 
			Eigen::Vector3d const c2 = (V.row(F(neigh, 0)) + V.row(F(neigh, 1)) + V.row(F(neigh, 2))) / 3.0;
			double dist = (c2 - c1).norm();
			double tempDist = distFromCenter + dist;
			
			/* updating the distance */
			if (tempDist < D(neigh)) {
				D(neigh) = tempDist;
				VertexPair vp2{ neigh,tempDist };
				DistPQueue.push(vp2);
			}
		}
	} while (!DistPQueue.empty());
}

void VectorFields::computeDijkstraDistanceFaceMultSource(const Eigen::VectorXi &source, Eigen::VectorXd &D)
{
	priority_queue<VertexPair, std::vector<VertexPair>, std::greater<VertexPair>> DistPQueue;

	// Computing distance for initial sample points S
	for (int i = 0; i < F.rows(); i++) {
		D(i) = numeric_limits<double>::infinity();
	}

	for (int i = 0; i < source.size(); i++) {
		D(source(i)) = 0.0f;
		VertexPair vp{ source(i),D(source(i)) };
		DistPQueue.push(vp);
	}

	double distFromCenter = numeric_limits<double>::infinity();

	// For other vertices in mesh
	do {
		if (DistPQueue.size() == 0) break;
		VertexPair vp1 = DistPQueue.top();
		distFromCenter = vp1.distance;
		DistPQueue.pop();

		// Updating the distance for neighbors of vertex of lowest distance in priority queue
		//auto const& elem = AdjMF3N_temp[vp1.vId];
		int const elem = vp1.vId;
		Eigen::Vector3d const c1 = (V.row(F(elem, 0)) + V.row(F(elem, 1)) + V.row(F(elem, 2))) / 3.0;
		for (auto it = 0; it != F.cols(); ++it) {
			/* Regular Dikjstra */
			const int neigh = AdjMF3N(elem, it);
			Eigen::Vector3d const c2 = (V.row(F(neigh, 0)) + V.row(F(neigh, 1)) + V.row(F(neigh, 2))) / 3.0;
			double dist = (c2 - c1).norm();
			double tempDist = distFromCenter + dist;

			/* updating the distance */
			if (tempDist < D(neigh)) {
				D(neigh) = tempDist;
				VertexPair vp2{ neigh,tempDist };
				DistPQueue.push(vp2);
			}
		}
	} while (!DistPQueue.empty());
}

void VectorFields::computeDijkstraDistanceFaceForSampling(const int &source, Eigen::VectorXd &D)
{
	priority_queue<VertexPair, std::vector<VertexPair>, std::greater<VertexPair>> DistPQueue;	

	D(source) = 0.0f;
	VertexPair vp{ source,D(source) };
	DistPQueue.push(vp);

	double distFromCenter = numeric_limits<double>::infinity();

	// For other vertices in mesh
	do {
		if (DistPQueue.size() == 0) break;
		VertexPair vp1 = DistPQueue.top();
		distFromCenter = vp1.distance;
		DistPQueue.pop();

		// Updating the distance for neighbors of vertex of lowest distance in priority queue
		//auto const& elem = AdjMF3N_temp[vp1.vId];
		int const elem = vp1.vId;
		Eigen::Vector3d const c1 = (V.row(F(elem, 0)) + V.row(F(elem, 1)) + V.row(F(elem, 2))) / 3.0;
		for (auto it = 0; it != F.cols(); ++it) {
			/* Regular Dikjstra */
			const int neigh = AdjMF3N(elem, it);
			Eigen::Vector3d const c2 = (V.row(F(neigh, 0)) + V.row(F(neigh, 1)) + V.row(F(neigh, 2))) / 3.0;
			double dist = (c2 - c1).norm();
			double tempDist = distFromCenter + dist;

			/* updating the distance */
			if (tempDist < D(neigh)) {
				D(neigh) = tempDist;
				VertexPair vp2{ neigh,tempDist };
				DistPQueue.push(vp2);
			}
		}
	} while (!DistPQueue.empty());
}

void VectorFields::computeEigenLaplace2D()
{
	//computeEigenMatlab(SF2D, MF2D, EigVects, eigVals);
	cout << "Eigenvalues: " << eigVals << endl << endl;
	printf("Dimension of eigenvectors %dx%d\n.", EigVects.rows(), EigVects.cols());
}

void VectorFields::computeEigenLaplace3D()
{
	//computeEigenMatlab(SF3D, MF3D, EigVects, eigVals);
	cout << "Eigenvalues: " << eigVals << endl; 
}

void VectorFields::computeEigenstructureGradient3D()
{
	//computeEigenMatlab(SV, MV, EigVectsDiv, eigValsDiv);
}

void VectorFields::constructLaplace2D()
{
	L2D = MF2Dinv * SF2D;
	laplaceGradArbField2D = L2D * gradArbField2D;
}

void VectorFields::constructLaplace3D()
{
	L3D = MF3Dinv * SF3D; 
	laplaceGradArbField3D = L3D * gradArbField3D; 
}

void VectorFields::computeEdges()
{
	igl::edges(F, E);

	printf("....E=%dx%d\n", E.rows(), E.cols());
}

void VectorFields::constructVFNeighbors()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Building \"Vertex-to-Face\" Adjacency (1)... ";

	VFNeighbors.resize(V.rows());
	int counter = 0; 

	for (int i = 0; i < F.rows(); i++) {
		for (int j = 0; j < F.cols(); j++) {
			VtoFPair vp1{counter++, F(i, j), i };
			for (int k = j + 1; k < F.cols(); k++) {
				VtoFPair vp2{counter++, F(i, k), i };
				if (vp2.vId > vp1.vId)	
					VFNeighbors[vp1.vId].insert(vp2);
				else						
					VFNeighbors[vp2.vId].insert(vp1);
			}
		}
	}

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}

void VectorFields::constructVFNeighborsFull()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Building \"Vertex-to-Face\" Adjacency (2)... ";

	VFNeighFull.resize(V.rows());
	int counter = 0;

	for (int i = 0; i < F.rows(); i++) {
		for (int j = 0; j < F.cols(); j++) {
			//VtoFPair vp1{ counter++, F(i, j), i };
			VtoFPair vp1{ i, F(i, j), i };
			for (int k = j + 1; k < F.cols(); k++) {
				//VtoFPair vp2{ counter++, F(i, k), i };
				VtoFPair vp2{ i, F(i, k), i };
				VFNeighFull[vp1.vId].insert(vp2);
				VFNeighFull[vp2.vId].insert(vp1);
			}
		}
	}

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}

void VectorFields::constructSubdomainSingle(const int &source)
{
	sample = source; 

	priority_queue<VertexPair, std::vector<VertexPair>, std::greater<VertexPair>> DistPQueue;
	Eigen::VectorXd D(F.rows());
	const double maxDist = 4.0 * avgEdgeLength; 

	// Computing distance for initial sample points S
	for (int i = 0; i < F.rows(); i++) {
		D(i) = numeric_limits<double>::infinity();
	}

	D(source) = 0.0f;
	VertexPair vp{ source,D(source) };
	DistPQueue.push(vp);
	SubDomain.insert(source);

	double distFromCenter = numeric_limits<double>::infinity();

	// For other vertices in mesh
	do {
		if (DistPQueue.size() == 0) break;
		VertexPair vp1 = DistPQueue.top();
		distFromCenter = vp1.distance;
		DistPQueue.pop();

		// Updating the distance for neighbors of vertex of lowest distance in priority queue
		//auto const& elem = AdjMF3N_temp[vp1.vId];
		int const elem = vp1.vId;
		Eigen::Vector3d const c1 = (V.row(F(elem, 0)) + V.row(F(elem, 1)) + V.row(F(elem, 2))) / 3.0;
		for (auto it = 0; it != F.cols(); ++it) {
			/* Regular Dikjstra */
			const int neigh = AdjMF3N(elem, it);
			Eigen::Vector3d const c2 = (V.row(F(neigh, 0)) + V.row(F(neigh, 1)) + V.row(F(neigh, 2))) / 3.0;
			double dist = (c2 - c1).norm();
			double tempDist = distFromCenter + dist;

			/* updating the distance */
			if (tempDist < D(neigh)) {
				D(neigh) = tempDist;
				VertexPair vp2{ neigh,tempDist };
				DistPQueue.push(vp2);
				SubDomain.insert(neigh);
			}
		}
	} while (distFromCenter < maxDist );
	//} while (!DistPQueue.empty());

	cout << "There are " << SubDomain.size() << " elements in this sub-domain." << endl; 
}

void VectorFields::constructBoundary()
{
	// Obtaining the OUTER-most ELEMENTS
	set<int> outerPart;
	set<int> innerBoundary, outerBoundary;
	for (std::set<int>::iterator it = SubDomain.begin(); it != SubDomain.end(); ++it) {
		for (int i = 0; i < F.cols(); i++) {
			if (SubDomain.find(AdjMF3N(*it, i)) == SubDomain.end()) {
				outerPart.insert(*it);
			}
		}
	}

	// Obtaining the BOUNDARY
	for (std::set<int>::iterator it = outerPart.begin(); it != outerPart.end(); ++it) {
		for (std::set<int>::iterator jt = AdjMF2Ring[*it].begin(); jt != AdjMF2Ring[*it].end(); ++jt) {
			if (SubDomain.find(*jt) == SubDomain.end()) {
				innerBoundary.insert(*jt);
				Boundary.insert(*jt);
			}
		}
	}

	for (std::set<int>::iterator it = innerBoundary.begin(); it != innerBoundary.end(); ++it) {
		for (std::set<int>::iterator jt = AdjMF2Ring[*it].begin(); jt != AdjMF2Ring[*it].end(); ++jt) {
			if (SubDomain.find(*jt) == SubDomain.end()) {
				outerBoundary.insert(*jt);
				Boundary.insert(*jt);
			}
		}
	}


	cout << "Such subdomain has " << Boundary.size() << " elements in its boundary." << endl;
}

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

void VectorFields::computeGradArbField3D()
{
	// Compute the Gradient
	gradArbField3D = GF3D * arbField;
}

void VectorFields::computeGradArbField2D()
{
	gradArbField2D = GF2D * arbField; 
}

void VectorFields::computeCoGradArbField2D()
{
	coGradArbField2D = J * gradArbField2D;
}

void VectorFields::computeCoGradArbField3D()
{
	if (coGradArbField2D.size() == 0) {
		coGradArbField2D = J * gradArbField2D;
	}
	printf("A^T=%dx%d, coGradField=%dx%d\n", A.rows(), A.cols(), coGradArbField2D.rows(), coGradArbField2D.cols());
	coGradArbField3D = A * coGradArbField2D; 
}

void VectorFields::computeCurlGradArbField3D()
{
	curlGradArbField3D = Curl3D * gradArbField3D;

	//cout << curlGradArbField3D << endl; 
}

void VectorFields::computeCurlGradArbField2D()
{
	curlGradArbField2D = Curl2D * gradArbField2D;
	cout << "Curl of Gradient field " << curlGradArbField2D << endl; 
}

void VectorFields::computeCurlCoGradArbField3D()
{
	curlCoGradArbField3D = Curl3D * coGradArbField3D;
}

void VectorFields::computeCurlCoGradArbField2D()
{
	curlCoGradArbField2D = Curl2D * coGradArbField2D;

	//cout << "CURL: " << endl << curlCoGradArbField2D << endl; 
}

void VectorFields::computeDivGradArbField3D()
{
	printf("Div=%dx%d, gradField=%dx%d\n", Div3D.rows(), Div3D.cols(), gradArbField3D.rows(), gradArbField3D.cols());
	divGradArbField3D = Div3D * gradArbField3D;

	//cout << divGradArbField3D << endl;
}

void VectorFields::computeDivGradArbField2D()
{
	printf("Div=%dx%d, gradField=%dx%d\n", Div2D.rows(), Div2D.cols(), gradArbField2D.rows(), gradArbField2D.cols());
	divGradArbField2D = Div2D * gradArbField2D;

	//cout << divGradArbField2D << endl;
}

void VectorFields::computeDivCoGradArbField3D()
{
	divCoGradArbField3D = Div3D * gradArbField3D;
}

void VectorFields::computeDivCoGradArbField2D()
{
	divCoGradArbField2D = Div2D * coGradArbField2D;

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
	const int idx = rand()%F.rows();
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
	double gradEnergy = gradArbField3D.transpose() * LapCurl3D * gradArbField3D;
	cout << "The energy is " << gradEnergy << endl; 
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
			if (it.row() == 200 && it.col()%2==0) {
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
		cout << i  << ", ";
	}
	cout << endl;
	
}

/* ====================== VISUALIZATION of TESTING ============================*/
void VectorFields::visualizeGradient3DArbField(igl::opengl::glfw::Viewer &viewer) 
{
	visualize3Dfields(viewer, gradArbField3D, Eigen::RowVector3d(0.9, 0.1, 0.2));

	/* TEST */
	//Eigen::Vector3d e;
	//Eigen::VectorXd gradLength(F.rows());
	//Eigen::MatrixXd vColor, GradVector(F.rows(), F.cols());
	//double totalGrad = 0.0, averageGrad;

	//for (int i = 0; i < F.rows(); i++){
	//	Eigen::RowVector3d c, g, v1, v2, v3;
	//	c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;		// center of each face
	//	//c = V.row(F(i, 0));													// first vertex of each face [NOT GOOD]
	//	g = gradArbField3D.block(3 * i, 0, 3, 1).transpose(); 
	//	GradVector.row(i) = g;
	//	gradLength(i) = g.norm();
	//	totalGrad += gradLength(i);
	//}

	//averageGrad = totalGrad / (double)F.rows();
	//double lengthScale = 1.0*avgEdgeLength / averageGrad;

	//for (int i = 0; i < F.rows(); i++)
	//{
	//	Eigen::RowVector3d c;
	//	c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;
	//	viewer.data().add_edges(c, c + GradVector.row(i)*lengthScale, Eigen::RowVector3d(0.9, 0.1, 0.2));
	//}

	/* END OF TEST */

	//igl::jet(arbField, true, vColor);
	//igl::parula(arbField, true, vColor);
	//viewer.data().set_colors(vColor);

	// Central point
	//Eigen::RowVectorXd centerV = (V.row(F(*(NeighRing[0].begin()), 0)) + V.row(F(*(NeighRing[0].begin()), 1)) + V.row(F(*(NeighRing[0].begin()), 2))) / 3.0;
	//viewer.data().add_points(centerV, Eigen::RowVector3d(1.0, 0.1, 0.0));
	//viewer.data().add_points(V.row(0), Eigen::RowVector3d(1.0, 0.1, 0.0));
}

void VectorFields::visualizeGradient2DArbField(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::Vector3d e;
	Eigen::VectorXd gradLength(F.rows());
	Eigen::MatrixXd vColor, GradVector(F.rows(), F.cols());	
	double totalGrad = 0.0, averageGrad;

	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::RowVector3d c, g, v1, v2, v3;
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;		// center of each face
		//c = V.row(F(i, 0));													// first vertex of each face [NOT GOOD]	
		g = (A.block(3*i,2*i,3,2) * gradArbField2D.block(2 * i, 0, 2, 1)).transpose();
		GradVector.row(i) = g;
		gradLength(i) = g.norm();
		totalGrad += gradLength(i);
	}

	averageGrad = totalGrad / (double)F.rows();
	double lengthScale = 1.0*avgEdgeLength / averageGrad;

	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::RowVector3d c;
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;
		viewer.data().add_edges(c, c + GradVector.row(i)*lengthScale, Eigen::RowVector3d(1.0, 0.1, 0.2));
	}

	//igl::jet(arbField, true, vColor);
	//viewer.data().set_colors(vColor);

	//viewer.data().add_points(V.row(0), Eigen::RowVector3d(1.0, 0.1, 0.0));
}

void VectorFields::visualizeCoGrad3DArbField(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::Vector3d e;
	Eigen::VectorXd gradLength(F.rows());
	Eigen::MatrixXd vColor, GradVector(F.rows(), F.cols());
	double totalGrad = 0.0, averageGrad;

	for (int i = 0; i < F.rows(); i++) {
		Eigen::RowVector3d c, g, v1, v2, v3;
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;		// center of each face
																			//c = V.row(F(i, 0));													// first vertex of each face [NOT GOOD]
		g = coGradArbField3D.block(3 * i, 0, 3, 1).transpose();
		GradVector.row(i) = g;
		gradLength(i) = g.norm();
		totalGrad += gradLength(i);
	}

	averageGrad = totalGrad / (double)F.rows();
	double lengthScale = 1.0*avgEdgeLength / averageGrad;

	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::RowVector3d c;
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;
		viewer.data().add_edges(c, c + GradVector.row(i)*lengthScale, Eigen::RowVector3d(0.2, 0.1, 1.0));
	}
}

void VectorFields::visualizeCoGrad2DArbField(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::Vector3d e;
	Eigen::VectorXd gradLength(F.rows());
	Eigen::MatrixXd vColor, GradVector(F.rows(), F.cols());
	double totalGrad = 0.0, averageGrad;

	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::RowVector3d c, g, v1, v2, v3;
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;		// center of each face
																			//c = V.row(F(i, 0));													// first vertex of each face [NOT GOOD]	
		g = (A.block(3 * i, 2 * i, 3, 2) * coGradArbField2D.block(2 * i, 0, 2, 1)).transpose();
		GradVector.row(i) = g;
		gradLength(i) = g.norm();
		totalGrad += gradLength(i);
	}

	averageGrad = totalGrad / (double)F.rows();
	double lengthScale = 1.0*avgEdgeLength / averageGrad;

	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::RowVector3d c;
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;
		viewer.data().add_edges(c, c + GradVector.row(i)*lengthScale, Eigen::RowVector3d(0.2, 0.1, 1.0));
	}

	//igl::jet(arbField, true, vColor);
	//viewer.data().set_colors(vColor);

	viewer.data().add_points(V.row(0), Eigen::RowVector3d(1.0, 0.1, 0.0));
}

void VectorFields::visualizeCurlGrad3DArbField(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::MatrixXd vColor;
	//igl::jet(curlGradArbField3D, true, vColor);
	igl::parula(curlGradArbField3D, false, vColor);
	viewer.data().set_colors(vColor);

	viewer.data().add_points(V.row(0), Eigen::RowVector3d(1.0, 0.1, 0.0));
}

void VectorFields::visualizeCurlGrad2DArbField(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::MatrixXd vColor;
	//igl::jet(curlGradArbField3D, true, vColor);
	igl::parula(curlGradArbField2D, false, vColor);
	viewer.data().set_colors(vColor);
}

void VectorFields::visualizeCurlCoGrad3DArbField(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::MatrixXd vColor;
	igl::parula(curlCoGradArbField3D, false, vColor);
	viewer.data().set_colors(vColor);
}

void VectorFields::visualizeCurlCoGrad2DArbField(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::MatrixXd vColor;
	igl::jet(curlCoGradArbField2D, true, vColor);
	viewer.data().set_colors(vColor);
}

void VectorFields::visualizeLaplaceGrad2DArbField(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::Vector3d e;
	Eigen::VectorXd gradLength(F.rows());
	Eigen::MatrixXd vColor, GradVector(F.rows(), F.cols());
	double totalGrad = 0.0, averageGrad;

	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::RowVector3d c, g, v1, v2, v3;
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;		// center of each face
																			//c = V.row(F(i, 0));													// first vertex of each face [NOT GOOD]	
		g = (A.block(3 * i, 2 * i, 3, 2) * laplaceGradArbField2D.block(2 * i, 0, 2, 1)).transpose();
		GradVector.row(i) = g;
		gradLength(i) = g.norm();
		totalGrad += gradLength(i);
	}

	averageGrad = totalGrad / (double)F.rows();
	double lengthScale = 1.0*avgEdgeLength / averageGrad;

	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::RowVector3d c;
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;
		viewer.data().add_edges(c, c + GradVector.row(i)*lengthScale, Eigen::RowVector3d(0.2, 0.1, 0.2));
	}

	igl::jet(arbField, true, vColor);
	viewer.data().set_colors(vColor);

	viewer.data().add_points(V.row(0), Eigen::RowVector3d(1.0, 0.1, 0.0));
}

void VectorFields::visualizeLaplaceGrad3DArbField(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::Vector3d e;
	Eigen::VectorXd gradLength(F.rows());
	Eigen::MatrixXd vColor, GradVector(F.rows(), F.cols());
	double totalGrad = 0.0, averageGrad;

	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::RowVector3d c, g, v1, v2, v3;
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;		// center of each face
																			//c = V.row(F(i, 0));													// first vertex of each face [NOT GOOD]	
		//g = (A.block(3 * i, 2 * i, 3, 2) * laplaceGradArbField2D.block(2 * i, 0, 2, 1)).transpose();
		g = laplaceGradArbField3D.block(3 * i, 0, 3, 1).transpose();
		GradVector.row(i) = g;
		gradLength(i) = g.norm();
		totalGrad += gradLength(i);
	}

	averageGrad = totalGrad / (double)F.rows();
	double lengthScale = 1.0*avgEdgeLength / averageGrad;

	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::RowVector3d c;
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;
		viewer.data().add_edges(c, c + GradVector.row(i)*lengthScale, Eigen::RowVector3d(0.2, 0.1, 0.2));
	}

	igl::jet(arbField, true, vColor);
	viewer.data().set_colors(vColor);

	//viewer.data().add_points(V.row(0), Eigen::RowVector3d(1.0, 0.1, 0.0));
}

void VectorFields::visualizeDivGrad3DArbField(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::MatrixXd vColor;
	//igl::jet(divGradArbField3D, true, vColor);
	igl::parula(divGradArbField3D, false, vColor);
	viewer.data().set_colors(vColor);
}

void VectorFields::visualizeDivCoGrad2DArbField(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::MatrixXd vColor;
	igl::jet(divCoGradArbField2D, false, vColor);
	viewer.data().set_colors(vColor);
}

void VectorFields::visualizeDivCoGrad3DArbField(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::MatrixXd vColor;
	igl::jet(divCoGradArbField3D, true, vColor);
	viewer.data().set_colors(vColor);
}

void VectorFields::visualizeDivGrad2DArbField(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::MatrixXd vColor;
	//igl::jet(divGradArbField2D, true, vColor);
	igl::parula(divGradArbField2D, true, vColor);
	viewer.data().set_colors(vColor);
}

void VectorFields::visualizeFaceNeighbors(igl::opengl::glfw::Viewer &viewer, const int &idx){
	Eigen::VectorXd z(F.rows());
	Eigen::MatrixXd vColor;

	for (int i = 0; i < F.rows(); i++) {
		z(i) = 0.0;
	}
	z(idx) = 1.0;		

	//for (int i = 0; i < AdjMF3N.cols(); i++) {
	//	if(i==0)
	//	z(AdjMF3N(idx, i)) = 0.5; 
	//}

	for (std::set<int, double>::iterator it = AdjMF2Ring[idx].begin(); it != AdjMF2Ring[idx].end(); ++it) {
		z(*it) = 0.5f;
	}

	//for (Eigen::SparseMatrix<double>::InnerIterator it(AjdMF_igl, idx); it; ++it) {
	//	z(it.row()) = 0.5;
	//	printf("col=%d, row=%d\n", it.col(), it.row());
	//}

	// Visualizing the COLOR
	igl::jet(z, true, vColor);
	printf("Size of color: %dx%d\n", vColor.rows(), vColor.cols());
	viewer.data().set_colors(vColor);

	// Visualilzing the EDGE
	for(int i=0; i<F.cols(); i++){
		if (i == 0)
		viewer.data().add_edges(V.row(EdgePairMatrix(idx, 2 * i)), V.row(EdgePairMatrix(idx, 2 * i + 1)), Eigen::RowVector3d(0.2, 0.1, 0.2));
	}
}

void VectorFields::visualizeVertexFacesNeighbors(igl::opengl::glfw::Viewer &viewer, const int &idx)
{
	Eigen::VectorXd z(F.rows());
	Eigen::MatrixXd vColor;

	for (int i = 0; i < F.rows(); i++) {
		z(i) = 0.0;
	}

	for (std::set<VtoFPair>::iterator it = VFNeighFull[idx].begin(); it != VFNeighFull[idx].end(); ++it) {
		z(it->fId) = 0.5; 
	}

	igl::jet(z, true, vColor);
	viewer.data().set_colors(vColor);

	viewer.data().add_points(V.row(idx), Eigen::RowVector3d(0.3, 0.3, 0.9));
}

void VectorFields::visualizeNeighboringRings(igl::opengl::glfw::Viewer &viewer) {
	// Data for coloring
	Eigen::VectorXd z(F.rows());
	Eigen::MatrixXd vColor;
	for (int i = 0; i < F.rows(); i++) {
		z(i) = 0.0;
	}

	Eigen::RowVectorXd centerV = (V.row(F(*(NeighRing[0].begin()), 0)) + V.row(F(*(NeighRing[0].begin()), 1)) + V.row(F(*(NeighRing[0].begin()), 2))) / 3.0;
	viewer.data().add_points(centerV, Eigen::RowVector3d(1.0, 0.1, 0.0));

	// Iterate through all elements
	int k = NeighRing.size();
	double delta = 1.0 / (double)k;

	for (int i = 0; i < k; i++) {
		//if (i % 2 == 1) continue; 
		for (std::set<int, double>::iterator it = NeighRing[i].begin(); it != NeighRing[i].end(); ++it) {
			z(*it) = 1.0 - (double)(i + 1)*delta;
			if (i % 2 == 1) z(*it) = 2.0 - (double)(i + 1)*delta;
		}
	}

	igl::jet(z, true, vColor);
	viewer.data().set_colors(vColor);
}

void VectorFields::visualizeDijkstra(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::MatrixXd vColor;
	igl::jet(arbField, true, vColor);
	viewer.data().set_colors(vColor);
}

void VectorFields::visualizeEigenfields(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::Vector3d e;
	Eigen::VectorXd eigfieldLength(F.rows());
	Eigen::MatrixXd vColor, EigfieldVector(F.rows(), F.cols());
	double totalEigfield = 0.0, averageEigfield;
	const int eigfieldID = 1;

	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::RowVector3d c, g, v1, v2, v3;
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;		// center of each face												// first vertex of each face [NOT GOOD]	
		//g = (A.block(3 * i, 2 * i, 3, 2) * EigVects.block(2 * i, eigfieldID, 2, 1)).transpose();
		g = EigVects.block(3 * i, eigfieldID, 3, 1).transpose();
		EigfieldVector.row(i) = g;
		eigfieldLength(i) = g.norm();
		totalEigfield += eigfieldLength(i);
	}

	averageEigfield = totalEigfield / (double)F.rows();
	double lengthScale = 1.0*avgEdgeLength / averageEigfield;

	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::RowVector3d c;
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;
		viewer.data().add_edges(c, c + EigfieldVector.row(i)*lengthScale, Eigen::RowVector3d(0.2, 0.1, 0.2));
		//viewer.data().add_edges(c, c + EigfieldVector.row(i).normalized()*avgEdgeLength, Eigen::RowVector3d(0.2, 0.1, 0.2));
	}

	//igl::jet(arbField, true, vColor);
	//viewer.data().set_colors(vColor);

	//viewer.data().add_points(V.row(0), Eigen::RowVector3d(1.0, 0.1, 0.0));
}

void VectorFields::visualizeArbField(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::MatrixXd vColor;
	igl::jet(arbField, true, vColor);
	viewer.data().set_colors(vColor);
}

void VectorFields::visualizeEigFieldsDiv(igl::opengl::glfw::Viewer &viewer, const int &eigID)
{
	Eigen::VectorXd eigField = GF3D * EigVectsDiv.col(eigID);
	visualize3Dfields(viewer, eigField, Eigen::RowVector3d(0.2, 0.1, 0.2));
}
void VectorFields::visualizeRandomFace(igl::opengl::glfw::Viewer &viewer, const int &faceID)
{
	viewer.data().clear();
	viewer.data().set_mesh(V, F);

	Eigen::RowVector3d c;
	c = (V.row(F(faceID, 0)) + V.row(F(faceID, 1)) + V.row(F(faceID, 2))) / 3.0;
	viewer.data().add_points(c, Eigen::RowVector3d(1.0, 0.1, 0.0));
}

void VectorFields::visualizeDijkstraFace(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::MatrixXd vColor;
	igl::jet(dijkstraFace, true, vColor);
	viewer.data().set_colors(vColor);
}

void VectorFields::visualizeSubdomain(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::VectorXd dom(F.rows());
	for (int i = 0; i < F.rows(); i++) dom(i) = 0.0; 

	for (std::set<int>::iterator it = SubDomain.begin(); it != SubDomain.end(); ++it) {
		dom(*it) = 0.5;
		if (*it == 542) dom(*it) = 1.0;
	}

	for (std::set<int>::iterator it = Boundary.begin(); it != Boundary.end(); ++it) {
		dom(*it) = 0.25;
	}
	
	
	Eigen::MatrixXd vColor;
	igl::jet(dom, true, vColor);
	viewer.data().set_colors(vColor);
}

void VectorFields::visualizeSamples(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::RowVector3d color(0.1, 0.1, 0.8);
	Eigen::RowVector3d c; 

	for (int i : Sample) {
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;
		viewer.data().add_points(c, color);
	}
}

/* ====================== MESH-RELATED FUNCTIONS ============================*/
void VectorFields::readMesh(const string &meshFile)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();	
	cout << "> Reading mesh... ";

	// For actual work of reading mesh object
	V.resize(0, 0);
	F.resize(0, 0);

	if (meshFile.substr(meshFile.find_last_of(".") + 1) == "off") {
		igl::readOFF(meshFile, V, F);
	}
	else if (meshFile.substr(meshFile.find_last_of(".") + 1) == "obj") {
		igl::readOBJ(meshFile, V, F);
	}
	else {
		cout << "Error! File type can be either .OFF or .OBJ only." << endl;
		cout << "Program will exit in 2 seconds." << endl;
		Sleep(2000);
		exit(10);
	}

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;

	// Printing Mesh-related information
	printf("....V=%dx%d\n", V.rows(), V.cols());
	printf("....F=%dx%d\n", F.rows(), F.cols());
}

void VectorFields::getVF(Eigen::MatrixXd &V, Eigen::MatrixXi &F)
{
	V = this->V;
	F = this->F;
}

void VectorFields::computeFaceCenter()
{
	FC.resize(F.rows(), 3);

	for (int i = 0; i < F.rows(); i++) {
		FC.row(i) = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;
	}

}
/* ====================== VISUALIZATION ============================*/
void VectorFields::visualizeMassMatrix(igl::opengl::glfw::Viewer &viewer, const MassMatrixToShow &type)
{
	Eigen::VectorXd z;
	Eigen::MatrixXd vColor;

	switch (type)
	{
	case MASS_MV:
		z.resize(V.rows());
		for (int i = 0; i < V.rows(); i++) {
			z(i) = MV.coeff(i, i);
		}
		break;

	case MASS_MVinv:
		z.resize(V.rows());
		for (int i = 0; i < V.rows(); i++) {
			z(i) = MVinv.coeff(i, i);
		}
		break;

	case MASS_MF2D:
		z.resize(F.rows());
		for (int i = 0; i < F.rows(); i++) {
			z(i) = MF2D.coeff(2 * i + 1, 2 * i + 1);
		}
		break; 

	case MASS_MF2Dinv:
		z.resize(F.rows());
		for (int i = 0; i < F.rows(); i++) {
			z(i) = MF2Dinv.coeff(2 * i + 1, 2 * i + 1);
		}
		break;

	case MASS_MF3D:
		z.resize(F.rows());
		for (int i = 0; i < F.rows(); i++) {
			z(i) = MF3D.coeff(3 * i + 2, 3 * i + 2);
		}
		break;

	case MASS_MF3Dinv:
		z.resize(F.rows());
		for (int i = 0; i < F.rows(); i++) {
			z(i) = MF3Dinv.coeff(3 * i + 2, 3 * i + 2);
		}
		break;

	default:
		break;
	}

	igl::jet(z, true, vColor);
	viewer.data().set_colors(vColor);
}

void VectorFields::visualizeGradient(igl::opengl::glfw::Viewer &viewer, const GradientToShow &type)
{
	switch (type)
	{
	case GRAD_2D:

		break;

	case GRAD_3D:

		break;

	default:
		break;
	}
}

void VectorFields::visualizeLocalFrames(igl::opengl::glfw::Viewer &viewer)
{
	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::RowVector3d c, e, f;
		e << A.coeffRef(3 * i, 2 * i),		A.coeffRef(3 * i + 1, 2 * i),		A.coeffRef(3 * i + 2, 2 * i);
		f << A.coeffRef(3 * i, 2 * i + 1),	A.coeffRef(3 * i + 1, 2 * i + 1),	A.coeffRef(3 * i + 2, 2 * i + 1);
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;

		// First basis
		viewer.data().add_edges(c, c + avgEdgeLength*e, Eigen::RowVector3d(1.0, 0.1, 0.2));
		// Second basis
		viewer.data().add_edges(c, c + avgEdgeLength*f, Eigen::RowVector3d(0.0, 0.1, 0.9));
	}
}

void VectorFields::visualizeApproximatedFields(igl::opengl::glfw::Viewer &viewer)
{
	viewer.data().clear();
	viewer.data().set_mesh(V, F);
	Eigen::RowVector3d color1 = Eigen::RowVector3d(1.0, 0.1, 0.2);
	Eigen::RowVector3d color2 = Eigen::RowVector3d(0.0, 0.1, 1.0);
	visualize2DfieldsScaled(viewer, Xf.col(0), color1);
	//visualize2Dfields(viewer, Xf.col(0), color1);
	//visualize2Dfields(viewer, Xf.col(1), color2);	
}

void VectorFields::visualize2Dfields(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color)
{
	Eigen::Vector3d e;
	Eigen::VectorXd blockLength(F.rows());
	Eigen::MatrixXd vColor, VectorBlock(F.rows(), F.cols());
	double totalGrad = 0.0, avgField;

	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::RowVector3d c, g, v1, v2, v3;
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;		// center of each face
																			//c = V.row(F(i, 0));													// first vertex of each face [NOT GOOD]	
		g = (A.block(3 * i, 2 * i, 3, 2) * field2D.block(2 * i, 0, 2, 1)).transpose();
		VectorBlock.row(i) = g;
		blockLength(i) = g.norm();
		totalGrad += blockLength(i);
	}


	avgField = totalGrad / (double)F.rows();
	//double lengthScale = 1.0*avgEdgeLength / avgField;
	double lengthScale = 1.0*avgEdgeLength;

	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::RowVector3d c;
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;	

		//if (i == *(NeighRing[0].begin())){
		//	viewer.data().add_edges(c, c + VectorBlock.row(i)*lengthScale*2.0, Eigen::RowVector3d(0.1, 0.1, 0.2));
		//	//viewer.data().add_edges(c, c + VectorBlock.row(i)*lengthScale*2.0, Eigen::RowVector3d(1.1, 0.1, 0.2));
		//}else 
		{
			viewer.data().add_edges(c, c + VectorBlock.row(i)*lengthScale, color);
			//viewer.data().add_edges(c, c + VectorBlock.row(i)*lengthScale, Eigen::RowVector3d(1.0, 0.1, 0.2));
			//viewer.data().add_edges(c, c + VectorBlock.row(i).normalized()*avgEdgeLength, Eigen::RowVector3d(1.0, 0.1, 0.2));
		}
		
	}
}

void VectorFields::visualize2DfieldsNormalized(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color)
{
	Eigen::Vector3d e;
	Eigen::VectorXd blockLength(F.rows());
	Eigen::MatrixXd vColor, VectorBlock(F.rows(), F.cols());
	double totalGrad = 0.0, avgField;

	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::RowVector3d c, g, v1, v2, v3;
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;		// center of each face
																			//c = V.row(F(i, 0));													// first vertex of each face [NOT GOOD]	
		g = (A.block(3 * i, 2 * i, 3, 2) * field2D.block(2 * i, 0, 2, 1)).transpose();
		VectorBlock.row(i) = g;
	}

	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::RowVector3d c;
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;
		viewer.data().add_edges(c, c + VectorBlock.row(i).normalized()*avgEdgeLength, color);
	}
}

void VectorFields::visualize2DfieldsScaled(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color)
{
	Eigen::Vector3d e;
	Eigen::MatrixXd vColor, VectorBlock(F.rows(), F.cols());

	for (int i = 0; i < F.rows(); i+=1)
	{
		Eigen::RowVector3d c, g, v1, v2, v3;
		//c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;		// center of each face
		c = FC.row(i);
																			//c = V.row(F(i, 0));													// first vertex of each face [NOT GOOD]	
		g = (A.block(3 * i, 2 * i, 3, 2) * field2D.block(2 * i, 0, 2, 1)).transpose();
		VectorBlock.row(i) = g;
	}

	double lengthScale = 1.0*avgEdgeLength;
	for (int i = 0; i < F.rows(); i+=1)
	{
		Eigen::RowVector3d c;
		//c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;
		c = FC.row(i);
		viewer.data().add_edges(c, c + VectorBlock.row(i)*lengthScale, color);
	}
}

void VectorFields::visualize2DfieldsRegular(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color)
{
	Eigen::Vector3d e;
	Eigen::VectorXd blockLength(F.rows());
	Eigen::MatrixXd vColor, VectorBlock(F.rows(), F.cols());
	double totalGrad = 0.0, avgField;

	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::RowVector3d c, g, v1, v2, v3;
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;		// center of each face
																			//c = V.row(F(i, 0));													// first vertex of each face [NOT GOOD]	
		g = (A.block(3 * i, 2 * i, 3, 2) * field2D.block(2 * i, 0, 2, 1)).transpose();
		VectorBlock.row(i) = g;
		blockLength(i) = g.norm();
		totalGrad += blockLength(i);
	}

	avgField = totalGrad / (double)F.rows();
	double lengthScale = 1.0*avgEdgeLength;

	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::RowVector3d c;
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;
		viewer.data().add_edges(c, c + VectorBlock.row(i)*lengthScale, color);
	}
}



void VectorFields::visualize3Dfields(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field3D, const Eigen::RowVector3d &color)
{
	// Reset
	viewer.data().clear();
	viewer.data().set_mesh(V, F);

	Eigen::Vector3d e;
	Eigen::VectorXd blockLength(F.rows());
	Eigen::MatrixXd vColor, VectorBlock(F.rows(), F.cols());
	double totalGrad = 0.0, averageGrad;

	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::RowVector3d c, g, v1, v2, v3;
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;		// center of each face
																			//c = V.row(F(i, 0));													// first vertex of each face [NOT GOOD]	
		g = field3D.block(3 * i, 0, 3, 1).transpose();
		VectorBlock.row(i) = g;
		blockLength(i) = g.norm();
		totalGrad += blockLength(i);
	}

	averageGrad = totalGrad / (double)F.rows();
	double lengthScale = 1.0*avgEdgeLength / averageGrad;

	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::RowVector3d c;
		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;
		//viewer.data().add_edges(c, c + VectorBlock.row(i)*lengthScale, Eigen::RowVector3d(1.0, 0.1, 0.2));
		viewer.data().add_edges(c, c + VectorBlock.row(i)*lengthScale, color);
		//viewer.data().add_edges(c, c + VectorBlock.row(i).normalized()*avgEdgeLength, Eigen::RowVector3d(1.0, 0.1, 0.2));
	}
}

//void VectorFields::visualize2Dfields(igl::opengl::glfw::Viewer &viewer, const Eigen::SparseVector<double> &field2D, const Eigen::RowVector3d &color)
//{
//	Eigen::Vector3d e;
//	Eigen::VectorXd blockLength(F.rows());
//	Eigen::MatrixXd vColor, VectorBlock(F.rows(), F.cols());
//	double totalGrad = 0.0, avgField;
//
//	for (int i = 0; i < F.rows(); i++)
//	{
//		Eigen::RowVector3d c, g, v1, v2, v3;
//		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;		// center of each face
//																			//c = V.row(F(i, 0));													// first vertex of each face [NOT GOOD]	
//		g = (A.block(3 * i, 2 * i, 3, 2) * field2D.block(2 * i, 0, 2, 1)).transpose();
//		VectorBlock.row(i) = g;
//		blockLength(i) = g.norm();
//		totalGrad += blockLength(i);
//	}
//
//	avgField = totalGrad / (double)F.rows();
//	double lengthScale = 1.0*avgEdgeLength; // / avgField;
//
//	for (int i = 0; i < F.rows(); i++)
//	{
//		Eigen::RowVector3d c;
//		c = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;	
//
//		//if (i == *(NeighRing[0].begin())){
//		//	viewer.data().add_edges(c, c + VectorBlock.row(i)*lengthScale*2.0, Eigen::RowVector3d(0.1, 0.1, 0.2));
//		//	//viewer.data().add_edges(c, c + VectorBlock.row(i)*lengthScale*2.0, Eigen::RowVector3d(1.1, 0.1, 0.2));
//		//}else 
//		{
//			viewer.data().add_edges(c, c + VectorBlock.row(i)*lengthScale, color);
//			//viewer.data().add_edges(c, c + VectorBlock.row(i)*lengthScale, Eigen::RowVector3d(1.0, 0.1, 0.2));
//			//viewer.data().add_edges(c, c + VectorBlock.row(i).normalized()*avgEdgeLength, Eigen::RowVector3d(1.0, 0.1, 0.2));
//		}
//		
//	}
//}

void VectorFields::visualizeBasis(igl::opengl::glfw::Viewer &viewer, const int &id)
{
	viewer.data().clear();
	viewer.data().set_mesh(V, F);

	int bId = id; 
	Eigen::RowVector3d color;
	if (id % 2 == 0) {
		color = Eigen::RowVector3d(1.0, 0.1, 0.2);
	}else {
		color = Eigen::RowVector3d(0.0, 0.1, 1.0);
	}

	if (id >= 2 * Sample.size()) {
		bId = 2 * Sample.size() - 1;
	}

	printf("Showing the %d BasisTemp field\n", bId);
	visualize2DfieldsScaled(viewer, BasisTemp.col(bId), color);

	Eigen::RowVector3d const c1 = (V.row(F(Sample[bId/2], 0)) + V.row(F(Sample[bId/2], 1)) + V.row(F(Sample[bId/2], 2))) / 3.0;
	viewer.data().add_points(c1, Eigen::RowVector3d(0.1, 0.1, 0.1));
}

void VectorFields::visualizeBasisNormalized(igl::opengl::glfw::Viewer &viewer, const int &id)
{
	viewer.data().clear();
	viewer.data().set_mesh(V, F);

	int bId = id;
	Eigen::RowVector3d color;
	if (id % 2 == 0) {
		color = Eigen::RowVector3d(1.0, 0.1, 0.2);
	}
	else {
		color = Eigen::RowVector3d(0.0, 0.1, 1.0);
	}

	if (id >= 2 * Sample.size()) {
		bId = 2 * Sample.size() - 1;
	}

	printf("Showing the %d BasisTemp field\n", bId);
	visualize2DfieldsScaled(viewer, Basis.col(bId), color);

	Eigen::RowVector3d const c1 = (V.row(F(Sample[bId / 2], 0)) + V.row(F(Sample[bId / 2], 1)) + V.row(F(Sample[bId / 2], 2))) / 3.0;
	viewer.data().add_points(c1, Eigen::RowVector3d(0.1, 0.1, 0.1));
}

void VectorFields::visualizeBasisSum(igl::opengl::glfw::Viewer &viewer, const int &id)
{
	viewer.data().clear();
	viewer.data().set_mesh(V, F);

	Eigen::RowVector3d color;
	if(id==0)
		color = Eigen::RowVector3d(1.0, 0.4, 0.4);
	else 
		color = Eigen::RowVector3d(0.0, 0.4, 0.9);
		
	visualize2DfieldsScaled(viewer, BasisSum.col(id), color);
	//visualize2DfieldsScaled(viewer, BasisSumN.col(id), color);
	//for (int i = 0; i < Sample.size(); i++) {
	//	Eigen::RowVector3d const c1 = (V.row(F(Sample[i], 0)) + V.row(F(Sample[i], 1)) + V.row(F(Sample[i], 2))) / 3.0;
	//	viewer.data().add_points(c1, Eigen::RowVector3d(0.1, 0.1, 0.1));
	//}
}

void VectorFields::visualizeApproxResult(igl::opengl::glfw::Viewer &viewer, const int &id)
{
	viewer.data().clear();
	viewer.data().set_mesh(V, F);

	Eigen::RowVector3d color;
	if (id == 0)
		color = Eigen::RowVector3d(1.0, 0.4, 0.4);
	else
		color = Eigen::RowVector3d(0.0, 0.4, 0.9);

	//cout << "Size of X_Lifted " << XFullDim.rows() << "x" << XFullDim.cols() << "." << endl; 
	//visualize2DfieldsNormalized(viewer, XFullDim.col(id), color);
	visualize2DfieldsScaled(viewer, XFullDim.col(id), color);
	//visualize2DfieldsRegular(viewer, XFullDim.col(id), color);
}

void VectorFields::visualizeUserConstraints(igl::opengl::glfw::Viewer &viewer)
{
	for (int i = 0; i < userConstraints.size(); i++) {
		Eigen::RowVector3d c, g, v1, v2, v3;
		//c = (V.row(F(userConstraints[i], 0)) + V.row(F(userConstraints[i], 1)) + V.row(F(userConstraints[i], 2))) / 3.0;														// first vertex of each face [NOT GOOD]	
		c = FC.row(userConstraints[i]);
		g = (A.block(3 * userConstraints[i], 2 * userConstraints[i], 3, 2) * cBar.block(2 * i, 0, 2, 1)).transpose();		
		viewer.data().add_edges(c, c + g/g.norm()*avgEdgeLength, Eigen::RowVector3d(0.1, 0.1, 0.2));
		viewer.data().add_points(c, Eigen::RowVector3d(0.1, 0.1, 0.2));
	}
}

void  VectorFields::visualizeGlobalConstraints(igl::opengl::glfw::Viewer &viewer)
{
	for (int i = 0; i < globalConstraints.size(); i++) {
		Eigen::RowVector3d cc, g, v1, v2, v3;
		cc = FC.row(globalConstraints[i]);
		g = (A.block(3 * globalConstraints[i], 2 * globalConstraints[i], 3, 2) * c.block(2 * i, 0, 2, 1)).transpose();
		viewer.data().add_edges(cc, cc + g / g.norm()*avgEdgeLength, Eigen::RowVector3d(0.1, 0.1, 0.2));
		viewer.data().add_points(cc, Eigen::RowVector3d(0.1, 0.1, 0.2));
	}
}

void VectorFields::visualizeSingularitiesConstraints(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::VectorXd z(F.rows());
	Eigen::MatrixXd vColor;

	for (int i = 0; i < F.rows(); i++) {
		z(i) = 0.0;
	}

	for (int id = 0; id < SingNeighCC.size(); id++) {
		const double diff = 0.6 / (double)SingNeighCC[id].size();
		for (int i = 0; i < SingNeighCC[id].size(); i++) {
			z(SingNeighCC[id][i]) = 0.3 + i*diff;
		}
	}
	igl::jet(z, false, vColor);
	//viewer.data().set_colors(vColor);

	for (int i : singularities) {
		viewer.data().add_points(V.row(i), Eigen::RowVector3d(0.1, 0.9, 0.3));
	}
}

//void VectorFields::visualizeSparseMatrixInMatlab(const Eigen::SparseMatrix<double> &M)
//{
//	printf("Size of M=%dx%d\n", M.rows(), M.cols());
//
//	using namespace matlab::engine;
//	Engine *ep;
//	mxArray *MM = NULL, *MS = NULL, *result = NULL, *eigVecResult, *nEigs;
//
//	const int NNZ_M = M.nonZeros();
//	int nnzMCounter = 0;
//
//	double	*srm = (double*)malloc(NNZ_M * sizeof(double));
//	mwIndex *irm = (mwIndex*)malloc(NNZ_M * sizeof(mwIndex));
//	mwIndex *jcm = (mwIndex*)malloc((M.cols() + 1) * sizeof(mwIndex));
//
//	MM = mxCreateSparse(M.rows(), M.cols(), NNZ_M, mxREAL);
//	srm = mxGetPr(MM);
//	irm = mxGetIr(MM);
//	jcm = mxGetJc(MM);
//
//	// Getting matrix M
//	jcm[0] = nnzMCounter;
//	for (int i = 0; i < M.outerSize(); i++) {
//		for (Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it) {
//			srm[nnzMCounter] = it.value();
//			irm[nnzMCounter] = it.row();
//			nnzMCounter++;
//		}
//		jcm[i + 1] = nnzMCounter;
//	}
//
//	// Start Matlab Engine
//	ep = engOpen(NULL);
//	if (!(ep = engOpen(""))) {
//		fprintf(stderr, "\nCan't start MATLAB engine\n");
//		cout << "CANNOT START MATLAB " << endl;
//	}
//	else {
//		cout << "MATLAB STARTS. OH YEAH!!!" << endl;
//	}
//
//	engPutVariable(ep, "M", MM);
//	engEvalString(ep, "spy(M)");
//}