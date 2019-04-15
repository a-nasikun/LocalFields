#include "NRoSyFields.h"

/* Reading data*/
void NRoSyFields::readMesh(const string &filename)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Reading mesh... ";

	// For actual work of reading mesh object
	V.resize(0, 0);
	F.resize(0, 0);

	if (filename.substr(filename.find_last_of(".") + 1) == "off") {
		igl::readOFF(filename, V, F);
	}
	else if (filename.substr(filename.find_last_of(".") + 1) == "obj") {
		igl::readOBJ(filename, V, F);
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

	igl::edges(F, E);
	const int genus = (2 - V.rows() + E.rows() - F.rows()) / 2;
	printf("This model is of genus %d\n", genus);
}

void NRoSyFields::readMesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
	this->V = V;
	this->F = F; 
}

void NRoSyFields::scaleMesh()
{
	cout << "NRoSyFields::Scaling the mesh\n";

	Eigen::RowVectorXd minV(V.rows(), 3);
	Eigen::RowVectorXd maxV(V.rows(), 3);
	Eigen::RowVector3d length;
	Eigen::MatrixXd MatSubs(V.rows(), V.cols());
	Eigen::MatrixXd MatAdd(V.rows(), V.cols());
	double scaleFactor;

	/* Get the min and max coefficients */
	for (int i = 0; i < V.cols(); i++)
	{
		minV(i) = V.col(i).minCoeff();
		maxV(i) = V.col(i).maxCoeff();
		length(i) = maxV(i) - minV(i);
	}

	/* Creating the substraction matrices */
	for (int i = 0; i < V.rows(); i++)
	{
		MatSubs.row(i) = minV;
	}

	scaleFactor = length.maxCoeff();

	/* Translate to the Origin */
	V = V - MatSubs;

	/* Scale w.r.t the longest axis */
	V = V * (1.0 / scaleFactor);


	for (int i = 0; i < V.cols(); i++)
	{
		maxV(i) = V.col(i).maxCoeff();
	}

	/* Creating the Addition matrices */
	for (int i = 0; i < V.rows(); i++)
	{
		MatAdd.row(i) = 0.5*maxV;
	}

	/* Translate s.t. the center is in the Origin */
	V = V - MatAdd;


	for (int i = 0; i < V.cols(); i++)
	{
		minV(i) = V.col(i).minCoeff();
		maxV(i) = V.col(i).maxCoeff();
		length(i) = maxV(i) - minV(i);
	}
	scaleFactor = length.maxCoeff();
}

void NRoSyFields::computeFaceCenter()
{
	FC.resize(F.rows(), 3);

	for (int i = 0; i < F.rows(); i++) {
		FC.row(i) = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;
	}
}



/* Setting up required structure */
void NRoSyFields::constructFrameBasis()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> NRoSyFields::Constructing the basis vector for local frame... ";

	Eigen::Vector3d e, f, n;
	frameBasis.resize(3 * F.rows());
	
	for (int i = 0; i < F.rows(); i++) {
		e = V.row(F(i, 1)) - V.row(F(i, 0));
		e.normalize();

		frameBasis(3 * i + 0) = e(0);
		frameBasis(3 * i + 1) = e(1);
		frameBasis(3 * i + 2) = e(2);
	}

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}

void NRoSyFields::constructMappingMatrix()
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

	Eigen::MatrixXd NF;
	igl::per_face_normals(V, F, NF);

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

/* Creating NRoSyFields */
void NRoSyFields::representingNRoSyFields(const Eigen::MatrixXd& NFields)
{

}

void NRoSyFields::constructNRoSyFields(const int& nRoSy, const Eigen::MatrixXd& NFields)
{
	cout << "> NRoSyFields::Construction (" << nRoSy << ") \n";
	NRoSy = nRoSy;
	scaleMesh();
	computeFaceCenter();
	constructFrameBasis();
	constructMappingMatrix();
	representingNRoSyFields(NFields);
}

void NRoSyFields::constructNRoSyFields(const int& nRoSy, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
	NRoSy = nRoSy;
	this->V = V; 
	this->F = F;
	scaleMesh();
	constructFrameBasis();
}

/* Rep. Vectors and N-RoSy Fields interface */
void NRoSyFields::convertNRoSyToRepVectors(Eigen::VectorXd& vectorFields)
{

}

void NRoSyFields::convertRepVectorsToNRoSy(const Eigen::VectorXd& vectorFields)
{
	theta.resize(F.rows());
	magnitude.resize(F.rows());

	/* Temp variables */
	Eigen::Vector2d v, b(1,0);

	for (int i = 0; i < F.rows(); i++)
	{
		v = vectorFields.block(2 * i, 0, 2, 1);

		magnitude(i) = v.norm();
		v.normalize();
		double angle;
		if(v(1)<0)
			theta(i) = 2*M_PI - acos(b.dot(v)); // / (double)NRoSy;
		else 
			theta(i) = acos(b.dot(v)); // / (double)NRoSy;
	}
}


/* Visualizing the NRoSyFields */
void NRoSyFields::visualizeNRoSyFields(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::Vector2d b(1, 0);
	Eigen::VectorXd col(NRoSy);
	Eigen::MatrixXd color;
	for (int i = 0; i < NRoSy; i++)
	{
		col(i) = i*(1.0 / NRoSy);
	}
	igl::jet(col, true, color);

	for (int i = 0; i < NRoSy; i++)
	{
		Eigen::VectorXd TempFields(2 * F.rows());
		cout << "Drawing the " << i << " fields \n";

		/* Construct rotation matrix*/
		for (int j = 0; j < F.rows(); j++)
		{
			double angle = theta(j) + ((double)i*2.0*M_PI / (double)NRoSy);
			if (j == 0)
			{
				printf("angle 0=%.5f, theta 0=%.5f\n", angle, theta(j));
			}
			Eigen::Matrix2d RotM;
			RotM(0, 0) =  cos(angle);
			RotM(0, 1) = -sin(angle);
			RotM(1, 0) =  sin(angle);
			RotM(1, 1) =  cos(angle);

			TempFields.block(2 * j, 0, 2, 1) = magnitude(j) * RotM * b; 
		}

		visualize2Dfields(viewer, TempFields, color.row(i), 0.3, false);
	}
}


void NRoSyFields::visualize2Dfields(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color, const double& scale, const bool& normalized)
{
	/* For Timing*/
	chrono::high_resolution_clock::time_point	t1, t2, te1, te2, ta1, ta2;
	chrono::duration<double>					duration, da, de;
	t1 = chrono::high_resolution_clock::now();
	//cout << "> Adding edges... ";

	/* Averaga Edge Length */
	Eigen::Vector3d e;
	double totalLength = 0.0;
	for (int i = 0; i < F.rows(); i++) {
		for (int j = 0; j < F.cols(); j++) {
			e = V.row(F(i, (j + 1) % F.cols())) - V.row(F(i, j));
			totalLength += e.norm();
		}
	}
	double avgEdgeLength = totalLength / (double)(F.rows()*F.cols());
	
	/* Some constants for arrow drawing */
	const double HEAD_RATIO = 3.0;
	const double EDGE_RATIO = scale;
	double lengthScale = EDGE_RATIO*avgEdgeLength;
	

	/* Computing the rotation angle for 1:3 ratio of arrow head */
	double rotAngle = M_PI - atan(1.0 / 3.0);
	Eigen::Matrix2d rotMat1, rotMat2;
	//Eigen::SparseMatrix<double> MRot1(2 * FaceToDraw.size(), 2 * FaceToDraw.size()), MRot2(2 * FaceToDraw.size(), 2 * FaceToDraw.size());
	Eigen::SparseMatrix<double> MRot1(2 * F.rows(), 2 * F.rows()), MRot2(2 * F.rows(), 2 * F.rows());
	vector<Eigen::Triplet<double>> R1Triplet, R2Triplet;
	R1Triplet.reserve(2 * 2 * F.rows());
	R2Triplet.reserve(2 * 2 * F.rows());
	Eigen::MatrixXd FCLoc(F.rows(), 3);

	/* Defining the rotation matrix (2-by-2) on the local frame */
	rotMat1 << cos(rotAngle), -sin(rotAngle), sin(rotAngle), cos(rotAngle);
	rotMat2 << cos(-rotAngle), -sin(-rotAngle), sin(-rotAngle), cos(-rotAngle);

	for (int i = 0; i < F.rows(); i++)
	{
		/* Rotation matrix for the first head */
		R1Triplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, rotMat1(0, 0)));
		R1Triplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 0, rotMat1(1, 0)));
		R1Triplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 1, rotMat1(0, 1)));
		R1Triplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, rotMat1(1, 1)));

		/* Rotation matrix for the second head */
		R2Triplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, rotMat2(0, 0)));
		R2Triplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 0, rotMat2(1, 0)));
		R2Triplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 1, rotMat2(0, 1)));
		R2Triplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, rotMat2(1, 1)));

		/* Getting the face center of selected faces */
		FCLoc.row(i) = FC.row(i);
	}
	MRot1.setFromTriplets(R1Triplet.begin(), R1Triplet.end());
	MRot2.setFromTriplets(R2Triplet.begin(), R2Triplet.end());

	/* Getting the local data from the population of data */
	Eigen::SparseMatrix<double> ALoc(3 * F.rows(), 2 * F.rows());
	vector<Eigen::Triplet<double>> ATriplet;
	ATriplet.reserve(6 * F.rows());
	Eigen::VectorXd fieldLoc(2 * F.rows()), fields3D(3 * F.rows()), rot1Field, rot2Field;
	Eigen::MatrixXd TFields(F.rows(), F.cols()), Head1Fields(F.rows(), F.cols()), Head2Fields(F.rows(), F.cols());

	for (int i = 0; i < F.rows(); i++)
	{
		/* Getting the selected ALoc from A */
		ATriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 2 * i + 0, A.coeff(3 * i + 0, 2 * i + 0)));
		ATriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 2 * i + 0, A.coeff(3 * i + 1, 2 * i + 0)));
		ATriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 2 * i + 0, A.coeff(3 * i + 2, 2 * i + 0)));
		ATriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 2 * i + 1, A.coeff(3 * i + 0, 2 * i + 1)));
		ATriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 2 * i + 1, A.coeff(3 * i + 1, 2 * i + 1)));
		ATriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 2 * i + 1, A.coeff(3 * i + 2, 2 * i + 1)));

		/* Getting the selected face */
		fieldLoc.block(2 * i, 0, 2, 1) = field2D.block(2 * i, 0, 2, 1);
	}
	ALoc.setFromTriplets(ATriplet.begin(), ATriplet.end());
	fields3D = ALoc * fieldLoc;

	/* The head of the arrows */
	rot1Field = MRot1*fieldLoc;
	rot1Field = ALoc * rot1Field;
	rot2Field = MRot2*fieldLoc;
	rot2Field = ALoc * rot2Field;

	/* Transform field to Matrix format */
	for (int i = 0; i < F.rows(); i++)
	{
		TFields.row(i) = (fields3D.block(3 * i, 0, 3, 1)).transpose();
		Head1Fields.row(i) = (rot1Field.block(3 * i, 0, 3, 1)).transpose();
		Head2Fields.row(i) = (rot2Field.block(3 * i, 0, 3, 1)).transpose();
	}

	/* If user wants normalized fields, then so do it */
	if (normalized)
	{
		TFields.rowwise().normalize();
		Head1Fields.rowwise().normalize();
		Head2Fields.rowwise().normalize();
	}

	/* Draw the fields */
	viewer.data().add_edges(FCLoc, FCLoc + TFields*lengthScale, color);
	viewer.data().add_edges(FCLoc + TFields*lengthScale, FCLoc + TFields*lengthScale + Head1Fields*lengthScale / HEAD_RATIO, color);
	viewer.data().add_edges(FCLoc + TFields*lengthScale, FCLoc + TFields*lengthScale + Head2Fields*lengthScale / HEAD_RATIO, color);

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	//cout << "in " << duration.count() << " seconds" << endl;
}
