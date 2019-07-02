#include "TensorFields.h"

#include <set>

#include <igl/per_vertex_normals.h>
#include <igl/edges.h>
#include <igl/readOBJ.h>
#include <igl/readOFF.h>

TensorFields::TensorFields()
{

}

TensorFields::~TensorFields()
{

}

/* ====================== MESH-RELATED FUNCTIONS ============================*/
void TensorFields::readMesh(const string &meshFile)
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
		cout << "Program will exit." << endl;		
		exit(10);
	}

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;

	// Printing Mesh-related information
	printf("....V=%dx%d\n", V.rows(), V.cols());
	printf("....F=%dx%d\n", F.rows(), F.cols());

	igl::edges(F, E);
	printf("....E=%dx%d\n", E.rows(), E.cols());

	igl::doublearea(V, F, doubleArea);

	const int genus = (2 - V.rows() + E.rows() - F.rows()) / 2;
	printf("This model is of genus %d\n", genus);
}

void TensorFields::scaleMesh()
{
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
	//cout << "BEFORE Scaling\n";
	//cout << "Max coeff \n" << maxV << endl;
	//cout << "Min coeff \n" << minV << endl;
	//cout << "Scale factor: " << scaleFactor << endl; 	


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

	//cout << "AFTER Scaling\n";
	//cout << "Max coeff \n" << maxV << endl;
	//cout << "Min coeff \n" << minV << endl;
	//cout << "Scale factor: " << scaleFactor << endl;
}

void TensorFields::computeFaceCenter()
{
	FC.resize(F.rows(), 3);

	for (int i = 0; i < F.rows(); i++) {
		FC.row(i) = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;
	}

}

void TensorFields::computeEdges()
{
	igl::edges(F, E);
	printf("....E=%dx%d\n", E.rows(), E.cols());
}


/* ====================== UTILITY FUNCTIONS ============================*/
void TensorFields::computeAverageEdgeLength()
{
	Eigen::Vector3d e;
	double totalLength = 0.0;
	for (int i = 0; i < F.rows(); i++) {
		for (int j = 0; j < F.cols(); j++) {
			e = V.row(F(i, (j + 1) % F.cols())) - V.row(F(i, j));
			totalLength += e.norm();
		}
	}
	avgEdgeLength = totalLength / (double)(F.rows()*F.cols());
}

void TensorFields::constructMappingMatrix()
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

void TensorFields::constructFaceAdjacency3NMatrix()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Building \"Face-to-Face\" Adjacency (3 Neigbors)... ";

	// TWO RING
	//AdjMF2Ring.clear();
	//AdjMF2Ring.resize(F.rows());

	/* _____________ VF Adjacency _________ */
	vector<set<VtoFPair>>			VFNeighbors;
	VFNeighbors.resize(V.rows());
	int counterAdj = 0;
	for (int i = 0; i < F.rows(); i++) {
		for (int j = 0; j < F.cols(); j++) {
			VtoFPair vp1{ counterAdj++, F(i, j), i };
			for (int k = j + 1; k < F.cols(); k++) {
				VtoFPair vp2{ counterAdj++, F(i, k), i };
				if (vp2.vId > vp1.vId)
					VFNeighbors[vp1.vId].insert(vp2);
				else
					VFNeighbors[vp2.vId].insert(vp1);
			}
		}
	}

	// ONE RING
	vector<set<FacePair>>			AdjMF3N_temp;
	AdjMF3N_temp.resize(F.rows());
	vector<set<Edge_VPair>>			EdgePairsList;
	EdgePairsList.resize(F.rows());
	int counter = 0, counter1 = 0;

	/* Populate the adjacency stuff */
	for (int i = 0; i < V.rows(); i++) {
		//printf("Vertex %d has %d neighbors, which arew: \n", i, VFNeighbors[i].size());
		for (std::set<VtoFPair>::iterator it1 = VFNeighbors[i].begin(); it1 != VFNeighbors[i].end(); ++it1) {
			//printf("%d (in F=%d) : ", it1->vId, it1->fId);
			for (std::set<VtoFPair>::iterator it2 = next(it1, 1); it2 != VFNeighbors[i].end(); ++it2) {
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
	// FREE MEMORY for VFNeighbors
	VFNeighbors.clear();
	VFNeighbors.shrink_to_fit();



	// MOVING THE ADJACENCY it to matrix format
	AdjMF3N.resize(F.rows(), F.cols());
	for (int i = 0; i < F.rows(); i++) {
		int counter = 0;
		//for (std::set<int>::iterator it = AdjMF3N_temp[i].begin(); it != AdjMF3N_temp[i].end(); ++it) {
		for (std::set<FacePair>::iterator it = AdjMF3N_temp[i].begin(); it != AdjMF3N_temp[i].end(); ++it) {
			AdjMF3N(i, counter++) = it->f1;
		}
	}

	// Clearing redundant data structure, used as temporary DS only
	AdjMF3N_temp.clear();
	AdjMF3N_temp.shrink_to_fit();

	// MOVING EDGE Adjacency to MATRIX format
	Eigen::MatrixXi					EdgePairMatrix;
	EdgePairMatrix.resize(F.rows(), 3 * F.cols());
	for (int i = 0; i < F.rows(); i++) {
		int eCounter = 0;
		for (std::set<Edge_VPair>::iterator it = EdgePairsList[i].begin(); it != EdgePairsList[i].end(); ++it) {
			EdgePairMatrix(i, eCounter++) = it->v1;
			EdgePairMatrix(i, eCounter++) = it->v2;
		}
	}
	// Save memory by free-ing space occupied by EdgePairList (temp data structure)
	EdgePairsList.clear();
	EdgePairsList.shrink_to_fit();

	printf("Size of EdgePairMatrix = %dx%d (F=%dx%d)\n", EdgePairMatrix.rows(), EdgePairMatrix.cols(), F.rows(), F.cols());

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;

}

/* ====================== MAIN METHODS OF THE CLASS ======================*/
void TensorFields::constructCurvatureTensor()
{
	cout << "Constructing curvature tensor \n";
	cout << "__Computing vertex normal\n";
	/* Getting normals on each vertex */
	Eigen::MatrixXd NV;
	igl::per_vertex_normals(V, F, NV);

	/* Declare local variable for the loop here, to avoid excessive allocation (constructor) + de-allocation (destructor) */
	double f2Form1, f2Form2, f2Form3;
	Eigen::Vector3d t1, t2, t3, e1, e2, e3, n1, n2, n3, nTemp, nT;
	Eigen::Matrix3d m1, m2, m3, mT;
	Eigen::Matrix2d mT2D;

	Tensor.resize(2 * F.rows(), 2);
	//CurvatureTensor2D.resize(2 * F.rows(), 2 * F.rows());
	//CurvatureTensor2D.reserve(2 * 2 * F.rows());			// 2*F rows, each with 2 non-zero entries
	vector<Eigen::Triplet<double>> CTriplet;
	CTriplet.reserve(2 * 2 * F.rows());						// 2*F rows, each with 2 non-zero entries
	const double scale = 0.2 * avgEdgeLength;

	cout << "__Computing curvature tensor\n";
	/* Loop over all faces */
	for (int i = 0; i < F.rows(); i++)
	{
		/* Getting the edges */
		e1 = (V.row(F(i, 1)) - V.row(F(i, 0))).transpose();
		e2 = (V.row(F(i, 2)) - V.row(F(i, 1))).transpose();
		e3 = (V.row(F(i, 0)) - V.row(F(i, 2))).transpose();

		/* Getting the rotated edge (CW, 90, on triangle plane->uses triangle normal)*/
		nT = (NF.row(i)).transpose();
		nT.normalize();
		t1 = e1.cross(nT);
		t2 = e2.cross(nT);
		t3 = e3.cross(nT);

		/* Getting normals for each edge center (average of two vertices) */
		/* NOT THE MOST Efficient implementation */
		/* Get the first neighbor's normal*/
		for (int j = 0; j < F.cols(); j++)		// Loop over 3 neighbors of i-th face
		{
			int neigh = AdjMF3N(i, j);
			for (int k = 0; k < F.cols(); k++)	// Loop over 3 vertices of the j-th neighbor
			{
				if (F(i, 0) == F(neigh, (k + 1) % F.cols()) && F(i, 1) == F(neigh, k))
				{
					nTemp = (NF.row(neigh)).transpose();
					nTemp.normalize();
					n1 = nT + nTemp;
					n1.normalize();
					//printf("____ The 1st edge: %d->%d to Triangle %d (%d, %d, %d)\n", F(i, 0), F(i, 1), neigh, F(neigh, 0), F(neigh, 1), F(neigh, 2));
				}
			}
		}
		/* Get the second neighbor's normal*/
		for (int j = 0; j < F.cols(); j++)		// Loop over 3 neighbors of i-th face
		{
			int neigh = AdjMF3N(i, j);
			for (int k = 0; k < F.cols(); k++)	// Loop over 3 vertices of the j-th neighbor
			{
				if (F(i, 1) == F(neigh, (k + 1) % F.cols()) && F(i, 2) == F(neigh, k))
				{
					nTemp = (NF.row(neigh)).transpose();
					nTemp.normalize();
					n2 = nT + nTemp;
					n2.normalize();
					//printf("____ The 2nd edge: %d->%d to Triangle %d (%d, %d, %d)\n", F(i, 1), F(i, 2), neigh, F(neigh, 0), F(neigh, 1), F(neigh, 2));
				}
			}
		}
		/* Get the third neighbor neighbor's normal */
		for (int j = 0; j < F.cols(); j++)		// Loop over 3 neighbors of i-th face
		{
			int neigh = AdjMF3N(i, j);
			for (int k = 0; k < F.cols(); k++)	// Loop over 3 vertices of the j-th neighbor
			{
				if (F(i, 2) == F(neigh, (k + 1) % F.cols()) && F(i, 0) == F(neigh, k))
				{
					nTemp = (NF.row(neigh)).transpose();
					nTemp.normalize();
					n3 = nT + nTemp;
					n3.normalize();
					//printf("____ The 3rd edge: %d->%d to Triangle %d (%d, %d, %d)\n", F(i, 2), F(i, 0), neigh, F(neigh, 0), F(neigh, 1), F(neigh, 2));
				}
			}
		}

		/* Computing 2nd Fundamental form */
		f2Form1 = 2.0 * (n2 - n3).dot(e1);
		f2Form2 = 2.0 * (n3 - n1).dot(e2);
		f2Form3 = 2.0 * (n1 - n2).dot(e3);

		/* Computing the outer product */
		m1 = (f2Form2 + f2Form3 - f2Form1) * (t1 * t1.transpose());
		m2 = (f2Form3 + f2Form1 - f2Form2) * (t2 * t2.transpose());
		m3 = (f2Form1 + f2Form2 - f2Form3) * (t3 * t3.transpose());

		/* Computing the curvature tensor on each face */
		mT = (m1 + m2 + m3) / (2.0*doubleArea(i)*doubleArea(i));

		/* Inserting the 2x2 matrix*/
		Eigen::MatrixXd ALoc = A.block(3 * i, 2 * i, 3, 2);
		mT2D = ALoc.transpose() * mT * ALoc;
		Tensor.block(2 * i, 0, 2, 2) = mT2D;

		//CTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, mT2D(0, 0)));
		//CTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 0, mT2D(1, 0)));
		//CTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 1, mT2D(0, 1)));
		//CTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, mT2D(1, 1)));

		/* For testing purpose only*/
		//if (i <= 20)
		//{
		//	// Showing edges
		//	viewer.data().add_edges(V.row(F(i, 0)), V.row(F(i, 0)) + e1.transpose(), Eigen::RowVector3d(0.9, 0.0, 0.0));
		//	viewer.data().add_edges(V.row(F(i, 1)), V.row(F(i, 1)) + e2.transpose(), Eigen::RowVector3d(0.0, 0.7, 0.0));
		//	viewer.data().add_edges(V.row(F(i, 2)), V.row(F(i, 2)) + e3.transpose(), Eigen::RowVector3d(0.0, 0.0, 1.0));
		//
		//	// Showing rotated edge => t
		//	viewer.data().add_edges(V.row(F(i, 0)) + e1.transpose() / 2.0, V.row(F(i, 0)) + e1.transpose() / 2.0 + scale*t1.transpose(), Eigen::RowVector3d(0.9, 0.0, 0.0));
		//	viewer.data().add_edges(V.row(F(i, 1)) + e2.transpose() / 2.0, V.row(F(i, 1)) + e2.transpose() / 2.0 + scale*t2.transpose(), Eigen::RowVector3d(0.0, 0.7, 0.0));
		//	viewer.data().add_edges(V.row(F(i, 2)) + e3.transpose() / 2.0, V.row(F(i, 2)) + e3.transpose() / 2.0 + scale*t3.transpose(), Eigen::RowVector3d(0.0, 0.0, 1.0));
		//
		//	// Showing the normals  ni
		//	viewer.data().add_edges(V.row(F(i, 0)) + e1.transpose() / 2.0, V.row(F(i, 0)) + e1.transpose() / 2.0 + scale*n1.transpose(), Eigen::RowVector3d(0.9, 0.0, 0.0));
		//	viewer.data().add_edges(V.row(F(i, 1)) + e2.transpose() / 2.0, V.row(F(i, 1)) + e2.transpose() / 2.0 + scale*n2.transpose(), Eigen::RowVector3d(0.0, 0.7, 0.0));
		//	viewer.data().add_edges(V.row(F(i, 2)) + e3.transpose() / 2.0, V.row(F(i, 2)) + e3.transpose() / 2.0 + scale*n3.transpose(), Eigen::RowVector3d(0.0, 0.0, 1.0));
		//
		//	cout << "MT=" << i << endl << ": " << mT << endl;
		//	cout << "MT2D=" << i << endl << ": " << mT2D << endl;
		//	//cout << "M1=" << f2Form1 << endl << ": " << m1 << endl;
		//	//cout << "M2=" << f2Form2 << endl << ": " << m2 << endl;
		//	//cout << "M3=" << f2Form3 << endl << ": " << m3 << endl;
		//}
	}
	//CurvatureTensor2D.setFromTriplets(CTriplet.begin(), CTriplet.end());
}