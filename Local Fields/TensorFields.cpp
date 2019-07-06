#include "TensorFields.h"
#include "EigenSolver.h"

#include <set>
#include <queue>

#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>
#include <igl/edges.h>
#include <igl/readOBJ.h>
#include <igl/readOFF.h>

#include <Eigen/PardisoSupport>

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

void TensorFields::getVF(Eigen::MatrixXd &V, Eigen::MatrixXi &F)
{
	V = this->V;
	F = this->F;
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

void TensorFields::computeFaceNormal()
{
	igl::per_face_normals(V, F, NF);
}

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

// For every vertex V, find where it belongs in edge E
void TensorFields::constructEVList()
{
	VENeighbors.resize(V.rows());

	for (int i = 0; i < E.rows(); i++) {
		VENeighbors[E(i, 0)].emplace(i);
		VENeighbors[E(i, 1)].emplace(i);

		//if (i < 100)
		//if (E(i, 0) == 0 || E(i, 1) == 0)
		//	printf("Edge %d has <%d, %d> vertices \n", i, E(i, 0), E(i, 1));
	}
}

void TensorFields::constructEFList()
{
	cout << "Constructing F-E lists\n";
	FE.resize(F.rows(), 3);
	EF.resize(E.rows(), 2);

	vector<set<int>> EFlist;
	EFlist.resize(E.rows());

	for (int ijk = 0; ijk < F.rows(); ijk++) {
		//printf("Test of F=%d: ", ijk);
		for (int i = 0; i < F.cols(); i++) {
			int ii = F(ijk, i);
			int in = F(ijk, (i + 1) % (int)F.cols());
			int ip;
			if (i < 1)
				ip = 2;
			else
				ip = (i - 1) % (int)F.cols();			

			/* Loop over all edges having element F(i,j) */
			for (set<int>::iterator ij = VENeighbors[ii].begin(); ij != VENeighbors[ii].end(); ++ij) {
				int edge = *ij;
				if ((ii == E(edge, 0) && in == E(edge, 1)) || (ii == E(edge, 1) && in == E(edge, 0)))
				{
					FE(ijk, ip) = edge;
					EFlist[edge].emplace(ijk);
				}
			}
		}
	}

	for (int i = 0; i < E.rows(); i++) {
		int counter = 0;
		for (set<int>::iterator ij = EFlist[i].begin(); ij != EFlist[i].end(); ++ij)
		{
			EF(i, counter) = *ij;
			counter++;
		}
	}



	/* For test */
	//int ft = 264;
	//printf("F(%d) has edges <%d, %d, %d>\n", ft, FE(ft, 0), FE(ft, 1), FE(ft, 2));
	//printf("__F(%d) has vertices (%d, %d, %d)\n", ft, F(ft, 0), F(ft, 1), F(ft, 2));
	//printf("__Edge(%d) has (%d,%d) vertices \n", FE(ft, 0), E(FE(ft, 0), 0), E(FE(ft, 0), 1));
	//printf("__Edge(%d) has (%d,%d) vertices \n", FE(ft, 1), E(FE(ft, 1), 0), E(FE(ft, 1), 1));
	//printf("__Edge(%d) has (%d,%d) vertices \n", FE(ft, 2), E(FE(ft, 2), 0), E(FE(ft, 2), 1));
	//
	//printf("Edge (%d) belongs to face <%d and %d>\n", FE(ft, 0), EF(FE(ft, 0), 0), EF(FE(ft, 0), 1));

}
/* ====================== UTILITY FUNCTIONS ============================*/
void TensorFields::constructMassMatrixMF3D()
{
	MF.resize(3 * F.rows(), 3 * F.rows());
	MFinv.resize(3 * F.rows(), 3 * F.rows());
	vector<Eigen::Triplet<double>> MFTriplet;
	MFTriplet.reserve(3 * F.rows());
	vector<Eigen::Triplet<double>> MFInvTriplet;
	MFInvTriplet.reserve(2 * F.rows());

	for (int i = 0; i < F.rows(); i++) {
		double area = doubleArea(i) / 2.0;
		MFTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * i + 0, area));
		MFTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * i + 1, area));
		MFTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * i + 2, area));
		MFInvTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * i + 0, 1.0/area));
		MFInvTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * i + 1, 1.0/area));
		MFInvTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * i + 2, 1.0/area));
	}
	MF.setFromTriplets(MFTriplet.begin(), MFTriplet.end());
	MFinv.setFromTriplets(MFInvTriplet.begin(), MFInvTriplet.end());

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

void TensorFields::selectFaceToDraw(const int& numFaces)
{
	/*Getting faces to draw, using farthest point sampling (could take some time, but still faster than drawing everything for huge mesh) */
	cout << "Selecting " << numFaces << " that will be drawn\n";

	if (numFaces < F.rows())
	{
		FaceToDraw.resize(numFaces);
		Eigen::VectorXd D(F.rows());

		/* Initialize the value of D */
		for (int i = 0; i < F.rows(); i++) {
			D(i) = numeric_limits<double>::infinity();
		}

		//srand(time(NULL));
		//FaceToDraw[0] = rand() % F.rows();
		FaceToDraw[0] = 0;

		for (int i = 1; i < numFaces; i++) {
			Eigen::VectorXi::Index maxIndex;
			computeDijkstraDistanceFaceForSampling(FaceToDraw[i - 1], D);
			D.maxCoeff(&maxIndex);
			FaceToDraw[i] = maxIndex;
		}
	}
	else
	{
		FaceToDraw.resize(F.rows());
		for (int i = 0; i < F.rows(); i++)
		{
			FaceToDraw[i] = i;
		}
	}
}

void TensorFields::computeDijkstraDistanceFaceForSampling(const int &source, Eigen::VectorXd &D)
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

void TensorFields::computeFrameRotation(igl::opengl::glfw::Viewer &viewer)
{
	cout << "Computing the angles between pairs of triangles \n";
	FrameRot.resize(E.rows(), 2);
	Eigen::Vector3d e_ij, e_ji;
	Eigen::Vector3d e_i, e_j;		// e_i: 1st edge vector of the 1st triangle     |  e_j: 1st edge vector of the 2nd triangle

	for (int ei = 0; ei < E.rows(); ei++)
	{
		/* Obtain two neighboring triangles TA and TB */
		int TA = EF(ei, 0);
		int TB = EF(ei, 1);

		/* Frame basis of each TA and TB */
		Eigen::VectorXd basisTA = V.row(F(TA, 1)) - V.row(F(TA, 0));
		Eigen::VectorXd basisTB = V.row(F(TB, 1)) - V.row(F(TB, 0));
		e_i = basisTA;
		e_j = basisTB;

		/* Finding the common edge direction + index of its first vertex */
		int eMatchA, eMatchB;
		for (int j = 0; j < 3; j++)
		{
			if (E(ei, 0) == F(TA, j))
			{				
				if (E(ei, 1) == F(TA, (j + 1) % 3))
				{
					e_ij = V.row(E(ei, 1)) - V.row(E(ei, 0));
					e_ji = V.row(E(ei, 0)) - V.row(E(ei, 1));
					eMatchA = j;
				}
				else
				{
					e_ij = V.row(E(ei, 0)) - V.row(E(ei, 1));
					e_ji = V.row(E(ei, 1)) - V.row(E(ei, 0));
					eMatchA = (3 + j - 1) % 3; 
				}
			}
		}

		/* for the 2nd triangle */
		for (int j = 0; j < 3; j++)
		{
			if (E(ei, 1) == F(TB, j))
			{
				if (E(ei, 0) == F(TB, (j + 1) % 3))
				{
					eMatchB = j;
				} 
				else
				{
					eMatchB = (3 + j - 1) % 3; 
				}
			} 
		}

		/* Computing angle for triangle A (the first one) */
		double dp_, angle_;
		switch (eMatchA)
		{
		case 0:
			FrameRot(ei, 0) = 0.0;
			break;
		case 1:
			dp_ = e_i.dot(e_ij)/ (e_i.norm() * e_ij.norm());
			angle_ = acos(dp_);
			FrameRot(ei, 0) = angle_;
			break;
		case 2:
			dp_ = e_ij.dot(e_i) / (e_i.norm() * e_ij.norm());
			angle_ = M_PI + acos(dp_);
			FrameRot(ei, 0) = angle_;
			break;
		default:
			break;
		}

		/* Computing angle for triangle B (the second one) */
		switch (eMatchB)
		{
		case 0:
			FrameRot(ei, 1) = 0.0;
			break;
		case 1:
			dp_ = e_j.dot(e_ij) / (e_j.norm() * e_ij.norm());
			angle_ = acos(dp_);
			FrameRot(ei, 1) = angle_;
			break;
		case 2:
			dp_ = e_ij.dot(e_j) / (e_j.norm() * e_ij.norm());
			angle_ = M_PI + acos(dp_);
			FrameRot(ei, 1) = angle_;
			break;
		default:
			break;
		}


		/** _____________________ DEBUG PURPOSE _____________________________*/
		//if (ei < 100 && ei%10==0)
		if(ei==1110)
		{
			// first basis of the triangle frame
			viewer.data().add_edges(V.row(F(TA, 0)), V.row(F(TA, 0)) + e_i.transpose(), Eigen::RowVector3d(0.9, 0.0, 0.0));
			// shared edge
			viewer.data().add_edges(V.row(F(TA, eMatchA)), V.row(F(TA, eMatchA)) + e_ij.transpose(), Eigen::RowVector3d(0.0, 0.0, 0.9));
			viewer.data().add_points(V.row(F(TA, eMatchA)), Eigen::RowVector3d(0.0, 0.1, 0.1));
			printf("Angle between them is %.5f degree \n", FrameRot(ei, 0)*180.0 / M_PI);
		}
	}
}

void TensorFields::obtainTransformationForLaplacian(double cT, double sT, double cF, double sF, Eigen::Matrix3d& G) 
{
	double t1 = cT;				// First rotation matrix T: map the first basis to basis in common edge
	double t2 = -sT;
	double t3 = sT;
	double t4 = cT;

	double f1 = cF;				// Second rotation matrix F: map the second basis to basis in common edge
	double f2 = -sF;
	double f3 = sF; 
	double f4 = cF; 

	double g1 = t1*f1 + t2*f2;	// T * F' (F.transpose())
	double g2 = t1*f3 + t2*f4;
	double g3 = t3*f1 + t4*f2;
	double g4 = t3*f3 + t4*f4;

	G(0, 0) = g1*g1;			// Transformation matrix in voigt notation
	G(0, 1) = g2*g2;
	G(0, 2) = 2 * g1*g2;
	G(1, 0) = g3*g3;
	G(1, 1) = g4*g4;
	G(1, 2) = 2 * g3*g4;
	G(2, 0) = g1*g3;
	G(2, 1) = g2*g4;
	G(2, 2) = g1*g4 + g2*g3;
}


void TensorFields::buildStiffnessMatrix()
{
	SF.resize(3 * F.rows(), 3 * F.rows());
	vector<Eigen::Triplet<double>> STriplet;
	STriplet.reserve(4 * 9 * F.rows());

	for (int ei = 0; ei < E.rows(); ei++)
	{
		/* Obtain two neighboring triangles TA and TB */
		int TA = EF(ei, 0);
		int TB = EF(ei, 1);

		/* Construct the rotation matrix RA and SB */
		double cosRA = cos(FrameRot(ei, 0));
		double sinRA = cos(FrameRot(ei, 0));
		double cosSB = cos(FrameRot(ei, 1));
		double sinSB = cos(FrameRot(ei, 1));


		/* Transformation T of entries of matrix B to basis of matrix A) => R*-Id*ST*B*S*-Id*RT 
		** having M = R*-Id*ST
		** then T = M*B*MT */
		

		Eigen::Matrix3d B1toB2, B2toB1;
		/* B1toB2 => parallel transport matrix A to B *
		** B2toB1 => parallel transport matrix B to A */
		obtainTransformationForLaplacian(cosRA, sinRA, cosSB, sinSB, B2toB1);
		obtainTransformationForLaplacian(cosSB, sinSB, cosRA, sinRA, B1toB2);

		/* (Combinatorial) Laplace matrix from A (first triangle) perspective */
		for (int i = 0; i < 3; i++) 
		{
			STriplet.push_back(Eigen::Triplet<double>(3 * TA + i, 3 * TA + i, 1.0 / 3.0));
			for (int j = 0; j < 3; j++) 
			{
				STriplet.push_back(Eigen::Triplet<double>(3 * TA + i, 3 * TB + j, -B2toB1(i, j) / 3.0));
			}
		}

		/* (Combinatorial) Laplace matrix from B (second triangle) perspective */
		for (int i = 0; i < 3; i++)
		{
			STriplet.push_back(Eigen::Triplet<double>(3 * TB + i, 3 * TB + i, 1.0 / 3.0));
			for (int j = 0; j < 3; j++)
			{
				STriplet.push_back(Eigen::Triplet<double>(3 * TB + i, 3 * TA + j, -B1toB2(i, j) / 3.0));
			}
		}
	}

	/* Populate the matrix with configured triplets */
	SF.setFromTriplets(STriplet.begin(), STriplet.end());
}


//void TensorFields::buildStiffnessMatrix_oldOne()
//{
//	SF.resize(3 * F.rows(), 3 * F.rows());
//	vector<Eigen::Triplet<double>> STriplet;
//	STriplet.reserve(4 * 9 * F.rows());
//
//	for (int ei = 0; ei < E.rows(); ei++)
//	{
//		/* Obtain two neighboring triangles TA and TB */
//		int TA = EF(ei, 0);
//		int TB = EF(ei, 1);
//
//		/* Construct the rotation matrix RA and SB */
//		double cosRA = cos(FrameRot(ei, 0));
//		double sinRA = cos(FrameRot(ei, 0));
//		double cosSB = cos(FrameRot(ei, 1));
//		double sinSB = cos(FrameRot(ei, 1));
//
//		/* Rbar = R*A*RT | Rtil = RT*A*R | Sbar=S*A*ST | Stil=ST*A*S */
//		Eigen::Matrix3d Rbar, Rtil, Sbar, Stil;
//
//		Rbar << cosRA*cosRA, -sinRA*-sinRA, 2 * cosRA*-sinRA,
//				sinRA*sinRA,  cosRA*cosRA,  2 * sinRA*cosRA,
//				cosRA*sinRA, -sinRA*cosRA,  cosRA*cosRA - sinRA*sinRA;
//		Rtil << cosRA*cosRA,  sinRA*sinRA, 2 * cosRA*sinRA,
//			   -sinRA*-sinRA, cosRA*cosRA, 2 * -sinRA*cosRA,
//				cosRA*-sinRA, sinRA*cosRA, cosRA*cosRA - sinRA*sinRA;
//
//		Sbar << cosSB*cosSB, -sinSB*-sinSB, 2 * cosSB*-sinSB,
//			    sinSB*sinSB,  cosSB*cosSB,  2 * sinSB*cosSB,
//			    cosSB*sinSB, -sinSB*cosSB,  cosSB*cosSB - sinSB*sinSB;
//		Stil << cosSB*cosSB,  sinSB*sinSB, 2 * cosSB*sinSB,
//			   -sinSB*-sinSB, cosSB*cosSB, 2 * -sinSB*cosSB,
//				cosSB*-sinSB, sinSB*cosSB, cosSB*cosSB - sinSB*sinSB;
//
//		/* Stiffness matrix from A (first triangle) perspective */
//		Eigen::Matrix3d RtRb = Rtil*Rbar;
//		Eigen::Matrix3d RtSb = -Rtil*Sbar;
//		for (int i = 0; i < 3; i++) {
//			for (int j = 0; j < 3; j++)
//			{
//				STriplet.push_back(Eigen::Triplet<double>(3 * TA + i, 3 * TA + j, RtRb(i, j)));
//				STriplet.push_back(Eigen::Triplet<double>(3 * TA + i, 3 * TB + j, RtSb(i, j)));
//			}
//		}
//
//		/* Stiffness matrix from B (second triangle) perspective */
//		Eigen::Matrix3d StSb = Stil*Sbar;
//		Eigen::Matrix3d StRb = -Stil*Rbar;
//		for (int i = 0; i < 3; i++) {
//			for (int j = 0; j < 3; j++)
//			{
//				STriplet.push_back(Eigen::Triplet<double>(3 * TB + i, 3 * TB + j, RtRb(i, j)));
//				STriplet.push_back(Eigen::Triplet<double>(3 * TB + i, 3 * TA + j, RtSb(i, j)));
//			}												   
//		}
//	}
//
//	/* Populate the matrix with configured triplets */
//	SF.setFromTriplets(STriplet.begin(), STriplet.end());
//}

/* Converting tensor fields (each of 2x2 size) to voigt's notation vector fields (3x1) 
** __input : tensor fields
** __output: voigt reprsentation vector fields */
void TensorFields::convertTensorToVoigt(const Eigen::MatrixXd& tensor, Eigen::VectorXd& voigt)
{
	voigt.resize(3 * F.rows());

	for (int i = 0; i < F.rows(); i++)
	{
		voigt(3 * i + 0) = tensor(2 * i + 0, 0);
		voigt(3 * i + 1) = tensor(2 * i + 1, 1);
		voigt(3 * i + 2) = tensor(2 * i + 1, 0);
	}
}

void TensorFields::convertTensorToVoigt_Elementary(const Eigen::Matrix2d& tensor, Eigen::Vector3d& voigt)
{	
	voigt(0) = tensor(0, 0);
	voigt(1) = tensor(1, 1);
	voigt(2) = tensor(1, 0);
}

/* Converting voigt's notation vector fields (3x1) back to tensor fields (each of 2x2 size)
** __input : voigt reprsentation vector fields
** __output: tensor fields */
void TensorFields::convertVoigtToTensor(const Eigen::VectorXd& voigt, Eigen::MatrixXd& tensor)
{
	tensor.resize(2 * F.rows(), 2);
	for (int i = 0; i < F.rows(); i++)
	{
		tensor(2 * i + 0, 0) = voigt(3 * i + 0);
		tensor(2 * i + 1, 0) = voigt(3 * i + 2);
		tensor(2 * i + 0, 1) = voigt(3 * i + 2);
		tensor(2 * i + 1, 1) = voigt(3 * i + 1);
	}
}

void TensorFields::convertVoigtToTensor_Elementary(const Eigen::Vector3d& voigt, Eigen::Matrix2d& tensor)
{
	tensor(0, 0) = voigt(0);
	tensor(1, 0) = voigt(2);
	tensor(0, 1) = voigt(2);
	tensor(1, 1) = voigt(1);
}

void TensorFields::constructTensorRepFields(const Eigen::MatrixXd& tensor, Eigen::MatrixXd& matrixRep)
{
	matrixRep.resize(2 * F.rows(), 2);
	for (int i = 0; i < F.rows(); i++)
	{
		Eigen::Matrix2d evect_;
		Eigen::Vector2d eval_;
		Eigen::Matrix2d block_ = tensor.block(2 * i, 0, 2, 2);
		computeEigenExplicit(block_, eval_, evect_);
		matrixRep.block(2 * i, 0, 2, 1) = eval_(0)*(evect_.col(0).normalized());
		matrixRep.block(2 * i, 1, 2, 1) = eval_(1)*(evect_.col(1).normalized());
	}
}

/* ====================== MAIN METHODS OF THE CLASS ======================*/
void TensorFields::computeTensorFields()
{
	cout << "Computing 2 vector fields to represent tensor (eigen problem). \n"; 
	constructTensorRepFields(Tensor, tensorFields);
	//tensorFields.resize(2 * F.rows(), 2);
	//
	//cout << Tensor.block(0, 0, 20, 2);
	//
	//for (int i = 0; i < F.rows(); i++)
	//{
	//	Eigen::Matrix2d evect_;
	//	Eigen::Vector2d eval_;
	//	Eigen::Matrix2d block_ = Tensor.block(2 * i, 0, 2, 2);
	//	computeEigenExplicit(block_, eval_, evect_);
	//	tensorFields.block(2 * i, 0, 2, 1) = eval_(0)*(evect_.col(0).normalized());
	//	tensorFields.block(2 * i, 1, 2, 1) = eval_(1)*(evect_.col(1).normalized());
	//}
}

void TensorFields::constructCurvatureTensor(igl::opengl::glfw::Viewer &viewer)
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
	//const double scale = avgEdgeLength;

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

		/* Computing the curvature tensor on each face 
		** 1/(8*area^2) * sum of (IIj+IIk-IIi) * outerproduct(ti,ti) */
		mT = (m1 + m2 + m3) / (2.0*doubleArea(i)*doubleArea(i));

		/* Inserting the 2x2 matrix*/
		Eigen::MatrixXd ALoc = A.block(3 * i, 2 * i, 3, 2);
		mT2D = ALoc.transpose() * mT * ALoc;
		Tensor.block(2 * i, 0, 2, 2) = mT2D;

		/* For testing purpose only*/
		//if (i <= 20)
		//{
		//	// Showing edges
		//	viewer.data().add_edges(V.row(F(i, 0)), V.row(F(i, 0)) + e1.transpose(), Eigen::RowVector3d(0.9, 0.0, 0.0));
		//	viewer.data().add_edges(V.row(F(i, 1)), V.row(F(i, 1)) + e2.transpose(), Eigen::RowVector3d(0.0, 0.7, 0.0));
		//	viewer.data().add_edges(V.row(F(i, 2)), V.row(F(i, 2)) + e3.transpose(), Eigen::RowVector3d(0.0, 0.0, 1.0));
		//
		//	// Showing rotated edge => t
		//	viewer.data().add_edges(V.row(F(i, 0)) + e1.transpose() / 2.0, V.row(F(i, 0)) + e1.transpose() / 2.0 + t1.transpose(), Eigen::RowVector3d(0.9, 0.0, 0.0));
		//	viewer.data().add_edges(V.row(F(i, 1)) + e2.transpose() / 2.0, V.row(F(i, 1)) + e2.transpose() / 2.0 + t2.transpose(), Eigen::RowVector3d(0.0, 0.7, 0.0));
		//	viewer.data().add_edges(V.row(F(i, 2)) + e3.transpose() / 2.0, V.row(F(i, 2)) + e3.transpose() / 2.0 + t3.transpose(), Eigen::RowVector3d(0.0, 0.0, 1.0));
		//	cout << i << " t1= " << t1.transpose() << " | n1" << nT.transpose() << " | e1" << e1.transpose() << endl;
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
}

void TensorFields::constructVoigtVector()
{
	cout << "Representing tensor in voigt notation\n";
	convertTensorToVoigt(Tensor, voigtReps);
	//voigtReps.resize(3 * F.rows());
	//
	//for (int i = 0; i < F.rows(); i++)
	//{
	//	voigtReps(3*i + 0) = Tensor(2*i+0, 0);
	//	voigtReps(3*i + 1) = Tensor(2*i+1, 1);
	//	voigtReps(3*i + 2) = Tensor(2*i+1, 0);
	//}

}


/* ==================== VISUALIZATION ================== */
void TensorFields::visualize2Dfields(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color, const double& scale, const bool& normalized)
{
	/* For Timing*/
	chrono::high_resolution_clock::time_point	t1, t2, te1, te2, ta1, ta2;
	chrono::duration<double>					duration, da, de;
	t1 = chrono::high_resolution_clock::now();
	//cout << "> Adding edges... ";

	//=======
	/* Some constants for arrow drawing */
	const double HEAD_RATIO = 3.0;
	const double EDGE_RATIO = scale;
	double lengthScale = EDGE_RATIO*avgEdgeLength;
	//>>>>>>> master

	/* Computing the rotation angle for 1:3 ratio of arrow head */
	double rotAngle = M_PI - atan(1.0 / 3.0);
	Eigen::Matrix2d rotMat1, rotMat2;
	Eigen::SparseMatrix<double> MRot1(2 * FaceToDraw.size(), 2 * FaceToDraw.size()), MRot2(2 * FaceToDraw.size(), 2 * FaceToDraw.size());
	vector<Eigen::Triplet<double>> R1Triplet, R2Triplet;
	R1Triplet.reserve(2 * 2 * FaceToDraw.size());
	R2Triplet.reserve(2 * 2 * FaceToDraw.size());
	Eigen::MatrixXd FCLoc(FaceToDraw.size(), 3);

	/* Defining the rotation matrix (2-by-2) on the local frame */
	rotMat1 << cos(rotAngle), -sin(rotAngle), sin(rotAngle), cos(rotAngle);
	rotMat2 << cos(-rotAngle), -sin(-rotAngle), sin(-rotAngle), cos(-rotAngle);

	for (int i = 0; i < FaceToDraw.size(); i++)
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
		FCLoc.row(i) = FC.row(FaceToDraw[i]);
	}
	MRot1.setFromTriplets(R1Triplet.begin(), R1Triplet.end());
	MRot2.setFromTriplets(R2Triplet.begin(), R2Triplet.end());

	/* Getting the local data from the population of data */
	Eigen::SparseMatrix<double> ALoc(3 * FaceToDraw.size(), 2 * FaceToDraw.size());
	vector<Eigen::Triplet<double>> ATriplet;
	ATriplet.reserve(6 * FaceToDraw.size());
	Eigen::VectorXd fieldLoc(2 * FaceToDraw.size()), fields3D(3 * FaceToDraw.size()), rot1Field, rot2Field;
	Eigen::MatrixXd TFields(FaceToDraw.size(), F.cols()), Head1Fields(FaceToDraw.size(), F.cols()), Head2Fields(FaceToDraw.size(), F.cols());

	for (int i = 0; i < FaceToDraw.size(); i++)
	{
		/* Getting the selected ALoc from A */
		ATriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 2 * i + 0, A.coeff(3 * FaceToDraw[i] + 0, 2 * FaceToDraw[i] + 0)));
		ATriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 2 * i + 0, A.coeff(3 * FaceToDraw[i] + 1, 2 * FaceToDraw[i] + 0)));
		ATriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 2 * i + 0, A.coeff(3 * FaceToDraw[i] + 2, 2 * FaceToDraw[i] + 0)));
		ATriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 2 * i + 1, A.coeff(3 * FaceToDraw[i] + 0, 2 * FaceToDraw[i] + 1)));
		ATriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 2 * i + 1, A.coeff(3 * FaceToDraw[i] + 1, 2 * FaceToDraw[i] + 1)));
		ATriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 2 * i + 1, A.coeff(3 * FaceToDraw[i] + 2, 2 * FaceToDraw[i] + 1)));

		/* Getting the selected face */
		fieldLoc.block(2 * i, 0, 2, 1) = field2D.block(2 * FaceToDraw[i], 0, 2, 1);
	}
	ALoc.setFromTriplets(ATriplet.begin(), ATriplet.end());
	fields3D = ALoc * fieldLoc;

	/* The head of the arrows */
	rot1Field = MRot1*fieldLoc;
	rot1Field = ALoc * rot1Field;
	rot2Field = MRot2*fieldLoc;
	rot2Field = ALoc * rot2Field;

	/* Transform field to Matrix format */
	for (int i = 0; i < FaceToDraw.size(); i++)
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

void TensorFields::visualizeTensorFields(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::RowVector3d color1(0.0, 0.0, 1.0);
	Eigen::RowVector3d color2(1.0, 0.0, 0.0);

	//double scale = 0.01;
	visualize2Dfields(viewer,  tensorFields.col(0), color1, scale);
	visualize2Dfields(viewer, -tensorFields.col(0), color1, scale);
	visualize2Dfields(viewer,  tensorFields.col(1), color2, scale);
	visualize2Dfields(viewer, -tensorFields.col(1), color2, scale);
}

void TensorFields::visualizeSmoothedTensorFields(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::RowVector3d color1(0.1, 0.5, 0.8);
	Eigen::RowVector3d color2(0.8, 0.5, 0.1);
	Eigen::MatrixXd smoothedFields;
	constructTensorRepFields(smoothedTensor, smoothedFields);

	//double scale = 0.01;
	visualize2Dfields(viewer,  smoothedFields.col(0), color1, scale);
	visualize2Dfields(viewer, -smoothedFields.col(0), color1, scale);
	visualize2Dfields(viewer,  smoothedFields.col(1), color2, scale);
	visualize2Dfields(viewer, -smoothedFields.col(1), color2, scale);
}


/* TESTING STUFF*/
void TensorFields::TEST_TENSOR(igl::opengl::glfw::Viewer &viewer, const string& meshFile)
{
	/* Read + construct utilities */
	readMesh(meshFile);
	scaleMesh();
	computeEdges();
	computeAverageEdgeLength();
	computeFaceCenter();
	computeFaceNormal();
	constructMappingMatrix();
	constructMassMatrixMF3D();
	constructFaceAdjacency3NMatrix();
	constructEVList();
	constructEFList();
	selectFaceToDraw(5000);	

	/* Construct necessary elements for tensor analysis */
	constructCurvatureTensor(viewer);
	computeTensorFields();
	constructVoigtVector();
	//visualizeTensorFields(viewer);

	computeFrameRotation(viewer);
	///testTransformation(viewer);
	buildStiffnessMatrix();
		
	/* Testing the result */
	testDirichletAndLaplace();
	testSmoothing(viewer);
}

void TensorFields::testDirichletAndLaplace()
{
	double dir = voigtReps.transpose()*SF*voigtReps;
	printf("The energy is %.10f\n", dir);
	double l2 = voigtReps.transpose()*MF*voigtReps;
	printf("The energy is %.10f\n", MF);
}

void TensorFields::testSmoothing(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::VectorXd id(MF.rows());
	id.setConstant(1.0);
	double factor1 = id.transpose()*MF*id;
	double factor2 = id.transpose()*SF*id;
	double lambda = 0.25;

	Eigen::VectorXd lap = (MFinv*SF)*voigtReps;

	cout << "Set and solve for smoothing \n";
	//Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> sparseSolver;
	//sparseSolver.compute(MF - factor1 / factor2*lambda*SF);
	//sparseSolver.compute(MF - lambda*SF);
	
	Eigen::MatrixXd smoothedFields;
	//Eigen::VectorXd smoothedVoigt = sparseSolver.solve(MF*voigtReps);
	Eigen::VectorXd smoothedVoigt = voigtReps - lambda*SF*voigtReps;
	convertVoigtToTensor(smoothedVoigt, smoothedTensor);
	
	visualizeSmoothedTensorFields(viewer);

	double dir = smoothedVoigt.transpose()*SF*smoothedVoigt;
	printf("The energy is %.10f\n", dir);

}

void TensorFields::testTransformation(igl::opengl::glfw::Viewer &viewer)
{
	cout << "Test on transformation on ";
	srand(time(NULL));
	const int ei = rand() % E.rows();
	cout << " edge " << ei << endl; 

	/* Obtain two neighboring triangles TA and TB */
	int TA = EF(ei, 0);
	int TB = EF(ei, 1);

	/* Drawing the edges: first edge of TA, common edge, first edge of TB */
	viewer.data().add_edges(V.row(F(TA, 0)), V.row(F(TA, 1)), Eigen::RowVector3d(0.8, 0.1, 0.1));
	viewer.data().add_edges(V.row(E(ei,0)),  V.row(E(ei,1)),  Eigen::RowVector3d(0.1, 0.1, 0.1));
	viewer.data().add_edges(V.row(F(TB, 0)), V.row(F(TB, 1)), Eigen::RowVector3d(0.1, 0.1, 0.8));
	cout << "___Angle between TA and the shared edge is: " << FrameRot(ei, 0)*180.0/M_PI << endl;
	cout << "___Angle between TB and the shared edge is: " << FrameRot(ei, 1)*180.0/M_PI << endl;

	/* Construct the rotation matrix RA and SB */
	double cosRA = cos(FrameRot(ei, 0));
	double sinRA = sin(FrameRot(ei, 0));
	double cosSB = cos(FrameRot(ei, 1));
	double sinSB = sin(FrameRot(ei, 1));

	/* Matrix representation */
	Eigen::Matrix2d RA; RA << cosRA, -sinRA, sinRA, cosRA;
	Eigen::Matrix2d SB; SB << cosSB, -sinSB, sinSB, cosSB;

	cout << "RA\n: " << RA << endl << endl << "SB:\n" << SB << endl << endl; 

	cout << "___Transforming back and forth\n";
	/* Transformation T of entries of matrix B to basis of matrix A) => R*-Id*ST*B*S*-Id*RT
	** having M = R*-Id*ST
	** then T = M*B*MT */
	Eigen::Matrix3d B1toB2, B2toB1;
	/* B1toB2 => parallel transport matrix A to B *
	** B2toB1 => parallel transport matrix B to A */
	obtainTransformationForLaplacian(cosRA, sinRA, cosSB, sinSB, B2toB1);
	obtainTransformationForLaplacian(cosSB, sinSB, cosRA, sinRA, B1toB2);

	cout << "___Check the transformation in voigt: \n";
	cout << "_____ B2 to B1: \n" << B2toB1 << endl;
	cout << "_____ (B2 to B1)-T \n" << B2toB1.inverse().transpose() << endl; 
	cout << "_____ B1 to B2: \n" << B1toB2 << endl << endl;

	cout << "___Testing on a random matrix \n";
	double aa = rand() % 10 - 5;
	double ab = rand() % 10 - 5;
	double ac = rand() % 10 - 5;
	Eigen::Matrix2d A1; A1 << aa, ab, ab, ac;
	Eigen::Matrix2d A1back;
	Eigen::Vector3d a1, a1inB2, a1inB1fromB2;
	cout << "__A1 was: \n" << A1 << endl << endl;
	cout << "__convert to voigt, becomes: \n";
	convertTensorToVoigt_Elementary(A1, a1);
	cout << a1 << endl;
	cout << "__transfer a1 to B2\n";
	a1inB2 = B1toB2*a1;
	cout << a1inB2 << endl; 
	cout << "____Operation in matrix: S * R' * A * R * S'\n";
	Eigen::Matrix2d A1inB2 = SB * RA.transpose() * A1 * RA * SB.transpose();
	cout << A1inB2 << endl; 
	cout << "__transfer a1 back to B1: \n";
	a1inB1fromB2 = B2toB1*a1inB2;
	cout << a1inB1fromB2 << endl;
	cout << "____Operation in matrix: R * S' * B * S * R'\n";
	cout << RA * SB.transpose() * A1inB2 * SB * RA.transpose() << endl;
	cout << "__convert a1 back to a matrix \n";
	convertVoigtToTensor_Elementary(a1inB1fromB2, A1back);
	cout << "A1 after going to B2 and back to B1 is now: \n" << A1back << endl << endl;
}