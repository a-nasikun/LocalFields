#include "NRoSyFields.h"

#include <igl/principal_curvature.h>

#include <random>

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

	igl::doublearea(V, F, doubleArea);
}

void NRoSyFields::getVF(Eigen::MatrixXd &V, Eigen::MatrixXi &F)
{
	V = this->V;
	F = this->F;
}

void NRoSyFields::computeFaceCenter()
{
	FC.resize(F.rows(), 3);

	for (int i = 0; i < F.rows(); i++) {
		FC.row(i) = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;
	}
}

void NRoSyFields::constructFaceAdjacency3NMatrix()
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

void NRoSyFields::constructFaceAdjacency2RingMatrix()
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

// For every vertex V, find where it belongs in edge E
void NRoSyFields::constructEVList()
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

void NRoSyFields::constructEFList()
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
}


void NRoSyFields::computeAverageEdgeLength()
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

	cout << "__Avg edgelength=" << avgEdgeLength << endl;
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

void NRoSyFields::selectFaceToDraw(const int& numFaces)
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

void NRoSyFields::computeDijkstraDistanceFaceForSampling(const int &source, Eigen::VectorXd &D)
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

void NRoSyFields::constructMassMatrixMF3D()
{
	MF.resize(2 * F.rows(), 2 * F.rows());
	MFinv.resize(2 * F.rows(), 2 * F.rows());
	MF2DhNeg.resize(2 * F.rows(), 2 * F.rows());
	MF2DhPos.resize(2 * F.rows(), 2 * F.rows());

	vector<Eigen::Triplet<double>> MFTriplet; MFTriplet.reserve(2 * F.rows());
	vector<Eigen::Triplet<double>> MFInvTriplet; MFInvTriplet.reserve(2 * F.rows());
	vector<Eigen::Triplet<double>> MhPosTriplet; MhPosTriplet.reserve(2 * F.rows());
	vector<Eigen::Triplet<double>> MhNegTriplet; MhNegTriplet.reserve(2 * F.rows());
		

	for (int i = 0; i < F.rows(); i++) {
		double area = doubleArea(i) / 2.0;
		MFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, area));
		MFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, area));
		MFInvTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, 1.0 / area));
		MFInvTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, 1.0 / area));
		double sqrt_area = sqrt(area);
		MhPosTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, sqrt_area));
		MhPosTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, sqrt_area));
		MhNegTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, 1.0 / sqrt_area));
		MhNegTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, 1.0 / sqrt_area));
	}
	MF.setFromTriplets(MFTriplet.begin(), MFTriplet.end());
	MFinv.setFromTriplets(MFInvTriplet.begin(), MFInvTriplet.end());
	MF2DhPos.setFromTriplets(MhPosTriplet.begin(), MhPosTriplet.end());
	MF2DhNeg.setFromTriplets(MhNegTriplet.begin(), MhNegTriplet.end());
}

//void NRoSyFields::buildStiffnessMatrix_Combinatorial()
//{
//	cout << "Try to build the harmonic energy \n";
//	SF.resize(2 * F.rows(), 2 * F.rows());
//	vector<Eigen::Triplet<double>> STriplet;
//	STriplet.reserve(4 * 4 * F.rows());
//
//	srand(time(NULL));
//	const int ee = rand() % E.rows();
//
//	for (int ei = 0; ei < E.rows(); ei++)
//	{
//		/* Obtain two neighboring triangles TA and TB */
//		//cout << "Obtain two neighboring triangles TA and TB \n";
//		int TA = EF(ei, 0);
//		int TB = EF(ei, 1);
//
//		if (ei == ee)
//		{
//			cout << "TA: " << TA << ", and TB: " << TB << endl;
//		}
//
//		/* Construct the rotation matrix RA and SB */
//		//cout << "Construct the rotation matrix RA and SB\n";
//		double cosRA = cos(FrameRot(ei, 0)/(double)nRot);
//		double sinRA = sin(FrameRot(ei, 0)/(double)nRot);
//		double cosSB = cos(FrameRot(ei, 1)/(double)nRot);
//		double sinSB = sin(FrameRot(ei, 1)/(double)nRot);
//
//		Eigen::Matrix2d R1; R1 << cosRA, -sinRA, sinRA, cosRA;
//		Eigen::Matrix2d R2; R2 << cosSB, -sinSB, sinSB, cosSB;
//
//		Eigen::Matrix2d B2toB1 = R1*R2.transpose();
//		Eigen::Matrix2d B1toB2 = R2*R1.transpose();
//
//		//double cosRA = cos(nRot*(FrameRot(ei, 0)-FrameRot(ei, 1)));
//		//double sinRA = sin(nRot*(FrameRot(ei, 0)-FrameRot(ei, 1)));
//		//double cosSB = cos(nRot*(FrameRot(ei, 1)-FrameRot(ei, 0)));
//		//double sinSB = sin(nRot*(FrameRot(ei, 1)-FrameRot(ei, 0)));
//		//
//		//Eigen::Matrix2d R1; R1 << cosRA, -sinRA, sinRA, cosRA;
//		//Eigen::Matrix2d R2; R2 << cosSB, -sinSB, sinSB, cosSB;
//		//
//		///* The transport matrix */
//		//Eigen::Matrix2d B2toB1 = R1;
//		//Eigen::Matrix2d B1toB2 = R2;
//
//		if (ei == ee)
//		{
//			cout << "Data of : " << ei << endl;
//			cout << "B2 to B1 : \n" << B2toB1 << endl;
//			cout << "B1 to B2 : \n" << B1toB2 << endl;
//		}
//
//		/* (Geometric) Laplace matrix from A (first triangle) perspective */
//		for (int i = 0; i < 2; i++)
//		{
//			STriplet.push_back(Eigen::Triplet<double>(2 * TA + i, 2 * TA + i, 1.0/3.0));
//			for (int j = 0; j < 2; j++)
//			{
//				STriplet.push_back(Eigen::Triplet<double>(2 * TA + i, 2 * TB + j, B2toB1(i, j) / 3.0));
//			}
//		}
//
//		/* (Geometric) Laplace matrix from B (first triangle) perspective */
//		for (int i = 0; i < 2; i++)
//		{
//			STriplet.push_back(Eigen::Triplet<double>(2 * TB + i, 2 * TB + i, 1.0 / 3.0));
//			for (int j = 0; j < 2; j++)
//			{
//				STriplet.push_back(Eigen::Triplet<double>(2 * TB + i, 2 * TA + j, B1toB2(i, j) / 3.0));
//			}
//		}
//
//	}
//	SF.setFromTriplets(STriplet.begin(), STriplet.end());
//	string filename = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/Arma_Lap_NRoSy_Comb";
//	WriteSparseMatrixToMatlab(SF, filename);
//}

void NRoSyFields::buildStiffnessMatrix_Combinatorial()
{
	cout << "Try to build the harmonic energy \n";
	SF.resize(2 * F.rows(), 2 * F.rows());
	vector<Eigen::Triplet<double>> STriplet;
	STriplet.reserve(4 * 4 * F.rows());

	srand(time(NULL));
	const int ee = rand() % E.rows();

	for (int ei = 0; ei < E.rows(); ei++)
	{
		/* Obtain two neighboring triangles TA and TB */
		//cout << "Obtain two neighboring triangles TA and TB \n";
		int TA = EF(ei, 0);
		int TB = EF(ei, 1);

		if (ei == ee)
		{
			cout << "TA: " << TA << ", and TB: " << TB << endl;
		}

		/* Construct the rotation matrix RA and SB */
		//cout << "Construct the rotation matrix RA and SB\n";
		double map_angle = nRot * (FrameRot(ei, 0) - FrameRot(ei, 1) + M_PI);
		Eigen::Matrix2d Rot; Rot << cos(map_angle), -sin(map_angle), sin(map_angle), cos(map_angle);


		if (ei == ee)
		{
			cout << "Data of : " << ei << endl;
			cout << "B2 to B1 : \n" << Rot << endl;
			cout << "B1 to B2 : \n" << Rot.transpose() << endl;
		}

		/* (Geometric) Laplace matrix from A (first triangle) perspective */
		for (int i = 0; i < 2; i++)
		{
			STriplet.push_back(Eigen::Triplet<double>(2 * TA + i, 2 * TA + i, 1.0 / 3.0));
			for (int j = 0; j < 2; j++)
			{
				STriplet.push_back(Eigen::Triplet<double>(2 * TA + i, 2 * TB + j, -Rot(i, j) / 3.0));
			}
		}

		/* (Geometric) Laplace matrix from B (first triangle) perspective */
		for (int i = 0; i < 2; i++)
		{
			STriplet.push_back(Eigen::Triplet<double>(2 * TB + i, 2 * TB + i, 1.0 / 3.0));
			for (int j = 0; j < 2; j++)
			{
				STriplet.push_back(Eigen::Triplet<double>(2 * TB + i, 2 * TA + j, -Rot(j, i) / 3.0));
			}
		}

	}
	SF.setFromTriplets(STriplet.begin(), STriplet.end());
	string filename = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/Arma_Lap_NRoSy_Comb";
	//WriteSparseMatrixToMatlab(SF, filename);
}
//void NRoSyFields::buildStiffnessMatrix_Geometric()
//{
//	cout << "Try to build the harmonic energy \n";
//	SF.resize(2 * F.rows(), 2 * F.rows());
//	vector<Eigen::Triplet<double>> STriplet;
//	STriplet.reserve(4 * 4 * F.rows());
//
//	srand(time(NULL));
//	const int ee = rand() % E.rows();
//
//	for (int ei = 0; ei < E.rows(); ei++)
//	{
//		/* Obtain two neighboring triangles TA and TB */
//		//cout << "Obtain two neighboring triangles TA and TB \n";
//		int TA = EF(ei, 0);
//		int TB = EF(ei, 1);
//		double weight;
//
//		if (ei == ee)
//		{
//			cout << "TA: " << TA << ", and TB: " << TB << endl; 
//		}
//
//		/* Construct the rotation matrix RA and SB */
//		//cout << "Construct the rotation matrix RA and SB\n";
//		double cosRA = cos(FrameRot(ei, 0));
//		double sinRA = sin(FrameRot(ei, 0));
//		double cosSB = cos(FrameRot(ei, 1));
//		double sinSB = sin(FrameRot(ei, 1));
//		//double cosRA = cos(FrameRot(ei, 0) / (double)nRot);
//		//double sinRA = sin(FrameRot(ei, 0) / (double)nRot);
//		//double cosSB = cos(FrameRot(ei, 1) / (double)nRot);
//		//double sinSB = sin(FrameRot(ei, 1) / (double)nRot);
//
//		Eigen::Matrix2d R1; R1 << cosRA, -sinRA, sinRA, cosRA;
//		Eigen::Matrix2d R2; R2 << cosSB, -sinSB, sinSB, cosSB;
//
//		//cout << "Copmute the weight \n";
//		/* Compute the weight on each edge */
//		Eigen::Vector3d edge_ = V.row(E(ei, 1)) - V.row(E(ei, 0));
//		double el = edge_.dot(edge_);
//		/* Should be multiplied by 3.0
//		** but because later will be divided by 3.0 again (due to 3 neighbors),
//		** I just dont multiply with anything */
//		//weight = el * el / (0.5*doubleArea(TA) + 0.5*doubleArea(TB));
//		weight = (3.0 * el * el) / (0.5*doubleArea(TA) + 0.5*doubleArea(TB));
//
//		/* The transport matrix */
//		Eigen::Matrix2d B2toB1 = R1*R2.transpose();
//		Eigen::Matrix2d B1toB2 = R2*R1.transpose();
//
//		if (ei == ee)
//		{
//			cout << "Data of : " << ei << endl;
//			cout << "B2 to B1 : \n" << B2toB1 << endl;
//			cout << "B1 to B2 : \n" << B1toB2 << endl; 
//		}
//
//		/* (Geometric) Laplace matrix from A (first triangle) perspective */
//		for (int i = 0; i < 2; i++)
//		{
//			STriplet.push_back(Eigen::Triplet<double>(2 * TA + i, 2 * TA + i, weight));
//			for (int j = 0; j < 2; j++)
//			{
//				STriplet.push_back(Eigen::Triplet<double>(2 * TA + i, 2 * TB + j, weight*B2toB1(i, j)));
//			}
//		}
//
//		/* (Geometric) Laplace matrix from B (first triangle) perspective */
//		for (int i = 0; i < 2; i++)
//		{
//			STriplet.push_back(Eigen::Triplet<double>(2 * TB + i, 2 * TB + i, weight));
//			for (int j = 0; j < 2; j++)
//			{
//				STriplet.push_back(Eigen::Triplet<double>(2 * TB + i, 2 * TA + j, weight*B1toB2(i, j)));
//			}
//		}
//	}
//	SF.setFromTriplets(STriplet.begin(), STriplet.end());
//}

void NRoSyFields::buildStiffnessMatrix_Geometric()
{
	cout << "Try to build the harmonic energy \n";
	SF.resize(2 * F.rows(), 2 * F.rows());
	vector<Eigen::Triplet<double>> STriplet;
	STriplet.reserve(4 * 4 * F.rows());

	srand(time(NULL));
	const int ee = rand() % E.rows();

	for (int ei = 0; ei < E.rows(); ei++)
	{
		/* Obtain two neighboring triangles TA and TB */
		//cout << "Obtain two neighboring triangles TA and TB \n";
		int TA = EF(ei, 0);
		int TB = EF(ei, 1);
		double weight;

		/* Construct the rotation matrix RA and SB */
		//cout << "Construct the rotation matrix RA and SB\n";
		double map_angle = nRot * (FrameRot(ei, 0) - FrameRot(ei, 1) + M_PI);
		Eigen::Matrix2d Rot; Rot << cos(map_angle), -sin(map_angle), sin(map_angle), cos(map_angle);


		if (ei == ee)
		{
			cout << "Data of : " << ei << endl;
			cout << "B2 to B1 : \n" << Rot << endl;
			cout << "B1 to B2 : \n" << Rot.transpose() << endl;
		}

		//cout << "Copmute the weight \n";
		/* Compute the weight on each edge */
		Eigen::Vector3d edge_ = V.row(E(ei, 1)) - V.row(E(ei, 0));
		double el = edge_.dot(edge_);
		/* Should be multiplied by 3.0
		** but because later will be divided by 3.0 again (due to 3 neighbors),
		** I just dont multiply with anything */
		//weight = el * el / (0.5*doubleArea(TA) + 0.5*doubleArea(TB));
		weight = (3.0 * el * el) / (0.5*doubleArea(TA) + 0.5*doubleArea(TB));
				

		/* (Geometric) Laplace matrix from A (first triangle) perspective */
		for (int i = 0; i < 2; i++)
		{
			STriplet.push_back(Eigen::Triplet<double>(2 * TA + i, 2 * TA + i, weight*1.0));
			for (int j = 0; j < 2; j++)
			{
				STriplet.push_back(Eigen::Triplet<double>(2 * TA + i, 2 * TB + j, -weight*Rot(i, j)));
			}
		}

		/* (Geometric) Laplace matrix from B (first triangle) perspective */
		for (int i = 0; i < 2; i++)
		{
			STriplet.push_back(Eigen::Triplet<double>(2 * TB + i, 2 * TB + i, weight));
			for (int j = 0; j < 2; j++)
			{
				STriplet.push_back(Eigen::Triplet<double>(2 * TB + i, 2 * TA + j, -weight*Rot(j, i) ));
			}
		}
	}
	SF.setFromTriplets(STriplet.begin(), STriplet.end());

	string model = "Brezel1920";
	string fileLaplace = "D:/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/" + model + "_SF_NRoSy_Geom3";
	//WriteSparseMatrixToMatlab(SF, fileLaplace);
	fileLaplace = "D:/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/" + model + "_M_NRoSy_Geom2";
	//WriteSparseMatrixToMatlab(MF, fileLaplace);
}

void NRoSyFields::computeFrameRotation(igl::opengl::glfw::Viewer &viewer)
{
	cout << "Computing the angles between pairs of triangles \n";
	FrameRot.resize(E.rows(), 2);
	Eigen::Vector3d e_ij, e_ji;
	Eigen::Vector3d e_i, e_j;		// e_i: 1st edge vector of the 1st triangle     |  e_j: 1st edge vector of the 2nd triangle

	srand(time(NULL));
	int testEdge = rand() % E.rows();

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
		double dp_1, angle_1;
		switch (eMatchA)
		{
		case 0:
			angle_1 = 0.0;
			FrameRot(ei, 0) = 0.0;
			break;
		case 1:
			//dp_1 = e_i.dot(e_ij)/ (e_i.norm() * e_ij.norm());
			dp_1 = (e_ij).dot(-e_i) / (e_i.norm()*e_ij.norm());
			angle_1 = acos(dp_1);					
			FrameRot(ei, 0) = M_PI - angle_1;
			break;
		case 2:
			//dp_1 = e_ij.dot(e_i) / (e_i.norm() * e_ij.norm());
			dp_1 = (e_i).dot(-e_ij) / (e_i.norm() * e_ij.norm());
			angle_1 = acos(dp_1);
			FrameRot(ei, 0) =  M_PI + angle_1;
			break;
		default:
			break;
		}

		/* Computing angle for triangle B (the second one) */
		double dp_2, angle_2;
		switch (eMatchB)
		{
		case 0:
			angle_2 = 0.0;
			FrameRot(ei, 1) = 0.0;
			break;
		case 1:
			//dp_2 = e_j.dot(e_ij) / (e_j.norm() * e_ij.norm());
			dp_2 = (e_ji).dot(-e_j) / (e_j.norm() * e_ji.norm());
			angle_2 = acos(dp_2);
			FrameRot(ei, 1) =  M_PI - angle_2;
			break;
		case 2:
			//dp_2 = e_ij.dot(e_j) / (e_j.norm() * e_ij.norm());
			dp_2 = (e_j).dot(-e_ji) / (e_j.norm() * e_ji.norm());
			angle_2 = acos(dp_2);
			FrameRot(ei, 1) = M_PI + angle_2;
			break;
		default:
			break;
		}


		/** _____________________ DEBUG PURPOSE _____________________________*/
		//if (ei < 100 && ei%10==0)
		if (ei == testEdge)
		{
			// first basis of the triangle frame (frist and second)
			viewer.data().add_edges(V.row(F(TA, 0)), V.row(F(TA, 0)) + e_i.transpose(), Eigen::RowVector3d(0.9, 0.0, 0.0));
			viewer.data().add_edges(V.row(F(TB, 0)), V.row(F(TB, 0)) + e_j.transpose(), Eigen::RowVector3d(0.0, 0.9, 0.0));
			// shared edge
			viewer.data().add_edges(V.row(E(ei, 0)), V.row(E(ei, 1)), Eigen::RowVector3d(0.0, 0.0, 0.9));
			viewer.data().add_points(V.row(E(ei, 0)), Eigen::RowVector3d(0.0, 0.1, 0.1));
			viewer.data().add_points(V.row(E(ei, 1)), Eigen::RowVector3d(0.9, 0.1, 0.1));
			//viewer.data().add_edges(V.row(F(TA, eMatchA)), V.row(F(TA, eMatchA)) + e_ij.transpose(), Eigen::RowVector3d(0.0, 0.0, 0.9));
			//viewer.data().add_points(V.row(F(TA, eMatchA)), Eigen::RowVector3d(0.0, 0.1, 0.1));
			printf("Angle between basis and shared edge (fram 1) is %.5f degree \n", FrameRot(ei, 0)*180.0 / M_PI);
			cout << "__case 1= " << eMatchA << endl;
			printf("__angle1 = %.10f \n", angle_1*180.0 / M_PI);
			printf("Angle between basis and shared edge (fram 2) is %.5f degree \n", FrameRot(ei, 1)*180.0 / M_PI);
			cout << "__case 2= " << eMatchB << endl;
			printf("__angle2 = %.10f \n", angle_2*180.0 / M_PI);
		}
	}
}

void NRoSyFields::computeEigenFields_generalized(const int &numEigs, const string& filename)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Computing reference generalized-eigenproblem (in Matlab)... ";

	printf("Size MF=%dx%d || size SF=%dx%d \n", MF.rows(), MF.cols(), SF.rows(), SF.cols());

	computeEigenMatlab(SF, MF, numEigs, eigFieldsNRoSyRef, eigValuesNRoSyRef, filename);
	//WriteSparseMatrixToMatlab(MF2D, "hello");

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}

void NRoSyFields::computeEigenFields_regular(const int &numEigs, const string& filename)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Computing reference regular-eigenproblem (in Matlab)... ";

	computeEigenMatlab(SF, numEigs, eigFieldsNRoSyRef, eigValuesNRoSyRef, filename);
	
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}

void NRoSyFields::storeBasis(const string& filename)
{
	writeEigenSparseMatrixToBinary(Basis, filename);
}

void NRoSyFields::retrieveBasis(const string& filename)
{
	readEigenSparseMatrixFromBinary(filename, Basis);

	//BasisTemp = Basis; 
	//normalizeBasisAbs();
	printf("Basis size=%dx%d (nnz=%.4f)\n", Basis.rows(), Basis.cols(), (double)Basis.nonZeros() / (double)Basis.rows());
}

/* Creating NRoSyFields */
void NRoSyFields::representingNRoSyFields(const Eigen::MatrixXd& NFields)
{

}

void NRoSyFields::constructNRoSyFields(const int& nRoSy, const Eigen::MatrixXd& NFields)
{
	cout << "> NRoSyFields::Construction (" << nRoSy << ") \n";
	nRot = nRoSy;
	scaleMesh();
	computeFaceCenter();
	constructFrameBasis();
	constructMappingMatrix();
	representingNRoSyFields(NFields);
}

void NRoSyFields::constructNRoSyFields(const int& nRoSy, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
	nRot = nRoSy;
	this->V = V; 
	this->F = F;
	scaleMesh();
	constructFrameBasis();
}

void NRoSyFields::computeMaximalPrincipalCurvature(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::VectorXd &PD, Eigen::VectorXd& PV)
{
	cout << "Computing principal curvature\n";
	Eigen::MatrixXd PD1, PD2, PDF;
	Eigen::VectorXd PD2D;
	Eigen::VectorXd PV1, PV2, PVF;

	igl::principal_curvature(V, F, PD1, PD2, PV1, PV2);

	cout << "Mapping to space of triangles \n";
	PDF.resize(F.rows(), F.cols());
	PVF.resize(F.rows());
	for (int i = 0; i < F.rows(); i++)
	{
		PDF.row(i) = (PD1.row(F(i, 0)) + PD1.row(F(i, 1)) + PD1.row(F(i, 2))) / 3.0;
		PVF(i) = (PV1(F(i, 0)) + PV1(F(i, 1)) + PV1(F(i, 2))) / 3.0;
	}

	PV = PVF;

	printf("Dim of PDF: %dx%d | A=%dx%d \n", PDF.rows(), PDF.cols(), A.rows(), A.cols());
	cout << "Converting to local coordinates \n";

	PDF.transposeInPlace();
	PD2D = Eigen::Map<Eigen::VectorXd>(PDF.data(), PDF.cols()*PDF.rows());
	//PD2D = Eigen::Map<Eigen::VectorXd>(PDF.data(), 1);

	printf("Dim of PD2D: %dx%d \n", PD2D.rows(), PD2D.cols());
	printf("Size of PV: %d \n", PV.size());

	PD = A.transpose()*PD2D;
	printf("Dim of PD (per face, to draw) : %dx%d \n", PD.rows(), PD.cols());

	cout << "Matrix rep: " << PDF.block(0, 0, 3, 3) << endl << endl;
	cout << "Vector rep: " << PD2D.block(0, 0, 1, 9) << endl << endl;
}

void NRoSyFields::smoothNRoSyFields(double lambda, Eigen::PardisoLDLT<Eigen::SparseMatrix<double>>& sparseSolver, const Eigen::VectorXd& inputFields, Eigen::VectorXd& outputFields)
{
	/* Re-Scale the Lambda */
	

	//Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> sparseSolver;
	//Eigen::SparseMatrix<double>						LHS = MF + mu*SF;
	Eigen::VectorXd									rhs = MF*inputFields;
	//sparseSolver.analyzePattern(LHS);
	//sparseSolver.factorize(LHS);
	outputFields = sparseSolver.solve(rhs);
}

/* Rep. Vectors and N-RoSy Fields interface */
void NRoSyFields::convertNRoSyToRepVectors(const NRoSy& nRoSyFields, Eigen::VectorXd& repVect)
{
	cout << "Converting N-Rosy to reprsentation vectors \n";
	double scale = 1.0;
	Eigen::Vector2d b(1, 0);	
	
	const int nSize = nRoSyFields.magnitude.size();
	repVect.resize(2 * nSize);

	/* Construct rotation matrix*/
	for (int j = 0; j < nSize; j++)
	{
		//double angle = nRoSyFields.theta(j) + (2.0*M_PI / (double)nRot);
		double angle = nRoSyFields.theta(j) * nRot;
		
		Eigen::Matrix2d RotM;
		RotM(0, 0) = cos(angle);
		RotM(0, 1) = -sin(angle);
		RotM(1, 0) = sin(angle);
		RotM(1, 1) = cos(angle);

		repVect.block(2 * j, 0, 2, 1) = nRoSyFields.magnitude(j) * RotM * b;
	}
}

void NRoSyFields::convertRepVectorsToNRoSy(const Eigen::VectorXd& repVect, NRoSy& nRoSyFields)
{
	const int nSize = repVect.size() / 2;
	nRoSyFields.theta.resize(nSize);
	nRoSyFields.magnitude.resize(nSize);

	/* Temp variables */
	Eigen::Vector2d v, b(1,0);

	for (int i = 0; i < nSize; i++)
	{
		v = repVect.block(2 * i, 0, 2, 1);

		nRoSyFields.magnitude(i) = v.norm();
		//nRoSyFields.magnitude(i) = 1.0;		// normalizing the n-fields
		v.normalize();
		
		double dotP = b.dot(v);
		double angle = acos(dotP);
		if (v(1) < 0)
		{
			angle = 2 * M_PI - angle;
			nRoSyFields.theta(i) = angle / (double) nRot;
		}
		else
		{
			nRoSyFields.theta(i) = angle / (double) nRot;
		}
	}
}

void NRoSyFields::createNRoSyFromVectors(const Eigen::VectorXd& vectorFields)
{
	cout << "Converting to nRoSy fields \n";
	const int nSize = vectorFields.size() / 2;
	this->nRoSy.theta.resize(nSize);
	this->nRoSy.magnitude.resize(nSize);

	/* Temp variables */
	Eigen::Vector2d v, b(1, 0);

	for (int i = 0; i < nSize; i++)
	{
		v = vectorFields.block(2 * i, 0, 2, 1);

		if (i % 50 == 0)
		{
			//printf("Data %d = (%.5f, %.5f) \n", i, v(0), v(1));
		}

		this->nRoSy.magnitude(i) = v.norm();
		v.normalize();
		double angle;
		if (v(1)<0)
			this->nRoSy.theta(i) = M_PI - acos(b.dot(v));
		else
			this->nRoSy.theta(i) = acos(b.dot(v));
	}
}

void NRoSyFields::createNRoSyFromVectors(const Eigen::VectorXd& vectorFields, NRoSy& nRoSyFields)
{
	cout << "Converting to nRoSy fields ( on " << vectorFields.size()/2 << " faces). \n";
	int nSize = vectorFields.size() / 2;
	nRoSyFields.theta.resize(nSize);
	nRoSyFields.magnitude.resize(nSize);

	/* Temp variables */
	Eigen::Vector2d v, b(1, 0);

	srand(time(NULL));
	const int iTest = rand() % nSize;

	for (int i = 0; i < nSize; i++)
	{
		v = vectorFields.block(2 * i, 0, 2, 1);

		
		if (i % 50 == 0)
		{
			//printf("Data %d = (%.5f, %.5f) \n", i, v(0), v(1));
		}

		nRoSyFields.magnitude(i) = v.norm();
		v.normalize();
		double dotP = b.dot(v);
		double angle = acos(dotP);
		if (v(1) < 0)
		{
			nRoSyFields.theta(i) = 2*M_PI - angle;
		}
		else
		{
			nRoSyFields.theta(i) = angle;
		}

		if (i == iTest)
		{
			//printf("data for %d: [%.2f, %.2f]\n", i, v(0), v(1));
			printf("[%d] v is <%.3f, %.3f>, with b is <%.3f, %.3f> \n", i, v(0), v(1), b(0), b(1));
			printf("__whose inner angle = %.3f yielding => %.3f\n", angle*180.0 / M_PI, nRoSyFields.theta(i)*180.0 / M_PI);
		}
	}
}

/* Visualizing the NRoSyFields */
void NRoSyFields::visualizeNRoSyFields(igl::opengl::glfw::Viewer &viewer, const NRoSy& nRoSyFields, const Eigen::RowVector3d& color)
{
	//double scale = 250.0;
	double scale = 1.0;
	//double scale = 2.5;
	//double scale = 0.25;
	//double scale = 0.1;
	//double scale = 50.0; 
	//double scale = 0.001;
	Eigen::Vector2d b(1, 0);
	

	for (int i = 0; i < nRot; i++)
	{
		Eigen::VectorXd TempFields(2 * F.rows());
		//cout << "Drawing the " << i << " fields \n";

		/* Construct rotation matrix*/
		//for (int j = 0; j < F.rows(); j++)
		for(int j:FaceToDraw)
		{
			double angle = nRoSyFields.theta(j) + ((double)i*2.0*M_PI / (double)nRot);
			//if (j == 0)
			//{
			//	printf("angle 0=%.5f, theta 0=%.5f\n", angle, nRoSyFields.theta(j));
			//}
			Eigen::Matrix2d RotM;
			RotM(0, 0) =  cos(angle);
			RotM(0, 1) = -sin(angle);
			RotM(1, 0) =  sin(angle);
			RotM(1, 1) =  cos(angle);

			TempFields.block(2 * j, 0, 2, 1) = nRoSyFields.magnitude(j) * RotM * b;
		}

		visualize2Dfields(viewer, TempFields, color, scale, false);
	}
}

void NRoSyFields::visualizeRepVectorFields(igl::opengl::glfw::Viewer &viewer, const NRoSy& nRoSyFields, const Eigen::RowVector3d& color)
{
	const double scale = 1.0;
	Eigen::Vector2d b(1, 0);	
	//Eigen::RowVector3d color (117.0/255.0, 107.0/255.0, 177.0/255.0);
	Eigen::VectorXd TempFields(2 * F.rows());

	/* Construct rotation matrix*/
	for (int j = 0; j < F.rows(); j++)
	{
		double angle = nRot * nRoSyFields.theta(j);
		TempFields(2 * j)		= nRoSyFields.magnitude(j)*cos(angle);
		TempFields(2 * j + 1)   = nRoSyFields.magnitude(j)*sin(angle);		
	}
	visualize2Dfields(viewer, TempFields, color, scale, false);	
}

void NRoSyFields::visualizeRepVectorFields(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd& repVector, const Eigen::RowVector3d& color)
{
	viewer.data().clear();
	viewer.data().set_mesh(V, F);

	//const double scale = 250.0;
	const double scale = 1.0;
	//Eigen::RowVector3d color(117.0 / 255.0, 107.0 / 255.0, 177.0 / 255.0);
	visualize2Dfields(viewer, repVector, color, scale, false);
}

void NRoSyFields::visualize2Dfields(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &field2D, const Eigen::RowVector3d &color, const double& scale, const bool& normalized)
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

void NRoSyFields::visualizeEigenFields(igl::opengl::glfw::Viewer &viewer, const int id)
{
	viewer.data().clear();
	viewer.data().set_mesh(V, F);

	Eigen::RowVector3d color(0.1, 0.1, 0.9);

	NRoSy nRoSy_eigenFields;
	convertRepVectorsToNRoSy(eigFieldsNRoSyRef.col(id), nRoSy_eigenFields);
	visualizeNRoSyFields(viewer, nRoSy_eigenFields, color);
}

void NRoSyFields::visualizeBasis(igl::opengl::glfw::Viewer &viewer, const int &id)
{
	viewer.data().clear();
	viewer.data().set_mesh(V, F);

	int bId = id;
	Eigen::RowVector3d color;
	if (id % 2 == 0) {
		color = Eigen::RowVector3d(1.0, 0.1, 0.2);
		//color = Eigen::RowVector3d(0.6, 0.2, 0.2);
	}
	else {
		color = Eigen::RowVector3d(0.0, 0.1, 1.0);
	}

	if (id >= 2 * Sample.size()) {
		bId = 2 * Sample.size() - 1;
	}

	printf("Showing the %d BasisTemp field (Sample=%d) \n", bId, Sample[id / 2]);
	Eigen::VectorXd basisRepVectors;
	NRoSy basisNRoSy; 

	Eigen::RowVector3d colorEven(0.1, 0.1, 0.9);
	Eigen::RowVector3d colorOdd(0.1, 0.1, 0.9);

	cout << "Getting the n-th basis as rep vctors\n";
	basisRepVectors = Basis.col(id);
	convertRepVectorsToNRoSy(basisRepVectors, basisNRoSy);

	cout << "Show is\n";
	if(id%2==0)
		visualizeNRoSyFields(viewer, basisNRoSy, colorEven);
	else 
		visualizeNRoSyFields(viewer, basisNRoSy, colorOdd);
	//visualizeRepVectorFields(viewer, basisRepVectors);
}

void NRoSyFields::visualizeConstrainedFields(igl::opengl::glfw::Viewer &viewer)
{
	cout << "Visualizing the fields \n";
	Eigen::RowVector3d color = Eigen::RowVector3d(0.1, 0.1, 0.9);
	//visualizeRepVectorFields(viewer, Xf);

	//viewer.data().clear();
	//viewer.data().set_mesh(V, F);
	//viewer.data().line_width = 1.0f;
	
	NRoSy nRoSy_;
	//visualizeRepVectorFields(viewer, Xf, color);
	convertRepVectorsToNRoSy(Xf, nRoSy_);
	visualizeNRoSyFields(viewer, nRoSy_, color);
}

void NRoSyFields::visualizeConstrainedFields_Reduced(igl::opengl::glfw::Viewer &viewer)
{
	cout << "Visualizing the fields \n";
	Eigen::RowVector3d color = Eigen::RowVector3d(0.9, 0.1, 0.1);
	//visualizeRepVectorFields(viewer, Xf);

	//viewer.data().clear();
	//viewer.data().set_mesh(V, F);
	//viewer.data().line_width = 1.0f;

	NRoSy nRoSy_;
	//visualizeRepVectorFields(viewer, XfBar, color);
	convertRepVectorsToNRoSy(XfBar, nRoSy_);
	visualizeNRoSyFields(viewer, nRoSy_, color);
}

void NRoSyFields::visualizeConstraints(igl::opengl::glfw::Viewer &viewer)
{
	cout << "visualizing the constraints \n";
	/* ORIGINAL + OVERLAY on 2nd Mesh */
	/* Some constants for arrow drawing */
	const double HEAD_RATIO = 3.0;
	const double ARRAW_RATIO = 4.0;
	const double EDGE_RATIO = 1.5;
	double lengthScale = EDGE_RATIO*avgEdgeLength;
	Eigen::RowVector3d color2(0.2, 0.8, 0.1);
	Eigen::RowVector3d color(0.1, 0.1, 0.1);

	cout << "Getting face normals \n";
	Eigen::MatrixXd NF;
	igl::per_face_normals(V, F, NF);

	cout << "Drawing on the overlay mesh \n";
	viewer.selected_data_index = 0;
	viewer.data().line_width = 5.0;
	viewer.data().point_size = 5.0;
	viewer.data().show_lines = false;

	cout << "Computing rotaiton angle\n";
	/* Computing the rotation angle for 1:3 ratio of arrow head */
	double rotAngle = M_PI - atan(1.0 / 3.0);
	Eigen::Matrix2d rotMat1, rotMat2;
	rotMat1 << cos(rotAngle), -sin(rotAngle), sin(rotAngle), cos(rotAngle);
	rotMat2 << cos(-rotAngle), -sin(-rotAngle), sin(-rotAngle), cos(-rotAngle);

	NRoSy nRoSy_;
	convertRepVectorsToNRoSy(c, nRoSy_);
	Eigen::VectorXd c2(c.size());


	Eigen::MatrixXd ALoc(3, 2);
	for (int j = 0; j < nRot; j++)
	{
		for (int i = 0; i < globalConstraints.size(); i++) {
			//cout << "Setting up the " << i << "-th constraint.\n";
			Eigen::RowVector3d cc, g, h1, h2, v1, v2, v3, e, f, n;
			Eigen::Vector2d v, id;


			id << 1.0, 0.0;
			double angle = nRoSy_.theta(i) + ((double)j*2.0*M_PI / (double)nRot);
			double cosA = cos(angle);
			double sinA = sin(angle);
			Eigen::Matrix2d RotA; RotA << cosA, -sinA, sinA, cosA;
			c2.block(2 * i, 0, 2, 1) = RotA*id;

			//cout << "Fetching the face center + nromal\n";
			cc = FC.row(globalConstraints[i]);
			n = NF.row(globalConstraints[i]);
			n *= (avgEdgeLength / 10.0);
			//cout << "Getting local c\n";
			ALoc = A.block(3 * globalConstraints[i], 2 * globalConstraints[i], 3, 2);
			v = c2.block(2 * i, 0, 2, 1);
			
			g = (ALoc * v).transpose();
			cc += n;
			f = cc + g*lengthScale;
			//cout << "Adding new entries on viewer \n";
			//cout << "cc: " << cc << endl;
			//cout << "f:" << f << endl;
			//cout << "colro: " << color << endl; 
			viewer.data().add_edges(cc, f, color);
			viewer.data().add_points(cc, color2);

			//cout << "Drawing the head \n";
			h1 = (ALoc * (rotMat1*c2.block(2 * i, 0, 2, 1))).transpose();
			h2 = (ALoc * (rotMat2*c2.block(2 * i, 0, 2, 1))).transpose();

			//cout << "h1=" << h1 << endl;
			//cout << "h2=" << h2 << endl;
			//cout << "scale: " << lengthScale << endl;
			//cout << "head ration: " << HEAD_RATIO << endl; 
			viewer.data().add_edges(f, f + h1*lengthScale / HEAD_RATIO, color);
			viewer.data().add_edges(f, f + h2*lengthScale / HEAD_RATIO, color);
		}
	}
	viewer.selected_data_index = 0;
	viewer.data().line_width = 1.0;
}

void NRoSyFields::visualizeSoftConstraints(igl::opengl::glfw::Viewer &viewer)
{
	/* Color */
	Eigen::RowVector3d color(0.1, 0.1, 0.1);

	/* Some constants for arrow drawing */
	const double HEAD_RATIO = 5.0;
	const double EDGE_RATIO = 2.0;

	/* Computing the rotation angle for 1:3 ratio of arrow head */
	double rotAngle = M_PI - atan(1.0 / 3.0);
	Eigen::Matrix2d rotMat1, rotMat2;
	rotMat1 << cos(rotAngle), -sin(rotAngle), sin(rotAngle), cos(rotAngle);
	rotMat2 << cos(-rotAngle), -sin(-rotAngle), sin(-rotAngle), cos(-rotAngle);

	double lengthScale = EDGE_RATIO*avgEdgeLength;
	Eigen::RowVector3d f, h1, h2, e, c;
	Eigen::Vector2d v;
	Eigen::MatrixXd ALoc(3, 2);
	int face;
	for (int i = 0; i < constraintVect2D.size(); i++)
	{
		for (int j = 0; j < constraintVect2D[i].size(); j++)
		{
			face = curvesConstraints[i][j];
			c = FC.row(face);
			//f = VectorBlock.row(i);
			v = constraintVect2D[i][j];
			ALoc = A.block(3 * face, 2 * face, 3, 2);
			f = (ALoc * v).transpose();
			h1 = (ALoc* (rotMat1*v)).transpose();
			h2 = (ALoc* (rotMat2*v)).transpose();
			e = c + f*lengthScale;
			//cout << "c: " << c << "e: " << e << endl; 
			viewer.data().add_edges(c, e, color);
			viewer.data().add_edges(e, e + h1*lengthScale / HEAD_RATIO, color);
			viewer.data().add_edges(e, e + h2*lengthScale / HEAD_RATIO, color);
		}
	}
}

/* ============================= N-FIELDS DESIGN ============================= */

void NRoSyFields::nRoSyFieldsDesignRef()
{
	//nRoSyFieldsDesignRef_HardConstraints();
	//nRoSyFieldsDesignRef_SoftConstraints();
	nRoSyFieldsDesignRef_Splines();
}

void NRoSyFields::nRoSyFieldsDesignRef_HardConstraints()
{
	Eigen::VectorXd					b, g, h, vEst;
	Eigen::SparseMatrix<double>		A_LHS;
	//Eigen::SparseMatrix<double>		B2F;
	if (MF.rows() < 1) constructMassMatrixMF3D();
	if (SF.rows() < 1) buildStiffnessMatrix_Geometric();
	//B2F = SF * MFinv * SF;

	constructRandomHardConstraints(C, c);
	setupRHSBiharmSystemRef(BF, C, c, g, h, vEst, b);
	setupLHSBiharmSystemRef(BF, C, c, A_LHS);
	solveBiharmSystemRef(vEst, A_LHS, b, Xf);
}

void NRoSyFields::createAlignmentField(Eigen::VectorXd& v)
{
	// Create alignment fields
	// based on maximal curvature direction
	Eigen::VectorXd PD, PV, v1;
	computeMaximalPrincipalCurvature(V, F, PD, PV);

	

	NRoSy nRoSy_;
	createNRoSyFromVectors(PD, nRoSy_);

	printf("PD size=%d | PV=%d |n-RoSy: %d and %d \n", PD.size(), PV.size(), nRoSy_.magnitude.size(), nRoSy_.theta.size());

	//nRoSy_.magnitude = PV;
	convertNRoSyToRepVectors(nRoSy_, v1);

	Eigen::VectorXd id(MF.rows());
	id.setConstant(1.0);
	double factor1 = id.transpose()*MF*id;
	double factor2 = id.transpose()*SF*id;
	double lambda = 50;
	double mu = lambda * factor1 / factor2;

	Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> sparseSolver;
	Eigen::SparseMatrix<double>						LHS = MF + mu*SF;
	sparseSolver.analyzePattern(LHS);
	sparseSolver.factorize(LHS);

	/* Smoothing the curvature fields */
	smoothNRoSyFields(mu, sparseSolver, v1, v);	
	

	for(int i=0; i<5; i++)
		smoothNRoSyFields(mu, sparseSolver, v, v);

	convertRepVectorsToNRoSy(v, nRoSy_);
	
}

void NRoSyFields::constructRandomHardConstraints(Eigen::SparseMatrix<double>& C, Eigen::VectorXd& c)
{
	// Define the constraints
	const bool readFromFile = false;			/// IMPORTANT!!!!!!!!!!!!!!!!!!!!
	bool lineNotFound = true;
	//string filename = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Constraints/Constraints_CDragon_Rand_20.txt";;
	string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Armadillo_randConstraints.txt";
	//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Fertility_randConstraints.txt";
	//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/CDragon_randConstraints.txt";
	//string filename = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Constraints/Constraints_Cube_Rand_25.txt";

	/* Reading the constraints from file */
	if (readFromFile)
	{
		/* Old type of reading */
		//cout << "____Loading constraints from file \n";
		//LoadSTDVectorFromTxtFile(filename, globalConstraints);

		/* New read: all constraints are in a single line */
		ifstream file(resultFile);
		string oneLine, oneWord;

		int lineNow = 0;
		int toRead = testID;
		vector<double> constraint;
		constraint.reserve(100);

		if (file.is_open())
		{
			//getline(file, oneLine);
			while (getline(file, oneLine) && lineNotFound)
			{
				if (toRead == lineNow)
				{
					istringstream iStream(oneLine);
					getline(iStream, oneWord, ',');

					while (oneWord != "") {
						constraint.push_back(stod(oneWord));
						cout << oneWord << "|";
						getline(iStream, oneWord, ',');
					}
					cout << endl;

					globalConstraints.resize(constraint.size());

					for (int i = 0; i < constraint.size(); i++) {
						globalConstraints[i] = constraint[i];
					}

					lineNotFound = false;
					//return;
				}

				cout << "line =" << lineNow << endl;
				lineNow++;
			}
		}
		file.close();
	}
	else
	{
		/* Random number generator */
		std::random_device rd;								// Will be used to obtain a seed for the random number engine
		std::mt19937 gen(rd());								// Standard mersenne_twister_engine seeded with rd()

		/* Creating random constraints */
		//std::uniform_int_distribution<> disConst(20, 50); // From 0 to F.rows()-1
		std::uniform_int_distribution<> disConst(5,10); // From 0 to F.rows()-1
		int numConstraints = disConst(gen);
		cout << "we gave you " << numConstraints << " number of constnraints. Good luck!\n";

		//int numConstraints = 20;
		set<int> constraints;
		globalConstraints.resize(numConstraints);

		do {
			std::uniform_int_distribution<> dis(0, F.rows() - 1); // From 0 to F.rows()-1
			int constraintFace = dis(gen);
			constraints.insert(constraintFace);
		} while (constraints.size() < numConstraints);

		int counter1 = 0;
		for (int i : constraints) {
			globalConstraints[counter1++] = i;
		}

		bool writeToFile = false;

		if (writeToFile) {
			std::ofstream ofs;
			ofs.open(resultFile, std::ofstream::out | std::ofstream::app);

			for (int i : globalConstraints) {
				ofs << i << ",";
			}
			ofs << "\n";																												// l2-norm

			ofs.close();
			//WriteSTDVectorToTxtFile(globalConstraints, filename);
		}
	}
	
	cout << "Setting up matrix C\n";
	// Setting up matrix C
	Eigen::SparseMatrix<double> CTemp;
	vector<Eigen::Triplet<double>> CTriplet;
	CTriplet.reserve(2 * globalConstraints.size());
	Eigen::VectorXd cTemp(2 * globalConstraints.size());
	//c.resize(2 * globalConstraints.size());
	Eigen::Vector2d cRand;

	int counter = 0;
	for (int i = 0; i < globalConstraints.size(); i++) {
		//cRand(0) = (double)(rand() % F.rows()) / (double) F.rows();
		//cRand(1) = (double)(rand() % F.rows()) / (double)F.rows();
		cRand << 1.0, 0.0;
		cRand.normalize();

		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * globalConstraints[i] + 0, 1.0));
		cTemp(counter, 0) = cRand(0);
		counter++;

		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * globalConstraints[i] + 1, 1.0));
		cTemp(counter, 0) = cRand(1);
		counter++;
	}
	cout << "__Setting the selection matrix for the constraints \n";
	C.resize(2 * globalConstraints.size(), SF.rows());
	C.setFromTriplets(CTriplet.begin(), CTriplet.end());

	/* representing c as representation vectors */
	NRoSy nRoSy_;
	createNRoSyFromVectors(cTemp, nRoSy_);
	convertNRoSyToRepVectors(nRoSy_, c);
}

void NRoSyFields::setupRHSBiharmSystemRef(const Eigen::SparseMatrix<double>& B2F, const Eigen::SparseMatrix<double>& C, const Eigen::VectorXd& c, Eigen::VectorXd& g, Eigen::VectorXd& h, Eigen::VectorXd& vEst, Eigen::VectorXd& b)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Constructing RHS... ";

	vEst.resize(B2F.cols());
	for (int i = 0; i < vEst.rows(); i++) {
		vEst(i) = 0.5;
	}

	// Create alignment fields based on maximal curvature direction
	Eigen::VectorXd repV;
	createAlignmentField(repV);

	//g = (B2F+MF) * vEst +MF*repV;
	g = (B2F) * vEst + MF*repV;
	b.resize(B2F.rows() + c.rows(), c.cols());

	// First column of b
	h = C * vEst - c;
	b << g, h;

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}
void NRoSyFields::setupLHSBiharmSystemRef(const Eigen::SparseMatrix<double>& B2F, const Eigen::SparseMatrix<double>& C, const Eigen::VectorXd& c, Eigen::SparseMatrix<double>& A_LHS)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Constructing LHS... ";

	A_LHS.resize(B2F.rows() + C.rows(), B2F.cols() + C.rows());

	vector<Eigen::Triplet<double>>	ATriplet;
	ATriplet.reserve(10 * B2F.rows());		// It should be #rows x 4 blocks @ 2 elements (8) + #constraints,
											// but made it 10 for safety + simplicity

	for (int k = 0; k < B2F.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(B2F, k); it; ++it) {
			ATriplet.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
		}
	}

	for (int k = 0; k < C.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(C, k); it; ++it) {
			ATriplet.push_back(Eigen::Triplet<double>(B2F.rows() + it.row(), it.col(), it.value()));
			ATriplet.push_back(Eigen::Triplet<double>(it.col(), B2F.cols() + it.row(), it.value()));
		}
	}
	A_LHS.setFromTriplets(ATriplet.begin(), ATriplet.end());

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}
void NRoSyFields::solveBiharmSystemRef(const Eigen::VectorXd& vEst, const Eigen::SparseMatrix<double>& A_LHS, const Eigen::VectorXd& b, Eigen::VectorXd& Xf)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Solving the global system (Pardiso LDLT)... \n";


	//cout << "Starting to solve problem." << endl;
	Xf.resize(SF.rows());

	// Setting up the solver
	//Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> sparseSolver(A_LHS);
	Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> sparseSolver(A_LHS);
	//Eigen::PastixLDLT<Eigen::SparseMatrix<double>,1> sparseSolver(A_LHS);

	// FIRST BASIS
	cout << "....Solving first problem (first frame)..." << endl;
	Eigen::VectorXd x = sparseSolver.solve(b);

	if (sparseSolver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system. " << endl;
		if (sparseSolver.info() == Eigen::NumericalIssue)
			cout << "NUMERICAL ISSUE. " << endl;
		if (sparseSolver.info() == Eigen::InvalidInput)
			cout << "Input is Invalid. " << endl;
		cout << sparseSolver.info() << endl;
		return;
	}

	Xf = -x.block(0, 0, SF.rows(), 1) + vEst;

	printf("____Xf size is %dx%d\n", Xf.rows(), Xf.cols());

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}

void NRoSyFields::nRoSyFieldsDesignRef_Splines()
{
	Eigen::VectorXd					b, g, h, vEst;
	Eigen::SparseMatrix<double>		A_LHS;
	//Eigen::SparseMatrix<double>		B2F;
	if (MF.rows() < 1) constructMassMatrixMF3D();
	if (SF.rows() < 1) buildStiffnessMatrix_Geometric();
	//B2F = SF * MFinv * SF;

	double weight = 0;
	for (int i = 0; i < userVisualConstraints.size(); i += 2)
	{
		weight += doubleArea(userVisualConstraints[i]);
	}

	vector<double> lambda(2);
	lambda[0] = 1.0;
	lambda[1] = 0.001/weight;

	//constructRandomHardConstraints(C, c);
	constructInteractiveConstraints();
	setupRHSBiharmSystemRef_Chris(BF, lambda, c, g, h, b);
	setupLHSBiharmSystemRef_Chris(BF, lambda, C, A_LHS);
	solveBiharmSystemRef_Chris(A_LHS, b, Xf);
}

void NRoSyFields::setupRHSBiharmSystemRef_Chris(const Eigen::SparseMatrix<double>& B2F, const vector<double>& lambda, const Eigen::VectorXd& c, Eigen::VectorXd& g, Eigen::VectorXd& h, Eigen::VectorXd& b)
{
	// Create alignment fields based on maximal curvature direction
	//Eigen::VectorXd repV;
	//createAlignmentField(repV);
	
	b.resize(B2F.nonZeros() + c.rows(), c.cols());

	g = lambda[1]*MF*alignFields;
	h = c;
	b << g, h; 

}
void NRoSyFields::setupLHSBiharmSystemRef_Chris(const Eigen::SparseMatrix<double>& B2F, const vector<double>& lambda, const Eigen::SparseMatrix<double>& C, Eigen::SparseMatrix<double>& A_LHS)
{
	A_LHS.resize(B2F.rows() + C.rows(), B2F.cols() + C.rows());

	vector<Eigen::Triplet<double>>	ATriplet;
	ATriplet.reserve(10 * B2F.rows());		// It should be #rows x 4 blocks @ 2 elements (8) + #constraints,
											// but made it 10 for safety + simplicity

	//Eigen::SparseMatrix<double> BM = lambda[0]*B2F + lambda[1]*MF;

	// From B and M matrices
	//for (int k = 0; k < B2F.outerSize(); ++k) {
	//	for (Eigen::SparseMatrix<double>::InnerIterator it(B2F, k); it; ++it) {
	//		ATriplet.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));			
	//	}
	//}

	BM = lambda[0] * BF + lambda[1] * MF;

	for (int k = 0; k < BM.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(BM, k); it; ++it) {
			ATriplet.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
		}
	}

	// From the constraint matrix C
	for (int k = 0; k < C.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(C, k); it; ++it) {
			ATriplet.push_back(Eigen::Triplet<double>(B2F.rows() + it.row(), it.col(), it.value()));
			ATriplet.push_back(Eigen::Triplet<double>(it.col(), B2F.cols() + it.row(), it.value()));
		}
	}
	A_LHS.setFromTriplets(ATriplet.begin(), ATriplet.end());

}
void NRoSyFields::solveBiharmSystemRef_Chris(const Eigen::SparseMatrix<double>& A_LHS, const Eigen::VectorXd& b, Eigen::VectorXd& Xf)
{
	Xf.resize(SF.rows());

	// Setting up the solver
	//Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> sparseSolver(A_LHS);
	Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> sparseSolver(A_LHS);

	// FIRST BASIS
	cout << "....Solving first problem (first frame)..." << endl;
	Eigen::VectorXd x = sparseSolver.solve(b);

	if (sparseSolver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system. " << endl;
		if (sparseSolver.info() == Eigen::NumericalIssue)
			cout << "NUMERICAL ISSUE. " << endl;
		if (sparseSolver.info() == Eigen::InvalidInput)
			cout << "Input is Invalid. " << endl;
		cout << sparseSolver.info() << endl;
		return;
	}

	Xf = x.block(0, 0, SF.rows(), 1);
}

/* ============================= SOFT CONSTRAINTS ============================= */

void NRoSyFields::nRoSyFieldsDesignRef_SoftConstraints()
{
	cout << "NFIelds design\n";
	Eigen::VectorXd					b, g, h, vEst;
	Eigen::SparseMatrix<double>		A_LHS;

	//constructSoftConstraints();

	Eigen::Vector3d lambda; 	
	lambda(0) = 1000 * MF.diagonal().sum() / SF.diagonal().sum();
	lambda(1) = 0;
	lambda(2) = 1.0;
	setupRHSGlobalProblemSoftConstraints(lambda, b);
	setupLHSGlobalProblemSoftConstraints(lambda, A_LHS);	
	solveGlobalSystemMappedLDLTSoftConstraints(vEst, A_LHS, b);

	cout << "Lambda 0 is " << lambda(0) << endl; 
}

void NRoSyFields::constructSoftConstraints()
{
	const int NUM_CURVES = 8;
	curvesConstraints.resize(NUM_CURVES);

	srand(time(NULL));
	int init_, end_;
	vector<int> aCurve;
	/* Automatic random set up */
	//for (int i = 0; i < NUM_CURVES; i++) {
	//	init_ = rand() % F.rows();
	//	//end_ = rand() % F.rows();
	//	end_ = init_ + 40;
	//	constructCurvesAsConstraints(init_, end_, aCurve);
	//	curvesConstraints[i] = aCurve; 
	//}


	int constCounter = 0;
	/* Manual set-up for Chinese Dragon */
	//// Face
	//constructCurvesAsConstraints(152474, 51474, aCurve);
	//curvesConstraints[constCounter++] = aCurve;
	//
	//// Back
	//constructCurvesAsConstraints(44109, 68907, aCurve);
	//curvesConstraints[constCounter++] = aCurve;
	//
	//// body - bottom
	//constructCurvesAsConstraints(13471, 195817, aCurve);
	//curvesConstraints[constCounter++] = aCurve;
	//
	//// body - right
	//constructCurvesAsConstraints(123036, 247143, aCurve);
	//curvesConstraints[constCounter++] = aCurve;
	//
	//// body - left
	//constructCurvesAsConstraints(234815, 232296, aCurve);
	//curvesConstraints[constCounter++] = aCurve;
	//
	//// front_right_leg
	//constructCurvesAsConstraints(75468, 7716, aCurve);
	//curvesConstraints[constCounter++] = aCurve;
	//
	//// front_left_leg
	//constructCurvesAsConstraints(231495, 77171, aCurve);
	//curvesConstraints[constCounter++] = aCurve;
	//
	//// tail
	//constructCurvesAsConstraints(230301, 113500, aCurve);
	//curvesConstraints[constCounter++] = aCurve;

	/* Manual set-up for Armadillo */
	// Head
	constructCurvesAsConstraints(68818,6278, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	// Stomach
	constructCurvesAsConstraints(56965, 41616, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	// Leg/Foot (R then L)
	constructCurvesAsConstraints(28590, 16119, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(25037, 571, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	// Arm/Hand
	constructCurvesAsConstraints(55454, 6877, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(49059, 36423, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	// Back
	constructCurvesAsConstraints(68331, 72522, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	// Tail
	constructCurvesAsConstraints(24056, 1075, aCurve);
	curvesConstraints[constCounter++] = aCurve;

	/* Project elements to local frame */
	projectCurvesToFrame();

	/* Get the number of constraints */
	int numConstraints = 0;
	for (int i = 0; i < curvesConstraints.size(); i++)
	{
		for (int j = 0; j < curvesConstraints[i].size() - 1; j++)
		{
			numConstraints++;
		}
	}

	/* Setup to constraint matrix */
	Eigen::VectorXd cTemp(2 * numConstraints);
	C.resize(2 * numConstraints, SF.cols());

	int counter = 0;
	int elem;
	vector<Eigen::Triplet<double>> CTriplet;
	for (int i = 0; i < curvesConstraints.size(); i++)
	{
		for (int j = 0; j < curvesConstraints[i].size() - 1; j++)
		{
			elem = curvesConstraints[i][j];
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * elem + 0, 1.0));
			cTemp(counter++) = constraintVect2D[i][j](0);
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * elem + 1, 1.0));
			cTemp(counter++) = constraintVect2D[i][j](1);
		}
	}

	C.setFromTriplets(CTriplet.begin(), CTriplet.end());

	/* representing c as representation vectors */
	NRoSy nRoSy_;
	createNRoSyFromVectors(cTemp, nRoSy_);
	convertNRoSyToRepVectors(nRoSy_, c);

	printf("Size of selector matrix C: %dx%d \n", C.rows(), C.cols());
	printf("Size of constraint vector c: %d\n", c.rows());
}

void NRoSyFields::constructCurvesAsConstraints(const int& init, const int& end, vector<int>& curve)
{
	priority_queue<VertexPair, std::vector<VertexPair>, std::greater<VertexPair>> DistPQueue;
	Eigen::VectorXd D(F.rows());
	Eigen::VectorXi prev(F.rows());

	// Computing distance for initial sample points S
	for (int i = 0; i < F.rows(); i++) {
		D(i) = numeric_limits<double>::infinity();
		prev(i) = -1;
	}

	D(end) = 0.0f;
	VertexPair vp{ end, D(end) };
	DistPQueue.push(vp);

	curve.resize(0);
	curve.shrink_to_fit();
	curve.reserve(F.rows() / 2);

	// For other vertices in mesh
	//double distFromCenter;
	int neigh;
	do {
		if (DistPQueue.size() == 0) break;
		VertexPair vp1 = DistPQueue.top();
		//distFromCenter = vp1.distance;
		DistPQueue.pop();

		// Updating the distance for neighbors of vertex of lowest distance in priority queue
		int const elem = vp1.vId;
		Eigen::Vector3d const c1 = (V.row(F(elem, 0)) + V.row(F(elem, 1)) + V.row(F(elem, 2))) / 3.0;
		for (auto it = 0; it != F.cols(); ++it) {
			/* Regular Dikjstra */
			neigh = AdjMF3N(elem, it);
			Eigen::Vector3d const c2 = FC.row(neigh);
			double dist = (c2 - c1).norm();
			//double tempDist = D(elem) + dist;
			double tempDist = (FC.row(end) - FC.row(neigh)).norm();

			/* updating the distance */
			if (tempDist < D(neigh)) {
				D(neigh) = tempDist;
				VertexPair vp2{ neigh,tempDist };
				DistPQueue.push(vp2);
				prev(neigh) = vp1.vId;
			}
		}
	} while (!DistPQueue.empty());

	// Obtaining the path <reverse>
	int u = init;
	while (prev[u] != -1 && u != end) {
		curve.push_back(u);
		u = prev(u);
	}

	printf("Path from %d to %d has %d elements.\n", init, end, curve.size());


	curve.shrink_to_fit();
}

void NRoSyFields::measureSoftConstraintError(const Eigen::Vector3d& lambda)
{
	///setupGlobalProblem(lambda);
	///setAndSolveUserSystem(lambda);

	Eigen::VectorXd diff = (Xf - XfBar);
	double error = diff.transpose()*MF*diff;
	double ref = Xf.transpose()*MF*Xf;
	double relError = sqrt(error / ref);
	cout << "The l2-norm error is " << relError << endl;
	cout << "__The rel error is: " << error / ref << endl;
	cout << "__The apprxo length is: " << XfBar.transpose()*MF*XfBar << endl;
	cout << "__The difference length is: " << error << endl;
	cout << "__The reference length is: " << ref << endl;

	error = diff.transpose()*SF*diff;
	ref = Xf.transpose()*SF*Xf;
	relError = sqrt(error / ref);
	cout << "The rel. energy error is " << relError << endl;
	cout << "__The rel energy is: " << error / ref << endl;
	cout << "__The apprxo energy is: " << XfBar.transpose()*SF*XfBar << endl;
	cout << "__The difference energy is: " << error << endl;
	cout << "__The reference energy is: " << ref << endl;
	cout << "__The apprxo energy is: " << XfBar.transpose()*SF*XfBar << endl;
}

void NRoSyFields::projectCurvesToFrame()
{
	Eigen::MatrixXd ALoc(3, 2);
	Eigen::Vector2d vec2D;
	Eigen::Vector3d vec3D;
	int face1, face2, face3;

	constraintVect2D.resize(curvesConstraints.size());
	for (int i = 0; i < curvesConstraints.size(); i++)
	{
		const int curveSize = curvesConstraints[i].size() - 1;
		constraintVect2D[i].resize(curveSize);
		for (int j = 0; j < curvesConstraints[i].size() - 1; j++)
		{
			face1 = curvesConstraints[i][j];
			face2 = curvesConstraints[i][j + 1];
			ALoc = A.block(3 * face1, 2 * face1, 3, 2);

			//if (j < curvesConstraints[i].size() - 2)
			//{
			//	face3 = curvesConstraints[i][j + 2];
			//	//face3 = curvesConstraints[i][curveSize-1];
			//	vec3D = (FC.row(face3) - FC.row(face1)).transpose();
			//}
			//else
			//{
			//	vec3D = (FC.row(face2) - FC.row(face1)).transpose();
			//}

			vec3D = FC.row(curvesConstraints[i][curvesConstraints[i].size() - 1]) - FC.row(curvesConstraints[i][0]);

			vec2D = ALoc.transpose() * vec3D;
			vec2D.normalize();
			constraintVect2D[i][j] = vec2D;
			//cout << "vec2D= " << vec2D << endl; 
		}
	}

	cout << "Fields are projected to 2d frame " << endl;
}

void NRoSyFields::setupRHSGlobalProblemSoftConstraints(const Eigen::Vector3d& lambda, Eigen::VectorXd& b)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Setting up the RHS of the system... ";

	Eigen::SparseMatrix<double> Mconst = C*MF*C.transpose();

	b = lambda(2) * C.transpose() * Mconst * c;
	//b = lambda(2) * C.transpose() * c;

	printf("Size of b: %d. \n", b.rows());

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}

void NRoSyFields::setupLHSGlobalProblemSoftConstraints(const Eigen::Vector3d& lambda, Eigen::SparseMatrix<double>& A_LHS)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Setting up the LHS of the system... ";

	//const double lambda_1 = 10000 / B2D.coeff(0, 0);

	Eigen::SparseMatrix<double> Mconst = C*MF*C.transpose();

	//A_LHS = lambda(0)*SF2DAsym + lambda(1)*B2D + lambda(2)*C.transpose()*C;
	A_LHS = lambda(0)*SF + lambda(2)*C.transpose()*Mconst*C;
	//A_LHS = lambda(0)*SF2DAsym + lambda(2)*C.transpose()*C;
	//A_LHS = lambda(0)*SF2DAsym + lambda_1*B2D + lambda(2)*C.transpose()*C;

	printf("Size of LHS: %dx%d. \n", A_LHS.rows(), A_LHS.cols());

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}

void NRoSyFields::solveGlobalSystemMappedLDLTSoftConstraints(const Eigen::VectorXd& vEst, Eigen::SparseMatrix<double>& A_LHS, Eigen::VectorXd& b)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Solving the global system (Pardiso LDLT)... \n";


	//cout << "Starting to solve problem." << endl;
	Xf.resize(SF.rows());

	// Setting up the solver
	//Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> sparseSolver(A_LHS);
	Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> sparseSolver(A_LHS);
	//Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>> sparseSolver(A_LHS);
	//Eigen::PastixLDLT<Eigen::SparseMatrix<double>,1> sparseSolver(A_LHS);

	// FIRST BASIS
	cout << "....Solving first problem (first frame)..." << endl;
	Xf = sparseSolver.solve(b);

	if (sparseSolver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system. " << endl;
		if (sparseSolver.info() == Eigen::NumericalIssue)
			cout << "NUMERICAL ISSUE. " << endl;
		if (sparseSolver.info() == Eigen::InvalidInput)
			cout << "Input is Invalid. " << endl;
		cout << sparseSolver.info() << endl;
		return;
	}

	//Xf = x;

	printf("____Xf size is %dx%d\n", Xf.rows(), Xf.cols());

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}

// user interactive constraints
void NRoSyFields::pushNewUserConstraints(const int& fInit, const int& fEnd)
{
	userVisualConstraints.push_back(fInit);
	userVisualConstraints.push_back(fEnd);
}

void NRoSyFields::constructInteractiveConstraints()
{
	/* Define the constraints */
	const int numConstraints = userVisualConstraints.size() / 2;
	globalConstraints.resize(numConstraints);
	vector<Eigen::Vector2d> constraintValues(numConstraints);

	/* Global constraints from user input */
	for (int i = 0; i < userVisualConstraints.size(); i += 2)
	{
		/* Location of constraints */
		globalConstraints[i / 2] = userVisualConstraints[i];

		/* Getting the constraints + making them into local coordinates */
		Eigen::RowVector3d dir = FC.row(userVisualConstraints[i + 1]) - FC.row(userVisualConstraints[i]);
		Eigen::MatrixXd ALoc(3, 2);
		ALoc = A.block(3 * userVisualConstraints[i], 2 * userVisualConstraints[i], 3, 2);
		Eigen::Vector2d normDir = ALoc.transpose() * dir.transpose();
		normDir.normalize();
		//cout << "___ constraint in 2D: " << normDir << endl;
		constraintValues[i / 2] = normDir;
	}

	/* Setting up matrix C and column vector c */
	Eigen::SparseMatrix<double> CTemp;
	vector<Eigen::Triplet<double>> CTriplet;
	CTriplet.reserve(2 * globalConstraints.size());
	//c.resize(2 * globalConstraints.size());
	Eigen::VectorXd cTemp(2 * globalConstraints.size());

	/* Putting the constraints into action */
	for (int i = 0; i < globalConstraints.size(); i++) {
		CTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * globalConstraints[i] + 0, 1.0));
		cTemp(2 * i + 0) = constraintValues[i](0);

		CTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * globalConstraints[i] + 1, 1.0));
		cTemp(2 * i + 1) = constraintValues[i](1);
	}
	C.resize(2 * globalConstraints.size(), SF.rows());
	C.setFromTriplets(CTriplet.begin(), CTriplet.end());

	NRoSy nRoSy_;
	createNRoSyFromVectors(cTemp, nRoSy_);
	convertNRoSyToRepVectors(nRoSy_, c);

}

void NRoSyFields::resetInteractiveConstraints()
{
	userVisualConstraints.clear();
	userVisualConstraints.shrink_to_fit();
}

/* ============================= SUBSPACE CONSTRUCTION ============================= */
void NRoSyFields::constructBasis()
{
	
	constructBasis_LocalEigenProblem();
	printf("Basis has size of %dx%d \n", Basis.rows(), Basis.cols());
	//cout << "Basis:: \n";
	//cout << Basis.block(0, 0, 20, 1) << endl << endl;
}

void NRoSyFields::constructSamples(const int &n)
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

void NRoSyFields::farthestPointSampling()
{
	Sample.resize(numSample);
	Eigen::VectorXd D;
	D.resize(F.rows());

	// Initialize the value of D
	for (int i = 0; i < F.rows(); i++) {
		D(i) = numeric_limits<double>::infinity();
	}

	srand(time(NULL));
	Sample[0] = 0;

	//computeDijkstraDistanceFaceForSampling(Sample[0], D);
	//Eigen::VectorXi::Index maxIndex1;
	//D.maxCoeff(&maxIndex1);
	//Sample[1] = maxIndex1;

	for (int i = 1; i < numSample; i++) {
		Eigen::VectorXi::Index maxIndex;
		computeDijkstraDistanceFaceForSampling(Sample[i - 1], D);
		D.maxCoeff(&maxIndex);
		Sample[i] = maxIndex;
	}

	sampleDistance = D;
}

void NRoSyFields::constructBasis_LocalEigenProblem()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Constructing Basis...\n";

	const int Num_fields = 2;
	// Setup sizes of each element to construct basis
	try {
		Basis.resize(1, 1);
		Basis.data().clear();
		Basis.data().squeeze();
		Basis.resize(Num_fields * F.rows(), 2 * Sample.size());
	}
	catch (string &msg) {
		cout << "Cannot allocate memory for basis.." << endl;
	}

	vector<vector<Eigen::Triplet<double>>> UiTriplet(Sample.size());
	Eigen::SparseMatrix<double> Sh = MF2DhNeg*SF*MF2DhNeg;

	cout << "....Constructing and solving local systems...";
	const int NUM_PROCESS = 4;
	vector<chrono::duration<double>>durations;
	durations.resize(NUM_PROCESS);

	for (int i = 0; i < NUM_PROCESS; i++) {
		durations[i] = t1 - t1;
		//cout << "Init dur " << i<< " = " << durations[i].count() << " seconds" << endl;
	}

	/* Default color for the domain selected */
	localSystem.resize(F.rows());
	for (int fid = 0; fid < F.rows(); fid++) {
		//localSystem(fid) = 1-0.3725;
		localSystem(fid) = 0;
	}

	int id, tid, ntids, ipts, istart, iproc;

	//const int NUM_THREADS = omp_get_num_procs();
	const int NUM_THREADS = omp_get_num_procs()/2;
	omp_set_num_threads(NUM_THREADS);
#pragma omp parallel private(tid,ntids,ipts,istart,id)	
	{
		iproc = omp_get_num_procs();
		//iproc = 1; 
		tid = omp_get_thread_num();
		ntids = omp_get_num_threads();
		//ntids = 2; 
		ipts = (int)ceil(1.00*(double)Sample.size() / (double)ntids);
		//ipts = (int)ceil(1.00*(double)ntids / (double)ntids);
		istart = tid * ipts;
		if (tid == ntids - 1) ipts = Sample.size() - istart;
		if (ipts <= 0) ipts = 0;

		std::vector<Engine*> ep;
		ep.resize(ntids);
		printf("num threads=%d, iproc=%d, ID=%d, start=%d, to end=%d, num els=%d\n", ntids, iproc, tid, istart, istart + ipts, ipts);

		Eigen::VectorXd				D(F.rows());
		for (int i = 0; i < F.rows(); i++) {
			D(i) = numeric_limits<double>::infinity();
		}

		//cout << "[" << tid << "] Number of processors " << iproc << ", with " << ntids << " threads." << endl;

		UiTriplet[tid].reserve(Num_fields * ((double)ipts / (double)Sample.size()) * 2 * 10.0 * F.rows());

		// Computing the values of each element
		for (id = istart; id < (istart + ipts); id++) {
			//for (id = istart; id < (istart + ipts) && id < 10; id++) {
			if (id >= Sample.size()) break;

			vector<Eigen::Triplet<double>> BTriplet, C1Triplet, C2Triplet;
			///cout << "[" << id << "] Constructing local eigen problem\n ";

			///cout << "__creating subdomain \n";
			LocalFields localField(id);
			t1 = chrono::high_resolution_clock::now();
			//localField.constructSubdomain(Sample[id], V, F, avgEdgeLength, AdjMF3N, distRatio);
			//localField.constructSubdomain(Sample[id], V, F, avgEdgeLength, AdjMF2Ring, distRatio);
			localField.constructSubdomain(Sample[id], V, F, D, AdjMF2Ring, Sample.size(), this->numSupport);
			t2 = chrono::high_resolution_clock::now();
			durations[0] += t2 - t1;

			///cout << "__creating boundary \n";
			t1 = chrono::high_resolution_clock::now();
			localField.constructBoundary(F, AdjMF3N, AdjMF2Ring);
			t2 = chrono::high_resolution_clock::now();
			durations[1] += t2 - t1;

			///cout << "__gathering local elements \n";
			t1 = chrono::high_resolution_clock::now();
			localField.constructLocalElements(Num_fields, F);
			t2 = chrono::high_resolution_clock::now();
			durations[2] += t2 - t1;

			///cout << "__solving local eigenvalue problem \n";
			t1 = chrono::high_resolution_clock::now();
			//ep[tid] = engOpen(NULL);
			//printf("Starting engine %d for element %d\n", tid, id);
			//if (id % ((int)(Sample.size() / 4)) == 0)
			//	cout << "[" << id << "] Constructing local eigen problem\n ";
			//localField.constructLocalEigenProblemWithSelector(ep[tid], tid, Num_fields, SF, MF, AdjMF2Ring, 2, doubleArea, UiTriplet[id]);		// Matlab
			localField.constructLocalEigenProblemWithSelector(Num_fields, Sh, MF2DhNeg, AdjMF2Ring, 2, doubleArea, UiTriplet[id]);						// Spectra
			//localField.constructLocalEigenProblemWithSelectorRotEig(ep[tid], tid, SF2DAsym, MF2D, AdjMF2Ring, 2, doubleArea, UiTriplet[id]);		// 2nd basis: 90 rotation of the first basis
			//engClose(ep[tid]);
			//localField.constructLocalEigenProblem(SF2D, AdjMF3N, doubleArea, UiTriplet[id]);
			t2 = chrono::high_resolution_clock::now();
			durations[3] += t2 - t1;


			if (id == 0)
			{
				SubDomain = localField.SubDomain;
				Boundary = localField.Boundary;
			}
		}
	}
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;

	cout << "....Gathering local elements as basis matrix... ";
	t1 = chrono::high_resolution_clock::now();

	bool writeBasisCompsToFile = false;
	//if (writeBasisCompsToFile)
	//	writeBasisElementsToFile(UiTriplet, 2);
	//else
	gatherBasisElements(UiTriplet, 2);

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;


	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "..in Total of " << duration.count() << " seconds" << endl;


	// Information about the process timing
	printf("> Basis Timing information \n");
	printf("....[0] Constructing internal elements: %.8f seconds.\n", durations[0].count());
	printf("....[1] Constructing boundary: %.8f seconds.\n", durations[1].count());
	printf("....[2] Constructing local subdomains: %.8f seconds.\n", durations[2].count());
	printf("....==> Total SET UP: %.8f\n", durations[0].count() + durations[1].count() + durations[2].count());
	printf("....[3] Solving local eigenvalue problems: %.8f seconds.\n", durations[3].count());

	// Information about Basis
	printf("> Basis Structure information \n");
	printf("....Size = %dx%d\n", Basis.rows(), Basis.cols());
	printf("....NNZ=%d, per row = %.4f\n", Basis.nonZeros(), (double)Basis.nonZeros() / (double)Basis.rows());
}

void NRoSyFields::gatherBasisElements(const vector<vector<Eigen::Triplet<double>>> &UiTriplet, const int& NUM_EIGEN)
{
	/* Set the basis sum to be zero on each column-group */
	vector<Eigen::Triplet<double>> BTriplet;

	/* Getting the number of total elements on the triplet*/
	int totalElements = 0;
	for (int j = 0; j < Sample.size(); j++) {
		totalElements += UiTriplet[j].size();
	}

	/* Constructing the basis matrix from the triplet */
	BTriplet.resize(totalElements);
	for (int j = 0; j < Sample.size(); j++) {
		int tripSize = 0;
		for (int k = 0; k < j; k++) {
			tripSize += UiTriplet[k].size();
		}
		//std::copy(UiTriplet[j].begin(), UiTriplet[j].end(), BTriplet.begin() + tripSize);
		std::move(UiTriplet[j].begin(), UiTriplet[j].end(), BTriplet.begin() + tripSize);
	}
	Basis.setFromTriplets(BTriplet.begin(), BTriplet.end());
}

/* ============================= REDUCED N-FIELDS DESIGN ============================= */
void NRoSyFields::nRoSyFieldsDesign_Reduced()
{
	//nRoSyFieldsDesign_Reduced_HardConstraints();
	nRoSyFieldsDesign_Reduced_SoftConstraints();
}

void NRoSyFields::nRoSyFieldsDesign_Reduced_HardConstraints()
{
	Eigen::SparseMatrix<double>		B2FBar = Basis.transpose()*(SF*MFinv*SF)*Basis;
	Eigen::VectorXd					bBar, gBar, hBar, vEstBar;
	Eigen::SparseMatrix<double>		A_LHSBar;

	constructRandomHardConstraints_Reduced();
	setupRHSBiharmSystem_Reduced(B2FBar, gBar, hBar, vEstBar, bBar);
	setupLHSBiharmSystem_Reduced(B2FBar, A_LHSBar);
	solveBiharmSystem_Reduced(vEstBar, A_LHSBar, bBar);
}

void NRoSyFields::constructRandomHardConstraints_Reduced()
{
	//userConstraints = globalConstraints; 
	CBar = C * Basis;
	cBar = c;	
}

void NRoSyFields::setupRHSBiharmSystem_Reduced(const Eigen::SparseMatrix<double>& B2FBar, Eigen::VectorXd& gBar, Eigen::VectorXd& hBar, Eigen::VectorXd& vEstBar, Eigen::VectorXd& bBar)
{
	vEstBar.resize(B2FBar.rows());
	for (int i = 0; i < vEstBar.rows(); i++) {
		vEstBar(i) = 0.5;
	}

	// Create alignment fields based on maximal curvature direction
	Eigen::VectorXd repV, repVbar;
	createAlignmentField(repV);
	repVbar = Basis*(MF*repV);

	gBar = B2FBar * vEstBar + repVbar;
	bBar.resize(B2FBar.rows() + cBar.rows());

	// Constructing b
	hBar = CBar * vEstBar - cBar;
	bBar << gBar, hBar;
}

void NRoSyFields::setupLHSBiharmSystem_Reduced(const Eigen::SparseMatrix<double>& B2FBar, Eigen::SparseMatrix<double>& A_LHSBar)
{
	A_LHSBar.resize(B2FBar.rows() + CBar.rows(), B2FBar.cols() + CBar.rows());
	vector<Eigen::Triplet<double>>	ATriplet;
	ATriplet.reserve(B2FBar.nonZeros() + 2 * CBar.nonZeros());

	for (int k = 0; k < B2FBar.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(B2FBar, k); it; ++it) {
			ATriplet.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
		}
	}

	for (int k = 0; k < CBar.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(CBar, k); it; ++it) {
			ATriplet.push_back(Eigen::Triplet<double>(B2FBar.rows() + it.row(), it.col(), it.value()));
			ATriplet.push_back(Eigen::Triplet<double>(it.col(), B2FBar.cols() + it.row(), it.value()));
		}
	}
	A_LHSBar.setFromTriplets(ATriplet.begin(), ATriplet.end());
}

void NRoSyFields::solveBiharmSystem_Reduced(const Eigen::VectorXd& vEstBar, const Eigen::SparseMatrix<double>& A_LHSBar, const Eigen::VectorXd& bBar)
{
	Eigen::VectorXd XLowDim(vEstBar.size());
	XfBar.resize(Basis.rows());
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> sparseSolver(A_LHSBar);

	cout << "....Solving for the first frame.\n";
	Eigen::VectorXd x = sparseSolver.solve(bBar);

	if (sparseSolver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system. " << endl;
		if (sparseSolver.info() == Eigen::NumericalIssue)
			cout << "NUMERICAL ISSUE. " << endl;
		cout << sparseSolver.info() << endl;
		return;
	}

	XLowDim = -x.block(0, 0, vEstBar.size(), 1) + vEstBar;

	XfBar = Basis * XLowDim;
}
void NRoSyFields::nRoSyFieldsDesign_Reduced_Splines()
{
	//Eigen::SparseMatrix<double>		B2FBar = Basis.transpose()*(SF*MFinv*SF)*Basis;
	Eigen::VectorXd					bBar, gBar, hBar, vEstBar;
	Eigen::SparseMatrix<double>		A_LHSBar;

	double weight = 0;
	for (int i = 0; i < userVisualConstraints.size(); i += 2)
	{
		weight += doubleArea(userVisualConstraints[i]);
	}

	vector<double> lambda(2);
	lambda[0] = 1.0;
	//lambda[1] = 0.001 / weight;
	lambda[1] = 0.000005; 
	
	cout << "Getting the constraints \n";
	getReducedConstraints();
	cout << "Getting the RHS \n";
	setupRHSSplines_Reduced(BFBar, lambda, gBar, hBar, bBar);
	cout << "Getting the LHS \n";
	setupLHSSplines_Reduced(BFBar, lambda, A_LHSBar);
	cout << "solving reduced system \n";
	solveSplines_Reduced(A_LHSBar, bBar);
}
void NRoSyFields::getReducedConstraints()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Obtaining user constraints ";

	//constructConstraints();

	//userConstraints = globalConstraints; 
	CBar = C * Basis;
	cBar = c;

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;

	printf(".... C_LoCal = %dx%d\n", CBar.rows(), CBar.cols());
	printf(".... c_LoCal = %dx%d\n", cBar.rows(), cBar.cols());
}
void NRoSyFields::setupRHSSplines_Reduced(const Eigen::SparseMatrix<double>& B2FBar, const vector<double>& lambda, Eigen::VectorXd& gBar, Eigen::VectorXd& hBar, Eigen::VectorXd& bBar)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Constructing RHS (mapped)...";
	
	// ================
	bBar.resize(BFBar.rows() + cBar.rows());
	cout << "__setting gbar:";
		t2 = chrono::high_resolution_clock::now();
		duration = t2 - t0;
		cout << "in " << duration.count() << " seconds." << endl;
	gBar = lambda[1] * Basis.transpose()*(MF*alignFields);
	cout << "__setting hbar:";
	hBar = cBar;
		t2 = chrono::high_resolution_clock::now();
		duration = t2 - t0;
		cout << "in " << duration.count() << " seconds." << endl;
	cout << "__setting bbar\n";
	bBar << gBar, hBar;
		t2 = chrono::high_resolution_clock::now();
		duration = t2 - t0;
		cout << "in " << duration.count() << " seconds." << endl;

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;

}
void NRoSyFields::setupLHSSplines_Reduced(const Eigen::SparseMatrix<double>& B2FBar, const vector<double>& lambda, Eigen::SparseMatrix<double>& A_LHSBar)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Setting up LHS: ";
		
	A_LHSBar.resize(BFBar.rows() + CBar.rows(), BFBar.cols() + CBar.rows());
	vector<Eigen::Triplet<double>>	ATriplet;
	ATriplet.reserve(BFBar.nonZeros() + 2 * CBar.nonZeros());

	//Eigen::SparseMatrix<double> BM = lambda[0] * B2FBar + lambda[1] * Basis.transpose()*MF*Basis;

	BMBar = BFBar + Basis.transpose()*(lambda[1] * MF)*Basis;

	cout << "Iterate over BM\n";
	for (int k = 0; k < BMBar.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(BMBar, k); it; ++it) {
			ATriplet.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
		}
	}

	cout << "Iterate over constraints \n";
	for (int k = 0; k < CBar.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(CBar, k); it; ++it) {
			ATriplet.push_back(Eigen::Triplet<double>(BMBar.rows() + it.row(), it.col(), it.value()));
			ATriplet.push_back(Eigen::Triplet<double>(it.col(), BMBar.cols() + it.row(), it.value()));
		}
	}

	cout << "Set from constraints. \n";
	A_LHSBar.setFromTriplets(ATriplet.begin(), ATriplet.end());

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;
}
void NRoSyFields::solveSplines_Reduced(const Eigen::SparseMatrix<double>& A_LHSBar, const Eigen::VectorXd& bBar)
{
	Eigen::VectorXd XLowDim(A_LHSBar.rows());
	XfBar.resize(Basis.rows());
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> sparseSolver(A_LHSBar);

	cout << "....Solving for the first frame.\n";
	Eigen::VectorXd x = sparseSolver.solve(bBar);

	if (sparseSolver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system. " << endl;
		if (sparseSolver.info() == Eigen::NumericalIssue)
			cout << "NUMERICAL ISSUE. " << endl;
		cout << sparseSolver.info() << endl;
		return;
	}

	XLowDim = x.block(0, 0, A_LHSBar.rows(), 1);

	XfBar = Basis * XLowDim;
}

// SOFT Constraints
void NRoSyFields::nRoSyFieldsDesign_Reduced_SoftConstraints()
{
	Eigen::VectorXd					bBar, gBar, hBar, vEstBar;
	Eigen::SparseMatrix<double>		A_LHSBar;

	Eigen::Vector3d lambda;
	lambda(0) = 1000 * MF.diagonal().sum() / SF.diagonal().sum();
	lambda(1) = 0;
	lambda(2) = 1.0;

	constructSoftConstraints_Reduced();
	setupRHSSoftConstraints_Reduced(lambda, bBar);
	setupLHSSoftConstraints_Reduced(lambda, A_LHSBar);
	solveUserSystemMappedLDLTSoftConstraints(A_LHSBar, bBar);
}

void NRoSyFields::constructSoftConstraints_Reduced()
{
	//userConstraints = globalConstraints; 
	CBar = C * Basis;
	cBar = c;
}

void NRoSyFields::setupRHSSoftConstraints_Reduced(const Eigen::Vector3d& lambda, Eigen::VectorXd& bBar)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Constructing RHS (mapped)...";

	/* Mass matrix of the selected faces (on constraints) */
	Eigen::SparseMatrix<double> Mconst = C*MF*C.transpose();

	//printf("Siz of Mconst: %dx%d\n", Mconst.rows(), Mconst.cols());
	//printf("Siz of Cbar: %dx%d\n", CBar.rows(), CBar.cols());

	//bBar = lambda(2)*CBar.transpose() * cBar; 
	bBar = lambda(2)*CBar.transpose() * Mconst * cBar;

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;
}

void NRoSyFields::setupLHSSoftConstraints_Reduced(const Eigen::Vector3d& lambda, Eigen::SparseMatrix<double>& A_LHSBar)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Constructing LHS (mapped)...";

	Eigen::SparseMatrix<double> SF2DBar = Basis.transpose() * SF* Basis;
	//Eigen::SparseMatrix<double> B2DBar = Basis.transpose() * B2D * Basis;
	//A_LHSBar = SF2DBar + lambda*CBar.transpose()*CBar; 


	//const double lambda_1 = 10000 / B2D.coeff(0, 0);
	cout << "lambda_2 " << lambda(2) << endl;

	/* Local matrix */
	Eigen::SparseMatrix<double> Mconst = C*MF*C.transpose();

	//printf("Siz of Mconst: %dx%d\n", Mconst.rows(), Mconst.cols());
	//printf("Siz of Cbar: %dx%d\n", CBar.rows(), CBar.cols());

	//A_LHSBar = lambda(0)*SF2DBar +  lambda(1)*B2DBar + lambda(2)*CBar.transpose()*CBar;
	A_LHSBar = lambda(0)*SF2DBar + lambda(2)*CBar.transpose()*Mconst*CBar;
	//A_LHSBar = lambda(0)*SF2DBar + lambda(2)*CBar.transpose()*CBar;
	//A_LHSBar = lambda(0)*SF2DBar + lambda_1*B2DBar + lambda(2)*CBar.transpose()*CBar;

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;
}

void NRoSyFields::solveUserSystemMappedLDLTSoftConstraints(Eigen::SparseMatrix<double>& A_LHSBar, Eigen::VectorXd& bBar)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Solving reduced system...\n";

	XfBar.resize(Basis.rows());
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> sparseSolver(A_LHSBar);

	cout << "....Solving for the first frame.\n";
	Eigen::VectorXd XLowDim = sparseSolver.solve(bBar);

	if (sparseSolver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system. " << endl;
		if (sparseSolver.info() == Eigen::NumericalIssue)
			cout << "NUMERICAL ISSUE. " << endl;
		cout << sparseSolver.info() << endl;
		return;
	}	

	XfBar = Basis * XLowDim;

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "..in Total of " << duration.count() << " seconds." << endl;
}

void NRoSyFields::measureAccuracy()
{
	if (Xf.size() < 1 || XfBar.size() < 1)
	{
		cout << "Either the reduced and the reference fields are not computed yet. Return \n";
		return;
	}

	cout << "____Computing L2-norm \n";
	Eigen::VectorXd diff = Xf - XfBar;
	double norm1 = diff.transpose()*MF*diff;
	double norm2 = Xf.transpose()*MF*Xf;
	double normL2 = sqrt(norm1 / norm2);
	double error = normL2;
	cout << "____The L2 Norm is << " << normL2 << ": sqrt(" << norm1 << "/" << norm2 << ")" << endl;

	cout << "____Harmonic energy \n";
	double harm_energy1_sym = Xf.transpose()*SF*Xf;
	double harm_energy2_sym = XfBar.transpose()*SF*XfBar;
	double harm_relEnergy_sym = abs(harm_energy1_sym - harm_energy2_sym) / harm_energy1_sym;
	cout << "____Harmonic Energy => Ref=" << harm_energy1_sym << ", Approx:" << harm_energy2_sym << endl;
	cout << "____Relative energy: " << harm_relEnergy_sym << endl;
}



/* ============================= Testing stuff ============================= */
void NRoSyFields::TEST_NROSY(igl::opengl::glfw::Viewer &viewer, const string& meshFile)
{
	nRot = 2;
	readMesh(meshFile);
	scaleMesh();
	igl::doublearea(V, F, doubleArea);
	string model = "Fertility_";
	NRoSy nRoSy_;

	viewer.data().set_mesh(V, F);
	viewer.append_mesh();
	viewer.data().set_mesh(V, F);
	viewer.data().show_lines = false;
	viewer.selected_data_index = 0;
	
	computeAverageEdgeLength();
	computeFaceCenter();
	constructFaceAdjacency3NMatrix();
	constructFaceAdjacency2RingMatrix();
	constructEVList();
	constructEFList();
	constructFrameBasis();
	constructMappingMatrix();

	constructMassMatrixMF3D();
	computeFrameRotation(viewer);
	buildStiffnessMatrix_Geometric();
	//buildStiffnessMatrix_Combinatorial();

	//selectFaceToDraw(20000);
	selectFaceToDraw(F.rows()/3.0);
	Eigen::VectorXd inputNFields;
	string fieldsfile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/"+model+"_InputFields";
	//ReadVectorFromMatlab(inputNFields,fieldsfile, 2*F.rows());
	//visualize2Dfields(viewer, inputNFields, Eigen::RowVector3d(0.1, 0.2, 0.2), 1.0);
	//createNRoSyFromVectors(inputNFields, nRoSy);
	//visualizeNRoSyFields(viewer, nRoSy );
	//convertNRoSyToRepVectors(nRoSy, repVector);
	//visualizeRepVectorFields(viewer, repVector);
	//convertRepVectorsToNRoSy(repVector, nRoSy);
	//visualizeNRoSyFields(viewer, nRoSy);

	/* Visualizing + Manipulating n-RoSy fields */
	//convertNRoSyToRepVectors(nRoSy, repVector);
	//visualizeRepVectorFields(viewer);
	//convertRepVectorsToNRoSy(repVector, nRoSy);
	//visualizeNRoSyFields(viewer, nRoSy);

	/* Working with eigenvectors of n-RoSy fields*/

	string fileEigen = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/" + model + to_string(nRot) + "-fields_25_Ref";
	//computeEigenFields_generalized(100, fieldsfile);
	//computeEigenFields_regular(50, fileEigen);
	NRoSy nRoSy_eigenFields;
	//convertRepVectorsToNRoSy(eigFieldsNRoSyRef.col(0), nRoSy_eigenFields);
	//visualizeNRoSyFields(viewer, nRoSy_eigenFields, Eigen::RowVector3d(0.0, 0.1, 0.9));
	////visualizeRepVectorFields(viewer, eigFieldsNRoSyRef.col(0));

	//Xf = eigFieldsNRoSyRef.col(0);

	

	/* Build reduced space */
	numSupport = 40.0;
	numSample = 5000;
	constructSamples(numSample);
	string basisFile = "D:/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_" + model + to_string(nRot) + "-fields_" + to_string(numSample*2) + "_Eigfields_"+ to_string((int)numSupport) + "sup";
	//constructBasis();
	//storeBasis(basisFile);
	retrieveBasis(basisFile);
	//visualizeBasis(viewer, 0);

	/* Projection */
	//Eigen::SparseMatrix<double>							B_NR = Basis.transpose() * MF * Basis;
	//Eigen::PardisoLDLT<Eigen::SparseMatrix<double>>		sparseSolver(B_NR);
	//Eigen::VectorXd inputFields = Basis.transpose()*MF*Xf;
	//double error;
	//testProjection_MyBasis_NoRegularizer(Basis, sparseSolver, MF, inputFields, Xf, error);

	/* Test on smoothing of the n-RoSy fields */
	//Eigen::VectorXd PD, PV;
	//Eigen::VectorXd inputF, outputF;
	//computeMaximalPrincipalCurvature(V, F, PD, PV);
	//createNRoSyFromVectors(PD, nRoSy_);
	//nRoSy_.magnitude = PV; 
	//convertNRoSyToRepVectors(nRoSy_, inputF);
	//smoothNRoSyFields(100, inputF, outputF);
	//Xf = inputF;
	//XfBar = outputF;

	/* Prepare for n-fields design */
	double weight = 0;
	for (int i = 0; i < userVisualConstraints.size(); i += 2){
		weight += doubleArea(userVisualConstraints[i]);
	}
	
	vector<double> lambda(2);
	lambda[0] = 1.0;
	if (userVisualConstraints.size() < 1)	lambda[1] = 0.0;
	else									lambda[1] = 0.001 / weight;
	
	createAlignmentField(alignFields);
	BF = SF*MFinv*SF;
	BM = lambda[0] * BF;
	BFBar = Basis.transpose()*BF*Basis;
	BMBar = BFBar;
	
	convertRepVectorsToNRoSy(alignFields, nRoSy_);
	XfBar = alignFields;
	visualizeNRoSyFields(viewer, nRoSy_, Eigen::RowVector3d(0.0, 0.8, 0.1));

	/* =========== N-ROSY FIELDS DESIGN ===============*/
	/* Constrained fields (biharmonic) */
	//nRoSyFieldsDesignRef();
	//visualizeConstraints(viewer);
	//visualizeConstrainedFields(viewer);
	
	///* Reduced Constrained fields (biharmonic)--hard constraints */
	//constructRandomHardConstraints(C, c);
	//nRoSyFieldsDesign_Reduced();
	//visualizeConstrainedFields_Reduced(viewer);
	//visualizeConstraints(viewer);
	//measureAccuracy();

	/* Constrained fields (SOFT constraints) */
	//constructSoftConstraints();
	//nRoSyFieldsDesignRef();
	//visualizeSoftConstraints(viewer);
	//visualizeConstrainedFields(viewer);

	
	//nRoSyFieldsDesign_Reduced();
	//measureAccuracy();
	//
	//NRoSy nRoSy_temp; 
	//convertRepVectorsToNRoSy(Xf, nRoSy_temp);
	//string fileNRoSyRef= "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/" + model + to_string(nRot) + "_reference.txt";
	//writeNRoSyFieldsToFile(nRoSy_temp, fileNRoSyRef);
	//
	//convertRepVectorsToNRoSy(XfBar, nRoSy_temp);
	//string fileNRoSyRed = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/" + model + to_string(nRot) + "_approximation.txt";
	//writeNRoSyFieldsToFile(nRoSy_temp, fileNRoSyRed);
}

void NRoSyFields::writeNRoSyFieldsToFile(const NRoSy& nRoSy, const string& filename)
{

	double scale = 1.0;
	Eigen::Vector2d b(1, 0);

	vector<Eigen::VectorXd> nFields(nRot);

	const int FIELD_TO_WRITE = 1;
	for (int i = 0; i < nRot; i++)
	{
		nFields[i].resize(2 * F.rows()); // = Eigen::VectorXd(2 * F.rows());
		cout << "Drawing the " << i << " fields \n";

		/* Construct rotation matrix*/
		//for (int j = 0; j < F.rows(); j++)
		for (int j = 0; j<F.rows(); j++)
		{
			double angle = nRoSy.theta(j) + ((double)i*2.0*M_PI / (double)nRot);
			if (j == 0)
			{
				printf("angle 0=%.5f, theta 0=%.5f\n", angle, nRoSy.theta(j));
			}
			Eigen::Matrix2d RotM;
			RotM(0, 0) = cos(angle);
			RotM(0, 1) = -sin(angle);
			RotM(1, 0) = sin(angle);
			RotM(1, 1) = cos(angle);

			nFields[i].block(2 * j, 0, 2, 1) = nRoSy.magnitude(j) * RotM * b;
			//nFields[i].block(2 * j, 0, 2, 1) = 1.0 * RotM * b;
		}		
	}

	vector<Eigen::VectorXd> nFields3d(nRot);
	for (int j = 0; j < nRot; j++)
	{
		nFields3d[j] = A*nFields[j];
	}

	ofstream myfile(filename.c_str());
	if (myfile.is_open())
	{
		cout << "Writing vector fields to file to text \n";
		printf("__|F|=%d  | vfields=%d\n", F.rows(), nFields3d[0].size());
		for (int i = 0; i < F.rows(); i++)
		{
			for (int j = 0; j < FIELD_TO_WRITE; j++)
			{
				myfile << nFields3d[j](3 * i) << " " << nFields3d[j](3 * i + 1) << " " << nFields3d[j](3 * i + 2) << "\n";
			}			
		}
		myfile.close();
	}
	else cout << "Unable to open file";

	cout << "Writing to file doen!!!!\n";
}

void NRoSyFields::writeConstraintsToFile(const string& filename)
{
	ofstream myfile(filename.c_str());

	Eigen::MatrixXd ALoc(3, 2);
	NRoSy nRoSy_;
	convertRepVectorsToNRoSy(c, nRoSy_);
	Eigen::VectorXd c2(c.size());
	Eigen::Vector2d id, v2;
	Eigen::Vector3d v3;

	if (myfile.is_open())
	{
		cout << "Writing constraints to file \n";
		//for (int j = 0; j < nRot; j++)
		for (int j = 0; j < 1; j++)
		{
			for (int i = 0; i < globalConstraints.size(); i++)
			{
				ALoc = A.block(3 * globalConstraints[i], 2 * globalConstraints[i], 3, 2);

				id << 1.0, 0.0;
				double angle = nRoSy_.theta(i) + ((double)j*2.0*M_PI / (double)nRot);
				double cosA = cos(angle);
				double sinA = sin(angle);
				Eigen::Matrix2d RotA; RotA << cosA, -sinA, sinA, cosA;
				//c2.block(2 * i, 0, 2, 1) = RotA*id;
				v2 = RotA*id;
				v3 = ALoc * v2;
				//printf("v2=[%.3f, %.3f] and v3=[%.3f, %.3f, %.3f] || Rot = [%.3f, %.3f, %.3f, %.3f. \n", v2(0), v2(1), v3(0), v3(1), v3(2), cosA, -sinA, sinA, cosA);
				cout << "ALoc: " << ALoc << endl; 

				myfile << globalConstraints[i] << " " << v3(0) << " " << v3(1) << " " << v3(2) << "\n";
			}
		}
		myfile.close();
	} else cout << "Unable to open file";

	cout << "Writing to file doen!!!!\n";
}

/* PROJECTION ON REDUCED FIELDS */
void NRoSyFields::testProjection_MyBasis_NoRegularizer(const Eigen::SparseMatrix<double>& Basis, const Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> &sparseSolver, const Eigen::SparseMatrix<double>& B, const Eigen::VectorXd& a, const Eigen::VectorXd& v, double &error)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> [Testing basis 2D of arbitrary field without regularizer...]\n";	

	Eigen::VectorXd w = sparseSolver.solve(a);

	if (sparseSolver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system. " << endl;
		if (sparseSolver.info() == Eigen::NumericalIssue)
			cout << "NUMERICAL ISSUE. " << endl;
		cout << sparseSolver.info() << endl;
		return;
	}	
	wb = Basis*w;

	// Compare their L2-Norm
	cout << "Computing L2-norm \n";
	Eigen::VectorXd diff = v - wb;
	double length1 = wb.transpose()*MF*wb;
	double norm1 = diff.transpose()*MF*diff;
	double norm2 = v.transpose()*MF*v;
	double normL2 = sqrt(norm1 / norm2);
	error = normL2;

	cout << "____The L2 Norm is << " << normL2 << ": sqrt(" << norm1 << "/" << norm2 << ")" << endl;
	
	/* Measuring the 'length' of each vector */
	cout << "ENERGY SYM \n";
	double harm_energy1_sym = v.transpose()*SF*v;
	double harm_energy2_sym = wb.transpose()*SF*wb;
	double harm_relEnergy_sym = abs(harm_energy1_sym - harm_energy2_sym) / harm_energy1_sym;
	cout << "____Harmonic Energy => Ref=" << harm_energy1_sym << ", Approx:" << harm_energy2_sym << endl;
	cout << "____Relative energy: " << harm_relEnergy_sym << endl;

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;

	bool writeToFile = false;
	if (writeToFile) {
		std::ofstream ofs;
		string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Arma43k_L2projection_eigenFields_" + to_string(nRot) + "_" + to_string(Basis.cols()) + "_40sup.txt";

		ofs.open(resultFile, std::ofstream::out | std::ofstream::app);

		ofs << globalConstraints.size() << "\t"
			<< harm_energy1_sym << "\t" << "\t" << norm2 << "\t"	// reference (harmAsym, harmSym, biharmAsym, biharmSym)
			<< harm_energy2_sym << "\t" << harm_relEnergy_sym << "\t"						// approx (harmAsym, harmRelEnergyAsym,harmSym, harmRelEnergySym)
			<< length1 << "\t" << normL2 << "\n";																												// l2-norm

		ofs.close();
	}
}
