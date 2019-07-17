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
	vector<Eigen::Triplet<double>> MFTriplet;
	MFTriplet.reserve(2 * F.rows());
	vector<Eigen::Triplet<double>> MFInvTriplet;
	MFInvTriplet.reserve(2 * F.rows());

	for (int i = 0; i < F.rows(); i++) {
		double area = doubleArea(i) / 2.0;
		MFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, area));
		MFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, area));
		MFInvTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, 1.0 / area));
		MFInvTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, 1.0 / area));
	}
	MF.setFromTriplets(MFTriplet.begin(), MFTriplet.end());
	MFinv.setFromTriplets(MFInvTriplet.begin(), MFInvTriplet.end());

}

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
		double cosRA = cos(FrameRot(ei, 0));
		double sinRA = sin(FrameRot(ei, 0));
		double cosSB = cos(FrameRot(ei, 1));
		double sinSB = sin(FrameRot(ei, 1));

		Eigen::Matrix2d R1; R1 << cosRA, -sinRA, sinRA, cosRA;
		Eigen::Matrix2d R2; R2 << cosSB, -sinSB, sinSB, cosSB;
		
		/* The transport matrix */
		Eigen::Matrix2d B2toB1 = R1*R2.transpose();
		Eigen::Matrix2d B1toB2 = R2*R1.transpose();

		if (ei == ee)
		{
			cout << "Data of : " << ei << endl;
			cout << "B2 to B1 : \n" << B2toB1 << endl;
			cout << "B1 to B2 : \n" << B1toB2 << endl;
		}

		/* (Geometric) Laplace matrix from A (first triangle) perspective */
		for (int i = 0; i < 2; i++)
		{
			STriplet.push_back(Eigen::Triplet<double>(2 * TA + i, 2 * TA + i, 1.0/3.0));
			for (int j = 0; j < 2; j++)
			{
				STriplet.push_back(Eigen::Triplet<double>(2 * TA + i, 2 * TB + j, B2toB1(i, j) / 3.0));
			}
		}

		/* (Geometric) Laplace matrix from B (first triangle) perspective */
		for (int i = 0; i < 2; i++)
		{
			STriplet.push_back(Eigen::Triplet<double>(2 * TB + i, 2 * TB + i, 1.0 / 3.0));
			for (int j = 0; j < 2; j++)
			{
				STriplet.push_back(Eigen::Triplet<double>(2 * TB + i, 2 * TA + j, B1toB2(i, j) / 3.0));
			}
		}

	}
	SF.setFromTriplets(STriplet.begin(), STriplet.end());
	string filename = "D:/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/Arma_Lap_NRoSy_Comb";
	WriteSparseMatrixToMatlab(SF, filename);
}

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

		if (ei == ee)
		{
			cout << "TA: " << TA << ", and TB: " << TB << endl; 
		}

		/* Construct the rotation matrix RA and SB */
		//cout << "Construct the rotation matrix RA and SB\n";
		double cosRA = cos(FrameRot(ei, 0));
		double sinRA = sin(FrameRot(ei, 0));
		double cosSB = cos(FrameRot(ei, 1));
		double sinSB = sin(FrameRot(ei, 1));

		Eigen::Matrix2d R1; R1 << cosRA, -sinRA, sinRA, cosRA;
		Eigen::Matrix2d R2; R2 << cosSB, -sinSB, sinSB, cosSB;

		//cout << "Copmute the weight \n";
		/* Compute the weight on each edge */
		Eigen::Vector3d edge_ = V.row(E(ei, 1)) - V.row(E(ei, 0));
		double el = edge_.dot(edge_);
		/* Should be multiplied by 3.0
		** but because later will be divided by 3.0 again (due to 3 neighbors),
		** I just dont multiply with anything */
		//weight = el * el / (0.5*doubleArea(TA) + 0.5*doubleArea(TB));
		weight = (3.0 * el * el) / (0.5*doubleArea(TA) + 0.5*doubleArea(TB));

		/* The transport matrix */
		Eigen::Matrix2d B2toB1 = R1*R2.transpose();
		Eigen::Matrix2d B1toB2 = R2*R1.transpose();

		if (ei == ee)
		{
			cout << "Data of : " << ei << endl;
			cout << "B2 to B1 : \n" << B2toB1 << endl;
			cout << "B1 to B2 : \n" << B1toB2 << endl; 
		}

		/* (Geometric) Laplace matrix from A (first triangle) perspective */
		for (int i = 0; i < 2; i++)
		{
			STriplet.push_back(Eigen::Triplet<double>(2 * TA + i, 2 * TA + i, weight));
			for (int j = 0; j < 2; j++)
			{
				STriplet.push_back(Eigen::Triplet<double>(2 * TA + i, 2 * TB + j, weight*B2toB1(i, j)));
			}
		}

		/* (Geometric) Laplace matrix from B (first triangle) perspective */
		for (int i = 0; i < 2; i++)
		{
			STriplet.push_back(Eigen::Triplet<double>(2 * TB + i, 2 * TB + i, weight));
			for (int j = 0; j < 2; j++)
			{
				STriplet.push_back(Eigen::Triplet<double>(2 * TB + i, 2 * TA + j, weight*B1toB2(i, j)));
			}
		}
	}
	SF.setFromTriplets(STriplet.begin(), STriplet.end());
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
			if (acos(dp_1) < 0.0)
			{
				angle_1 = M_PI + acos(dp_1);
			}
			else {
				angle_1 = acos(dp_1);
			}
			
			FrameRot(ei, 0) = M_PI - angle_1;
			break;
		case 2:
			//dp_1 = e_ij.dot(e_i) / (e_i.norm() * e_ij.norm());
			dp_1 = (e_i).dot(-e_ij) / (e_i.norm() * e_ij.norm());
			//angle_1 =acos(dp_1);
			if (acos(dp_1) < 0.0)
			{
				angle_1 = M_PI + acos(dp_1);
			}
			else {
				angle_1 = acos(dp_1);
			}
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
			//angle_2 =acos(dp_2);
			if (acos(dp_2) < 0.0)
			{
				angle_2 = M_PI + acos(dp_2);
			}
			else {
				angle_2 = acos(dp_2);
			}
			FrameRot(ei, 1) =  M_PI - angle_2;
			break;
		case 2:
			//dp_2 = e_ij.dot(e_j) / (e_j.norm() * e_ij.norm());
			dp_2 = (e_j).dot(-e_ji) / (e_j.norm() * e_ji.norm());
			//angle_2 = acos(dp_2);
			if (acos(dp_2) < 0.0)
			{
				angle_2 = M_PI + acos(dp_2);
			}
			else {
				angle_2 = acos(dp_2);
			}
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

/* Rep. Vectors and N-RoSy Fields interface */
void NRoSyFields::convertNRoSyToRepVectors(const NRoSy& nRoSyFields, Eigen::VectorXd& repVect)
{
	double scale = 1.0;
	Eigen::Vector2d b(1, 0);	
	
	repVect.resize(2 * F.rows());

	/* Construct rotation matrix*/
	for (int j = 0; j < F.rows(); j++)
	{
		double angle = nRoSyFields.theta(j) + (2.0*M_PI / (double)nRot);
		
		Eigen::Matrix2d RotM;
		RotM(0, 0) = cos(angle);
		RotM(0, 1) = -sin(angle);
		RotM(1, 0) = sin(angle);
		RotM(1, 1) = cos(angle);

		repVector.block(2 * j, 0, 2, 1) = nRoSyFields.magnitude(j) * RotM * b;
	}
}

void NRoSyFields::convertRepVectorsToNRoSy(const Eigen::VectorXd& repVect, NRoSy& nRoSyFields)
{
	nRoSyFields.theta.resize(F.rows());
	nRoSyFields.magnitude.resize(F.rows());

	/* Temp variables */
	Eigen::Vector2d v, b(1,0);

	for (int i = 0; i < F.rows(); i++)
	{
		v = repVect.block(2 * i, 0, 2, 1);

		nRoSyFields.magnitude(i) = v.norm();
		v.normalize();
		double angle;
		if(v(1)<0)
			nRoSyFields.theta(i) = M_PI - acos(b.dot(v));
		else 
			nRoSyFields.theta(i) = acos(b.dot(v));
	}
}

void NRoSyFields::createNRoSyFromVectors(const Eigen::VectorXd& vectorFields)
{
	cout << "Converting to nRoSy fields \n";
	this->nRoSy.theta.resize(F.rows());
	this->nRoSy.magnitude.resize(F.rows());

	/* Temp variables */
	Eigen::Vector2d v, b(1, 0);

	for (int i = 0; i < F.rows(); i++)
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

/* Visualizing the NRoSyFields */
void NRoSyFields::visualizeNRoSyFields(igl::opengl::glfw::Viewer &viewer, const NRoSy& nRoSyFields)
{
	//double scale = 250.0;
	double scale = 1.0;
	Eigen::Vector2d b(1, 0);
	Eigen::RowVector3d color(0.1, 0.1, 0.9);
	//Eigen::VectorXd col(nRot);
	//Eigen::MatrixXd color;
	//for (int i = 0; i < nRot; i++)
	//{
	//	col(i) = i*(1.0 / nRot);
	//}
	//igl::jet(col, true, color);

	for (int i = 0; i < nRot; i++)
	{
		Eigen::VectorXd TempFields(2 * F.rows());
		cout << "Drawing the " << i << " fields \n";

		/* Construct rotation matrix*/
		for (int j = 0; j < F.rows(); j++)
		{
			double angle = nRoSyFields.theta(j) + ((double)i*2.0*M_PI / (double)nRot);
			if (j == 0)
			{
				printf("angle 0=%.5f, theta 0=%.5f\n", angle, nRoSyFields.theta(j));
			}
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

void NRoSyFields::visualizeRepVectorFields(igl::opengl::glfw::Viewer &viewer, const NRoSy& nRoSyFields)
{
	const double scale = 1.0;
	Eigen::Vector2d b(1, 0);	
	Eigen::RowVector3d color (117.0/255.0, 107.0/255.0, 177.0/255.0);
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

void NRoSyFields::visualizeRepVectorFields(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd& repVector)
{
	viewer.data().clear();
	viewer.data().set_mesh(V, F);

	//const double scale = 250.0;
	const double scale = 1.0;
	Eigen::RowVector3d color(117.0 / 255.0, 107.0 / 255.0, 177.0 / 255.0);
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

	NRoSy nRoSy_eigenFields;
	convertRepVectorsToNRoSy(eigFieldsNRoSyRef.col(id), nRoSy_eigenFields);
	visualizeNRoSyFields(viewer, nRoSy_eigenFields);
}

/* ============================= Testing stuff ============================= */
void NRoSyFields::TEST_NROSY(igl::opengl::glfw::Viewer &viewer, const string& meshFile)
{
	nRot = 4;
	readMesh(meshFile);
	scaleMesh();
	
	computeAverageEdgeLength();
	computeFaceCenter();
	constructFaceAdjacency3NMatrix();
	constructEVList();
	constructEFList();
	constructFrameBasis();
	constructMappingMatrix();
	constructMassMatrixMF3D();

	selectFaceToDraw(2500);
	Eigen::VectorXd inputNFields;
	string fieldsfile = "D:/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/Arma_InputFields";
	ReadVectorFromMatlab(inputNFields,fieldsfile, 2*F.rows());
	createNRoSyFromVectors(inputNFields);
	//visualizeNRoSyFields(viewer);

	//convertNRoSyToRepVectors(nRoSy, repVector);
	//visualizeRepVectorFields(viewer);
	//convertRepVectorsToNRoSy(repVector, nRoSy);
	//visualizeNRoSyFields(viewer, nRoSy);

	computeFrameRotation(viewer);
	buildStiffnessMatrix_Geometric();
	computeEigenFields_generalized(25, fieldsfile);
	//buildStiffnessMatrix_Combinatorial();
	//computeEigenFields_regular(50, fieldsfile);
	NRoSy nRoSy_eigenFields;
	//convertRepVectorsToNRoSy(eigFieldsNRoSyRef.col(1), nRoSy_eigenFields);
	//visualizeNRoSyFields(viewer, nRoSy_eigenFields);
	visualizeRepVectorFields(viewer, eigFieldsNRoSyRef.col(0));
	

}


