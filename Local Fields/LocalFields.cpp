#include "LocalFields.h"

LocalFields::LocalFields(const int &sampleID)
{
	id = sampleID;
}

void LocalFields::constructSubdomain(const int &sampleID, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const double &avgEdgeLength, const Eigen::MatrixXi &AdjMF3N)
{
	int center = sampleID;
	this->sampleID = sampleID; 

	priority_queue<VertexPair, std::vector<VertexPair>, std::greater<VertexPair>> DistPQueue;
	Eigen::VectorXd D(F.rows());
	const double maxDist = 12.5 * avgEdgeLength;

	// Computing distance for initial sample points S
	for (int i = 0; i < F.rows(); i++) {
		D(i) = numeric_limits<double>::infinity();
	}

	D(center) = 0.0f;
	VertexPair vp{ center,D(center) };
	DistPQueue.push(vp);
	SubDomain.insert(center);

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
	} while (distFromCenter < maxDist);
	//cout << "Sample ID [" << id << "] = (" << sampleID << ") has " << SubDomain.size() << " elements." << endl; 
}

void LocalFields::constructSubdomain(const int &sampleID, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const double &avgEdgeLength, const Eigen::MatrixXi &AdjMF3N, const double& distRatio)
{
	int center = sampleID;
	this->sampleID = sampleID;

	priority_queue<VertexPair, std::vector<VertexPair>, std::greater<VertexPair>> DistPQueue;
	Eigen::VectorXd D(F.rows());
	const double maxDist = distRatio * avgEdgeLength;

	// Computing distance for initial sample points S
	for (int i = 0; i < F.rows(); i++) {
		D(i) = numeric_limits<double>::infinity();
	}

	D(center) = 0.0f;
	VertexPair vp{ center,D(center) };
	DistPQueue.push(vp);
	SubDomain.insert(center);

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
	} while (distFromCenter < maxDist);
}

void LocalFields::constructBoundary(const Eigen::MatrixXi& F, const Eigen::MatrixXi &AdjMF3N, const vector<set<int>> &AdjMF2Ring)
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

	// to get the Beyond Boundary regions
	for (std::set<int>::iterator it = outerBoundary.begin(); it != outerBoundary.end(); ++it) {
		for (std::set<int>::iterator jt = AdjMF2Ring[*it].begin(); jt != AdjMF2Ring[*it].end(); ++jt) {
			if (innerBoundary.find(*jt) == innerBoundary.end()) {
				BeyondBoundary.insert(*jt);
			}
		}
	}

	//cout << "Sample[" << id << "] has " << Boundary.size() << " elements in its boundary." << endl;
}

void LocalFields::constructLocalElements(const Eigen::MatrixXi &F)
{
	LocalElements.resize(SubDomain.size() + Boundary.size());
	//LocToGlobMap.resize(LocalElements.size());
	GlobToLocMap.resize(F.rows());
	for (int i = 0; i < F.rows(); i++) GlobToLocMap[i] = -1;

	int counter = 0;
	for (int face : SubDomain) {
		LocalElements[counter] = face;
		//LocToGlobMap[counter] = face;
		GlobToLocMap[face] = counter;
		counter++;
	}

	for (int face : Boundary) {
		LocalElements[counter] = face;
		GlobToLocMap[face] = counter;
		counter++;
	}

	//for (int face : BeyondBoundary) {
	//	GlobToLocMap[face] = counter++;
	//}

	//printf("Domain %d has %d elements (%d inner + %d boundary) \n", id, LocalElements.size(), SubDomain.size(), Boundary.size());
	//cout << "Sample ID [" << id << "] has " << LocalElements.size() << " elements." << endl;
}

void LocalFields::computeDijkstraFaceDistance(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &FC, const Eigen::MatrixXi &AdjMF3N)
{
	priority_queue<VertexPair, std::vector<VertexPair>, std::greater<VertexPair>> DistPQueue;
	
	/* Computing distance for initial sample points S */
	dijksFaceDist.resize(LocalElements.size());
	for (int i = 0; i < LocalElements.size(); i++) {
		dijksFaceDist(i) = numeric_limits<double>::infinity();
	}
	
	/* Initializing values for the sample/center element */
	const int center = GlobToLocMap[sampleID];
	dijksFaceDist(center) = 0.0f;
	VertexPair vp{ center, dijksFaceDist(center) };
	DistPQueue.push(vp);
	Eigen::RowVector3d const c0 = FC.row(sampleID);
	
	double distFromCenter = numeric_limits<double>::infinity();
	double maxDist = 0; 

	/* For other vertices in mesh */
	do {
		if (DistPQueue.size() == 0) break;
		VertexPair vp1 = DistPQueue.top();
		distFromCenter = vp1.distance;
		DistPQueue.pop();
	
		/* Updating the distance for neighbors of vertex of lowest distance in priority queue */
		int const elem = vp1.vId;
		//Eigen::Vector3d const c1 = (V.row(F(elem, 0)) + V.row(F(elem, 1)) + V.row(F(elem, 2))) / 3.0;
		Eigen::RowVector3d const c1 = FC.row(LocalElements[elem]);
		for (auto it = 0; it != F.cols(); ++it) {
			const int neigh = GlobToLocMap[AdjMF3N(LocalElements[elem], it)];
			if (neigh < 0) continue; // only work with those in local elements			
			Eigen::RowVector3d const c2 = FC.row(LocalElements[neigh]);
			
			/* Regular Dikjstra */
			//double dist = (c2 - c1).norm();
			//double tempDist = distFromCenter + dist;
			/* Dijkstra with distance correction */
			double tempDist = (c2 - c0).norm();
	
			/* updating the distance */
			if (tempDist <dijksFaceDist(neigh)) {
				dijksFaceDist(neigh) = tempDist;
				VertexPair vp2{ neigh,tempDist };
				DistPQueue.push(vp2);

				if (tempDist > maxDist) maxDist = tempDist; 
			}
		}
	} while (!DistPQueue.empty());

	/* Taking care of components in boundary */
	for (int i : Boundary)
	{
		//dijksFaceDist(GlobToLocMap[i]) = numeric_limits<double>::infinity();
		dijksFaceDist(GlobToLocMap[i]) = maxDist;
	}

	/* Conversion to the scaling factor */
	scalingFactor.setZero(LocalElements.size());
	const double pCube = 1 / (pow(maxDist, 3));
	const double pSquare = 1 / (pow(maxDist, 2));
	for (int i = 0; i < LocalElements.size(); i++)
	{
		scalingFactor(i) = 2 * dijksFaceDist(i)*dijksFaceDist(i)*dijksFaceDist(i)*pCube - 3 * dijksFaceDist(i)*dijksFaceDist(i)*pSquare + 1; 
	}
	
	/* Getting values for visualization */
	if (id == 0)
	{
		dijksFaceDistMapped.setOnes(F.rows());
		dijksFaceDistMapped *= numeric_limits<double>::infinity();
	
		for (int i = 0; i < LocalElements.size(); i++)
		{
			//dijksFaceDistMapped(LocalElements[i]) = dijksFaceDist(GlobToLocMap[LocalElements[i]]);
			//dijksFaceDistMapped(LocalElements[i]) = dijksFaceDist(i);
			dijksFaceDistMapped(LocalElements[i]) = scalingFactor(i);
		}
	}

}

void LocalFields::constructMatrixBLocal(const Eigen::SparseMatrix<double>& B2D)
{
	
	//======================== BLoc from B2D =========================
	BLoc.resize(2 * LocalElements.size(), 2 * LocalElements.size());
	vector<Eigen::Triplet<double>> BTriplet;
	BTriplet.reserve(4 * 2 * BLoc.rows());

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
	BLoc.prune(0.0);
	
	//cout << "Sample ID [" << id << "] has " << BLoc.rows() << "x" << BLoc.cols() << " size, with " << BLoc.nonZeros() << " non-zeros." << endl;

}

void LocalFields::constructMatrixBLocal(const Eigen::SparseMatrix<double>& B2D, const vector<set<int>>& AdjMF2Ring, vector<Eigen::Triplet<double>>& BTriplet)
{
	//======================== BLoc from B2D =========================
	BLoc.resize(2 * LocalElements.size(), 2 * LocalElements.size());

	//vector<Eigen::Triplet<double>> BTriplet;
	//BTriplet.reserve(100 * BLoc.rows());
	BTriplet.reserve(20 * 2 * LocalElements.size());
	
	


	for (int i = 0; i < LocalElements.size(); i++) {
		int li = LocalElements[i];
		// Get the DIAGONAL Elements from B2D Matrix
		BTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, B2D.coeff(2 * li + 0, 2 * li + 0)));
		BTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 1, B2D.coeff(2 * li + 0, 2 * li + 1)));
		BTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 0, B2D.coeff(2 * li + 1, 2 * li + 0)));
		BTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, B2D.coeff(2 * li + 1, 2 * li + 1)));

		// Get the NEIGHBORING elements
		//for (int j = 0; j < AdjMF3N.cols(); j++) {
		for(int j:AdjMF2Ring[LocalElements[i]]){
			const int neigh = j;
			if (GlobToLocMap[neigh] >= 0) {
				int neighLoc = GlobToLocMap[neigh];
				BTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * neighLoc + 0, B2D.coeff(2 * li + 0, 2 * neigh + 0)));
				BTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * neighLoc + 1, B2D.coeff(2 * li + 0, 2 * neigh + 1)));
				BTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * neighLoc + 0, B2D.coeff(2 * li + 1, 2 * neigh + 0)));
				BTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * neighLoc + 1, B2D.coeff(2 * li + 1, 2 * neigh + 1)));
			}
		}
	}
	BLoc.setFromTriplets(BTriplet.begin(), BTriplet.end());
}

void LocalFields::constructMatrixBLocalDirectInsert(const Eigen::SparseMatrix<double>& B2D, const vector<set<int>>& AdjMF2Ring)
{
	//======================== BLoc from B2D =========================
	BLoc.resize(2 * LocalElements.size(), 2 * LocalElements.size());
	BLoc.reserve(Eigen::VectorXd::Constant(2 * LocalElements.size(),20));
	vector<Eigen::Triplet<double>> BTriplet;


	for (int i = 0; i < LocalElements.size(); i++) {
		int li = LocalElements[i];
		// Get the DIAGONAL Elements from B2D Matrix
		BLoc.insert(2 * i + 0, 2 * i + 0) = B2D.coeff(2 * li + 0, 2 * li + 0);
		BLoc.insert(2 * i + 0, 2 * i + 1) = B2D.coeff(2 * li + 0, 2 * li + 1);
		BLoc.insert(2 * i + 1, 2 * i + 0) = B2D.coeff(2 * li + 1, 2 * li + 0);
		BLoc.insert(2 * i + 1, 2 * i + 1) = B2D.coeff(2 * li + 1, 2 * li + 1);
		//BTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, B2D.coeff(2 * li + 0, 2 * li + 0)));
		//BTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 1, B2D.coeff(2 * li + 0, 2 * li + 1)));
		//BTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 0, B2D.coeff(2 * li + 1, 2 * li + 0)));
		//BTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, B2D.coeff(2 * li + 1, 2 * li + 1)));

		// Get the NEIGHBORING elements
		//for (int j = 0; j < AdjMF3N.cols(); j++) {
		for (int j : AdjMF2Ring[LocalElements[i]]) {
			const int neigh = j;
			if (GlobToLocMap[neigh] >= 0) {
				int neighLoc = GlobToLocMap[neigh];
				//BTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * neighLoc + 0, B2D.coeff(2 * li + 0, 2 * neigh + 0)));
				//BTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * neighLoc + 1, B2D.coeff(2 * li + 0, 2 * neigh + 1)));
				//BTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * neighLoc + 0, B2D.coeff(2 * li + 1, 2 * neigh + 0)));
				//BTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * neighLoc + 1, B2D.coeff(2 * li + 1, 2 * neigh + 1)));
				BLoc.insert(2 * i + 0, 2 * neighLoc + 0) = B2D.coeff(2 * li + 0, 2 * neigh + 0);
				BLoc.insert(2 * i + 0, 2 * neighLoc + 1) = B2D.coeff(2 * li + 0, 2 * neigh + 1);
				BLoc.insert(2 * i + 1, 2 * neighLoc + 0) = B2D.coeff(2 * li + 1, 2 * neigh + 0);
				BLoc.insert(2 * i + 1, 2 * neighLoc + 1) = B2D.coeff(2 * li + 1, 2 * neigh + 1);
			}
		}
	}
	BLoc.setFromTriplets(BTriplet.begin(), BTriplet.end());
}

//void LocalFields::constructLocalConstraints()
void LocalFields::constructLocalConstraints(vector<Eigen::Triplet<double>>& C1Triplet, vector<Eigen::Triplet<double>>& C2Triplet)
{
	// Setting up matrix C
	vector<Eigen::Triplet<double>> CTriplet;
	CTriplet.reserve(2 + 2 * Boundary.size());
	//C1Triplet.reserve(2 * Boundary.size());
	//C2Triplet.reserve(2 * Boundary.size());
	int counter = 0;
	const int BCols = BLoc.cols();
	const int BRows = BLoc.rows();

	// Setting up vector c (There are 2 vector c)
	cLoc.resize(2 + 2 * Boundary.size(), 2);

	/* Set-up the constraint matrix C */
	CTriplet.push_back(Eigen::Triplet<double>(2*counter,     2 * GlobToLocMap[sampleID] + 0, 1.0));
	CTriplet.push_back(Eigen::Triplet<double>(2*counter + 1, 2 * GlobToLocMap[sampleID] + 1, 1.0));
	//C1Triplet.push_back(Eigen::Triplet<double>(BRows + counter,     2 * GlobToLocMap[sampleID] + 0, 1.0));
	//C1Triplet.push_back(Eigen::Triplet<double>(BRows + counter + 1, 2 * GlobToLocMap[sampleID] + 1, 1.0));
	//C2Triplet.push_back(Eigen::Triplet<double>(2 * GlobToLocMap[sampleID],     BCols + counter, 1.0));
	//C2Triplet.push_back(Eigen::Triplet<double>(2 * GlobToLocMap[sampleID] + 1, BCols + counter+1, 1.0));
	cLoc(2*counter, 0) = 1.0; cLoc(2*counter+1, 0) = 0.0; 
	cLoc(2*counter, 1) = 0.0; cLoc(2*counter+1, 1) = 1.0;
	counter++;
	for (int bound : Boundary) {
		CTriplet.push_back(Eigen::Triplet<double>(2*counter, 2 * GlobToLocMap[bound] + 0, 1.0));
		CTriplet.push_back(Eigen::Triplet<double>(2*counter + 1, 2 * GlobToLocMap[bound] + 1, 1.0));
		//C1Triplet.push_back(Eigen::Triplet<double>(BRows + counter, 2 * GlobToLocMap[bound], 1.0));
		//C1Triplet.push_back(Eigen::Triplet<double>(BRows + counter + 1, 2 * GlobToLocMap[bound] + 1, 1.0));
		//C2Triplet.push_back(Eigen::Triplet<double>(2 * GlobToLocMap[bound], BCols + counter, 1.0));
		//C2Triplet.push_back(Eigen::Triplet<double>(2 * GlobToLocMap[bound] + 1, BCols + counter + 1, 1.0));
		cLoc(2*counter, 0) = 0.0; cLoc(2*counter + 1, 0) = 0.0;
		cLoc(2*counter, 1) = 0.0; cLoc(2*counter + 1, 1) = 0.0;
		counter++;
	}
	CLoc.resize(2 + 2 * Boundary.size(), BLoc.cols());
	CLoc.setFromTriplets(CTriplet.begin(), CTriplet.end());	
}

void LocalFields::constructLocalConstraintsWithLaplacian(const Eigen::VectorXd& doubleArea, const Eigen::MatrixXi &AdjMF3N, const Eigen::SparseMatrix<double>& SF2D, vector<Eigen::Triplet<double>>& C1Triplet, vector<Eigen::Triplet<double>>& C2Triplet)
{
	// Setting up matrix C
	vector<Eigen::Triplet<double>> CTriplet;
	CTriplet.reserve(4 + 2 * Boundary.size());
	//C1Triplet.reserve(2 * Boundary.size());
	//C2Triplet.reserve(2 * Boundary.size());
	int counter = 0;
	const int BCols = BLoc.cols();
	const int BRows = BLoc.rows();

	// Setting up vector c (There are 2 vector c)
	cLoc.resize(4 + 2 * Boundary.size(), 2);
	

	/* Getting local Laplacian */
	//Eigen::Matrix2d SF2DLoc;
	//SF2DLoc(0, 0) = SF2D.coeff(2 * sampleID + 0, 2 * sampleID + 0);
	//SF2DLoc(0, 1) = SF2D.coeff(2 * sampleID + 0, 2 * sampleID + 1);
	//SF2DLoc(1, 0) = SF2D.coeff(2 * sampleID + 1, 2 * sampleID + 0);
	//SF2DLoc(1, 1) = SF2D.coeff(2 * sampleID + 1, 2 * sampleID + 1);
	//const double mInvLoc = 2.0/doubleArea(sampleID);				// 1/area of a triangle = 2 / doubleArea of a triangle
	//SF2DLoc *= mInvLoc; 

	


	/* Set-up the constraint matrix C 
	*  ==> Hard constraint on direction 	*/
	CTriplet.push_back(Eigen::Triplet<double>(2*counter, 2 * GlobToLocMap[sampleID] + 0, 1.0));
	CTriplet.push_back(Eigen::Triplet<double>(2*counter + 1, 2 * GlobToLocMap[sampleID] + 1, 1.0));
	cLoc(2*counter, 0) = 1.0; cLoc(2*counter+1, 0) = 0.0;
	cLoc(2*counter, 1) = 0.0; cLoc(2*counter+1, 1) = 1.0;
	counter++;
	
	
	/* Hard constraint -> Laplacian of the sample face */
	Eigen::MatrixXd SF2DBlock(2, 8);
	double mInv = 2.0 / doubleArea(sampleID);
	SF2DBlock(0, 0) = mInv * SF2D.coeff(2 * sampleID + 0, 2 * sampleID + 0);
	SF2DBlock(0, 1) = mInv * SF2D.coeff(2 * sampleID + 0, 2 * sampleID + 1);
	SF2DBlock(1, 0) = mInv * SF2D.coeff(2 * sampleID + 1, 2 * sampleID + 0);
	SF2DBlock(1, 1) = mInv * SF2D.coeff(2 * sampleID + 1, 2 * sampleID + 1);	
	CTriplet.push_back(Eigen::Triplet<double>(2 * counter,     2 * GlobToLocMap[sampleID] + 0,     SF2DBlock(0, 0)));
	CTriplet.push_back(Eigen::Triplet<double>(2 * counter,     2 * GlobToLocMap[sampleID] + 1,     SF2DBlock(0, 1)));
	CTriplet.push_back(Eigen::Triplet<double>(2 * counter + 1, 2 * GlobToLocMap[sampleID] + 0, SF2DBlock(1, 0)));
	CTriplet.push_back(Eigen::Triplet<double>(2 * counter + 1, 2 * GlobToLocMap[sampleID] + 1, SF2DBlock(1, 1)));
	for (int i = 0; i < 3; i++)
	{
		const int neigh = AdjMF3N(sampleID, i);
		mInv = 2.0 / doubleArea(neigh);
		SF2DBlock(0, 2+2*i+0) = mInv * SF2D.coeff(2 * neigh + 0, 2 * neigh + 0);
		SF2DBlock(0, 2+2*i+1) = mInv * SF2D.coeff(2 * neigh + 0, 2 * neigh + 1);
		SF2DBlock(1, 2+2*i+0) = mInv * SF2D.coeff(2 * neigh + 1, 2 * neigh + 0);
		SF2DBlock(1, 2+2*i+1) = mInv * SF2D.coeff(2 * neigh + 1, 2 * neigh + 1);
		CTriplet.push_back(Eigen::Triplet<double>(2 * counter,     2 * GlobToLocMap[neigh] + 0, SF2DBlock(0, 2 + 2 * i + 0)));
		CTriplet.push_back(Eigen::Triplet<double>(2 * counter,     2 * GlobToLocMap[neigh] + 1, SF2DBlock(0, 2 + 2 * i + 1)));
		CTriplet.push_back(Eigen::Triplet<double>(2 * counter + 1, 2 * GlobToLocMap[neigh] + 0, SF2DBlock(1, 2 + 2 * i + 0)));
		CTriplet.push_back(Eigen::Triplet<double>(2 * counter + 1, 2 * GlobToLocMap[neigh] + 1, SF2DBlock(1, 2 + 2 * i + 1)));
	}
	cLoc(2 * counter, 0) = 0.0;  cLoc(2 * counter + 1, 0) = 1.0;
	cLoc(2 * counter, 1) = -1.0; cLoc(2 * counter + 1, 1) = 0.0;
	counter++;

	/* == > Hard constraint on Vector Laplacian 	*/
	//CTriplet.push_back(Eigen::Triplet<double>(2*counter, ))
	//CTriplet.push_back(Eigen::Triplet<double>(2*counter,     2 * GlobToLocMap[sampleID] + 0, SF2DLoc(0,0)));
	//CTriplet.push_back(Eigen::Triplet<double>(2*counter,     2 * GlobToLocMap[sampleID] + 1, SF2DLoc(0,1)));
	//CTriplet.push_back(Eigen::Triplet<double>(2*counter + 1, 2 * GlobToLocMap[sampleID] + 0, SF2DLoc(1,0)));
	//CTriplet.push_back(Eigen::Triplet<double>(2*counter + 1, 2 * GlobToLocMap[sampleID] + 1, SF2DLoc(1,1)));
	//cLoc(2*counter, 0) = 0.0;  cLoc(2*counter + 1, 0) = 1.0;
	//cLoc(2*counter, 1) = -1.0; cLoc(2*counter + 1, 1) = 0.0;
	//counter++;


	//C1Triplet.push_back(Eigen::Triplet<double>(BRows + counter,     2 * GlobToLocMap[sampleID] + 0, 1.0));
	//C1Triplet.push_back(Eigen::Triplet<double>(BRows + counter + 1, 2 * GlobToLocMap[sampleID] + 1, 1.0));
	//C2Triplet.push_back(Eigen::Triplet<double>(2 * GlobToLocMap[sampleID],     BCols + counter, 1.0));
	//C2Triplet.push_back(Eigen::Triplet<double>(2 * GlobToLocMap[sampleID] + 1, BCols + counter+1, 1.0));
	for (int bound : Boundary) {
		CTriplet.push_back(Eigen::Triplet<double>(2*counter,     2 * GlobToLocMap[bound] + 0, 1.0));
		CTriplet.push_back(Eigen::Triplet<double>(2*counter + 1, 2 * GlobToLocMap[bound] + 1, 1.0));
		//C1Triplet.push_back(Eigen::Triplet<double>(BRows + counter, 2 * GlobToLocMap[bound], 1.0));
		//C1Triplet.push_back(Eigen::Triplet<double>(BRows + counter + 1, 2 * GlobToLocMap[bound] + 1, 1.0));
		//C2Triplet.push_back(Eigen::Triplet<double>(2 * GlobToLocMap[bound], BCols + counter, 1.0));
		//C2Triplet.push_back(Eigen::Triplet<double>(2 * GlobToLocMap[bound] + 1, BCols + counter + 1, 1.0));
		cLoc(2*counter, 0) = 0.0; cLoc(2*counter + 1, 0) = 0.0;
		cLoc(2*counter, 1) = 0.0; cLoc(2*counter + 1, 1) = 0.0;
		counter++;
	}
	CLoc.resize(4 + 2 * Boundary.size(), BLoc.cols());
	CLoc.setFromTriplets(CTriplet.begin(), CTriplet.end());

	// Setting up vector c (There are 2 vector c)
	//cLoc.resize(4 + 2 * Boundary.size(), 2);
	//Eigen::VectorXd zeroElements;
	//zeroElements.setZero(2 * Boundary.size());
	//cLoc.col(0) << 1.0, 0.0, 0.0, 1.0, zeroElements;
	//cLoc.col(1) << 0.0, 1.0, -1.0, 0.0, zeroElements;
	//cout << "Sample ID [" << id << "] has CLoc=" << CLoc.rows() << "x" << CLoc.cols() << ". cLoc " << cLoc.size() << "." << endl;
}


void LocalFields::setupRHSLocalProblemMapped()
{
	vEstimateLoc.resize(BLoc.cols());
	for (int i = 0; i < vEstimateLoc.rows(); i++) {
		vEstimateLoc(i) = 0.5;
	}

	gLoc = BLoc * vEstimateLoc;
	hLoc = CLoc * vEstimateLoc - cLoc.col(0);
	
	// First column of b
	bLoc.resize(BLoc.rows() + cLoc.rows(), 2);
	bLoc.col(0) << gLoc, hLoc;

	// Second Column of b
	hLoc = CLoc * vEstimateLoc - cLoc.col(1);
	bLoc.col(1) << gLoc, hLoc;
}

void LocalFields::setupLHSLocalProblemMapped(const vector<Eigen::Triplet<double>>& BTriplet, const vector<Eigen::Triplet<double>>& C1Triplet, const vector<Eigen::Triplet<double>>& C2Triplet)
{
	ALoc.resize(BLoc.rows() + CLoc.rows(), BLoc.cols() + CLoc.rows());

	vector<Eigen::Triplet<double>>				ATriplet;
	ATriplet.reserve(20 * BLoc.rows() + 2*CLoc.rows()); // In practice it's only 16, but allocate 20 for safety
	//ATriplet.reserve(BTriplet.size() + C1Triplet.size() + C2Triplet.size());
	//ATriplet.insert(ATriplet.end(), BTriplet.begin(), BTriplet.end());
	//ATriplet.insert(ATriplet.end(), C1Triplet.begin(), C1Triplet.end());
	//ATriplet.insert(ATriplet.end(), C2Triplet.begin(), C2Triplet.end());

	for (int k = 0; k < BLoc.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(BLoc, k); it; ++it) {
			ATriplet.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));		
		}
	}
	
	for (int k = 0; k < CLoc.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(CLoc, k); it; ++it) {
			/* PUSH BACK */
			ATriplet.push_back(Eigen::Triplet<double>(BLoc.rows() + it.row(), it.col(), it.value()));
			ATriplet.push_back(Eigen::Triplet<double>(it.col(), BLoc.cols() + it.row(), it.value()));			
		}
	}
	ALoc.setFromTriplets(ATriplet.begin(), ATriplet.end());
}

void LocalFields::solveLocalSystemMappedLDLT(vector<Eigen::Triplet<double>> &BTriplet)
{
	XfLoc.resize(BLoc.rows(), 2);
	//Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> sparseSolver(ALoc);
	Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> sparseSolver(ALoc);

	Eigen::VectorXd x = sparseSolver.solve(bLoc.col(0));

	if (sparseSolver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system. " << endl;
		if (sparseSolver.info() == Eigen::NumericalIssue)
			cout << "NUMERICAL ISSUE. " << endl;
		cout << sparseSolver.info() << endl;
		return;
	}

	XfLoc.col(0) = -x.block(0, 0, BLoc.rows(), 1) + vEstimateLoc;

	// SECOND BASIS
	x = sparseSolver.solve(bLoc.col(1));
	if (sparseSolver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system." << endl;
		return;
	}

	XfLoc.col(1) = -x.block(0, 0, BLoc.rows(), 1) + vEstimateLoc;

	/* Scale the vector fields accordingly */
	Eigen::MatrixXd XfLScaled(XfLoc.rows(), XfLoc.cols());
	Eigen::Vector2d vTemp;
	for (int i = 0; i < LocalElements.size(); i++)
	{
		if (scalingFactor(i) < 1e-8)
		{
			/* The boundary elements */
			XfLScaled(2 * i + 0, 0) = 0;
			XfLScaled(2 * i + 1, 0) = 0;
			XfLScaled(2 * i + 0, 1) = 0;
			XfLScaled(2 * i + 1, 1) = 0;
		}
		else
		{
			/* First components */
			vTemp = XfLoc.block(2 * i, 0, 2, 1);
			vTemp.normalize();
			XfLScaled.block(2 * i, 0, 2, 1) = scalingFactor(i)*vTemp;

			/* Second components */
			vTemp = XfLoc.block(2 * i, 1, 2, 1);
			vTemp.normalize();
			XfLScaled.block(2 * i, 1, 2, 1) = scalingFactor(i)*vTemp;
		}
	}

	/* With Scaling
	for (int i = 0; i < LocalElements.size(); i++) {
		// FIRST Option => Construct the BASIS first, then NORMALIZE it.
		// SECOND Option => NORMALIZE the element for each vector then Construct the BASIS first..
		// Let's do the FIRST option first, my gut feeling tells me this is a better opion..^^
		// First column ==> First basis (2 elements per-local frame)
		BTriplet.push_back(Eigen::Triplet<double>(2 * LocalElements[i] + 0, 2 * id + 0, XfLScaled(2 * i + 0, 0)));
		BTriplet.push_back(Eigen::Triplet<double>(2 * LocalElements[i] + 1, 2 * id + 0, XfLScaled(2 * i + 1, 0)));

		// Second column ==> Second basis (2 elements per-local frame)
		BTriplet.push_back(Eigen::Triplet<double>(2 * LocalElements[i] + 0, 2 * id + 1, XfLScaled(2 * i + 0, 1)));
		BTriplet.push_back(Eigen::Triplet<double>(2 * LocalElements[i] + 1, 2 * id + 1, XfLScaled(2 * i + 1, 1)));
	}
	*/ 


	/* No Scaling */
	for (int i = 0; i < LocalElements.size(); i++) {
		// FIRST Option => Construct the BASIS first, then NORMALIZE it.
		// SECOND Option => NORMALIZE the element for each vector then Construct the BASIS first..
		// Let's do the FIRST option first, my gut feeling tells me this is a better opion..^^
		// First column ==> First basis (2 elements per-local frame)
		BTriplet.push_back(Eigen::Triplet<double>(2 * LocalElements[i] + 0, 2 * id + 0, XfLoc(2 * i + 0, 0)));
		BTriplet.push_back(Eigen::Triplet<double>(2 * LocalElements[i] + 1, 2 * id + 0, XfLoc(2 * i + 1, 0)));

		// Second column ==> Second basis (2 elements per-local frame)
		BTriplet.push_back(Eigen::Triplet<double>(2 * LocalElements[i] + 0, 2 * id + 1, XfLoc(2 * i + 0, 1)));
		BTriplet.push_back(Eigen::Triplet<double>(2 * LocalElements[i] + 1, 2 * id + 1, XfLoc(2 * i + 1, 1)));
	}
	
}

void LocalFields::measureXF(const Eigen::VectorXd& doubleArea, const Eigen::SparseMatrix<double>& J)
{
	// Construct local mass matrix
	Eigen::SparseMatrix<double> MLocal;
	MLocal.resize(2 * LocalElements.size(), 2 * LocalElements.size());
	vector<Eigen::Triplet<double>> MTriplet;
	MTriplet.reserve(2 * LocalElements.size());

	for (int i = 0; i < LocalElements.size(); i++) {
		MTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, doubleArea(LocalElements[i])));
		MTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, doubleArea(LocalElements[i])));
	}
	MLocal.setFromTriplets(MTriplet.begin(), MTriplet.end());

	// Construct local Rotation matrix
	Eigen::SparseMatrix<double> JLocal;
	JLocal.resize(2 * LocalElements.size(), 2 * LocalElements.size());
	vector<Eigen::Triplet<double>> JTriplet;
	JTriplet.reserve(2 * LocalElements.size());

	for (int i = 0; i < LocalElements.size(); i++) {
		JTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, 0.0));
		JTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 0, 1.0));
		JTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 1, -1.0));
		JTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, 0.0));
	}
	JLocal.setFromTriplets(JTriplet.begin(), JTriplet.end());

	// Measure elements
	Eigen::VectorXd JXf0 = JLocal * XfLoc.col(0);
	Eigen::VectorXd diff = XfLoc.col(1) - JXf0; 

	double norm1 = diff.transpose() * MLocal * diff;
	double norm2 = XfLoc.col(0).transpose() * MLocal * XfLoc.col(0);
	double normXf = sqrt(norm1 / norm2);
	//printf("_______Error of %d is %.16f \n", id, normXf);
}