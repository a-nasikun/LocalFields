#include "LocalFields.h"
#include "EigenSolver.h"
#include "Utility.h"

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
	/* Defining the size of respective matrices */
	LocalElements.resize(SubDomain.size() + Boundary.size());
	InnerElements.resize(SubDomain.size());
	//LocToGlobMap.resize(LocalElements.size());
	GlobToLocMap.resize(F.rows());
	GlobToLocInnerMap.resize(F.rows());
	SelectorA.resize(2 * SubDomain.size(), 2 * LocalElements.size());
	vector<Eigen::Triplet<double>> SATriplet;
	SATriplet.reserve(2 * SubDomain.size());

	for (int i = 0; i < F.rows(); i++) 
	{ 
		GlobToLocMap[i] = -1; 
		GlobToLocInnerMap[i] = -1; 
	}

	int counter = 0;
	for (int face : SubDomain) {
		LocalElements[counter] = face;
		InnerElements[counter] = face; 
		//LocToGlobMap[counter] = face;
		GlobToLocMap[face] = counter;
		GlobToLocInnerMap[face] = counter; 
		SATriplet.push_back(Eigen::Triplet<double>(2 * counter + 0, 2 * counter + 0, 1.0));			/* Selector matrix */
		SATriplet.push_back(Eigen::Triplet<double>(2 * counter + 1, 2 * counter + 1, 1.0));
		counter++;
	}
	SelectorA.setFromTriplets(SATriplet.begin(), SATriplet.end());

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

void LocalFields::obtainLocalMatrixPatch2D(const Eigen::SparseMatrix<double>& MGlob, Eigen::SparseMatrix<double>& MPatch)
{
	//cout << id <<  "Setting up local patch matrix\n";
	const int localSize = LocalElements.size();
	MPatch.resize(2 * localSize, 2 * localSize);
	Eigen::SparseMatrix<double> MTemp;
	vector<Eigen::Triplet<double>> MTriplet, MTempTriplet;
	MTriplet.reserve(20 * 2 * localSize);
	MTempTriplet.reserve(20 * MGlob.rows());
	//Eigen::MatrixXd MTempDense;
	//MTempDense.setZero(MGlob.rows(), 2 * localSize);
	
	// Setup the big chunk of large matrix in dense format
	//cout << id << "Creating the big chunk matrix\n";
	for (int i = 0; i <localSize; i++)										// iterate over all local elements
	{
		const int localEl = LocalElements[i];
		for (int k = 2 * localEl; k <= 2 * localEl + 1; k++)			// there are 2 columns for each local elements
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(MGlob, k); it; ++it) {	// obtain elements on each column
				//int nEl = GlobToLocMap[floor(it.row() / 2)];
				//if (nEl < 0) continue;
				int elRow = 2 * i + (k-2*localEl);
				MTempTriplet.push_back(Eigen::Triplet<double>(elRow, it.row(), it.value()));
			}
		}
	}
	MTemp.resize(2 * localSize, MGlob.rows());
	MTemp.setFromTriplets(MTempTriplet.begin(), MTempTriplet.end());

	// Create local patch of the global matrix M (MGlob -> MPatch);
	//cout << id << "Creating the smaller chunk matrix\n";
	for (int i = 0; i <localSize; i++)										// iterate over all local elements
	{
		const int localEl = LocalElements[i];
		for (int k = 2 * localEl; k <= 2 * localEl + 1; k++)			// there are 2 columns for each local elements
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(MTemp, k); it; ++it) {	// obtain elements on each column
				int elCol = 2 * i + (k - 2 * localEl);
				MTriplet.push_back(Eigen::Triplet<double>(it.row(), elCol, it.value()));
			}
		}
	}	
	MPatch.setFromTriplets(MTriplet.begin(), MTriplet.end());
	//visualizeSparseMatrixInMatlab(MPatch);
	//cout << "[DONE!] Setting up local patch matrix\n";
}

/* Naive way, getting all values on the patch, then discarding the 0 elements/entries */
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

/* Look at the 1-ring elements -> correspond to the Laplace matrix neighbor */
void LocalFields::constructMatrixBLocal(const Eigen::SparseMatrix<double>& B2D, const vector<set<int>>& AdjMF2Ring)
{
	//======================== BLoc from B2D =========================
	BLoc.resize(2 * LocalElements.size(), 2 * LocalElements.size());
	
	vector<Eigen::Triplet<double>> BTriplet;
	BTriplet.reserve(20 * 2 * LocalElements.size());
	
	for (int i = 0; i < LocalElements.size(); i++) {
		int li = LocalElements[i];
		// Get the DIAGONAL Elements from B2D Matrix
		BTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, B2D.coeff(2 * li + 0, 2 * li + 0)));
		BTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 1, B2D.coeff(2 * li + 0, 2 * li + 1)));
		BTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 0, B2D.coeff(2 * li + 1, 2 * li + 0)));
		BTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, B2D.coeff(2 * li + 1, 2 * li + 1)));
	
		// Get the NEIGHBORING elements
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

/* Perhaps a bit smarter construction, visiting the non-zeros only 
*  currenlty work for symmetric matrix only */
void LocalFields::constructMatrixBLocal(const Eigen::SparseMatrix<double>& B2D, const vector<set<int>>& AdjMF2Ring, vector<Eigen::Triplet<double>>& BTriplet)
{
	obtainLocalMatrixPatch2D(B2D, BLoc);	
}

/* Direct insertion, slow as h*ll */
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

void LocalFields::constructLocalConstraintsWithLaplacian(const Eigen::VectorXd& doubleArea, const vector<set<int>>& AdjMF2Ring, const Eigen::SparseMatrix<double>& SF2D, vector<Eigen::Triplet<double>>& C1Triplet, vector<Eigen::Triplet<double>>& C2Triplet)
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
	CTriplet.push_back(Eigen::Triplet<double>(2 * counter + 0, 2 * GlobToLocMap[sampleID] + 0, mInv * SF2D.coeff(2 * sampleID + 0, 2 * sampleID + 0)));
	CTriplet.push_back(Eigen::Triplet<double>(2 * counter + 1, 2 * GlobToLocMap[sampleID] + 1, mInv * SF2D.coeff(2 * sampleID + 0, 2 * sampleID + 1)));
	CTriplet.push_back(Eigen::Triplet<double>(2 * counter + 1, 2 * GlobToLocMap[sampleID] + 0, mInv * SF2D.coeff(2 * sampleID + 1, 2 * sampleID + 0)));
	CTriplet.push_back(Eigen::Triplet<double>(2 * counter + 1, 2 * GlobToLocMap[sampleID] + 1, mInv * SF2D.coeff(2 * sampleID + 1, 2 * sampleID + 1)));
	
	/* Laplacian constraints on the neighboring faces */
	for (int neigh : AdjMF2Ring[sampleID]) {
		int neighLoc = GlobToLocMap[neigh];
		if (neighLoc >= 0) {			

			mInv = 2.0 / doubleArea(neigh);
			CTriplet.push_back(Eigen::Triplet<double>(2 * counter + 0, 2 * neighLoc + 0, mInv*SF2D.coeff(2*sampleID+0, 2 + 2 * neigh + 0)));
			CTriplet.push_back(Eigen::Triplet<double>(2 * counter + 0, 2 * neighLoc + 1, mInv*SF2D.coeff(2*sampleID+0, 2 + 2 * neigh + 1)));
			CTriplet.push_back(Eigen::Triplet<double>(2 * counter + 1, 2 * neighLoc + 0, mInv*SF2D.coeff(2*sampleID+1, 2 + 2 * neigh + 0)));
			CTriplet.push_back(Eigen::Triplet<double>(2 * counter + 1, 2 * neighLoc + 1, mInv*SF2D.coeff(2*sampleID+1, 2 + 2 * neigh + 1)));

			//BTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * neighLoc + 0, B2D.coeff(2 * li + 0, 2 * neigh + 0)));
			//BTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * neighLoc + 1, B2D.coeff(2 * li + 0, 2 * neigh + 1)));
			//BTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * neighLoc + 0, B2D.coeff(2 * li + 1, 2 * neigh + 0)));
			//BTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * neighLoc + 1, B2D.coeff(2 * li + 1, 2 * neigh + 1)));
		}
	}

	cout << "Sample " << id << " is done\n";

	cLoc(2 * counter, 0) = 0;  cLoc(2 * counter + 1, 0) = 0.01;
	cLoc(2 * counter, 1) = 0; cLoc(2 * counter + 1, 1) = 0;
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

//void LocalFields::constructLocalConstraintsWithLaplacian(const Eigen::VectorXd& doubleArea, const Eigen::SparseMatrix<double>& SF2D, vector<Eigen::Triplet<double>>& C1Triplet, vector<Eigen::Triplet<double>>& C2Triplet)
//{
//	// Setting up matrix C
//	vector<Eigen::Triplet<double>> CTriplet;
//	CTriplet.reserve(4 + 2 * Boundary.size());
//	//C1Triplet.reserve(2 * Boundary.size());
//	//C2Triplet.reserve(2 * Boundary.size());
//	int counter = 0;
//	const int BCols = BLoc.cols();
//	const int BRows = BLoc.rows();
//
//	// Setting up vector c (There are 2 vector c)
//	cLoc.resize(4 + 2 * Boundary.size(), 2);
//
//	/* Set-up the constraint matrix C
//	*  ==> Hard constraint on direction 	*/
//	CTriplet.push_back(Eigen::Triplet<double>(2 * counter, 2 * GlobToLocMap[sampleID] + 0, 1.0));
//	CTriplet.push_back(Eigen::Triplet<double>(2 * counter + 1, 2 * GlobToLocMap[sampleID] + 1, 1.0));
//	cLoc(2 * counter, 0) = 1.0; cLoc(2 * counter + 1, 0) = 0.0;
//	cLoc(2 * counter, 1) = 0.0; cLoc(2 * counter + 1, 1) = 1.0;
//	counter++;
//
//	int idB = GlobToLocMap[sampleID];
//	for (int k = 2*sampleID; k < 2*sampleID + 1; k++) {
//		for (Eigen::SparseMatrix<double>::InnerIterator it(SF2D, k); it; ++it) 
//		{	
//			int rowB = GlobToLocMap[floor(it.row()/2)];
//			//double mInv = 1.0 / doubleArea(floor(it.row() / 2));
//			if (rowB > 0)
//			{
//				CTriplet.push_back(Eigen::Triplet<double>(counter + (k-2*sampleID), rowB, it.value()));
//			}
//		}
//	}
//	cLoc(2 * counter, 0) = 0; cLoc(2 * counter + 1, 0) = 0;		// first contraint, i.e. on the first basis
//	cLoc(2 * counter, 1) = 0; cLoc(2 * counter + 1, 1) = 0;		// second constraint, i.e. on the second basis
//	counter++;
//
//
//	for (int bound : Boundary) {
//		CTriplet.push_back(Eigen::Triplet<double>(2 * counter, 2 * GlobToLocMap[bound] + 0, 1.0));
//		CTriplet.push_back(Eigen::Triplet<double>(2 * counter + 1, 2 * GlobToLocMap[bound] + 1, 1.0));
//		cLoc(2 * counter, 0) = 0.0; cLoc(2 * counter + 1, 0) = 0.0;
//		cLoc(2 * counter, 1) = 0.0; cLoc(2 * counter + 1, 1) = 0.0;
//		counter++;
//	}
//
//
//	/* Setting up the size of the [selection] constraint matrix */
//	CLoc.resize(4 + 2 * Boundary.size(), BLoc.cols());
//	CLoc.setFromTriplets(CTriplet.begin(), CTriplet.end());
//
//	cout << "ID " << id << " is done!\n";
//}

void LocalFields::constructLocalConstraintsWithLaplacian(const Eigen::VectorXd& doubleArea, const Eigen::SparseMatrix<double>& SF2D, vector<Eigen::Triplet<double>>& C1Triplet, vector<Eigen::Triplet<double>>& C2Triplet)
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

	/* Set-up the constraint matrix C
	*  ==> Hard constraint on direction 	*/
	CTriplet.push_back(Eigen::Triplet<double>(2 * counter, 2 * GlobToLocMap[sampleID] + 0, 1.0));
	CTriplet.push_back(Eigen::Triplet<double>(2 * counter + 1, 2 * GlobToLocMap[sampleID] + 1, 1.0));
	cLoc(2 * counter, 0) = 1.0; cLoc(2 * counter + 1, 0) = 0.0;
	cLoc(2 * counter, 1) = 0.0; cLoc(2 * counter + 1, 1) = 1.0;
	counter++;

	Eigen::SparseMatrix<double> SF2DLoc;
	obtainLocalMatrixPatch2D(SF2D, SF2DLoc);

	for (int i = 0; i < 2; i++)
	{
		for (Eigen::SparseMatrix<double>::InnerIterator it(SF2DLoc, i); it; ++it)
		{
			CTriplet.push_back(Eigen::Triplet<double>(2*counter+it.col(), it.row(), it.value()));
		}
	}
	cLoc(2 * counter, 0) = 0; cLoc(2 * counter + 1, 0) = 0;		// first contraint, i.e. on the first basis
	cLoc(2 * counter, 1) = 0; cLoc(2 * counter + 1, 1) = 0;		// second constraint, i.e. on the second basis
	counter++;

	//int idB = GlobToLocMap[sampleID];
	//for (int k = 2 * sampleID; k < 2 * sampleID + 1; k++) {
	//	for (Eigen::SparseMatrix<double>::InnerIterator it(SF2D, k); it; ++it)
	//	{
	//		int rowB = GlobToLocMap[floor(it.row() / 2)];
	//		//double mInv = 1.0 / doubleArea(floor(it.row() / 2));
	//		//if (rowB > 0)
	//		{
	//			CTriplet.push_back(Eigen::Triplet<double>(counter + (k - 2 * sampleID), rowB, it.value()));
	//		}
	//	}
	//}	


	for (int bound : Boundary) {
		CTriplet.push_back(Eigen::Triplet<double>(2 * counter, 2 * GlobToLocMap[bound] + 0, 1.0));
		CTriplet.push_back(Eigen::Triplet<double>(2 * counter + 1, 2 * GlobToLocMap[bound] + 1, 1.0));
		cLoc(2 * counter, 0) = 0.0; cLoc(2 * counter + 1, 0) = 0.0;
		cLoc(2 * counter, 1) = 0.0; cLoc(2 * counter + 1, 1) = 0.0;
		counter++;
	}


	/* Setting up the size of the [selection] constraint matrix */
	CLoc.resize(4 + 2 * Boundary.size(), BLoc.cols());
	CLoc.setFromTriplets(CTriplet.begin(), CTriplet.end());

	cout << "ID " << id << " is done!\n";
}

void LocalFields::constructLocalEigenProblem(const Eigen::SparseMatrix<double>& SF2D, const vector<set<int>>& AdjMF2Ring, Eigen::VectorXd& doubleArea, Eigen::MatrixXd &EigVectLocal)
{
	cout << "[MOD] Constructing local eigen problem\n ";
	//======================== BLoc from B2D =========================
	//BLoc.resize(2 * LocalElements.size(), 2 * LocalElements.size());
	Eigen::SparseMatrix<double> SF2DLoc, MF2DLoc;
	SF2DLoc.resize(2 * InnerElements.size(), 2 * InnerElements.size());
	MF2DLoc.resize(2 * InnerElements.size(), 2 * InnerElements.size());
	Eigen::VectorXd eigValsLoc;
	Eigen::MatrixXd EigVectLoc; 

	vector<Eigen::Triplet<double>> SFTriplet;
	SFTriplet.reserve(20 * 2 * InnerElements.size());
	vector<Eigen::Triplet<double>> MFTriplet;
	MFTriplet.reserve(2*InnerElements.size());

	cout << "Collecting inner elements\n";
	for (int i = 0; i < InnerElements.size(); i++) {
		int li = InnerElements[i];
		// Get the DIAGONAL Elements from B2D Matrix
		MFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, doubleArea(li)/2.0));
		MFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, doubleArea(li)/2.0));

		SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, SF2D.coeff(2 * li + 0, 2 * li + 0)));
		SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 1, SF2D.coeff(2 * li + 0, 2 * li + 1)));
		SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 0, SF2D.coeff(2 * li + 1, 2 * li + 0)));
		SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, SF2D.coeff(2 * li + 1, 2 * li + 1)));

		// Get the NEIGHBORING elements
		//for (int j = 0; j < AdjMF3N.cols(); j++) {
		for (int j : AdjMF2Ring[InnerElements[i]]) {
			const int neigh = j;
			if (GlobToLocInnerMap[neigh] >= 0) {
				int neighLoc = GlobToLocInnerMap[neigh];
				/* For non-diagonal elements */
				SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * neighLoc + 0, SF2D.coeff(2 * li + 0, 2 * neigh + 0)));
				SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * neighLoc + 1, SF2D.coeff(2 * li + 0, 2 * neigh + 1)));
				SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * neighLoc + 0, SF2D.coeff(2 * li + 1, 2 * neigh + 0)));
				SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * neighLoc + 1, SF2D.coeff(2 * li + 1, 2 * neigh + 1)));

				/* Diagonal elements */
				//SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, -SF2D.coeff(2 * li + 0, 2 * neigh + 0)));
				//SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 1, -SF2D.coeff(2 * li + 0, 2 * neigh + 1)));
				//SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 0, -SF2D.coeff(2 * li + 1, 2 * neigh + 0)));
				//SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, -SF2D.coeff(2 * li + 1, 2 * neigh + 1)));
			}
		}
	}
	cout << "Makign the matrix\n";
	SF2DLoc.setFromTriplets(SFTriplet.begin(), SFTriplet.end());
	cout << "SF2d loc done\n";
	MF2DLoc.setFromTriplets(MFTriplet.begin(), MFTriplet.end());
	cout << "MF2d loc done\n";

	computeEigenMatlab(SF2DLoc, MF2DLoc, 2, EigVectLoc, eigValsLoc, "hello");

	cout << "Eigenvectors \n";
	cout << EigVectLoc.block(0, 0, 100, 2) << endl << endl; 

	EigVectLocal.resize(SF2D.rows(), 2);
	for (int i = 0; i < InnerElements.size(); i++) 
	{
		EigVectLocal(2 * InnerElements[i] + 0, 0) = EigVectLoc(2 * i + 0, 0);
		EigVectLocal(2 * InnerElements[i] + 1, 0) = EigVectLoc(2 * i + 1, 0);
		EigVectLocal(2 * InnerElements[i] + 0, 1) = EigVectLoc(2 * i + 0, 1);
		EigVectLocal(2 * InnerElements[i] + 1, 1) = EigVectLoc(2 * i + 1, 1);		
	}
}

void LocalFields::constructLocalEigenProblemWithSelector(const Eigen::SparseMatrix<double>& SF2D, const vector<set<int>>& AdjMF2Ring, Eigen::VectorXd& doubleArea, Eigen::MatrixXd &EigLocal)
{
	cout << "[MOD] Constructing local eigen problem\n ";
	//======================== BLoc from B2D =========================
	//BLoc.resize(2 * LocalElements.size(), 2 * LocalElements.size());
	/* Setting up the matrices */
	Eigen::SparseMatrix<double> SF2DLoc, MF2DLoc, MF2DRed, SF2DRed;
	SF2DLoc.resize(2 * LocalElements.size(), 2 * LocalElements.size());
	MF2DLoc.resize(2 * LocalElements.size(), 2 * LocalElements.size());
	Eigen::VectorXd eigValsLoc;
	Eigen::MatrixXd EigVectLoc, eigTemp;
	vector<Eigen::Triplet<double>> SFTriplet;
	SFTriplet.reserve(20 * 2 * LocalElements.size());
	vector<Eigen::Triplet<double>> MFTriplet;
	MFTriplet.reserve(2 * LocalElements.size());

	/* Collecting the elements from the large matrices (in the whole surface) */
	cout << "Collecting inner elements\n";
	for (int i = 0; i < LocalElements.size(); i++) {
		int li = LocalElements[i];
		/*Get the DIAGONAL Elements for local patch mass Matrix*/
		MFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, doubleArea(li)/2.0));
		MFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, doubleArea(li)/2.0));

		/*Get the DIAGONAL Elements for local patch stiffness Matrix*/
		SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, SF2D.coeff(2 * li + 0, 2 * li + 0)));
		SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 1, SF2D.coeff(2 * li + 0, 2 * li + 1)));
		SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 0, SF2D.coeff(2 * li + 1, 2 * li + 0)));
		SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, SF2D.coeff(2 * li + 1, 2 * li + 1)));

		/*Get the NEIGHBORING elements 
		* i.e. Get the non-DIAGONAL Elements for local patch stiffness Matrix*/
		//for (int j = 0; j < AdjMF3N.cols(); j++) {
		for (int j : AdjMF2Ring[LocalElements[i]]) {
			const int neigh = j;
			if (GlobToLocMap[neigh] >= 0) {
				int neighLoc = GlobToLocMap[neigh];
				SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * neighLoc + 0, SF2D.coeff(2 * li + 0, 2 * neigh + 0)));
				SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * neighLoc + 1, SF2D.coeff(2 * li + 0, 2 * neigh + 1)));
				SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * neighLoc + 0, SF2D.coeff(2 * li + 1, 2 * neigh + 0)));
				SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * neighLoc + 1, SF2D.coeff(2 * li + 1, 2 * neigh + 1)));
			}
		}
	}

	SF2DLoc.setFromTriplets(SFTriplet.begin(), SFTriplet.end());
	MF2DLoc.setFromTriplets(MFTriplet.begin(), MFTriplet.end());

	/* Reduced matrices */
	MF2DRed = SelectorA * MF2DLoc * SelectorA.transpose();
	SF2DRed = SelectorA * SF2DLoc * SelectorA.transpose();

	/* Getting the eigenfields*/
	computeEigenMatlab(SF2DRed, MF2DRed, 2, eigTemp, eigValsLoc, "hello");
	EigVectLoc = SelectorA.transpose() * eigTemp;

	/* Mapping to larger matrix */
	EigLocal.setZero(SF2D.rows(), 2);
	for (int i = 0; i < InnerElements.size(); i++)
	{
		EigLocal(2 * InnerElements[i] + 0, 0) = EigVectLoc(2 * i + 0, 0);
		EigLocal(2 * InnerElements[i] + 1, 0) = EigVectLoc(2 * i + 1, 0);
		EigLocal(2 * InnerElements[i] + 0, 1) = EigVectLoc(2 * i + 0, 1);
		EigLocal(2 * InnerElements[i] + 1, 1) = EigVectLoc(2 * i + 1, 1);
	}
}

void LocalFields::constructLocalEigenProblem(const Eigen::SparseMatrix<double>& SF2D, const vector<set<int>>& AdjMF2Ring, Eigen::VectorXd& doubleArea, vector<Eigen::Triplet<double>>& BTriplet)
{
	//cout << "Constructing local eigen problem\n ";
	//======================== BLoc from B2D =========================
	//BLoc.resize(2 * LocalElements.size(), 2 * LocalElements.size());
	Eigen::SparseMatrix<double> SF2DLoc, MF2DLoc;
	SF2DLoc.resize(2 * InnerElements.size(), 2 * InnerElements.size());
	MF2DLoc.resize(2 * InnerElements.size(), 2 * InnerElements.size());
	Eigen::VectorXd eigValsLoc;
	Eigen::MatrixXd EigVectLoc;

	vector<Eigen::Triplet<double>> SFTriplet;
	SFTriplet.reserve(20 * 2 * InnerElements.size());
	vector<Eigen::Triplet<double>> MFTriplet;
	MFTriplet.reserve(2 * InnerElements.size());

	//cout << "MOD == Collecting inner elements\n";
	printf("Setting up %d [%d] components\n", id, sampleID);
	for (int i = 0; i < InnerElements.size(); i++) {
		int li = InnerElements[i];
		// Get the DIAGONAL Elements from B2D Matrix
		MFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, doubleArea(li)/2.0));
		MFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, doubleArea(li)/2.0));

		//SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, SF2D.coeff(2 * li + 0, 2 * li + 0)));
		//SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 1, SF2D.coeff(2 * li + 0, 2 * li + 1)));
		//SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 0, SF2D.coeff(2 * li + 1, 2 * li + 0)));
		//SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, SF2D.coeff(2 * li + 1, 2 * li + 1)));

		// Get the NEIGHBORING elements		
		//for (int j = 0; j < AdjMF3N.cols(); j++) {
		for (int j : AdjMF2Ring[InnerElements[i]]) {
			const int neigh = j;
			if (GlobToLocInnerMap[neigh] >= 0) {
				int neighLoc = GlobToLocInnerMap[neigh];
				/* Non diagonal elements*/
				SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * neighLoc + 0, SF2D.coeff(2 * li + 0, 2 * neigh + 0)));
				SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * neighLoc + 1, SF2D.coeff(2 * li + 0, 2 * neigh + 1)));
				SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * neighLoc + 0, SF2D.coeff(2 * li + 1, 2 * neigh + 0)));
				SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * neighLoc + 1, SF2D.coeff(2 * li + 1, 2 * neigh + 1)));

				/* Diagonal elements */
				SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, -1.0 * SF2D.coeff(2 * li + 0, 2 * neigh + 0)));
				SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 1, -1.0 * SF2D.coeff(2 * li + 0, 2 * neigh + 1)));
				SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 0, -1.0 * SF2D.coeff(2 * li + 1, 2 * neigh + 0)));
				SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, -1.0 * SF2D.coeff(2 * li + 1, 2 * neigh + 1)));
			}
		}
	}
	//cout << "Makign the matrix\n";
	SF2DLoc.setFromTriplets(SFTriplet.begin(), SFTriplet.end());
	//cout << "SF2d loc done\n";
	MF2DLoc.setFromTriplets(MFTriplet.begin(), MFTriplet.end());
	//cout << "MF2d loc done\n";
	//if (id == 0) visualizeSparseMatrixInMatlab(SF2DLoc);


	computeEigenMatlab(SF2DLoc, MF2DLoc, 4, EigVectLoc, eigValsLoc, "hello");
	//printf("EigVec Size=%dx%d \t\t elements=%d\n", EigVectLoc.rows(), EigVectLoc.cols(), InnerElements.size());

	printf("[%d] Eigenvalues: %.5f, %.5f, %.5f, %.5f\n", id, eigValsLoc(0), eigValsLoc(1), eigValsLoc(2), eigValsLoc(3));
	Eigen::VectorXd vect1; vect1.resize(EigVectLoc.col(0).size());
	Eigen::VectorXd vect2; vect2.resize(EigVectLoc.col(1).size());
	Eigen::VectorXd EigVectTemp1, EigVectTemp2;
	EigVectTemp1.resize(EigVectLoc.rows());
	EigVectTemp2.resize(EigVectLoc.rows());

	//if (id == 0)
	//{
	//	for (int k = 0; k < InnerElements.size(); k++)
	//	{
	//		Eigen::Vector2d v1;
	//		v1 << EigVectLoc(2 * k, 0), EigVectLoc(2 * k + 1, 0); //;;, .block(2 * k, 0, 2, 1);
	//		v1.normalize();
	//		EigVectTemp1(2 * k+0) = v1(0);
	//		EigVectTemp1(2 * k+1) = v1(1);
	//
	//		Eigen::Vector2d v2;
	//		v2 << EigVectLoc(2 * k, 1), EigVectLoc(2 * k + 1, 1); //;;, .block(2 * k, 0, 2, 1);
	//		v2.normalize();
	//		EigVectTemp2(2 * k + 0) = v2(0);
	//		EigVectTemp2(2 * k + 1) = v2(1);
	//		//cout << k << ": " << v1 << endl; 
	//	}
	//}

	//double energy1 = EigVectTemp1.transpose()*SF2DLoc*EigVectTemp1;
	//double energy2 = EigVectTemp2.transpose()*SF2DLoc*EigVectTemp2;
	//double energy1 = EigVectLoc.col(0).transpose()*SF2DLoc* EigVectLoc.col(0);
	//double energy2 = EigVectLoc.col(1).transpose()*SF2DLoc* EigVectLoc.col(1);
	//printf("[%d] Energy: %.11f \t\t %.11f\n", id, energy1, energy2);
	//cout << "Eigenvectors \n";
	//cout << EigVectLoc.block(0, 0, 100, 2) << endl << endl;
	//printf("Local elements: %d, eigveacts col=%d \n", InnerElements.size(), EigVectLoc.rows());
	//for (int k = 0; k < InnerElements.size(); k++)
	//{
	//	//cout << "k= " << k << endl; 
	//	Eigen::Vector2d v1; 
	//	v1 = EigVectLoc.block(2 * k, 0, 2, 1);
	//	//v1 << EigVectLoc( 2 * k + 0,0), EigVectLoc(2 * k + 1,0);
	//	//cout << "\t" << v1;
	//	v1.normalize();
	//	vect1.block(2 * k, 0, 2, 1) = v1; 
	//	//cout << "\t" << v1;
	//
	//	Eigen::Vector2d v2;
	//	//v2 << EigVectLoc(2 * k + 0,1), EigVectLoc(2 * k + 1,1);
	//	v2 = EigVectLoc.block(2 * k, 1, 2, 1);
	//	v2.normalize();
	//	vect2.block(2 * k, 1, 2, 1) = v2;
	//	//cout << vect2 << endl; 
	//}
	//printf("Size: Vec1=%d, SF2D=%dx%d; Vec2=%d\n", vect1.rows(), SF2DLoc.rows(), SF2DLoc.cols(), vect2.rows());
	
	for (int i = 0; i < InnerElements.size(); i++)
	{
		// First column ==> First basis (2 elements per-local frame)
		BTriplet.push_back(Eigen::Triplet<double>(2 * InnerElements[i] + 0, 2 * id + 0, EigVectLoc(2 * i + 0, 0)));
		BTriplet.push_back(Eigen::Triplet<double>(2 * InnerElements[i] + 1, 2 * id + 0, EigVectLoc(2 * i + 1, 0)));

		// Second column ==> Second basis (2 elements per-local frame)
		BTriplet.push_back(Eigen::Triplet<double>(2 * InnerElements[i] + 0, 2 * id + 1, EigVectLoc(2 * i + 0, 1)));
		BTriplet.push_back(Eigen::Triplet<double>(2 * InnerElements[i] + 1, 2 * id + 1, EigVectLoc(2 * i + 1, 1)));
	}	
}

void LocalFields:: constructLocalEigenProblem(const Eigen::SparseMatrix<double>& SF2D, const Eigen::MatrixXd &AdjMF3N, Eigen::VectorXd& doubleArea, vector<Eigen::Triplet<double>>& BTriplet)
{
	//cout << "Constructing local eigen problem\n ";
	//======================== BLoc from B2D =========================
	//BLoc.resize(2 * LocalElements.size(), 2 * LocalElements.size());
	Eigen::SparseMatrix<double> SF2DLoc, MF2DLoc;
	SF2DLoc.resize(2 * InnerElements.size(), 2 * InnerElements.size());
	MF2DLoc.resize(2 * InnerElements.size(), 2 * InnerElements.size());
	Eigen::VectorXd eigValsLoc;
	Eigen::MatrixXd EigVectLoc;

	vector<Eigen::Triplet<double>> SFTriplet;
	SFTriplet.reserve(20 * 2 * InnerElements.size());
	vector<Eigen::Triplet<double>> MFTriplet;
	MFTriplet.reserve(2 * InnerElements.size());

	//cout << "MOD == Collecting inner elements\n";
	for (int i = 0; i < InnerElements.size(); i++) {
		int li = InnerElements[i];
		// Get the DIAGONAL Elements from B2D Matrix
		MFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, doubleArea(li)));
		MFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, doubleArea(li)));

		SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, SF2D.coeff(2 * li + 0, 2 * li + 0)));
		SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 1, SF2D.coeff(2 * li + 0, 2 * li + 1)));
		SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 0, SF2D.coeff(2 * li + 1, 2 * li + 0)));
		SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, SF2D.coeff(2 * li + 1, 2 * li + 1)));

		// Get the NEIGHBORING elements
		for (int j = 0; j < AdjMF3N.cols(); j++) {
			//for (int j : AdjMF2Ring[InnerElements[i]]) {
			const int neigh = j;
			if (GlobToLocInnerMap[neigh] >= 0) {
				int neighLoc = GlobToLocInnerMap[neigh];
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
	//cout << "Makign the matrix\n";
	SF2DLoc.setFromTriplets(SFTriplet.begin(), SFTriplet.end());
	//cout << "SF2d loc done\n";
	MF2DLoc.setFromTriplets(MFTriplet.begin(), MFTriplet.end());
	//cout << "MF2d loc done\n";
	//if (id == 0) visualizeSparseMatrixInMatlab(SF2DLoc);


	computeEigenMatlab(SF2DLoc, MF2DLoc, 4, EigVectLoc, eigValsLoc, "hello");
	//printf("EigVec Size=%dx%d \t\t elements=%d\n", EigVectLoc.rows(), EigVectLoc.cols(), InnerElements.size());

	printf("[%d] Eigenvalues: %.5f, %.5f, %.5f, %.5f\n", id, eigValsLoc(0), eigValsLoc(1), eigValsLoc(2), eigValsLoc(3));
	Eigen::VectorXd vect1; vect1.resize(EigVectLoc.col(0).size());
	Eigen::VectorXd vect2; vect2.resize(EigVectLoc.col(1).size());
	Eigen::VectorXd EigVectTemp1, EigVectTemp2;
	EigVectTemp1.resize(EigVectLoc.rows());
	EigVectTemp2.resize(EigVectLoc.rows());

	//if (id == 0)
	{
		for (int k = 0; k < InnerElements.size(); k++)
		{
			Eigen::Vector2d v1;
			v1 << EigVectLoc(2 * k, 0), EigVectLoc(2 * k + 1, 0); //;;, .block(2 * k, 0, 2, 1);
			v1.normalize();
			EigVectTemp1(2 * k + 0) = v1(0);
			EigVectTemp1(2 * k + 1) = v1(1);

			Eigen::Vector2d v2;
			v2 << EigVectLoc(2 * k, 1), EigVectLoc(2 * k + 1, 1); //;;, .block(2 * k, 0, 2, 1);
			v2.normalize();
			EigVectTemp2(2 * k + 0) = v2(0);
			EigVectTemp2(2 * k + 1) = v2(1);
			//cout << k << ": " << v1 << endl; 
		}
	}

	//double energy1 = EigVectTemp1.transpose()*SF2DLoc*EigVectTemp1;
	//double energy2 = EigVectTemp2.transpose()*SF2DLoc*EigVectTemp2;
	double energy1 = EigVectLoc.col(0).transpose()*SF2DLoc* EigVectLoc.col(0);
	double energy2 = EigVectLoc.col(1).transpose()*SF2DLoc* EigVectLoc.col(1);
	printf("[%d] Energy: %.11f \t\t %.11f\n", id, energy1, energy2);
	//cout << "Eigenvectors \n";
	//cout << EigVectLoc.block(0, 0, 100, 2) << endl << endl;
	//printf("Local elements: %d, eigveacts col=%d \n", InnerElements.size(), EigVectLoc.rows());
	//for (int k = 0; k < InnerElements.size(); k++)
	//{
	//	//cout << "k= " << k << endl; 
	//	Eigen::Vector2d v1; 
	//	v1 = EigVectLoc.block(2 * k, 0, 2, 1);
	//	//v1 << EigVectLoc( 2 * k + 0,0), EigVectLoc(2 * k + 1,0);
	//	//cout << "\t" << v1;
	//	v1.normalize();
	//	vect1.block(2 * k, 0, 2, 1) = v1; 
	//	//cout << "\t" << v1;
	//
	//	Eigen::Vector2d v2;
	//	//v2 << EigVectLoc(2 * k + 0,1), EigVectLoc(2 * k + 1,1);
	//	v2 = EigVectLoc.block(2 * k, 1, 2, 1);
	//	v2.normalize();
	//	vect2.block(2 * k, 1, 2, 1) = v2;
	//	//cout << vect2 << endl; 
	//}
	//printf("Size: Vec1=%d, SF2D=%dx%d; Vec2=%d\n", vect1.rows(), SF2DLoc.rows(), SF2DLoc.cols(), vect2.rows());

	
	for (int i = 0; i < InnerElements.size(); i++)
	{
		// First column ==> First basis (2 elements per-local frame)
		BTriplet.push_back(Eigen::Triplet<double>(2 * InnerElements[i] + 0, 2 * id + 0, EigVectLoc(2 * i + 0, 0)));
		BTriplet.push_back(Eigen::Triplet<double>(2 * InnerElements[i] + 1, 2 * id + 0, EigVectLoc(2 * i + 1, 0)));

		// Second column ==> Second basis (2 elements per-local frame)
		BTriplet.push_back(Eigen::Triplet<double>(2 * InnerElements[i] + 0, 2 * id + 1, EigVectLoc(2 * i + 0, 1)));
		BTriplet.push_back(Eigen::Triplet<double>(2 * InnerElements[i] + 1, 2 * id + 1, EigVectLoc(2 * i + 1, 1)));
	}
}

//void LocalFields::constructLocalEigenProblemWithSelector(const Eigen::SparseMatrix<double>& SF2D, const Eigen::SparseMatrix<double>& MF2D, const vector<set<int>>& AdjMF2Ring, Eigen::VectorXd& doubleArea, vector<Eigen::Triplet<double>>& BTriplet)
//{
//	cout << "[" << id << "] Constructing local eigen problem\n ";
//	Eigen::SparseMatrix<double> SF2DLoc, MF2DLoc, MF2DRed, SF2DRed;
//	Eigen::VectorXd eigValsLoc;
//	Eigen::MatrixXd EigVectLoc, eigTemp;
//	//SF2DLoc.resize(2 * LocalElements.size(), 2 * LocalElements.size());
//	//MF2DLoc.resize(2 * LocalElements.size(), 2 * LocalElements.size());
//	//
//	//vector<Eigen::Triplet<double>> SFTriplet;
//	//SFTriplet.reserve(20 * 2 * LocalElements.size());
//	//vector<Eigen::Triplet<double>> MFTriplet;
//	//MFTriplet.reserve(2 * LocalElements.size());
//	//
//	
//	obtainLocalMatrixPatch2D(MF2D, MF2DLoc);
//	obtainLocalMatrixPatch2D(SF2D, SF2DLoc);
//
//	/* Reduced matrices */
//	MF2DRed = SelectorA * MF2DLoc * SelectorA.transpose();
//	SF2DRed = SelectorA * SF2DLoc * SelectorA.transpose();
//
//	/* Getting the eigenfields*/
//	computeEigenMatlab(SF2DRed, MF2DRed, 2, eigTemp, eigValsLoc, "hello");
//	EigVectLoc = SelectorA.transpose() * eigTemp;
//
//	/* Mapping to larger matrix */
//	for (int i = 0; i < InnerElements.size(); i++)
//	{
//		// First column ==> First basis (2 elements per-local frame)
//		BTriplet.push_back(Eigen::Triplet<double>(2 * InnerElements[i] + 0, 2 * id + 0, EigVectLoc(2 * i + 0, 0)));
//		BTriplet.push_back(Eigen::Triplet<double>(2 * InnerElements[i] + 1, 2 * id + 0, EigVectLoc(2 * i + 1, 0)));
//
//		// Second column ==> Second basis (2 elements per-local frame)
//		BTriplet.push_back(Eigen::Triplet<double>(2 * InnerElements[i] + 0, 2 * id + 1, EigVectLoc(2 * i + 0, 1)));
//		BTriplet.push_back(Eigen::Triplet<double>(2 * InnerElements[i] + 1, 2 * id + 1, EigVectLoc(2 * i + 1, 1)));
//	}	
//}

void LocalFields::constructLocalEigenProblemWithSelector(const Eigen::SparseMatrix<double>& SF2D, const Eigen::SparseMatrix<double>& MF2D, const vector<set<int>>& AdjMF2Ring, const int& NUM_EIG, const Eigen::VectorXd& doubleArea, vector<Eigen::Triplet<double>>& BTriplet)
{
	cout << "[" << id << "] Constructing local eigen problem\n ";
	Eigen::SparseMatrix<double> SF2DLoc, MF2DLoc, MF2DRed, SF2DRed;
	Eigen::VectorXd eigValsLoc;
	Eigen::MatrixXd EigVectLoc, eigTemp;
	MF2DLoc.resize(2 * LocalElements.size(), 2 * LocalElements.size());		

	obtainLocalMatrixPatch2D(MF2D, MF2DLoc);
	obtainLocalMatrixPatch2D(SF2D, SF2DLoc);

	//if (id == 0)
	//{
	//	visualizeSparseMatrixInMatlab(MF2DLoc);
	//}

	/* Reduced matrices */
	MF2DRed = SelectorA * MF2DLoc * SelectorA.transpose();
	SF2DRed = SelectorA * SF2DLoc * SelectorA.transpose();

	/* Getting the eigenfields*/
	computeEigenMatlab(SF2DRed, MF2DRed, NUM_EIG, eigTemp, eigValsLoc, "hello");
	//cusolverDnHandle_t	cusolverH;
	//computeEigenGPU(SF2DRed, MF2DRed, eigTemp, eigValsLoc);
	//cout << "ID=" << id << ", eig vals=\n" << eigValsLoc << endl; 
	EigVectLoc = SelectorA.transpose() * eigTemp;

	/* Mapping to larger matrix */
	for (int i = 0; i < InnerElements.size(); i++)
	{
		// First column ==> First basis (2 elements per-local frame)
		for (int j = 0; j < NUM_EIG; j++)
		{
			BTriplet.push_back(Eigen::Triplet<double>(2 * InnerElements[i] + 0, NUM_EIG * id + j, EigVectLoc(2 * i + 0, j)));
			BTriplet.push_back(Eigen::Triplet<double>(2 * InnerElements[i] + 1, NUM_EIG * id + j, EigVectLoc(2 * i + 1, j)));
		}

		//BTriplet.push_back(Eigen::Triplet<double>(NUM_EIG * InnerElements[i] + 0, NUM_EIG * id + 0, EigVectLoc(2 * i + 0, 0)));
		//BTriplet.push_back(Eigen::Triplet<double>(NUM_EIG * InnerElements[i] + 1, NUM_EIG * id + 0, EigVectLoc(2 * i + 1, 0)));

		// Second column ==> Second basis (2 elements per-local frame)
		//BTriplet.push_back(Eigen::Triplet<double>(NUM_EIG * InnerElements[i] + 0, NUM_EIG * id + 1, EigVectLoc(2 * i + 0, 1)));
		//BTriplet.push_back(Eigen::Triplet<double>(NUM_EIG * InnerElements[i] + 1, NUM_EIG * id + 1, EigVectLoc(2 * i + 1, 1)));
	}
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
	if (false)
	{
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
	}

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