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
	LocalElements.resize(SubDomain.size() + Boundary.size());
	InnerElements.resize(SubDomain.size());
	//LocToGlobMap.resize(LocalElements.size());
	GlobToLocMap.resize(F.rows());
	GlobToLocInnerMap.resize(F.rows());
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
	CTriplet.reserve(2 * Boundary.size());
	C1Triplet.reserve(2 * Boundary.size());
	C2Triplet.reserve(2 * Boundary.size());
	int counter = 0;
	const int BCols = BLoc.cols();
	const int BRows = BLoc.rows();
	CTriplet.push_back(Eigen::Triplet<double>(counter,     2 * GlobToLocMap[sampleID] + 0, 1.0));
	CTriplet.push_back(Eigen::Triplet<double>(counter + 1, 2 * GlobToLocMap[sampleID] + 1, 1.0));
	C1Triplet.push_back(Eigen::Triplet<double>(BRows + counter,     2 * GlobToLocMap[sampleID] + 0, 1.0));
	C1Triplet.push_back(Eigen::Triplet<double>(BRows + counter + 1, 2 * GlobToLocMap[sampleID] + 1, 1.0));
	C2Triplet.push_back(Eigen::Triplet<double>(2 * GlobToLocMap[sampleID],     BCols + counter, 1.0));
	C2Triplet.push_back(Eigen::Triplet<double>(2 * GlobToLocMap[sampleID] + 1, BCols + counter+1, 1.0));
	counter++;
	for (int bound : Boundary) {
		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * GlobToLocMap[bound] + 0, 1.0));
		CTriplet.push_back(Eigen::Triplet<double>(counter + 1, 2 * GlobToLocMap[bound] + 1, 1.0));
		C1Triplet.push_back(Eigen::Triplet<double>(BRows + counter, 2 * GlobToLocMap[bound], 1.0));
		C1Triplet.push_back(Eigen::Triplet<double>(BRows + counter + 1, 2 * GlobToLocMap[bound] + 1, 1.0));
		C2Triplet.push_back(Eigen::Triplet<double>(2 * GlobToLocMap[bound], BCols + counter, 1.0));
		C2Triplet.push_back(Eigen::Triplet<double>(2 * GlobToLocMap[bound] + 1, BCols + counter + 1, 1.0));
		counter++;
	}
	CLoc.resize(2 + 2 * Boundary.size(), BLoc.cols());
	CLoc.setFromTriplets(CTriplet.begin(), CTriplet.end());

	// Setting up vector c (There are 2 vector c)
	cLoc.resize(2 + 2 * Boundary.size(), 2);
	Eigen::VectorXd zeroElements(2 * Boundary.size());
	for (int i = 0; i < zeroElements.size(); i++) zeroElements(i) = 0.0;
	cLoc.col(0) << 1.0, 0.0, zeroElements;
	cLoc.col(1) << 0.0, 1.0, zeroElements;
		
	//cout << "Sample ID [" << id << "] has CLoc=" << CLoc.rows() << "x" << CLoc.cols() << ". cLoc " << cLoc.size() << "." << endl;
}

void LocalFields::constructLocalEigenProblem(const Eigen::SparseMatrix<double>& SF2D, const vector<set<int>>& AdjMF2Ring, Eigen::VectorXd& doubleArea, Eigen::MatrixXd &EigVectLocal)
{
	cout << "Constructing local eigen problem\n ";
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
		MFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, doubleArea(li)));
		MFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, doubleArea(li)));

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
				SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * neighLoc + 0, SF2D.coeff(2 * li + 0, 2 * neigh + 0)));
				SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * neighLoc + 1, SF2D.coeff(2 * li + 0, 2 * neigh + 1)));
				SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * neighLoc + 0, SF2D.coeff(2 * li + 1, 2 * neigh + 0)));
				SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * neighLoc + 1, SF2D.coeff(2 * li + 1, 2 * neigh + 1)));
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

	//cout << "Collecting inner elements\n";
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
		//for (int j = 0; j < AdjMF3N.cols(); j++) {
		for (int j : AdjMF2Ring[InnerElements[i]]) {
			const int neigh = j;
			if (GlobToLocInnerMap[neigh] >= 0) {
				int neighLoc = GlobToLocInnerMap[neigh];
				SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * neighLoc + 0, SF2D.coeff(2 * li + 0, 2 * neigh + 0)));
				SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * neighLoc + 1, SF2D.coeff(2 * li + 0, 2 * neigh + 1)));
				SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * neighLoc + 0, SF2D.coeff(2 * li + 1, 2 * neigh + 0)));
				SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * neighLoc + 1, SF2D.coeff(2 * li + 1, 2 * neigh + 1)));
			}
		}
	}
	//cout << "Makign the matrix\n";
	SF2DLoc.setFromTriplets(SFTriplet.begin(), SFTriplet.end());
	//cout << "SF2d loc done\n";
	MF2DLoc.setFromTriplets(MFTriplet.begin(), MFTriplet.end());
	//cout << "MF2d loc done\n";
	if (id == 0) visualizeSparseMatrixInMatlab(SF2DLoc);


	computeEigenMatlab(SF2DLoc, MF2DLoc, 4, EigVectLoc, eigValsLoc, "hello");
	//printf("EigVec Size=%dx%d \t\t elements=%d\n", EigVectLoc.rows(), EigVectLoc.cols(), InnerElements.size());

	//printf("[%d] Eigenvalues: %.5f, %.5f, %.5f, %.5f\n", id, eigValsLoc(0), eigValsLoc(1), eigValsLoc(2), eigValsLoc(3));
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
			EigVectTemp1(2 * k+0) = v1(0);
			EigVectTemp1(2 * k+1) = v1(1);

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

	//printf("System %d (%d) is solved.\n", id, XfLoc.rows());
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