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

void LocalFields::constructMatrixBLocal(const Eigen::SparseMatrix<double>& B2D)
{
	
	//======================== BLoc from B2D =========================
	BLoc.resize(2 * LocalElements.size(), 2 * LocalElements.size());
	vector<Eigen::Triplet<double>> BTriplet;

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

void LocalFields::constructMatrixBLocal(const Eigen::SparseMatrix<double>& B2D, const vector<set<int>>& AdjMF2Ring)
{
	//======================== BLoc from B2D =========================
	BLoc.resize(2 * LocalElements.size(), 2 * LocalElements.size());
	vector<Eigen::Triplet<double>> BTriplet;

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

void LocalFields::constructLocalConstraints()
{
	// Setting up matrix C
	vector<Eigen::Triplet<double>> CTriplet;
	int counter = 0;
	CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * GlobToLocMap[sampleID] + 0, 1.0));
	CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * GlobToLocMap[sampleID] + 1, 1.0));
	for (int bound : Boundary) {
		CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * GlobToLocMap[bound] + 0, 1.0));
		CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * GlobToLocMap[bound] + 1, 1.0));
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

void LocalFields::setupLHSLocalProblemMapped()
{
	ALoc.resize(BLoc.rows() + CLoc.rows(), BLoc.cols() + CLoc.rows());

	vector<Eigen::Triplet<double>>				ATriplet;
	ATriplet.reserve(10 * BLoc.rows() + CLoc.rows()); // In practice it's only 8, but allocate 10 for safety
	//vector<Eigen::Triplet<double>>::iterator	it2;
	//it2 = ATriplet.begin();
	//it2 = std::prev(ATriplet.end());
	for (int k = 0; k < BLoc.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(BLoc, k); it; ++it) {
			ATriplet.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
			//it2 = ATriplet.insert(it2, Eigen::Triplet<double>(it.row(), it.col(), it.value()));
			//it2 = std::next(it2, 1);

			//ATriplet.insert(ATriplet.end(), Eigen::Triplet<double>(it.row(), it.col(), it.value()));
		}
	}

	for (int k = 0; k < CLoc.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(CLoc, k); it; ++it) {
			/* PUSH BACK */
			ATriplet.push_back(Eigen::Triplet<double>(BLoc.rows() + it.row(), it.col(), it.value()));
			ATriplet.push_back(Eigen::Triplet<double>(it.col(), BLoc.cols() + it.row(), it.value()));
			
			/* INSERT with NEXT Iterator*/
			//it2 = ATriplet.insert(it2, Eigen::Triplet<double>(BLoc.rows() + it.row(), it.col(), it.value()));
			//it2 = std::next(it2, 1);
			//it2 = ATriplet.insert(it2, Eigen::Triplet<double>(it.col(), BLoc.cols() + it.row(), it.value()));
			//it2 = std::next(it2, 1);

			/* INSERT at the END */
			//ATriplet.insert(ATriplet.end(), Eigen::Triplet<double>(BLoc.rows() + it.row(), it.col(), it.value()));
			//ATriplet.insert(ATriplet.end(), Eigen::Triplet<double>(it.col(), BLoc.cols() + it.row(), it.value()));
		}
	}
	ALoc.setFromTriplets(ATriplet.begin(), ATriplet.end());

	//printf("Size of A Triplet = %d\n", ATriplet.size());
	//cout << "Sample ID [" << id << "]'s LHS is constructed, with " << ALoc.rows() << "x" << ALoc.cols() << "with nonzeros= " << ALoc.nonZeros() << endl;

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
	printf("_______Error of %d is %.16f \n", id, normXf);
}