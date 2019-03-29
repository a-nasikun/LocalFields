#include "VectorFields.h"

void VectorFields::computeAverageEdgeLength()
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
	int counter = 0, counter1 = 0;

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

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;

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
			VtoFPair vp1{ counter++, F(i, j), i };
			for (int k = j + 1; k < F.cols(); k++) {
				VtoFPair vp2{ counter++, F(i, k), i };
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

void VectorFields::constructVFAdjacency()
{
	VFAdjacency.resize(F.rows(), V.rows());
	vector<Eigen::Triplet<bool>> VFTriplet;
	VFTriplet.reserve(7 * V.rows());

	for (int i = 0; i < F.rows(); i++) {
		for (int j = 0; j < F.cols(); j++) {
			VFTriplet.push_back(Eigen::Triplet<bool>(i, F(i, j), true));
		}
	}

	VFAdjacency.setFromTriplets(VFTriplet.begin(), VFTriplet.end());
}

void VectorFields::testAdjacency()
{
	for (int i = 0; i < V.rows(); i += 100) {
		// VF_FULL
		//cout << "VF_Neigbor Full\n";
		printf("FULL_[%d] => ");
		for (std::set<VtoFPair>::iterator it1 = VFNeighFull[i].begin(); it1 != VFNeighFull[i].end(); ++it1) {
			printf("%d, ", it1->fId);
		} 
		cout << "\n";

		// VFAdjacency (boolean)
		//cout << "VF_Neigbor Full\n";
		printf("ADJ_[%d] => ");		
		for (Eigen::SparseMatrix<bool>::InnerIterator it(VFAdjacency, i); it; ++it) {
			printf("%d, ", it.row());
		}		
		cout << "\n";
	}
}

//void VectorFields::constructVFNeighborsFull()
//{
//	// For Timing
//	chrono::high_resolution_clock::time_point	t1, t2;
//	chrono::duration<double>					duration;
//	t1 = chrono::high_resolution_clock::now();
//	cout << "> Building \"Vertex-to-Face\" Adjacency (2)... ";
//
//	VFNeighFull.resize(singularities.size());
//	int counter = 0;
//
//	// Form singularities in "set" for easier member "finding/comparison"
//	set<int> singSet;
//	map<int, int> singGlobToLocMap;
//	for (int i = 0; i < singularities.size(); i++) {
//		singSet.insert(singularities[i]);
//	}
//
//	for (int i = 0; i < F.rows(); i++) {
//		for (int j = 0; j < F.cols(); j++) {
//			//VtoFPair vp1{ counter++, F(i, j), i };
//			if (singSet.find(F(i, j)) == singSet.end()) continue;
//			VtoFPair vp1{ i, F(i, j), i };
//			for (int k = j + 1; k < F.cols(); k++) {
//				//VtoFPair vp2{ counter++, F(i, k), i };
//				VtoFPair vp2{ i, F(i, k), i };
//				VFNeighFull[vp1.vId].insert(vp2);
//				VFNeighFull[vp2.vId].insert(vp1);
//			}
//		}
//	}
//
//	t2 = chrono::high_resolution_clock::now();
//	duration = t2 - t1;
//	cout << "in " << duration.count() << " seconds" << endl;
//}

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
	//constructMassMatrixMV();
	//constructMassMatrixMVinv();

	// Face-based Mass Matrices
	constructMassMatrixMF2D();
	constructMassMatrixMF2Dinv();
	//constructMassMatrixMF3D();
	//constructMassMatrixMF3Dinv();


	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "..in Total of " << duration.count() << " seconds" << endl;

	/* For the purpose of TESTING */
	//cout << "Size of mass matrices " << endl;
	//printf("____MV=%dx%d \n", MV.rows(), MV.cols());
	//printf("____MF2D=%dx%d \n", MF2D.rows(), MF2D.cols());
	//printf("____MF3D=%dx%d \n", MF3D.rows(), MF3D.cols());
	//cout << "MV \n" << MV.block(0, 0, 30, 30) << endl;
	//cout << "MF2D \n" << MF2D.block(0, 0, 30, 30) << endl;
}

void VectorFields::constructMassMatrixMV()
{
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "....Constructing vertex-based Mass matrix... ";

	//igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, MV);
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, MV);

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
		double area = doubleArea(i) / 2.0;
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
		double area = doubleArea(i) / 2.0;
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

	// Declare variables, only required as temporary structure for creating SF2D (Stiffness matrix in local frame)
	Eigen::SparseMatrix<double> LapCurl3D_NonConform, LapCurl2D_NonConform, LapDiv3D_NonConform, LapDiv2D_NonConform;
	Eigen::SparseMatrix<double> LapCurl3D_Conform, LapCurl2D_Conform, LapDiv3D_Conform, LapDiv2D_Conform;

	constructStiffnessMatrixSF3D(LapCurl3D_Conform, LapDiv3D_Conform, LapCurl3D_NonConform, LapDiv3D_NonConform);

	constructStiffnessMatrixSF2D(LapCurl3D_NonConform, LapCurl2D_NonConform, LapDiv3D_NonConform, LapDiv2D_NonConform);
	//constructStiffnessMatrixSF2D(LapCurl3D_Conform, LapCurl2D_Conform, LapDiv3D_Conform, LapDiv2D_Conform);	
	//constructStiffnessMatrixSF2D(LapCurl3D_NonConform, LapCurl2D_NonConform, LapDiv3D_Conform, LapDiv2D_Conform);
	//constructStiffnessMatrixSF2D(LapCurl3D_Conform, LapCurl2D_Conform, LapDiv3D_NonConform, LapDiv2D_NonConform);

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "..in total of" << duration.count() << " seconds" << endl;

}

void VectorFields::loadStiffnessMatrices()
{
	Eigen::SparseMatrix<double> MVerts, MEdges, MFields, G, GStar, JMat; 
	Eigen::SparseMatrix<double> MVertsInv, MEdgesInv, MFieldsInv;
	printf("> Loading matrices from Christopher's JavaView\n");

	/* File locations */
	string folder = "D:/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Local Fields/Matrices/Torus_73k/";
	//string folder = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Local Fields/Matrices/Torus_73k/";
	string file_MVerts	= folder + "M_Verts.txt";
	string file_MEdges	= folder + "M_Edges.txt";
	string file_MFields = folder + "M_Fields.txt";
	string file_G		= folder + "M_G.txt";
	string file_GStar	= folder + "M_GStar.txt";
	string file_J		= folder + "M_J.txt";

	/* Loading the files */
	cout << "MVertices \n";
	LoadSparseMatrixFromTxtFile(file_MVerts, MVerts);
	cout << "MEdges \n";
	LoadSparseMatrixFromTxtFile(file_MEdges, MEdges);
	cout << "MFields \n";
	LoadSparseMatrixFromTxtFile(file_MFields, MFields);
	cout << "Grad \n";
	LoadSparseMatrixFromTxtFile(file_G, G);
	cout << "GStar \n";
	LoadSparseMatrixFromTxtFile(file_GStar, GStar);
	cout << "JRot \n";
	LoadSparseMatrixFromTxtFile(file_J, JMat);

	/* Construct Mass Matrices */
	ConstructInverseMassMatrix(MVerts, MVertsInv);
	ConstructInverseMassMatrix(MEdges, MEdgesInv);
	ConstructInverseMassMatrix(MFields, MFieldsInv);
	MV = MVerts; 
	MVinv = MVertsInv;
	MF3D = MFields;
	MF3Dinv = MFieldsInv;


	/* For the purpose of testing */
	cout << "MVerts \n" << MVerts.block(0, 0, 30, 30) << endl;
	cout << "MEdges \n" << MEdges.block(0, 0, 30, 30) << endl;
	cout << "MFields \n" << MFields.block(0, 0, 30, 30) << endl; 
	

	/* Converting the matrices */
	cout << "Computing the SF3D\n";	
	/* Non-conforming + Non-conforming */
	//SF3D = MFields * (GStar*MEdgesInv*GStar.transpose() - JMat * GStar * MEdgesInv * GStar.transpose() * JMat) * MFields;
	SF3D = MF3D * (GStar * MEdgesInv * GStar.transpose() - JMat*GStar*MEdgesInv*GStar.transpose()*JMat) * MF3D;
	SF2D = A.transpose() * SF3D * A;	
}

void VectorFields::constructStiffnessMatrixSF2D(Eigen::SparseMatrix<double>& Matrix3D, Eigen::SparseMatrix<double>& Matrix2D)
{
	Matrix2D = A.transpose() * Matrix3D * A;
}

void VectorFields::constructStiffnessMatrixSF2D(Eigen::SparseMatrix<double>& LapCurl3D, Eigen::SparseMatrix<double>& LapCurl2D, Eigen::SparseMatrix<double>& LapDiv3D, Eigen::SparseMatrix<double>& LapDiv2D)
{
	// Explicit construction
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;

	t1 = chrono::high_resolution_clock::now();
	cout << "....Constructing Stiffness Matrix (2D) Curl part ";
	constructStiffnessMatrixCurlPart2D(LapCurl3D, LapCurl2D);
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
	

	t1 = chrono::high_resolution_clock::now();
	cout << "....Constructing Stiffness Matrix (2D) Divergent part ";
		constructStiffnessMatrixDivPart2D(LapDiv3D, LapDiv2D);
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;

	t1 = chrono::high_resolution_clock::now();
	cout << "....Divergent Part (2D) - Curl Part (2D) ";
		SF2D = LapCurl2D + LapDiv2D;
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}

void VectorFields::constructStiffnessMatrixSF3D(Eigen::SparseMatrix<double>& LapCurl3D_Conform, Eigen::SparseMatrix<double>& LapDiv3D_Conform, Eigen::SparseMatrix<double>& LapCurl3D_NonConform, Eigen::SparseMatrix<double>& LapDiv3D_NonConform)
{
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;

	t1 = chrono::high_resolution_clock::now();
	cout << "....Constructing Stiffness Matrix (3D) Curl part ";
		constructStiffnessMatrixCurlPart3D_Conform(LapCurl3D_Conform);
		constructStiffnessMatrixCurlPart3D_NonConform(LapCurl3D_NonConform);		
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;


	t1 = chrono::high_resolution_clock::now();
	cout << "....Constructing Stiffness Matrix (3D) Divergent part ";
		constructStiffnessMatrixDivPart3D_Conform(LapDiv3D_Conform);
		constructStiffnessMatrixDivPart3D_NonConform(LapDiv3D_NonConform);
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;


	//t1 = chrono::high_resolution_clock::now();
	//cout << "....Divergent Part - Curl Part ";
	//	//SF3D = LapDiv3D + LapCurl3D;
	//t2 = chrono::high_resolution_clock::now();
	//duration = t2 - t1;
	//cout << "in " << duration.count() << " seconds" << endl;
}

void VectorFields::constructStiffnessMatrixCurlPart3D_Conform(Eigen::SparseMatrix<double>& LapCurl3D_Conform)
{
	LapCurl3D_Conform = ((-MF3D * J3D) * (GF3D * MVinv)) * (GF3D.transpose() * (J3D * MF3D));
	//LapDiv3D = MF3D*(GF3D*MVinv*GF3D.transpose())*MF3D;
}

void VectorFields::constructStiffnessMatrixCurlPart3D_NonConform(Eigen::SparseMatrix<double>& LapCurl3D_NonConform)
{

	LapCurl3D_NonConform.resize(3 * F.rows(), 3 * F.rows());
	LapCurl3D_NonConform.reserve(3 * F.rows() * 4);
	vector<Eigen::Triplet<double>> LTriplet;
	LTriplet.reserve(12 * LapCurl3D_NonConform.rows());

	for (int i = 0; i < F.rows(); i++) {
		double				area1 = doubleArea(i)/2.0;
		for (int j = 0; j < F.cols(); j++) {
			int				neigh = AdjMF3N(i, j);

			if (neigh > i) {
				double				area2 = doubleArea(neigh)/2.0;
				double				area = area1 + area2;

				Eigen::Vector3d		edge = V.row(EdgePairMatrix(i, 2 * j + 1)) - V.row(EdgePairMatrix(i, 2 * j));
				Eigen::Matrix3d		block = (-3.0 / area) * edge * edge.transpose();

				/* For testing only*/
				if (i == 0) cout << "****First block: \n" << block << endl; 

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

				// ITS BLOCK => DIAGONAL
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * i + 0, -block(0, 0)));	// row 1
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * i + 1, -block(0, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * i + 2, -block(0, 2)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * i + 0, -block(1, 0)));	// row 2
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * i + 1, -block(1, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * i + 2, -block(1, 2)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * i + 0, -block(2, 0)));	// row 3
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * i + 1, -block(2, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * i + 2, -block(2, 2)));

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

				// TRANSPOSE => DIAGONAL
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 0, 3 * neigh + 0, -block(0, 0)));	// row 1
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 0, 3 * neigh + 1, -block(0, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 0, 3 * neigh + 2, -block(0, 2)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 1, 3 * neigh + 0, -block(1, 0)));	// row 2
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 1, 3 * neigh + 1, -block(1, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 1, 3 * neigh + 2, -block(1, 2)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 2, 3 * neigh + 0, -block(2, 0)));	// row 3
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 2, 3 * neigh + 1, -block(2, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 2, 3 * neigh + 2, -block(2, 2)));
			}
		}
	}
	LapCurl3D_NonConform.setFromTriplets(LTriplet.begin(), LTriplet.end());
}

void VectorFields::constructStiffnessMatrixCurlPart2D(Eigen::SparseMatrix<double>& LapCurl3D, Eigen::SparseMatrix<double>& LapCurl2D)
{
	//LapCurl2D = A.transpose() * LapCurl3D * A;
	constructStiffnessMatrixSF2D(LapCurl3D, LapCurl2D);
}

void VectorFields::constructStiffnessMatrixDivPart3D(Eigen::SparseMatrix<double>& LapDiv3D)
{
	// IMPLICIT Construction for Divergent Part
	//constructStiffnessMatrixDivPart3D_Conform(LapDiv3D);

	// EXPLICIT Construction
	constructStiffnessMatrixDivPart3D_NonConform(LapDiv3D);
}

void VectorFields::constructStiffnessMatrixDivPart3DFromCurl3D(Eigen::SparseMatrix<double>& LapCurl3D, Eigen::SparseMatrix<double>& LapDiv3D)
{
	/* Constructing the rotation matrix J on 3D/world coordinate */
	Eigen::SparseMatrix<double> J3D(3 * F.rows(), 3 * F.rows());
	vector<Eigen::Triplet<double>> JTriplet;
	JTriplet.reserve(2 * 3 * F.rows());
		
	Eigen::RowVector3d n3d; 
	for (int i = 0; i < F.rows(); i++)
	{
		n3d = NF.row(i);
		n3d.normalize();
		JTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * i + 1, -n3d(2)));
		JTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * i + 2, n3d(1)));
		JTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * i + 0, n3d(2)));
		JTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * i + 2, -n3d(0)));
		JTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * i + 0, -n3d(1)));
		JTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * i + 1, n3d(0)));
	}

	/* Building the Divergent operator from Curl Operator */
	J3D.setFromTriplets(JTriplet.begin(), JTriplet.end());
	LapDiv3D = (-J3D * LapCurl3D) *J3D;
}

void VectorFields::constructStiffnessMatrixDivPart3D_Conform(Eigen::SparseMatrix<double>& LapDiv3D)
{
	printf("GF3D=%dx%d | MF=%dx%d | MVinv=%dx%d\n", GF3D.rows(), GF3D.cols(), MF3D.rows(), MF3D.cols(), MVinv.rows(), MVinv.cols());
	LapDiv3D = MF3D*(GF3D*MVinv*GF3D.transpose())*MF3D;
}

void VectorFields::constructStiffnessMatrixDivPart3D_NonConform(Eigen::SparseMatrix<double>& LapDiv3D)
{
	LapDiv3D.resize(3 * F.rows(), 3 * F.rows());
	LapDiv3D.reserve(3 * F.rows() * 4);
	vector<Eigen::Triplet<double>> LTriplet;
	LTriplet.reserve(12 * 3 * F.rows());

	for (int i = 0; i < F.rows(); i++) {
		double				area1 = doubleArea(i)/2.0;
		Eigen::RowVector3d	n1 = NF.row(i);
		for (int j = 0; j < F.cols(); j++) {
			int				neigh = AdjMF3N(i, j);

			if (neigh > i) {
				double				area2 = doubleArea(neigh)/2.0;
				double				area = area1 + area2;
				Eigen::RowVector3d	n2 = NF.row(neigh);
				Eigen::RowVector3d	n = (n1 + n2) / 2.0;
				Eigen::Vector3d		edge = V.row(EdgePairMatrix(i, 2 * j + 1)) - V.row(EdgePairMatrix(i, 2 * j));

				edge = n.cross(edge);
				//edge.normalize();
				Eigen::Matrix3d		block = (-3.0 / area) * edge * edge.transpose();


				// ITS BLOCK ==> OFF_DIAGONAL
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * neigh + 0, block(0, 0)));	// row 1
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * neigh + 1, block(0, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * neigh + 2, block(0, 2)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * neigh + 0, block(1, 0)));	// row 2
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * neigh + 1, block(1, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * neigh + 2, block(1, 2)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * neigh + 0, block(2, 0)));	// row 3
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * neigh + 1, block(2, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * neigh + 2, block(2, 2)));

				// ITS BLOCK ==> DIAGONAL
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * i + 0, -block(0, 0)));	// row 1
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * i + 1, -block(0, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * i + 2, -block(0, 2)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * i + 0, -block(1, 0)));	// row 2
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * i + 1, -block(1, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * i + 2, -block(1, 2)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * i + 0, -block(2, 0)));	// row 3
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * i + 1, -block(2, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * i + 2, -block(2, 2)));

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

				// THE TRANSPOSE BLOCK ==> DIAGONAL
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 0, 3 * neigh + 0, -block(0, 0)));	// row 1
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 0, 3 * neigh + 1, -block(0, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 0, 3 * neigh + 2, -block(0, 2)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 1, 3 * neigh + 0, -block(1, 0)));	// row 2
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 1, 3 * neigh + 1, -block(1, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 1, 3 * neigh + 2, -block(1, 2)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 2, 3 * neigh + 0, -block(2, 0)));	// row 3
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 2, 3 * neigh + 1, -block(2, 1)));
				LTriplet.push_back(Eigen::Triplet<double>(3 * neigh + 2, 3 * neigh + 2, -block(2, 2)));
			}
		}
	}
	LapDiv3D.setFromTriplets(LTriplet.begin(), LTriplet.end());
}

void VectorFields::constructStiffnessMatrixDivPart2D(Eigen::SparseMatrix<double>& LapDiv3D, Eigen::SparseMatrix<double>& LapDiv2D)
{
	//LapDiv2D = A.transpose() * LapDiv3D * A;
	constructStiffnessMatrixSF2D(LapDiv3D, LapDiv2D);
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

void VectorFields::rearrangeGradient3D(Eigen::SparseMatrix<double> &Grad3D)
{
	//MFinv.resize(MF.rows(), MF.cols());
	vector<Eigen::Triplet<double>> GTriplet;
	GTriplet.reserve(Grad3D.nonZeros());
	int nRows = F.rows(), nCols = Grad3D.cols();

	/* Getting the sum of every non-zero elements in a row */
	for (int k = 0; k < Grad3D.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(Grad3D, k); it; ++it) {
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
	Grad3D.resize(0, 0);
	Grad3D.resize(3 * nRows, nCols);
	Grad3D.setFromTriplets(GTriplet.begin(), GTriplet.end());
}

void VectorFields::constructRotationMatrix()
{
	/* Rotation in 2D/Local Space */
	J.resize(2 * F.rows(), 2 * F.rows());
	vector<Eigen::Triplet<double>> JTriplet;
	JTriplet.reserve(2 * 2 * F.rows());
	const double cosT = 0.0, sinT = 1.0;

	for (int i = 0; i < F.rows(); i++) 
	{
		// Constructing the triplet
		JTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, cosT));		// column 1
		JTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 0, sinT));
		JTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 1, -sinT));	// column 2
		JTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, cosT));
	}

	J.setFromTriplets(JTriplet.begin(), JTriplet.end());
	//cout << "Matrix J" << endl << J.block(0, 0, 10, 10) << endl << endl;
	//visualizeSparseMatrixInMatlab(J);

	/* Rotation in 3D/Global Space */
	J3D.resize(3 * F.rows(), 3 * F.rows());
	vector<Eigen::Triplet<double>> J3DTriplet;
	J3DTriplet.reserve(3 * 3 * F.rows());

	for (int i = 0; i < F.rows(); i++) 
	{
		J3DTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * i + 1, -NF(i, 2)));
		J3DTriplet.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * i + 2,  NF(i, 1)));
		J3DTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * i + 0,  NF(i, 2)));
		J3DTriplet.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * i + 2, -NF(i, 0)));
		J3DTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * i + 0, -NF(i, 1)));
		J3DTriplet.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * i + 1,  NF(i, 0)));
	}
	J3D.setFromTriplets(J3DTriplet.begin(), J3DTriplet.end());
	//visualizeSparseMatrixInMatlab(J3D);
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

	//printf("____B2D=%dx%d [%.5f per-row (%d) filled]\n", B2D.rows(), B2D.cols(), (double)B2D.nonZeros() / (double)B2D.rows(), B2D.nonZeros());
	
	// FREE-ing Memory, not the best practice but used to save space for large mesh
	//SF2D.resize(0, 0);
}

// GETTER AND SETTER of IMPORTANT ELEMENTS
Eigen::VectorXd VectorFields::getRefFields() const
{
	return Xf; 
}