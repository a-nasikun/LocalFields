#include "VectorFields.h"

/* ====================== SETTING UP MATRICES ============================*/
void VectorFields::constructConstraints()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Constructing constraints... ";

	//construct1CentralConstraint();
	//constructRingConstraints();
	constructSpecifiedConstraints();
	
	//constructSingularities();
	//constructSpecifiedConstraintsWithSingularities();

	

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;

	// Information about constraints
	printf("....Num of Constraints = %d\n", globalConstraints.size());
	printf("....Matrix C size = %dx%d \n", C.rows(), C.cols());
	printf("....Matrix c size = %dx%d \n", c.rows(), c.cols());
}

void VectorFields::construct1CentralConstraint()
{
	vector<Eigen::Triplet<double>>	CTriplet;
	CTriplet.reserve(2);
	const int constNum = 1;
	//srand(time(NULL));
	//const int centralID = rand()%F.rows(); 
	const int centralID = *(NeighRing[0].begin());
	cout << "Face ID is : " << centralID << endl;

	// Setting up matrix C
	C.resize(2 * constNum, B2D.cols());

	CTriplet.push_back(Eigen::Triplet<double>(0, 2 * centralID + 0, 1.0));
	CTriplet.push_back(Eigen::Triplet<double>(1, 2 * centralID + 1, 1.0));

	C.setFromTriplets(CTriplet.begin(), CTriplet.end());

	// Setting up vector c (There are 2 vector c)
	c.resize(2 * constNum, 2);
	c.col(0) << 1.0, 0.0;
	c.col(1) << 0.0, 1.0;
}

void VectorFields::constructRingConstraints()
{
	// Define necessary data/variables
	vector<Eigen::Triplet<double>>	CTriplet;
	
	const int outerRingID = 9;
	const int outerBoundaryID = min(outerRingID, (int) NeighRing.size()-3);
	const int constNum = 1 + (int) NeighRing[outerBoundaryID+1].size() + (int) NeighRing[outerBoundaryID + 2].size();
	const int centralID = *(NeighRing[0].begin());
	int counter = 0;

	// Setting up matrix C
	C.resize(2 * constNum, B2D.cols());
	CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * centralID + 0, 1.0));
	CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * centralID + 1, 1.0));
	for (int i = outerBoundaryID + 1; i <= outerBoundaryID + 2; i++) {
		for (std::set<int, double>::iterator it = NeighRing[i].begin(); it != NeighRing[i].end(); ++it) {
			CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * (*it) + 0, 1.0));
			CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * (*it) + 1, 1.0));
		}
	}
	C.setFromTriplets(CTriplet.begin(), CTriplet.end());

	// Setting up vector c (There are 2 vector c)
	c.resize(2 * constNum, 2);
	Eigen::VectorXd zeroElements(2 * constNum - 2);
	for (int i = 0; i < zeroElements.size(); i++) zeroElements(i) = 0.0;
	c.col(0) << 1.0, 0.0, zeroElements;
	c.col(1) << 0.0, 1.0, zeroElements;
}

void VectorFields::constructSpecifiedConstraints()
{
	// Define the constraints
	const int numConstraints = 50;
	set<int> constraints;
	//vector<int> globalConstraints(numConstraints);
	globalConstraints.resize(numConstraints);
	Eigen::VectorXd D;
	D.resize(F.rows());

	// Initialize the value of D
	for (int i = 0; i < F.rows(); i++) {
		D(i) = numeric_limits<double>::infinity();
	}

	srand(time(NULL));
	//int curPoint = rand() % F.rows();
	int curPoint = 0; 
	constraints.insert(curPoint);

	// Creating constraints using farthest point sampling
	do {
		Eigen::VectorXi::Index maxIndex;
		computeDijkstraDistanceFaceForSampling(curPoint, D);
		D.maxCoeff(&maxIndex);
		constraints.insert(maxIndex);
		curPoint = maxIndex;
	} while (constraints.size() < numConstraints);

	int counter1 = 0;
	for (int i : constraints) {
		globalConstraints[counter1++] = i;
	}
	//printf("Constraints = %d\n", globalConstraints.size());

	// Setting up matrix C
	Eigen::SparseMatrix<double> CTemp;
	vector<Eigen::Triplet<double>> CTriplet;
	CTriplet.reserve(2 * globalConstraints.size());

	int counter = 0;
	for (int i = 0; i < globalConstraints.size(); i++) {
		CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * globalConstraints[i] + 0, 1.0));
		CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * globalConstraints[i] + 1, 1.0));
	}
	C.resize(2 * globalConstraints.size(), B2D.rows());
	C.setFromTriplets(CTriplet.begin(), CTriplet.end());
	//printf("Cp=%dx%d\n", C.rows(), C.cols());


	// Setting up vector c (There are 1 vector c)
	srand(time(NULL));
	c.resize(2 * globalConstraints.size());
	for (int i = 0; i < globalConstraints.size(); i++) {
		c(2 * i + 0, 0) = sqrt(2.0);
		c(2 * i + 1, 0) = sqrt(2.0);
	}
	//printf("cBar=%dx%d\n", c.rows(), c.cols());
}

void VectorFields::constructSingularities()
{
	const int NUM_SINGS = 2;

	if (NUM_SINGS > 0)
		constructVFAdjacency();
		//constructVFNeighborsFull();

	singularities.resize(NUM_SINGS);
	SingNeighCC.resize(NUM_SINGS);

	// For testing
	sharedEdgesVect.resize(NUM_SINGS);

	//time_t t;
	//srand((unsigned)time(&t));
	srand(time(NULL));
	for (int id = 0; id < NUM_SINGS; id++) {
		// Defining varaibles for singularities
		const int SingLocation	= rand() % V.rows();
		singularities[id]		= SingLocation;
		const int SingNeighNum	= VFAdjacency.col(SingLocation).nonZeros(); 
		Eigen::SparseMatrix<bool>::InnerIterator it0(VFAdjacency, SingLocation);
		const int firstNeigh	= it0.row(); 


		// Inserting the first neighbor (the one with lowest index/row number)
		SingNeighCC[id].resize(SingNeighNum);
		SingNeighCC[id][0]		= firstNeigh;
		int curNeigh			= firstNeigh;
		int vertex1				= SingLocation;

		// Getting the neighboring valence triangles in order
		for (int i2 = 1; i2<SingNeighNum; i2++) {
			int vertex2;
			// Setting the vertex on the edge pointing to vertex1 as vertex2 (edge = v2->v1)
			for (int i = 0; i < F.cols(); i++) {
				if (F(curNeigh, i%F.cols()) == vertex1) {
					vertex2 = F(curNeigh, (i + F.cols() - 1) % F.cols());
				}
			}
			//for (std::set<VtoFPair>::iterator it1 = next(VFNeighFull[SingLocation].begin(), 1); it1 != VFNeighFull[SingLocation].end(); ++it1) {
			// Getting the neighboring triangles in order (CCW) 
			for (Eigen::SparseMatrix<bool>::InnerIterator it1(VFAdjacency, SingLocation); it1; ++it1) {
				for (int i = 0; i < F.cols(); i++) {
					if (F(it1.row(), i) == vertex1 && F(it1.row(), (i + 1) % F.cols()) == vertex2) {
						SingNeighCC[id][i2] = it1.row();
						curNeigh = it1.row();
					}
				}
			}
		}
	}

	// Free Memory for VFNeighborsFull
	VFAdjacency.resize(0, 0);
	//VFNeighFull.clear();
	//VFNeighbors.shrink_to_fit();
}

void VectorFields::constructSpecifiedConstraintsWithSingularities()
{
	// Define the constraints
	const int numConstraints = 4;
	set<int> constraints;

	globalConstraints.resize(numConstraints);
	Eigen::VectorXd D;
	D.resize(F.rows());

	// Initialize the value of D
	for (int i = 0; i < F.rows(); i++) {
		D(i) = numeric_limits<double>::infinity();
	}

	srand(time(NULL));
	int curPoint = rand() % F.rows();
	constraints.insert(curPoint);

	// Creating constraints using farthest point sampling
	do {
		Eigen::VectorXi::Index maxIndex;
		computeDijkstraDistanceFaceForSampling(curPoint, D);
		D.maxCoeff(&maxIndex);
		constraints.insert(maxIndex);
		curPoint = maxIndex;
	} while (constraints.size() < numConstraints);

	int counter1 = 0;
	for (int i : constraints) {
		globalConstraints[counter1++] = i;
	}
	
	int numSingConstraints = 0;
	for (int i = 0; i < SingNeighCC.size(); i++) {
		// Use only n-1 neighboring faces as constraints
		for (int j = 0; j < (SingNeighCC[i].size()-1); j++) {
			numSingConstraints++;
		}
	}

	// Setting up matrix C and vector c
	c.resize(2 * (globalConstraints.size()+numSingConstraints));


	// HARD CONSTRAINTS
	Eigen::SparseMatrix<double> CTemp;
	vector<Eigen::Triplet<double>> CTriplet;
	CTriplet.reserve(2 * globalConstraints.size() +  2 * 4 * 7 * SingNeighCC.size());
	int counter = 0;
	for (int i = 0; i < globalConstraints.size(); i++) {
		// Matrix C
		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * globalConstraints[i] + 0, 1.0));
		c(counter++, 0) = sqrt(2.0);
		//c(counter++, 1) = sqrt(2.0);
		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * globalConstraints[i] + 1, 1.0));
		c(counter++, 0) = sqrt(2.0);
		//c(counter++, 1) = sqrt(2.0);
	}

	
	// SINGULARITIES CONSTRAINTS
	for (int id = 0; id < SingNeighCC.size(); id++) {
		// For testing
		sharedEdgesVect[id].resize(2*SingNeighCC[id].size()-2);

		// 4. Compute rotation of among its valence
		const double rotAngle = 2 * M_PI / (double)SingNeighCC[id].size();
		const double cosA = cos(rotAngle);
		const double sinA = sin(rotAngle);

		// Which case? => determining which edge is the common edge
		enum class SharedEdgeCase {Case1, Case2, Case3};
		SharedEdgeCase edgeCase1, edgeCase2; 
		for (int i = 0; i < (SingNeighCC[id].size() - 1); i++) {
			// 1. Find shared edge (naively)
			Eigen::RowVector3d es;
			for (int f1 = 0; f1 < F.cols(); f1++) {
				for (int f2 = 0; f2 < F.cols(); f2++) {
					bool b1 = F(SingNeighCC[id][i], (f1+1)%F.cols()) == F(SingNeighCC[id][i + 1], f2);
					bool b2 = F(SingNeighCC[id][i], f1) == F(SingNeighCC[id][i + 1], (f2 + 1) % F.cols());
					if (b1 && b2) {
						sharedEdgesVect[id][2 * i + 0] = F(SingNeighCC[id][i], f1);
						sharedEdgesVect[id][2 * i + 1] = F(SingNeighCC[id][i], (f1 + 1) % F.cols());
						es = V.row(F(SingNeighCC[id][i], (f1 + 1) % F.cols())) - V.row(F(SingNeighCC[id][i], f1));
						printf("Shared edge=%d->%d\n", F(SingNeighCC[id][i], f1), F(SingNeighCC[id][i], (f1 + 1) % F.cols()));
						
						if (f1 == 0)		edgeCase1 = SharedEdgeCase::Case1;	// => edge V0->V1 is the shared edge => it takes 0 step to reach v0
						else if(f1==1)		edgeCase1 = SharedEdgeCase::Case3;	// => edge V1->V2 is the shared edge => it takes 2 step to reach v0
						else if(f1==2)		edgeCase1 = SharedEdgeCase::Case2;	// => edge V2->V0 is the shared edge => it takes 1 step to reach v0

						if (f2 == 0)		edgeCase2 = SharedEdgeCase::Case1;
						else if (f2 == 1)	edgeCase2 = SharedEdgeCase::Case3;
						else if (f2 == 2)	edgeCase2 = SharedEdgeCase::Case2;
					}
				}
			}
			// 2. Find angles between basis1 and es
			Eigen::VectorXd eVect;
			Eigen::RowVector3d b11, b12;
			eVect = A.block(3 * SingNeighCC[id][i], 2 * SingNeighCC[id][i]+0, 3, 1);
			b11 << eVect(0), eVect(1), eVect(2);
			eVect = A.block(3 * SingNeighCC[id][i], 2 * SingNeighCC[id][i]+1, 3, 1);
			b12 << eVect(0), eVect(1), eVect(2);
			cout << "______B11: " << b11 << ", B12: " << b12 << endl;

			// Basis 1, Frame 1
			double cosR12 = (b11.dot(es)) / (b11.norm()*es.norm());
			if (cosR12 > 1.0) cosR12 = 1.0; 
			if (cosR12 <-1.0) cosR12 = -1.0;
			const double angleR12_1 = (edgeCase1==SharedEdgeCase::Case2 ? 2*M_PI - acos(cosR12) : acos(cosR12));
			const double cosR12_1 = cos(angleR12_1);
			const double sinR12_1 = sin(angleR12_1);
			printf("______[%.2f] Rotation matrix R12_1 = [%.3f,%.3f; %.3f, %.3f]\n", angleR12_1*180.0/M_PI, cosR12_1, -sinR12_1, sinR12_1, cosR12_1);

			// Basis 1, Frame 2
			cosR12 = (b12.dot(es)) / (b12.norm()*es.norm());
			if (cosR12 > 1.0) cosR12 = 1.0;
			if (cosR12 <-1.0) cosR12 = -1.0;
			const double angleR12_2 = (edgeCase1 == SharedEdgeCase::Case1 ? 2 * M_PI - acos(cosR12) : acos(cosR12));
			const double cosR12_2 = cos(angleR12_2);
			const double sinR12_2 = sin(angleR12_2);
			printf("______[%.2f] Rotation matrix R12_2 = [%.2f,%.2f; %.2f, %.2f]\n", angleR12_2*180.0 / M_PI, cosR12_2, -sinR12_2, sinR12_2, cosR12_2);

			// 3. Find angles between basis2 and es
			es = -es; 
			Eigen::RowVector3d b21, b22;
			eVect = A.block(3 * SingNeighCC[id][i+1], 2 * SingNeighCC[id][i+1] + 0, 3, 1);
			b21 << eVect(0), eVect(1), eVect(2);
			eVect = A.block(3 * SingNeighCC[id][i+1], 2 * SingNeighCC[id][i+1] + 1, 3, 1);
			b22 << eVect(0), eVect(1), eVect(2);
			cout << "______B21: " << b21 << ", B22: " << b22 << endl;
			
			// Basis 2, Frame 1
			double cosR21 = (b21.dot(es)) / (b21.norm()*es.norm());
			if (cosR21 > 1.0) cosR21 = 1.0;
			if (cosR21 < -1.0) cosR21 = -1.0;
			const double angleR21_1 = (edgeCase2 == SharedEdgeCase::Case2 ? 2 * M_PI - acos(cosR21) : acos(cosR21));
			const double cosR21_1 = cos(angleR21_1);
			const double sinR21_1 = sin(angleR21_1);
			printf("______[%.2f] Rotation matrix R22_1 = [%.2f,%.2f; %.2f, %.2f]\n", angleR21_1*180.0/M_PI, cosR21_1, -sinR21_1, sinR21_1, cosR21_1);

			// Basis 2, Frame 2
			cosR21 = (b22.dot(es)) / (b22.norm()*es.norm());
			if (cosR21 > 1.0) cosR21 = 1.0;
			if (cosR21 < -1.0) cosR21 = -1.0;
			const double angleR21_2 = (edgeCase2 == SharedEdgeCase::Case1 ? 2 * M_PI - acos(cosR21) : acos(cosR21));
			const double cosR21_2 = cos(angleR21_2);
			const double sinR21_2 = sin(angleR21_2);
			printf("______[%.2f] Rotation matrix R22_2 = [%.2f,%.2f; %.2f, %.2f]\n", angleR21_2*180.0/M_PI, cosR21_2, -sinR21_2, sinR21_2, cosR21_2);

			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, cosR12_1));
			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, -sinR12_1));
			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, -cosR21_1));
			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, sinR21_1));
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, cosR21_1*cosR12_1+sinR21_1*sinR12_1));
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, -cosR21_1*sinR12_1+sinR21_1*cosR12_1));
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, -1.0));
			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, 0.0));
			c(counter) = 0.0;
			counter++;
			
			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, sinR12_1));
			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, cosR12_1 ));
			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, -sinR21_1));
			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, -cosR21_1 ));
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, -sinR21_1*cosR12_1+cosR21_1*sinR12_1));
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, sinR21_1*sinR12_1+cosR21_1*cosR21_1));
			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, 0.0));
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, -1.0));
			c(counter) = 0.0;
			counter++;
		}
	}

	C.resize(2 * (globalConstraints.size()+numSingConstraints), B2D.rows());
	C.setFromTriplets(CTriplet.begin(), CTriplet.end());
	//printf("Cp=%dx%d\n", C.rows(), C.cols());	
}

//void VectorFields::constructSpecifiedConstraintsWithSingularities() ==> VERSION 1.1
//{
//	// Define the constraints
//	const int numConstraints = 4;
//	set<int> constraints;
//
//	globalConstraints.resize(numConstraints);
//	Eigen::VectorXd D;
//	D.resize(F.rows());
//
//	// Initialize the value of D
//	for (int i = 0; i < F.rows(); i++) {
//		D(i) = numeric_limits<double>::infinity();
//	}
//
//	srand(time(NULL));
//	int curPoint = rand() % F.rows();
//	constraints.insert(curPoint);
//
//	// Creating constraints using farthest point sampling
//	do {
//		Eigen::VectorXi::Index maxIndex;
//		computeDijkstraDistanceFaceForSampling(curPoint, D);
//		D.maxCoeff(&maxIndex);
//		constraints.insert(maxIndex);
//		curPoint = maxIndex;
//	} while (constraints.size() < numConstraints);
//
//	int counter1 = 0;
//	for (int i : constraints) {
//		globalConstraints[counter1++] = i;
//	}
//	//printf("Constraints = %d\n", globalConstraints.size());
//
//	//constructSingularities();
//
//	int numSingConstraints = 0;
//	for (int i = 0; i < SingNeighCC.size(); i++) {
//		// Use only n-1 neighboring faces as constraints
//		for (int j = 0; j < (SingNeighCC[i].size() - 1); j++) {
//			numSingConstraints++;
//		}
//	}
//
//	// Setting up matrix C and vector c
//	c.resize(2 * (globalConstraints.size() + 2 * numSingConstraints), 2);
//	//printf("cBar=%dx%d\n", c.rows(), c.cols());
//
//
//	// HARD CONSTRAINTS
//	Eigen::SparseMatrix<double> CTemp;
//	vector<Eigen::Triplet<double>> CTriplet;
//	CTriplet.reserve(2 * globalConstraints.size() + 2 * 2 * 4 * 7 * SingNeighCC.size());
//	int counter = 0;
//	for (int i = 0; i < globalConstraints.size(); i++) {
//		// Matrix C
//		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * globalConstraints[i] + 0, 1.0));
//		c(counter, 0) = sqrt(2.0);
//		c(counter++, 1) = sqrt(2.0);
//		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * globalConstraints[i] + 1, 1.0));
//		c(counter, 0) = sqrt(2.0);
//		c(counter++, 1) = sqrt(2.0);
//	}
//
//
//	// SINGULARITIES CONSTRAINTS
//	for (int id = 0; id < SingNeighCC.size(); id++) {
//		// For testing
//		sharedEdgesVect[id].resize(2 * SingNeighCC[id].size() - 2);
//
//		// 4. Compute rotation of among its valence
//		const double rotAngle = 2 * M_PI / (double)SingNeighCC[id].size();
//		const double cosA = cos(rotAngle);
//		const double sinA = sin(rotAngle);
//		// Which case?
//		enum class SharedEdgeCase { Case1, Case2, Case3 };
//		SharedEdgeCase edgeCase1, edgeCase2;
//		for (int i = 0; i < (SingNeighCC[id].size() - 1); i++) {
//			// 1. Find shared edge (naively)
//			Eigen::RowVector3d es;
//			for (int f1 = 0; f1 < F.cols(); f1++) {
//				for (int f2 = 0; f2 < F.cols(); f2++) {
//					bool b1 = F(SingNeighCC[id][i], (f1 + 1) % F.cols()) == F(SingNeighCC[id][i + 1], f2);
//					bool b2 = F(SingNeighCC[id][i], f1) == F(SingNeighCC[id][i + 1], (f2 + 1) % F.cols());
//					if (b1 && b2) {
//						sharedEdgesVect[id][2 * i + 0] = F(SingNeighCC[id][i], f1);
//						sharedEdgesVect[id][2 * i + 1] = F(SingNeighCC[id][i], (f1 + 1) % F.cols());
//						es = V.row(F(SingNeighCC[id][i], (f1 + 1) % F.cols())) - V.row(F(SingNeighCC[id][i], f1));
//						printf("Shared edge=%d->%d\n", F(SingNeighCC[id][i], f1), F(SingNeighCC[id][i], (f1 + 1) % F.cols()));
//
//						if (f1 == 0) edgeCase1 = SharedEdgeCase::Case1;
//						else if (f1 == 1) edgeCase1 = SharedEdgeCase::Case3;
//						else if (f1 == 2) edgeCase1 = SharedEdgeCase::Case2;
//
//						if (f2 == 0) edgeCase2 = SharedEdgeCase::Case1;
//						else if (f2 == 1) edgeCase2 = SharedEdgeCase::Case3;
//						else if (f2 == 2) edgeCase2 = SharedEdgeCase::Case2;
//					}
//				}
//			}
//			// 2. Find angles between basis1 and es
//			Eigen::VectorXd eVect;
//			Eigen::RowVector3d b11, b12;
//			eVect = A.block(3 * SingNeighCC[id][i], 2 * SingNeighCC[id][i] + 0, 3, 1);
//			b11 << eVect(0), eVect(1), eVect(2);
//			eVect = A.block(3 * SingNeighCC[id][i], 2 * SingNeighCC[id][i] + 1, 3, 1);
//			b12 << eVect(0), eVect(1), eVect(2);
//			cout << "______B11: " << b11 << ", B12: " << b12 << endl;
//
//			// Basis 1, Frame 1
//			double cosR12 = (b11.dot(es)) / (b11.norm()*es.norm());
//			if (cosR12 > 1.0) cosR12 = 1.0;
//			if (cosR12 <-1.0) cosR12 = -1.0;
//			const double angleR12_1 = (edgeCase1 == SharedEdgeCase::Case2 ? 2 * M_PI - acos(cosR12) : acos(cosR12));
//			const double cosR12_1 = cos(angleR12_1);
//			const double sinR12_1 = sin(angleR12_1);
//			printf("______[%.2f] Rotation matrix R12_1 = [%.3f,%.3f; %.3f, %.3f]\n", angleR12_1*180.0 / M_PI, cosR12_1, -sinR12_1, sinR12_1, cosR12_1);
//
//			// Basis 1, Frame 2
//			cosR12 = (b12.dot(es)) / (b12.norm()*es.norm());
//			if (cosR12 > 1.0) cosR12 = 1.0;
//			if (cosR12 <-1.0) cosR12 = -1.0;
//			const double angleR12_2 = (edgeCase1 == SharedEdgeCase::Case1 ? 2 * M_PI - acos(cosR12) : acos(cosR12));
//			const double cosR12_2 = cos(angleR12_2);
//			const double sinR12_2 = sin(angleR12_2);
//			printf("______[%.2f] Rotation matrix R12_2 = [%.2f,%.2f; %.2f, %.2f]\n", angleR12_2*180.0 / M_PI, cosR12_2, -sinR12_2, sinR12_2, cosR12_2);
//
//			// 3. Find angles between basis2 and es
//			es = -es;
//			Eigen::RowVector3d b21, b22;
//			eVect = A.block(3 * SingNeighCC[id][i + 1], 2 * SingNeighCC[id][i + 1] + 0, 3, 1);
//			b21 << eVect(0), eVect(1), eVect(2);
//			eVect = A.block(3 * SingNeighCC[id][i + 1], 2 * SingNeighCC[id][i + 1] + 1, 3, 1);
//			b22 << eVect(0), eVect(1), eVect(2);
//			cout << "______B21: " << b21 << ", B22: " << b22 << endl;
//
//			// Basis 2, Frame 1
//			double cosR21 = (b21.dot(es)) / (b21.norm()*es.norm());
//			if (cosR21 > 1.0) cosR21 = 1.0;
//			if (cosR21 < -1.0) cosR21 = -1.0;
//			const double angleR21_1 = (edgeCase2 == SharedEdgeCase::Case2 ? 2 * M_PI - acos(cosR21) : acos(cosR21));
//			const double cosR21_1 = cos(angleR21_1);
//			const double sinR21_1 = sin(angleR21_1);
//			printf("______[%.2f] Rotation matrix R22_1 = [%.2f,%.2f; %.2f, %.2f]\n", angleR21_1*180.0 / M_PI, cosR21_1, -sinR21_1, sinR21_1, cosR21_1);
//
//			// Basis 2, Frame 2
//			cosR21 = (b22.dot(es)) / (b22.norm()*es.norm());
//			if (cosR21 > 1.0) cosR21 = 1.0;
//			if (cosR21 < -1.0) cosR21 = -1.0;
//			const double angleR21_2 = (edgeCase2 == SharedEdgeCase::Case1 ? 2 * M_PI - acos(cosR21) : acos(cosR21));
//			const double cosR21_2 = cos(angleR21_2);
//			const double sinR21_2 = sin(angleR21_2);
//			printf("______[%.2f] Rotation matrix R22_2 = [%.2f,%.2f; %.2f, %.2f]\n", angleR21_2*180.0 / M_PI, cosR21_2, -sinR21_2, sinR21_2, cosR21_2);
//
//
//			// 5. Assigning singularities constraints
//			// Basis 1 => Frame 1
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, cosA*cosR12_1 - sinA*sinR12_1));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, -cosA*sinR12_1 - sinA*cosR12_1));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, -cosR21_1));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, sinR21_1));
//			c(counter, 0) = 0.0;
//			c(counter, 1) = 0.0;
//			counter++;
//			// Basis 1 => Frame 2
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, cosA*cosR12_2 - sinA*sinR12_2));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, -cosA*sinR12_2 - sinA*cosR12_2));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, -cosR21_2));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, sinR21_2));
//			c(counter, 0) = 0.0;
//			c(counter, 1) = 0.0;
//			counter++;
//
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, sinA*cosR12_1 + cosA*sinR12_1));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, -sinA*sinR12_1 + cosA*cosR12_1));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, -sinR21_1));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, -cosR21_1));
//			c(counter, 0) = 0.0;
//			c(counter, 1) = 0.0;
//			counter++;
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, sinA*cosR12_2 + cosA*sinR12_2));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, -sinA*sinR12_2 + cosA*cosR12_2));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, -sinR21_2));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, -cosR21_2));
//			c(counter, 0) = 0.0;
//			c(counter, 1) = 0.0;
//			counter++;
//
//			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, cosA));
//			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, -sinA));
//			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, -1.0));
//			//c(counter, 0) = 0.0;
//			//c(counter, 1) = 0.0;
//			//counter++;
//			//
//			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, sinA));
//			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, cosA));
//			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, -1.0));
//			//c(counter, 0) = 0.0;
//			//c(counter, 1) = 0.0;
//			//counter++;
//		}
//	}
//
//	C.resize(2 * (globalConstraints.size() + 2 * numSingConstraints), B2D.rows());
//	C.setFromTriplets(CTriplet.begin(), CTriplet.end());
//	//printf("Cp=%dx%d\n", C.rows(), C.cols());	
//}

//
//void VectorFields::constructSpecifiedConstraintsWithSingularities() ==> VERSION 1.0
//{
//	// Define the constraints
//	const int numConstraints = 5;
//	set<int> constraints;
//
//	globalConstraints.resize(numConstraints);
//	Eigen::VectorXd D;
//	D.resize(F.rows());
//
//	// Initialize the value of D
//	for (int i = 0; i < F.rows(); i++) {
//		D(i) = numeric_limits<double>::infinity();
//	}
//
//	srand(time(NULL));
//	int curPoint = rand() % F.rows();
//	constraints.insert(curPoint);
//
//	// Creating constraints using farthest point sampling
//	do {
//		Eigen::VectorXi::Index maxIndex;
//		computeDijkstraDistanceFaceForSampling(curPoint, D);
//		D.maxCoeff(&maxIndex);
//		constraints.insert(maxIndex);
//		curPoint = maxIndex;
//	} while (constraints.size() <= numConstraints);
//
//	int counter1 = 0;
//	for (int i : constraints) {
//		globalConstraints[counter1++] = i;
//	}
//	//printf("Constraints = %d\n", globalConstraints.size());
//
//	//constructSingularities();
//
//	int numSingConstraints = 0;
//	for (int i = 0; i < SingNeighCC.size(); i++) {
//		// Use only n-1 neighboring faces as constraints
//		for (int j = 0; j < (SingNeighCC[i].size() - 1); j++) {
//			numSingConstraints++;
//		}
//	}
//
//	// Setting up matrix C and vector c
//	c.resize(2 * (globalConstraints.size() + 2 * numSingConstraints), 2);
//	//printf("cBar=%dx%d\n", c.rows(), c.cols());
//
//
//	// HARD CONSTRAINTS
//	Eigen::SparseMatrix<double> CTemp;
//	vector<Eigen::Triplet<double>> CTriplet;
//	CTriplet.reserve(2 * globalConstraints.size() + 2 * 2 * 4 * 7 * SingNeighCC.size());
//	int counter = 0;
//	for (int i = 0; i < globalConstraints.size(); i++) {
//		// Matrix C
//		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * globalConstraints[i] + 0, 1.0));
//		c(counter, 0) = sqrt(2.0);
//		c(counter++, 1) = sqrt(2.0);
//		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * globalConstraints[i] + 1, 1.0));
//		c(counter, 0) = sqrt(2.0);
//		c(counter++, 1) = sqrt(2.0);
//	}
//
//
//	// SINGULARITIES CONSTRAINTS
//	for (int id = 0; id < SingNeighCC.size(); id++) {
//		// For testing
//		sharedEdgesVect[id].resize(2 * SingNeighCC[id].size() - 2);
//
//		// 4. Compute rotation of among its valence
//		const double rotAngle = 2 * M_PI / (double)SingNeighCC[id].size();
//		const double cosA = cos(rotAngle);
//		const double sinA = sin(rotAngle);
//		// Which case?
//		enum class SharedEdgeCase { Case1, Case2, Case3 };
//		SharedEdgeCase edgeCase1, edgeCase2;
//		for (int i = 0; i < (SingNeighCC[id].size() - 1); i++) {
//			// 1. Find shared edge (naively)
//			Eigen::RowVector3d es;
//			for (int f1 = 0; f1 < F.cols(); f1++) {
//				for (int f2 = 0; f2 < F.cols(); f2++) {
//					bool b1 = F(SingNeighCC[id][i], (f1 + 1) % F.cols()) == F(SingNeighCC[id][i + 1], f2);
//					bool b2 = F(SingNeighCC[id][i], f1) == F(SingNeighCC[id][i + 1], (f2 + 1) % F.cols());
//					if (b1 && b2) {
//						sharedEdgesVect[id][2 * i + 0] = F(SingNeighCC[id][i], f1);
//						sharedEdgesVect[id][2 * i + 1] = F(SingNeighCC[id][i], (f1 + 1) % F.cols());
//						es = V.row(F(SingNeighCC[id][i], (f1 + 1) % F.cols())) - V.row(F(SingNeighCC[id][i], f1));
//						printf("Shared edge=%d->%d\n", F(SingNeighCC[id][i], f1), F(SingNeighCC[id][i], (f1 + 1) % F.cols()));
//
//						if (f1 == 0) edgeCase1 = SharedEdgeCase::Case1;
//						else if (f1 == 1) edgeCase1 = SharedEdgeCase::Case3;
//						else if (f1 == 2) edgeCase1 = SharedEdgeCase::Case2;
//
//						if (f2 == 0) edgeCase2 = SharedEdgeCase::Case1;
//						else if (f2 == 1) edgeCase2 = SharedEdgeCase::Case3;
//						else if (f2 == 2) edgeCase2 = SharedEdgeCase::Case2;
//					}
//				}
//			}
//			// 2. Find angles between basis1 and es
//			Eigen::VectorXd eVect;
//			Eigen::RowVector3d b11, b12;
//			eVect = A.block(3 * SingNeighCC[id][i], 2 * SingNeighCC[id][i] + 0, 3, 1);
//			b11 << eVect(0), eVect(1), eVect(2);
//			eVect = A.block(3 * SingNeighCC[id][i], 2 * SingNeighCC[id][i] + 1, 3, 1);
//			b12 << eVect(0), eVect(1), eVect(2);
//			cout << "______B11: " << b11 << ", B12: " << b12 << endl;
//
//			// Basis 1, Frame 1
//			double cosR12 = (b11.dot(es)) / (b11.norm()*es.norm());
//			if (cosR12 > 1.0) cosR12 = 1.0;
//			if (cosR12 <-1.0) cosR12 = -1.0;
//			const double angleR12_1 = (edgeCase1 == SharedEdgeCase::Case2 ? 2 * M_PI - acos(cosR12) : acos(cosR12));
//			const double cosR12_1 = cos(angleR12_1);
//			const double sinR12_1 = sin(angleR12_1);
//			printf("______[%.2f] Rotation matrix R12_1 = [%.3f,%.3f; %.3f, %.3f]\n", angleR12_1*180.0 / M_PI, cosR12_1, -sinR12_1, sinR12_1, cosR12_1);
//
//			// Basis 1, Frame 2
//			cosR12 = (b12.dot(es)) / (b12.norm()*es.norm());
//			if (cosR12 > 1.0) cosR12 = 1.0;
//			if (cosR12 <-1.0) cosR12 = -1.0;
//			const double angleR12_2 = (edgeCase1 == SharedEdgeCase::Case1 ? 2 * M_PI - acos(cosR12) : acos(cosR12));
//			const double cosR12_2 = cos(angleR12_2);
//			const double sinR12_2 = sin(angleR12_2);
//			printf("______[%.2f] Rotation matrix R12_2 = [%.2f,%.2f; %.2f, %.2f]\n", angleR12_2*180.0 / M_PI, cosR12_2, -sinR12_2, sinR12_2, cosR12_2);
//
//			// 3. Find angles between basis2 and es
//			es = -es;
//			Eigen::RowVector3d b21, b22;
//			eVect = A.block(3 * SingNeighCC[id][i + 1], 2 * SingNeighCC[id][i + 1] + 0, 3, 1);
//			b21 << eVect(0), eVect(1), eVect(2);
//			eVect = A.block(3 * SingNeighCC[id][i + 1], 2 * SingNeighCC[id][i + 1] + 1, 3, 1);
//			b22 << eVect(0), eVect(1), eVect(2);
//			cout << "______B21: " << b21 << ", B22: " << b22 << endl;
//
//			// Basis 2, Frame 1
//			double cosR21 = (b21.dot(es)) / (b21.norm()*es.norm());
//			if (cosR21 > 1.0) cosR21 = 1.0;
//			if (cosR21 < -1.0) cosR21 = -1.0;
//			const double angleR21_1 = (edgeCase2 == SharedEdgeCase::Case2 ? 2 * M_PI - acos(cosR21) : acos(cosR21));
//			const double cosR21_1 = cos(angleR21_1);
//			const double sinR21_1 = sin(angleR21_1);
//			printf("______[%.2f] Rotation matrix R22_1 = [%.2f,%.2f; %.2f, %.2f]\n", angleR21_1*180.0 / M_PI, cosR21_1, -sinR21_1, sinR21_1, cosR21_1);
//
//			// Basis 2, Frame 2
//			cosR21 = (b22.dot(es)) / (b22.norm()*es.norm());
//			if (cosR21 > 1.0) cosR21 = 1.0;
//			if (cosR21 < -1.0) cosR21 = -1.0;
//			const double angleR21_2 = (edgeCase2 == SharedEdgeCase::Case1 ? 2 * M_PI - acos(cosR21) : acos(cosR21));
//			const double cosR21_2 = cos(angleR21_2);
//			const double sinR21_2 = sin(angleR21_2);
//			printf("______[%.2f] Rotation matrix R22_2 = [%.2f,%.2f; %.2f, %.2f]\n", angleR21_2*180.0 / M_PI, cosR21_2, -sinR21_2, sinR21_2, cosR21_2);
//
//
//			// 5. Assigning singularities constraints
//			// Basis 1 => Frame 1
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, cosA*cosR12_1 - sinA*sinR12_1));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, -cosA*sinR12_1 - sinA*cosR12_1));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, -cosR21_1));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, sinR21_1));
//			c(counter, 0) = 0.0;
//			c(counter, 1) = 0.0;
//			counter++;
//			// Basis 1 => Frame 2
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, cosA*cosR12_2 - sinA*sinR12_2));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, -cosA*sinR12_2 - sinA*cosR12_2));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, -cosR21_2));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, sinR21_2));
//			c(counter, 0) = 0.0;
//			c(counter, 1) = 0.0;
//			counter++;
//
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, sinA*cosR12_1 + cosA*sinR12_1));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, -sinA*sinR12_1 + cosA*cosR12_1));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, -sinR21_1));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, -cosR21_1));
//			c(counter, 0) = 0.0;
//			c(counter, 1) = 0.0;
//			counter++;
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, sinA*cosR12_2 + cosA*sinR12_2));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, -sinA*sinR12_2 + cosA*cosR12_2));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, -sinR21_2));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, -cosR21_2));
//			c(counter, 0) = 0.0;
//			c(counter, 1) = 0.0;
//			counter++;
//
//			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, cosA));
//			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, -sinA));
//			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, -1.0));
//			//c(counter, 0) = 0.0;
//			//c(counter, 1) = 0.0;
//			//counter++;
//			//
//			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, sinA));
//			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, cosA));
//			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, -1.0));
//			//c(counter, 0) = 0.0;
//			//c(counter, 1) = 0.0;
//			//counter++;
//		}
//	}
//
//	C.resize(2 * (globalConstraints.size() + 2 * numSingConstraints), B2D.rows());
//	C.setFromTriplets(CTriplet.begin(), CTriplet.end());
//	//printf("Cp=%dx%d\n", C.rows(), C.cols());	
//}
//
void VectorFields::setupGlobalProblem()
{	
	Eigen::VectorXd					b, g, h, vEst;
	Eigen::SparseMatrix<double>		A_LHS;
	//Eigen::VectorXd					vEst;

	//constructConstraints();	
	setupRHSGlobalProblemMapped(g, h, vEst, b);
	setupLHSGlobalProblemMapped(A_LHS);
	//setupRHSGlobalProblem();
	//setupLHSGlobalProblem();	

	//solveGlobalSystem();
	solveGlobalSystemMappedLDLT(vEst, A_LHS, b);
	//solveGlobalSystemMappedLU_GPU();
}

void VectorFields::setupRHSGlobalProblemMapped(Eigen::VectorXd& g, Eigen::VectorXd& h, Eigen::VectorXd& vEst, Eigen::VectorXd& b)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Constructing RHS... ";

	vEst.resize(B2D.cols());
	for (int i = 0; i < vEst.rows(); i++) {
		vEst(i) = 0.5;
	}

	g = B2D * vEst;
	b.resize(B2D.rows() + c.rows(), c.cols());

	// First column of b
	h = C * vEst - c;
	b << g, h;

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}

void VectorFields::setupLHSGlobalProblemMapped(Eigen::SparseMatrix<double>& A_LHS)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Constructing LHS... ";

	A_LHS.resize(B2D.rows() + C.rows(), B2D.cols() + C.rows());

	vector<Eigen::Triplet<double>>	ATriplet;
	ATriplet.reserve(10 * B2D.rows());		// It should be #rows x 4 blocks @ 2 elements (8) + #constraints,
											// but made it 10 for safety + simplicity

	for (int k = 0; k < B2D.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(B2D, k); it; ++it) {
			ATriplet.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
		}
	}

	for (int k = 0; k < C.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(C, k); it; ++it) {
			ATriplet.push_back(Eigen::Triplet<double>(B2D.rows() + it.row(), it.col(), it.value()));
			ATriplet.push_back(Eigen::Triplet<double>(it.col(), B2D.cols() + it.row(), it.value()));
		}
	}
	A_LHS.setFromTriplets(ATriplet.begin(), ATriplet.end());

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}


void VectorFields::solveGlobalSystemMappedLDLT(Eigen::VectorXd& vEst, Eigen::SparseMatrix<double>& A_LHS, Eigen::VectorXd& b)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Solving the global system (Pardiso LDLT)... \n";



	//cout << "Starting to solve problem." << endl;
	Xf.resize(B2D.rows());
	
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
		if(sparseSolver.info()==Eigen::InvalidInput)
			cout << "Input is Invalid. " << endl;
		cout << sparseSolver.info() << endl;
		return;
	}
	
	Xf = -x.block(0, 0, B2D.rows(), 1) + vEst;

	printf("____Xf size is %dx%d\n", Xf.rows(), Xf.cols());	
}

void VectorFields::solveGlobalSystemMappedLU_GPU(Eigen::VectorXd& vEst, Eigen::SparseMatrix<double>& A_LHS, Eigen::VectorXd& b)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Solving the global system (LU in GPU)... \n";

	//cout << "Starting to solve problem." << endl;
	Xf.resize(B2D.rows());

	Eigen::MatrixXd X;
	solveLUinCUDA(A_LHS, b, X);


	Xf.col(0) = -X.block(0, 0, B2D.rows(), 1) + vEst;
	Xf.col(1) = -X.block(0, 1, B2D.rows(), 1) + vEst;
	//cout << Xf.block(0, 0, 100, 2) << endl; 

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "..in " << duration.count() << " seconds" << endl;
}

void VectorFields::constructSamples(const int &n)
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

void VectorFields::farthestPointSampling()
{
	Sample.resize(numSample);
	Eigen::VectorXd D;
	D.resize(F.rows());

	// Initialize the value of D
	for (int i = 0; i < F.rows(); i++) {
		D(i) = numeric_limits<double>::infinity();
	}

	srand(time(NULL));
	Sample[0] = rand() % F.rows();
	//Sample[0] = 0;

	//computeDijkstraDistanceFaceForSampling(Sample[0], D);
	//Eigen::VectorXi::Index maxIndex1;
	//D.maxCoeff(&maxIndex1);
	//Sample[1] = maxIndex1;

	for (int i = 1; i < numSample; i++) {
		Eigen::VectorXi::Index maxIndex;
		computeDijkstraDistanceFaceForSampling(Sample[i-1], D);
		D.maxCoeff(&maxIndex);
		Sample[i] = maxIndex;
	}
}

void VectorFields::constructBasis()
{	
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Construting Basis...\n";

	double	coef = sqrt(pow(1.1, 2) + pow(1.1, 2));
	double distRatio = coef * sqrt((double)V.rows() / (double) Sample.size());

	// Setup sizes of each element to construct basis
	try {
		BasisTemp.resize(2 * F.rows(), 2 * Sample.size());
	}
	catch(string &msg) {
		cout << "Cannot allocate memory for basis.." << endl;
	}
	
	Basis.resize(BasisTemp.rows(), BasisTemp.cols());
	vector<vector<Eigen::Triplet<double>>> UiTriplet(Sample.size());

	
	cout << "....Constructing and solving local systems...";
	const int NUM_PROCESS = 8;
	durations.resize(NUM_PROCESS);

	for (int i = 0; i < NUM_PROCESS; i++) {
		durations[i] = t1 - t1;
		//cout << "Init dur " << i<< " = " << durations[i].count() << " seconds" << endl;
	}
	
	int id, tid, ntids, ipts, istart, iproc;
	

#pragma omp parallel private(tid,ntids,ipts,istart,id)	
	{		
		iproc = omp_get_num_procs();
		//iproc = 1; 
		tid		= omp_get_thread_num();
		ntids	= omp_get_num_threads();
		ipts	= (int)ceil(1.00*(double)Sample.size() / (double)ntids);
		istart	= tid * ipts;
		if (tid == ntids - 1) ipts = Sample.size() - istart;
		if (ipts <= 0) ipts = 0;

		Eigen::VectorXd				D(F.rows());
		for (int i = 0; i < F.rows(); i++) {
			D(i) = numeric_limits<double>::infinity();
		}
		
		//cout << "[" << tid << "] Number of processors " << iproc << ", with " << ntids << " threads." << endl;

		UiTriplet[tid].reserve(2.0 * ((double)ipts / (double)Sample.size()) * 2 * 10.0 * F.rows());

		// Computing the values of each element
		for (id = istart; id < (istart + ipts); id++) {
			if (id >= Sample.size()) break;

			LocalFields localField(id);
				t1 = chrono::high_resolution_clock::now();
			localField.constructSubdomain(Sample[id], V, F, avgEdgeLength, AdjMF3N, distRatio);
				t2 = chrono::high_resolution_clock::now();
				durations[0] += t2 - t1;

				t1 = chrono::high_resolution_clock::now();
			localField.constructBoundary(F, AdjMF3N, AdjMF2Ring);
				t2 = chrono::high_resolution_clock::now();
				durations[1] += t2 - t1;

				t1 = chrono::high_resolution_clock::now();
			localField.constructLocalElements(F);
				t2 = chrono::high_resolution_clock::now();
				durations[2] += t2 - t1;

				t1 = chrono::high_resolution_clock::now();
			//localField.constructMatrixBLocal(B2D);
			localField.constructMatrixBLocal(B2D, AdjMF2Ring);
				t2 = chrono::high_resolution_clock::now();
				durations[3] += t2 - t1;

				t1 = chrono::high_resolution_clock::now();
			localField.constructLocalConstraints();
				t2 = chrono::high_resolution_clock::now();
				durations[4] += t2 - t1;

				t1 = chrono::high_resolution_clock::now();
			localField.setupRHSLocalProblemMapped();
				t2 = chrono::high_resolution_clock::now();
				durations[5] += t2 - t1;

				t1 = chrono::high_resolution_clock::now();
			localField.setupLHSLocalProblemMapped();
				t2 = chrono::high_resolution_clock::now();
				durations[6] += t2 - t1;

				t1 = chrono::high_resolution_clock::now();
			localField.solveLocalSystemMappedLDLT(UiTriplet[id]);
				t2 = chrono::high_resolution_clock::now();
				durations[7] += t2 - t1;
			//cout << "System " << id << " ( " << XfLoc.rows() << ") is solved." << endl; 
			//printf("System %d (%d) is solved.\n", id, XfLoc.rows());

				localField.measureXF(doubleArea, J);
		}

	}
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl; 

	cout << "....Gathering local elements as basis matrix... ";
	t1 = chrono::high_resolution_clock::now();
	gatherBasisElements(UiTriplet);
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
	
	cout << "....Partition of unity of the basis matrix... ";
	t1 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	//normalizeBasis();
	normalizeBasisAbs();

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "..in Total of " << duration.count() << " seconds" << endl;

	//for (int i = 0; i < NUM_PROCESS; i++) {
		//printf("Process %d takes %.8f seconds.\n", i, durations[i].count());
		//cout << "BasisTemp matrix is normalized in " << duration.count() << " seconds" << endl;
	//}

	// Information about the process timing
	printf("> Basis Timing information \n");
	printf("....[0] Constructing internal elements: %.8f seconds.\n", durations[0].count());
	printf("....[1] Constructing boundary: %.8f seconds.\n", durations[1].count());
	printf("....[2] Constructing local subdomains: %.8f seconds.\n", durations[2].count());
	printf("....[3] Constructing matrix B local: %.8f seconds.\n", durations[3].count());
	printf("....[4] Constructing local constraints: %.8f seconds.\n", durations[4].count());
	printf("....[5] Constructing RHS (mapped): %.8f seconds.\n", durations[5].count());
	printf("....[6] Constructing LHS (mapped): %.8f seconds.\n", durations[6].count());
	printf("....[7] Solving local systems (mapped, Pardiso LDLT): %.8f seconds.\n", durations[7].count());

	// Information about Basis
	printf("> Basis Structure information \n");
	printf("....Size = %dx%d\n", Basis.rows(), Basis.cols());
	printf("....NNZ per row = %.2f\n", (double) Basis.nonZeros() / (double) Basis.rows());
}

void VectorFields::gatherBasisElements(const vector<vector<Eigen::Triplet<double>>> &UiTriplet)
{
	vector<Eigen::Triplet<double>> BTriplet;
	BasisSum.resize(2 * F.rows(), 2);
	for (int i = 0; i < BasisSum.rows(); i++) {
		BasisSum(i, 0) = 0.0;
		BasisSum(i, 1) = 0.0;
	}

	int totalElements = 0;
	for (int j = 0; j < Sample.size(); j++) {
		totalElements += UiTriplet[j].size();
	}

	BTriplet.resize(totalElements);
	for (int j = 0; j < Sample.size(); j++) {
		int tripSize = 0;
		for (int k = 0; k < j; k++) {
			tripSize += UiTriplet[k].size();
		}
		std::copy(UiTriplet[j].begin(), UiTriplet[j].end(), BTriplet.begin() + tripSize);
	}
	BasisTemp.setFromTriplets(BTriplet.begin(), BTriplet.end());

	//printf("A basis matrix (%dx%d) is constructed.\n", BasisTemp.rows(), BasisTemp.cols());

	for (int k = 0; k < BasisTemp.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(BasisTemp, k); it; ++it) {
			BasisSum(it.row(), it.col() % 2) += it.value();
		}
	}

}

void VectorFields::normalizeBasis()
{
	Eigen::MatrixXd BasisSum(BasisTemp.rows(), 2);
	BasisSumN.resize(BasisTemp.rows(), 2);
	Eigen::MatrixXd normSum(BasisTemp.rows(), 2);
	Eigen::MatrixXd BasisNorm(F.rows(), 2);
	Eigen::MatrixXi nonZeros(BasisTemp.rows(), 2);
		for (int i = 0; i < nonZeros.rows(); i++) {
			for (int j = 0; j < nonZeros.cols(); j++) {
				nonZeros(i, j) = 0;
				BasisSum(i, j) = 0.0;
				BasisSumN(i, j) = 0.0;
			}
		}
	vector<Eigen::Triplet<double>> BNTriplet;
	BNTriplet.reserve(BasisTemp.nonZeros());

	// Getting the sum of each pair on each frame AND
	// Counting the non-zeros per rows
	for (int k = 0; k < BasisTemp.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(BasisTemp, k); it; ++it) {
			BasisSum(it.row(), it.col() % 2) += it.value();
			nonZeros(it.row(), it.col() % 2) += 1;
		}
	}
	
	// Computing the norm
	for (int i = 0; i < F.rows(); i++) {
		double frame1Norm, frame2Norm, a,b;
		a = BasisSum(2 * i + 0, 0);
		b = BasisSum(2 * i + 1, 0);
		frame1Norm = sqrt(a*a + b*b);

		a = BasisSum(2 * i + 0, 1);
		b = BasisSum(2 * i + 1, 1);
		frame2Norm = sqrt(a*a + b*b);

		BasisNorm(i, 0) = frame1Norm;
		BasisNorm(i, 1) = frame2Norm;		
	}

	// Constructing normalized basis each element of basis
	for (int k = 0; k < BasisTemp.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(BasisTemp, k); it; ++it) {
			//double newValue = it.value() / BasisNorm(it.row() / 2, it.col() % 2);
			double newValue = it.value();
			BNTriplet.push_back(Eigen::Triplet<double>(it.row(), it.col(), newValue));
		}
	}
	Basis.setFromTriplets(BNTriplet.begin(), BNTriplet.end());	

	// Check Normalization
	// Getting the sum of each pair on each frame
	for (int k = 0; k < Basis.outerSize(); ++k) {		
		for (Eigen::SparseMatrix<double>::InnerIterator it(Basis, k); it; ++it) {
			BasisSumN(it.row(), it.col() % 2) += it.value();
		}
	}

	// Computing the norm
	for (int i = 0; i < F.rows(); i++) {
		double frame1Norm, frame2Norm, a, b;
		a = BasisSumN(2 * i + 0, 0);
		b = BasisSumN(2 * i + 1, 0);
		frame1Norm = sqrt(a*a + b*b);

		a = BasisSumN(2 * i + 0, 1);
		b = BasisSumN(2 * i + 1, 1);
		frame2Norm = sqrt(a*a + b*b);

		BasisNorm(i, 0) = frame1Norm;
		BasisNorm(i, 1) = frame2Norm;
	}
	
	// Show result (The sum=> should all be equal to 1
	//for (int i = 0; i < F.rows(); i++) {
	//	printf("[%.6f][%.6f]\n", BasisSumN.block(2*i,0,2,1).norm(), BasisSumN.block(2 * i, 1, 2, 1).norm());
	//}

	int numNonZeroes = Basis.nonZeros();
	int numElements = Basis.rows();
	//cout << "Average non-zeros is " << (double)numNonZeroes / (double)numElements << endl; 
}

//void VectorFields::normalizeBasisAbs()
//{
//	Eigen::MatrixXd BasisSum(BasisTemp.rows(), 2);
//	BasisSumN.resize(BasisTemp.rows(), 2);
//	Eigen::MatrixXd normSum(BasisTemp.rows(), 2);
//	Eigen::MatrixXd BasisNorm(F.rows(), 2);
//	Eigen::MatrixXi nonZeros(BasisTemp.rows(), 2);
//	for (int i = 0; i < nonZeros.rows(); i++) {
//		for (int j = 0; j < nonZeros.cols(); j++) {
//			nonZeros(i, j) = 0;
//			BasisSum(i, j) = 0.0;
//			BasisSumN(i, j) = 0.0;
//		}
//	}
//	vector<Eigen::Triplet<double>> BNTriplet;
//
//
//	// Getting the sum of each pair on each frame AND
//	// Counting the non-zeros per rows
//	for (int k = 0; k < BasisTemp.outerSize(); ++k) {
//		for (Eigen::SparseMatrix<double>::InnerIterator it(BasisTemp, k); it; ++it) {
//			BasisSum(it.row(), it.col() % 2) += abs(it.value());
//			nonZeros(it.row(), it.col() % 2) += 1;
//		}
//	}
//
//	// Computing the norm
//	for (int i = 0; i < F.rows(); i++) {
//		double frame1Norm, frame2Norm, a, b;
//		a = BasisSum(2 * i + 0, 0);
//		b = BasisSum(2 * i + 1, 0);
//		frame1Norm = sqrt(a*a + b*b);
//
//		a = BasisSum(2 * i + 0, 1);
//		b = BasisSum(2 * i + 1, 1);
//		frame2Norm = sqrt(a*a + b*b);
//
//		BasisNorm(i, 0) = frame1Norm;
//		BasisNorm(i, 1) = frame2Norm;
//	}
//
//	// Constructing normalized basis each element of basis
//	for (int k = 0; k < BasisTemp.outerSize(); ++k) {
//		for (Eigen::SparseMatrix<double>::InnerIterator it(BasisTemp, k); it; ++it) {
//			double newValue = it.value() / BasisNorm(it.row() / 2, it.col() % 2);
//			//double newValue = it.value();
//			BNTriplet.push_back(Eigen::Triplet<double>(it.row(), it.col(), newValue));
//		}
//	}
//	Basis.setFromTriplets(BNTriplet.begin(), BNTriplet.end());
//
//	// Check Normalization
//	// Getting the sum of each pair on each frame
//	for (int k = 0; k < Basis.outerSize(); ++k) {
//		for (Eigen::SparseMatrix<double>::InnerIterator it(Basis, k); it; ++it) {
//			BasisSumN(it.row(), it.col() % 2) += it.value();
//		}
//	}
//
//	// Computing the norm
//	for (int i = 0; i < F.rows(); i++) {
//		double frame1Norm, frame2Norm, a, b;
//		a = BasisSumN(2 * i + 0, 0);
//		b = BasisSumN(2 * i + 1, 0);
//		frame1Norm = sqrt(a*a + b*b);
//
//		a = BasisSumN(2 * i + 0, 1);
//		b = BasisSumN(2 * i + 1, 1);
//		frame2Norm = sqrt(a*a + b*b);
//
//		BasisNorm(i, 0) = frame1Norm;
//		BasisNorm(i, 1) = frame2Norm;
//	}
//
//	// Show result (The sum=> should all be equal to 1
//	//for (int i = 0; i < F.rows(); i++) {
//	//	printf("[%.6f][%.6f]\n", BasisSumN.block(2*i,0,2,1).norm(), BasisSumN.block(2 * i, 1, 2, 1).norm());
//	//}
//
//	int numNonZeroes = Basis.nonZeros();
//	int numElements = Basis.rows();
//	cout << "Average non-zeros is " << (double)numNonZeroes / (double)numElements << endl;
//}

void VectorFields::normalizeBasisAbs()
{
	Eigen::MatrixXd normSum(F.rows(), 2);
	BasisSumN.resize(BasisTemp.rows(), 2);
	vector<Eigen::Triplet<double>> BNTriplet;
	BNTriplet.reserve(BasisTemp.nonZeros());

	Eigen::MatrixXd BasisNorm(F.rows(), 2);

	for (int i = 0; i < normSum.rows(); i++) {
		for (int j = 0; j < normSum.cols(); j++) {
			normSum(i, j) = 0.0;
			BasisSumN(2 * i + 0, j) = 0.0;
			BasisSumN(2 * i + 1, j) = 0.0;
		}
	}

	// Getting the sum of norm on each pair on each frame 
	for (int k = 0; k < BasisTemp.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(BasisTemp, k); it; ++it) {
			if (it.row() % 2 == 1) continue;

			double a = it.value();
			double b = BasisTemp.coeff(it.row() + 1, it.col());
			double norm = sqrt(a*a + b*b);
			normSum(it.row() / 2, it.col() % 2) += norm;
		}
	}

	// Normalize the system
	for (int k = 0; k < BasisTemp.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(BasisTemp, k); it; ++it) {
			double newValue = it.value() / normSum(it.row() / 2, it.col() % 2);
			// To have the basis with norm 2.0
			BNTriplet.push_back(Eigen::Triplet<double>(it.row(), it.col(), newValue/2.0));
			BasisSumN(it.row(), it.col() % 2) += newValue;
		}
	}

	Basis.setFromTriplets(BNTriplet.begin(), BNTriplet.end());	
}
void VectorFields::setAndSolveUserSystem()
{
	// Declare function-scoped variables
	Eigen::VectorXd					bBar, gBar, hBar, vEstBar;
	Eigen::SparseMatrix<double>		A_LHSBar;

	setupUserBasis();
	getUserConstraints();
	setupRHSUserProblemMapped(gBar, hBar, vEstBar, bBar);
	setupLHSUserProblemMapped(A_LHSBar);
	solveUserSystemMappedLDLT(vEstBar, A_LHSBar, bBar);
	mapSolutionToFullRes();
}

void VectorFields::setupUserBasis()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Constructing Local Basis...";

	B2DBar = Basis.transpose() * B2D * Basis; 

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;
	printf(".... Local Basis = %dx%d\n", B2DBar.rows(), B2DBar.cols());
}
void VectorFields::getUserConstraints()
{	
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Obtaining user constraints ";

	userConstraints = globalConstraints; 
	CBar			= C * Basis;
	cBar			= c;

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;

	printf(".... C_Lobal = %dx%d\n", CBar.rows(), CBar.cols());
	printf(".... c_Lobal = %dx%d\n", cBar.rows(), cBar.cols());
}

void VectorFields::setupRHSUserProblemMapped(Eigen::VectorXd& gBar, Eigen::VectorXd& hBar, Eigen::VectorXd& vEstBar, Eigen::VectorXd& bBar)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Constructing RHS (mapped)...";
	
	vEstBar.resize(B2DBar.rows());
	for (int i = 0; i < vEstBar.rows(); i++) {
		vEstBar(i) = 0.5;
	}

	gBar = B2DBar * vEstBar;
	bBar.resize(B2DBar.rows() + cBar.rows());

	// Constructing b
	hBar = CBar * vEstBar - cBar;
	bBar<< gBar, hBar;

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;
}

void VectorFields::setupLHSUserProblemMapped(Eigen::SparseMatrix<double>& A_LHSBar)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Constructing LHS (mapped)...";


	A_LHSBar.resize(B2DBar.rows() + CBar.rows(), B2DBar.cols() + CBar.rows());
	vector<Eigen::Triplet<double>>	ATriplet;
	ATriplet.reserve(B2DBar.nonZeros() + 2 * CBar.nonZeros());

	for (int k = 0; k < B2DBar.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(B2DBar, k); it; ++it) {
			ATriplet.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
		}
	}

	for (int k = 0; k < CBar.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(CBar, k); it; ++it) {
			ATriplet.push_back(Eigen::Triplet<double>(B2DBar.rows() + it.row(), it.col(), it.value()));
			ATriplet.push_back(Eigen::Triplet<double>(it.col(), B2DBar.cols() + it.row(), it.value()));
		}
	}
	A_LHSBar.setFromTriplets(ATriplet.begin(), ATriplet.end());

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;
	printf("....Local LHS = %dx%d\n", A_LHSBar.rows(), A_LHSBar.cols());
}

void VectorFields::solveUserSystemMappedLDLT(Eigen::VectorXd& vEstBar, Eigen::SparseMatrix<double>& A_LHSBar, Eigen::VectorXd& bBar)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Solving reduced system...\n";

	XLowDim.resize(B2DBar.rows());
	XFullDim.resize(Basis.rows());
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

	XLowDim.col(0) = -x.block(0, 0, B2DBar.rows(), 1) + vEstBar;	
	
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "..in Total of " << duration.count() << " seconds." << endl;

	//cout << "Solution (LowDim) \n" << XLowDim << endl; 	
}

void VectorFields::mapSolutionToFullRes()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Mapping to full-resolution...";

	XFullDim = Basis * XLowDim;


	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << " in " << duration.count() << " seconds." << endl;

	printf("....XFull (%dx%d) =  Basis (%dx%d) * XLowDim (%dx%d) \n", XFullDim.rows(), XFullDim.cols(), Basis.rows(), Basis.cols(), XLowDim.rows(), XLowDim.cols());
}

void VectorFields::obtainUserVectorFields()
{
	cout << "Hello there \n" << endl; 
}

void VectorFields::measureApproxAccuracyL2Norm()
{
	Eigen::VectorXd diff = Xf - XFullDim;
	double xf = Xf.transpose() * MF2D * Xf;
	double L2norm = diff.transpose() * MF2D * diff; 
	printf("Diff 0 = %.10f\n", sqrt(L2norm / xf)); 	
}

void VectorFields::measureU1andJU0()
{
	//Eigen::VectorXd JXf0 = J*BasisTemp.col(0);
	//Eigen::VectorXd diff = Eigen::VectorXd(BasisTemp.col(1)) - JXf0;
	//
	//double norm1 = JXf0.transpose() * MF2D * JXf0;
	//double norm2 = diff.transpose()* MF2D * diff; 
	//double diffNorm = sqrt(norm2 / norm1);
	//cout << "_____|Xf(1)-J*Xf(0)|M = " << diffNorm << endl; 
}

/* ====================== MESH-RELATED FUNCTIONS ============================*/
void VectorFields::readMesh(const string &meshFile)
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
}

void VectorFields::getVF(Eigen::MatrixXd &V, Eigen::MatrixXi &F)
{
	V = this->V;
	F = this->F;
}

void VectorFields::computeFaceCenter()
{
	FC.resize(F.rows(), 3);

	for (int i = 0; i < F.rows(); i++) {
		FC.row(i) = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.0;
	}

}


//void VectorFields::visualizeSparseMatrixInMatlab(const Eigen::SparseMatrix<double> &M)
//{
//	printf("Size of M=%dx%d\n", M.rows(), M.cols());
//
//	using namespace matlab::engine;
//	Engine *ep;
//	mxArray *MM = NULL, *MS = NULL, *result = NULL, *eigVecResult, *nEigs;
//
//	const int NNZ_M = M.nonZeros();
//	int nnzMCounter = 0;
//
//	double	*srm = (double*)malloc(NNZ_M * sizeof(double));
//	mwIndex *irm = (mwIndex*)malloc(NNZ_M * sizeof(mwIndex));
//	mwIndex *jcm = (mwIndex*)malloc((M.cols() + 1) * sizeof(mwIndex));
//
//	MM = mxCreateSparse(M.rows(), M.cols(), NNZ_M, mxREAL);
//	srm = mxGetPr(MM);
//	irm = mxGetIr(MM);
//	jcm = mxGetJc(MM);
//
//	// Getting matrix M
//	jcm[0] = nnzMCounter;
//	for (int i = 0; i < M.outerSize(); i++) {
//		for (Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it) {
//			srm[nnzMCounter] = it.value();
//			irm[nnzMCounter] = it.row();
//			nnzMCounter++;
//		}
//		jcm[i + 1] = nnzMCounter;
//	}
//
//	// Start Matlab Engine
//	ep = engOpen(NULL);
//	if (!(ep = engOpen(""))) {
//		fprintf(stderr, "\nCan't start MATLAB engine\n");
//		cout << "CANNOT START MATLAB " << endl;
//	}
//	else {
//		cout << "MATLAB STARTS. OH YEAH!!!" << endl;
//	}
//
//	engPutVariable(ep, "M", MM);
//	engEvalString(ep, "spy(M)");
//}