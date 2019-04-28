#include "VectorFields.h"

#include <igl/per_vertex_normals.h>
#include <Eigen/Eigenvalues>
#include <random>
#include <Eigen/OrderingMethods>

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
	//constructSpecifiedHardConstraints();
	//constructRandomHardConstraints();
	constructSoftConstraints();
	//constructInteractiveConstraints();
	//constructInteractiveConstraintsWithLaplacian();

	//constructSingularities();
	//constructHardConstraintsWithSingularities();
	//constructHardConstraintsWithSingularities_Cheat();
	//constructHardConstraintsWithSingularitiesWithGauss();

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;

	// Information about constraints
	//printf("....Num of Constraints = %d\n", globalConstraints.size());
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

void VectorFields::constructSpecifiedHardConstraints()
{
	// Define the constraints
	const int numConstraints = 20;
	set<int> constraints;
	//vector<int> globalConstraints(numConstraints);
	globalConstraints.resize(numConstraints);
	Eigen::VectorXd D;
	D.resize(F.rows());

	// Initialize the value of D
	for (int i = 0; i < F.rows(); i++) {
		D(i) = numeric_limits<double>::infinity();
	}

	/* Random number generator */
	std::random_device rd;								// Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd());								// Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<> dis(0, F.rows() - 1); // From 0 to F.rows()-1

	srand(time(NULL));
	//int curPoint = rand() % F.rows();
	//int curPoint = 0; 
	int curPoint = dis(gen);
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
		
	// For testing only
	//computeDijkstraDistanceFaceForSampling(curPoint, D);
	//Eigen::VectorXi::Index counterPart;
	//D.maxCoeff(&counterPart);
	//const int counterPart = AdjMF3N(curPoint, 0);
	//globalConstraints[1] = AdjMF3N(0, 0);
	//globalConstraints[2] = AdjMF3N(0, 1);
	//globalConstraints[3] = AdjMF3N(0, 2);
	//printf("Constraints = %d\n", globalConstraints.size());

	// Setting up matrix C
	Eigen::SparseMatrix<double> CTemp;
	vector<Eigen::Triplet<double>> CTriplet;
	CTriplet.reserve(2 * globalConstraints.size());
	c.resize(2 * globalConstraints.size());
	Eigen::Vector2d cRand;

	int counter = 0;
	for (int i = 0; i < globalConstraints.size(); i++) {
		cRand(0) = (double)(rand() % F.rows()) / (double)F.rows();
		cRand(1) = (double)(rand() % F.rows()) / (double)F.rows();
		cRand.normalize();

		//const double alpha = M_PI / 2.0; 
		//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * globalConstraints[i] + 0, cos(alpha)));
		//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * globalConstraints[i] + 1, -sin(alpha)));
		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * globalConstraints[i] + 0, 1.0));
		//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * counterPart + 0, 1.0));
		//c(2 * counter + 0, 0) = 1.0;
		c(counter, 0) = cRand(0);
		counter++;

		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * globalConstraints[i] + 1, 1.0));
		//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * globalConstraints[i] + 0, sin(alpha)));
		//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * globalConstraints[i] + 1, cos(alpha)));
		//c(2 * counter + 1, 0) = 1.0;
		//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * counterPart + 1, 1.0));
		c(counter, 0) = cRand(1);
		counter++;
	}
	C.resize(2 * globalConstraints.size(), B2D.rows());
	C.setFromTriplets(CTriplet.begin(), CTriplet.end());
	//printf("Cp=%dx%d\n", C.rows(), C.cols());
}

void VectorFields::constructRandomHardConstraints()
{
	// Define the constraints
	const bool readFromFile = true; 
	string filename = "D:/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Constraints/Constraints_CDragon_Rand_25.txt";;
	//string filename = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Constraints/Constraints_Cube_Rand_25.txt";

	if (readFromFile)
	{
		LoadSTDVectorFromTxtFile(filename, globalConstraints);
	} 
	else
	{
		const int numConstraints = 25;
		set<int> constraints;
		globalConstraints.resize(numConstraints);

		/* Random number generator */
		std::random_device rd;								// Will be used to obtain a seed for the random number engine
		std::mt19937 gen(rd());								// Standard mersenne_twister_engine seeded with rd()
		std::uniform_int_distribution<> dis(0, F.rows() - 1); // From 0 to F.rows()-1

															  /* Creating random constraints */
		do {
			int constraintFace = dis(gen);
			constraints.insert(constraintFace);
		} while (constraints.size() < numConstraints);

		int counter1 = 0;
		for (int i : constraints) {
			globalConstraints[counter1++] = i;
		}

		WriteSTDVectorToTxtFile(globalConstraints, filename);
	}

	
	cout << "Setting up matrix C\n";
	// Setting up matrix C
	Eigen::SparseMatrix<double> CTemp;
	vector<Eigen::Triplet<double>> CTriplet;
	CTriplet.reserve(2 * globalConstraints.size());
	c.resize(2 * globalConstraints.size());
	Eigen::Vector2d cRand;

	int counter = 0;
	for (int i = 0; i < globalConstraints.size(); i++) {
		cRand(0) = (double)(rand() % F.rows()) / (double) F.rows();
		cRand(1) = (double)(rand() % F.rows()) / (double)F.rows();
		cRand.normalize();

		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * globalConstraints[i] + 0, 1.0));
		c(counter, 0) = cRand(0);
		counter++;

		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * globalConstraints[i] + 1, 1.0));
		c(counter, 0) = cRand(1);
		counter++;
	}
	C.resize(2 * globalConstraints.size(), B2D.rows());
	C.setFromTriplets(CTriplet.begin(), CTriplet.end());
}

void VectorFields::pushNewUserConstraints(const int& fInit, const int& fEnd)
{
	userVisualConstraints.push_back(fInit);
	userVisualConstraints.push_back(fEnd);
}
void VectorFields::constructInteractiveConstraints()
{
	/* Define the constraints */
	const int numConstraints = userVisualConstraints.size() / 2; 
	globalConstraints.resize(numConstraints);
	vector<Eigen::Vector2d> constraintValues(numConstraints);

	/* Global constraints from user input */
	for (int i = 0; i < userVisualConstraints.size(); i += 2)
	{
		/* Location of constraints */
		globalConstraints[i/2] = userVisualConstraints[i];

		/* Getting the constraints + making them into local coordinates */
		Eigen::RowVector3d dir = FC.row(userVisualConstraints[i+1]) - FC.row(userVisualConstraints[i]);
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
	c.resize(2 * globalConstraints.size());

	/* Putting the constraints into action */
	for (int i = 0; i < globalConstraints.size(); i++) {
		CTriplet.push_back(Eigen::Triplet<double>(2*i+0, 2 * globalConstraints[i] + 0, 1.0));
		c(2*i+0) = constraintValues[i](0);		
	
		CTriplet.push_back(Eigen::Triplet<double>(2*i+1, 2 * globalConstraints[i] + 1, 1.0));
		c(2*i+1) = constraintValues[i](1);
	}	
	C.resize(2 * globalConstraints.size(), B2D.rows());
	C.setFromTriplets(CTriplet.begin(), CTriplet.end());
}

void VectorFields::constructInteractiveConstraintsWithLaplacian()
{
	/* Define the constraints */
	const int numConstraints = userVisualConstraints.size() / 2;
	globalConstraints.clear();
	globalConstraints.shrink_to_fit();
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
		constraintValues[i / 2] = normDir;
	}

	/* Setting up matrix C and column vector c */
	Eigen::SparseMatrix<double> CTemp;
	vector<Eigen::Triplet<double>> CTriplet;
	CTriplet.reserve(40 * globalConstraints.size());
	c.resize(4 * globalConstraints.size());

	/* Putting the constraints into action */
	for (int i = 0; i < globalConstraints.size(); i++) {
		CTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * globalConstraints[i] + 0, 1.0));
		c(2 * i + 0) = constraintValues[i](0);

		CTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * globalConstraints[i] + 1, 1.0));
		c(2 * i + 1) = constraintValues[i](1);
	}

	/* Putting the constraints into action with LAPLACIAN CONSTRAINTS */
	Eigen::SparseMatrix<double> LapVFields = /*MF2Dinv * */SF2DAsym;

	//int counter = 0;
	for (int i = 0; i < globalConstraints.size(); i++) {
		const int gC = 2 * globalConstraints[i];
		for (int k = gC; k <= (gC + 1); k++)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(LapVFields, k); it; ++it)
			{
				printf("[%d, %d] = %.4f\n", it.row(), it.col(), it.value());
				//const double mInv = 2 / doubleArea(floor(it.row() / 2));
				CTriplet.push_back(Eigen::Triplet<double>(2*numConstraints+2*i+(k - gC), it.row(), it.value()));
			}
		}
		c(2*numConstraints+2*i, 0) = 0.0;
		c(2*numConstraints+2*i+ 1, 0) = -1;
	}

	C.resize(4 * globalConstraints.size(), B2D.rows());
	C.setFromTriplets(CTriplet.begin(), CTriplet.end());
	//visualizeSparseMatrixInMatlab(C);
}

void VectorFields::resetInteractiveConstraints()
{
	userVisualConstraints.clear();
	userVisualConstraints.shrink_to_fit();
}

void VectorFields::constructSingularities()
{
	const int NUM_SINGS = 2;

	if (NUM_SINGS > 0)
		constructVFAdjacency();
		//constructVFNeighborsFull();

	singularities.resize(NUM_SINGS);
	SingNeighCC.resize(NUM_SINGS);
	mappedBasis.resize(NUM_SINGS);
	mappedBasis2.resize(NUM_SINGS);

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
		mappedBasis[id].resize(SingNeighNum);
		mappedBasis2[id].resize(SingNeighNum);

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

void VectorFields::constructHardConstraintsWithSingularities() 
{
	// Define the constraints
	const int numConstraints = 2;
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
		for (int j = 0; j < (SingNeighCC[i].size() - 1); j++) {
			numSingConstraints++;
		}
	}

	// Setting up matrix C and vector c
	c.resize(2 * (globalConstraints.size() + numSingConstraints));


	// HARD CONSTRAINTS
	Eigen::SparseMatrix<double> CTemp;
	vector<Eigen::Triplet<double>> CTriplet;
	CTriplet.reserve(2 * globalConstraints.size() + 2 * 4 * 7 * SingNeighCC.size());
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
		// Getting the shared-edges of two neighboring faces For testing
		sharedEdgesVect[id].resize(2 * SingNeighCC[id].size() - 2);

		// 4. Compute rotation of among its valence
		const double rotAngle = 2 * M_PI / (double)SingNeighCC[id].size();
		const double cosConst = cos(rotAngle);
		const double sinConst = sin(rotAngle);

		// Which case? => determining which edge is the common edge
		enum class SharedEdgeCase { Case1, Case2, Case3 };
		SharedEdgeCase edgeCase1, edgeCase2;
		for (int i = 0; i < (SingNeighCC[id].size() - 1); i++) {
			// 1. Find shared edge (naively)
			Eigen::RowVector3d es;
			for (int f1 = 0; f1 < F.cols(); f1++) {
				for (int f2 = 0; f2 < F.cols(); f2++) {
					bool b1 = F(SingNeighCC[id][i], (f1 + 1) % F.cols()) == F(SingNeighCC[id][i + 1], f2);
					bool b2 = F(SingNeighCC[id][i], f1) == F(SingNeighCC[id][i + 1], (f2 + 1) % F.cols());
					if (b1 && b2) {
						sharedEdgesVect[id][2 * i + 0] = F(SingNeighCC[id][i], f1);
						sharedEdgesVect[id][2 * i + 1] = F(SingNeighCC[id][i], (f1 + 1) % F.cols());
						es = V.row(F(SingNeighCC[id][i], (f1 + 1) % F.cols())) - V.row(F(SingNeighCC[id][i], f1));
						printf("Shared edge=%d->%d\n", F(SingNeighCC[id][i], f1), F(SingNeighCC[id][i], (f1 + 1) % F.cols()));

						if (f1 == 0)		edgeCase1 = SharedEdgeCase::Case1;	// => edge V0->V1 is the shared edge => it takes 0 step to reach v0
						else if (f1 == 1)	edgeCase1 = SharedEdgeCase::Case3;	// => edge V1->V2 is the shared edge => it takes 2 step to reach v0
						else if (f1 == 2)	edgeCase1 = SharedEdgeCase::Case2;	// => edge V2->V0 is the shared edge => it takes 1 step to reach v0

						if (f2 == 0)		edgeCase2 = SharedEdgeCase::Case1;
						else if (f2 == 1)	edgeCase2 = SharedEdgeCase::Case3;
						else if (f2 == 2)	edgeCase2 = SharedEdgeCase::Case2;
					}
				}
			}
			// 2. Find angles between basis1 and shared_edges es
			Eigen::VectorXd eVect;
			Eigen::RowVector3d b11, b12;
			//b11 = (A.block(3 * SingNeighCC[id][i], 2 * SingNeighCC[id][i] + 0, 3, 1)).transpose();
			//b12 = (A.block(3 * SingNeighCC[id][i], 2 * SingNeighCC[id][i] + 1, 3, 1)).transpose();
			eVect = A.block(3 * SingNeighCC[id][i], 2 * SingNeighCC[id][i] + 0, 3, 1);
			b11 << eVect(0), eVect(1), eVect(2);
			eVect = A.block(3 * SingNeighCC[id][i], 2 * SingNeighCC[id][i] + 1, 3, 1);
			b12 << eVect(0), eVect(1), eVect(2);
			//cout << "______B11: " << b11 << ", B12: " << b12 << endl;

			// Basis 1, Frame 1
			double cosR12 = (b11.dot(es)) / (b11.norm()*es.norm());
			if (cosR12 > 1.0) cosR12 = 1.0;
			if (cosR12 <-1.0) cosR12 = -1.0;
			const double angleR12_1 = (edgeCase1 == SharedEdgeCase::Case2 ? 2 * M_PI - acos(cosR12) : acos(cosR12));
			printf("______[%.2f] Rotation matrix R12_1\n", angleR12_1*180.0 / M_PI);
			//const double cosR12_1 = cos(angleR12_1);
			//const double sinR12_1 = sin(angleR12_1);
			//printf("______[%.2f] Rotation matrix R12_1 = [%.3f,%.3f; %.3f, %.3f]\n", angleR12_1*180.0 / M_PI, cosR12_1, -sinR12_1, sinR12_1, cosR12_1);
			//printf("______[%.2f] Rotation matrix R12_1 = [%.3f,%.3f; %.3f, %.3f]\n", angleR12_1*180.0 / M_PI, cosR12_1, -sinR12_1, sinR12_1, cosR12_1);

			
			// 3. Find angles between basis2 and es
			//es = -es;
			Eigen::RowVector3d b21, b22;
			//b21 = (A.block(3 * SingNeighCC[id][i + 1], 2 * SingNeighCC[id][i + 1] + 0, 3, 1)).transpose();
			//b22 = (A.block(3 * SingNeighCC[id][i + 1], 2 * SingNeighCC[id][i + 1] + 1, 3, 1)).transpose();
			eVect = A.block(3 * SingNeighCC[id][i + 1], 2 * SingNeighCC[id][i + 1] + 0, 3, 1);
			b21 << eVect(0), eVect(1), eVect(2);
			eVect = A.block(3 * SingNeighCC[id][i + 1], 2 * SingNeighCC[id][i + 1] + 1, 3, 1);
			b22 << eVect(0), eVect(1), eVect(2);
			//cout << "______B21: " << b21 << ", B22: " << b22 << endl;

			// Basis 2, Frame 1
			double cosR21 = (b21.dot(es)) / (b21.norm()*es.norm());
			if (cosR21 > 1.0) cosR21 = 1.0;
			if (cosR21 < -1.0) cosR21 = -1.0;
			double angleR21_1 = (edgeCase2 == SharedEdgeCase::Case3 ? 2 * M_PI - acos(cosR21) : acos(cosR21));
			angleR21_1 = 2 * M_PI - angleR21_1; 
			printf("______[%.2f] Rotation matrix R22_1 = [%.2f]\n", angleR21_1*180.0 / M_PI);
			//const double cosR21_1 = cos(angleR21_1);
			//const double sinR21_1 = sin(angleR21_1);
			//printf("______[%.2f] Rotation matrix R22_1 = [%.2f,%.2f; %.2f, %.2f]\n", angleR21_1*180.0 / M_PI, cosR21_1, -sinR21_1, sinR21_1, cosR21_1);
			//printf("____ To map basis1 -> basis2: rotate by %.2f degree\n", (angleR12_1 + angleR21_1)*180.0 / M_PI);

			const double RotAngle = (angleR12_1 + angleR21_1 > 2 * M_PI ? (angleR12_1 + angleR21_1) - 2 * M_PI : angleR12_1 + angleR21_1);
			printf("____ To map basis1 -> basis2: rotate by %.2f degree\n", (RotAngle)*180.0 / M_PI);
			const double cosBasis = cos(RotAngle);
			const double sinBasis = sin(RotAngle);

			// Obtain the mapped basis
			Eigen::Matrix2d RotMat2D;
			RotMat2D << cosBasis, -sinBasis, sinBasis, cosBasis; 
			Eigen::Vector2d initBasis;
			initBasis << 1.0, 0.0;
			Eigen::Vector2d newBasis = RotMat2D * initBasis; 
			mappedBasis[id][i] = newBasis; 
			
			initBasis(0) = 0.0; initBasis(1) = 1.0;
			newBasis = RotMat2D * initBasis;
			mappedBasis2[id][i] = newBasis;


			// Basis 2, Frame 2
			//cosR21 = (b22.dot(es)) / (b22.norm()*es.norm());
			//if (cosR21 > 1.0) cosR21 = 1.0;
			//if (cosR21 < -1.0) cosR21 = -1.0;
			//const double angleR21_2 = (edgeCase2 == SharedEdgeCase::Case1 ? 2 * M_PI - acos(cosR21) : acos(cosR21));
			//const double cosR21_2 = cos(angleR21_2);
			//const double sinR21_2 = sin(angleR21_2);
			//printf("______[%.2f] Rotation matrix R22_2 = [%.2f,%.2f; %.2f, %.2f]\n", angleR21_2*180.0 / M_PI, cosR21_2, -sinR21_2, sinR21_2, cosR21_2);


			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, cosR12_1));
			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, -sinR12_1));
			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, -cosR21_1));
			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, sinR21_1));
			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, cosR21_1*cosR12_1 + sinR21_1*sinR12_1));
			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, -cosR21_1*sinR12_1 + sinR21_1*cosR12_1));
			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, -1.0));
			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, 0.0));
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, cosBasis));
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, -sinBasis));
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, -1.0));
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, 0));
			c(counter) = 0.0;
			counter++;

			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, sinR12_1));
			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, cosR12_1 ));
			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, -sinR21_1));
			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, -cosR21_1 ));
			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, -sinR21_1*cosR12_1 + cosR21_1*sinR12_1));
			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, sinR21_1*sinR12_1 + cosR21_1*cosR21_1));
			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, 0.0));
			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, -1.0));

			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, sinBasis));
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, cosBasis ));
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, 0.0));
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, -1.0));
			c(counter) = 0.0;
			counter++;
		}
	}

	C.resize(2 * (globalConstraints.size() + numSingConstraints), B2D.rows());
	C.setFromTriplets(CTriplet.begin(), CTriplet.end());
	//printf("Cp=%dx%d\n", C.rows(), C.cols());	
}

void VectorFields::constructHardConstraintsWithSingularities_Cheat()
{
	// Define the constraints
	const int numConstraints = 100;
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
		//for (int j = 0; j < (SingNeighCC[i].size()-1); j++) {
		for (int j = 0; j < (SingNeighCC[i].size()); j++) {
			numSingConstraints++;
		}
	}

	// Setting up matrix C and vector c
	c.resize(2 * (globalConstraints.size() + numSingConstraints));


	// HARD CONSTRAINTS
	Eigen::SparseMatrix<double> CTemp;
	vector<Eigen::Triplet<double>> CTriplet;
	CTriplet.reserve(2 * globalConstraints.size() + 2 * 4 * 7 * SingNeighCC.size());
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

	// Setting up hard constraints for neighboring faces
	Eigen::MatrixXd ALoc(3, 2);
	Eigen::RowVector3d c1, c2, field3D;
	Eigen::Vector2d field2D; 
	for (int id = 0; id < SingNeighCC.size(); id++) {	
		//printf("This sing has %d neighbors....", SingNeighCC[id].size());
		for (int i = 0; i < (SingNeighCC[id].size()); i++) {
			int i2 = (i < (SingNeighCC[id].size() - 1) ? i+1 : 0);

			// Computing the field from one barycenter pointing to another's barycenter
			ALoc = A.block(3 * SingNeighCC[id][i], 2 * SingNeighCC[id][i], 3, 2);
			c1 = FC.row(SingNeighCC[id][i]);
			c2 = FC.row(SingNeighCC[id][i2]); 
			field3D = c2 - c1;
			//printf("<%.15f, %.15f, %.15f>\n", field3D(0), field3D(1), field3D(2));
			field2D = ALoc.transpose() * field3D.transpose();
			field2D.normalize();
			field2D /= 4.0; 

			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, 1.0));
			c(counter) = field2D(0);
			counter++;

			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, 1.0));
			c(counter) = field2D(1);
			counter++;
			//printf("Writing to %d(%.3f) and %d(%.3f) \n", SingNeighCC[id][i], field2D(0), SingNeighCC[id][i2], field2D(1));
		}
	}

	C.resize(2 * (globalConstraints.size() + numSingConstraints), B2D.rows());
	C.setFromTriplets(CTriplet.begin(), CTriplet.end());
}

void VectorFields::constructHardConstraintsWithSingularitiesWithGauss()
{
	// Define the constraints
	const int numConstraints = 20;
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
		//for (int j = 0; j < (SingNeighCC[i].size()-1); j++) {
		for (int j = 0; j < (SingNeighCC[i].size()); j++) {
			numSingConstraints++;
		}
	}

	// Setting up matrix C and vector c
	c.resize(2 * (globalConstraints.size() + numSingConstraints));


	// HARD CONSTRAINTS
	Eigen::SparseMatrix<double> CTemp;
	vector<Eigen::Triplet<double>> CTriplet;
	CTriplet.reserve(2 * globalConstraints.size() + 2 * 4 * 7 * SingNeighCC.size());
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

	// Setting up hard constraints for neighboring faces
	Eigen::MatrixXd		ALoc(3, 2);
	Eigen::RowVector3d	c1, c2, field3D;
	Eigen::Vector2d		field2D;
	vector<vector<double>>		internalAngle(SingNeighCC.size()); 
	vector<double>				gaussAngle(SingNeighCC.size());

	// Local variables to compute angle on each triangle
	Eigen::Vector3d		edge1, edge2;
	double				angle; 
	int					singLoc;
	for (int id = 0; id < SingNeighCC.size(); id++)
	{
		internalAngle[id].resize(SingNeighCC[id].size());
		gaussAngle[id] = 0.0;

		for (int i = 0; i < (SingNeighCC[id].size()); i++) 
		{
			int i2 = (i < (SingNeighCC[id].size() - 1) ? i + 1 : 0);

			// [a] obtain shared edges
			for (int f = 0; f < F.cols(); f++) {
				if (F(SingNeighCC[id][i],f) == singularities[id])
				{
					// [b] get the two edges
					edge1 = V.row(F(SingNeighCC[id][i], (f == 0 ? 2 : f - 1))) - V.row(F(SingNeighCC[id][i], f));
					edge2 = V.row(F(SingNeighCC[id][i], (f == 2 ? 0 : f + 1))) - V.row(F(SingNeighCC[id][i], f)) ;
					angle = edge2.dot(edge1) / (edge1.norm()*edge2.norm());
					angle = acos(angle);

					// [c] get the angle
					internalAngle[id][i] = angle;
					gaussAngle[id]		+= angle; 
				}
			}
		}
	}

	// Show the angles
	for (int id = 0; id < SingNeighCC.size(); id++)
	{
		printf("__Gauss angle = %.4f\n", gaussAngle[id]*180.0/M_PI);
		for (int i = 0; i < (SingNeighCC[id].size()); i++) {
			printf("______ angle %d = %.3f \n", i, internalAngle[id][i] * 180.0 / M_PI);
		}
	}

	// SINGULARITIES CONSTRAINTS
	for (int id = 0; id < SingNeighCC.size(); id++) 
	{
		// Getting the shared-edges of two neighboring faces For testing
		sharedEdgesVect[id].resize(2 * SingNeighCC[id].size() - 2);

		// 4. Compute rotation of among its valence
		const double rotAngle = gaussAngle[id] / (double)SingNeighCC[id].size();	// all edges with similar angle => next iter: relative angle
		const double cosConst = cos(/*2*M_PI - */rotAngle);
		const double sinConst = sin(/*2*M_PI - */rotAngle);
		printf("Angle=%.3f, sin=%.3f, cos=%.3f\n", rotAngle*180.0 / M_PI, sinConst, cosConst);

		/* Give hard constraint on face 1 */
		Eigen::MatrixXd ALoc(3, 2);
		for (int f = 0; f < F.cols(); f++) {
			if (F(SingNeighCC[id][0], f) == singularities[id])
			{
				Eigen::Vector3d edge = V.row(F(SingNeighCC[id][0], (f==0 ? 2 : f-1))) - V.row(F(SingNeighCC[id][0], (f == 2 ? 0 : f + 1)));
				ALoc = A.block(3 * SingNeighCC[id][0], 2 * SingNeighCC[id][0], 3, 2);
				Eigen::Vector2d edge2D = ALoc.transpose() * edge; 
				edge2D = edge2D.normalized() / 4.0; 
				CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][0]+0, 1.0));
				c(counter) = edge2D(0);
				counter++;
				CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][0]+1, 1.0));
				c(counter) = edge2D(1);
				counter++;
			}
		}

		// Which case? => determining which edge is the common edge
		enum class SharedEdgeCase { Case1, Case2, Case3 };
		SharedEdgeCase edgeCase1, edgeCase2;
		for (int i = 0; i < (SingNeighCC[id].size() - 1); i++) 
		{
			// 1. Find shared edge (naively)
			Eigen::RowVector3d es;
			for (int f1 = 0; f1 < F.cols(); f1++) {
				for (int f2 = 0; f2 < F.cols(); f2++) {
					bool b1 = F(SingNeighCC[id][i], (f1 + 1) % F.cols()) == F(SingNeighCC[id][i + 1], f2);
					bool b2 = F(SingNeighCC[id][i], f1) == F(SingNeighCC[id][i + 1], (f2 + 1) % F.cols());
					if (b1 && b2) {
						sharedEdgesVect[id][2 * i + 0] = F(SingNeighCC[id][i], f1);
						sharedEdgesVect[id][2 * i + 1] = F(SingNeighCC[id][i], (f1 + 1) % F.cols());
						es = V.row(F(SingNeighCC[id][i], (f1 + 1) % F.cols())) - V.row(F(SingNeighCC[id][i], f1));
						printf("Shared edge=%d->%d\n", F(SingNeighCC[id][i], f1), F(SingNeighCC[id][i], (f1 + 1) % F.cols()));

						if (f1 == 0)		edgeCase1 = SharedEdgeCase::Case1;	// => edge V0->V1 is the shared edge => it takes 0 step to reach v0
						else if (f1 == 1)	edgeCase1 = SharedEdgeCase::Case3;	// => edge V1->V2 is the shared edge => it takes 2 step to reach v0
						else if (f1 == 2)	edgeCase1 = SharedEdgeCase::Case2;	// => edge V2->V0 is the shared edge => it takes 1 step to reach v0

						if (f2 == 0)		edgeCase2 = SharedEdgeCase::Case1;
						else if (f2 == 1)	edgeCase2 = SharedEdgeCase::Case3;
						else if (f2 == 2)	edgeCase2 = SharedEdgeCase::Case2;
					}
				}
			}

			// 2. Find angles between basis1 and shared_edges es
			Eigen::VectorXd eVect;
			Eigen::RowVector3d b11, b12;
			eVect = A.block(3 * SingNeighCC[id][i], 2 * SingNeighCC[id][i] + 0, 3, 1);
			b11 << eVect(0), eVect(1), eVect(2);
			eVect = A.block(3 * SingNeighCC[id][i], 2 * SingNeighCC[id][i] + 1, 3, 1);
			b12 << eVect(0), eVect(1), eVect(2);
			//cout << "______B11: " << b11 << ", B12: " << b12 << endl;

			// Basis 1, Frame 1
			double cosR12 = (b11.dot(es)) / (b11.norm()*es.norm());
			if (cosR12 > 1.0) cosR12 = 1.0;
			if (cosR12 <-1.0) cosR12 = -1.0;
			const double angleR12_1 = (edgeCase1 == SharedEdgeCase::Case2 ? 2 * M_PI - acos(cosR12) : acos(cosR12));
			printf("______[%.2f] Rotation matrix R12_1\n", angleR12_1*180.0 / M_PI);

			// 3. Find angles between basis2 and es
			Eigen::RowVector3d b21, b22;
			eVect = A.block(3 * SingNeighCC[id][i + 1], 2 * SingNeighCC[id][i + 1] + 0, 3, 1);
			b21 << eVect(0), eVect(1), eVect(2);
			eVect = A.block(3 * SingNeighCC[id][i + 1], 2 * SingNeighCC[id][i + 1] + 1, 3, 1);
			b22 << eVect(0), eVect(1), eVect(2);
			//cout << "______B21: " << b21 << ", B22: " << b22 << endl;

			// Basis 2, Frame 1
			double cosR21 = (b21.dot(es)) / (b21.norm()*es.norm());
			if (cosR21 > 1.0) cosR21 = 1.0;
			if (cosR21 < -1.0) cosR21 = -1.0;
			double angleR21_1 = (edgeCase2 == SharedEdgeCase::Case3 ? 2 * M_PI - acos(cosR21) : acos(cosR21));
			angleR21_1 = 2 * M_PI - angleR21_1;
			printf("______[%.2f] Rotation matrix R22_1 = [%.2f]\n", angleR21_1*180.0 / M_PI);
			
			const double RotAngle = (angleR12_1 + angleR21_1 > 2 * M_PI ? (angleR12_1 + angleR21_1) - 2 * M_PI : angleR12_1 + angleR21_1);
			const double cosBasis = cos(RotAngle);
			const double sinBasis = sin(RotAngle);
			printf("____ To map basis1 -> basis2: rotate by %.2f degree (cos=%.2f, cin=%.2f))\n", (RotAngle)*180.0 / M_PI, cosBasis, sinBasis);


			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, cosBasis));
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, -sinBasis));
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, -cosConst));
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, sinConst));
			c(counter) = 0.0;
			counter++;

			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, sinBasis));
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, cosBasis));
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, -sinConst));
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, -cosConst));
			c(counter) = 0.0;
			counter++;
		}
	}

	C.resize(2 * (globalConstraints.size() + numSingConstraints), B2D.rows());
	C.setFromTriplets(CTriplet.begin(), CTriplet.end());
}

void VectorFields::constructSoftConstraints()
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
	// Face
	constructCurvesAsConstraints(152474, 51474, aCurve);
	curvesConstraints[constCounter++] = aCurve;

	// Back
	constructCurvesAsConstraints(44109, 68907, aCurve);
	curvesConstraints[constCounter++] = aCurve;

	// body - bottom
	constructCurvesAsConstraints(13471, 195817, aCurve);
	curvesConstraints[constCounter++] = aCurve;

	// body - right
	constructCurvesAsConstraints(123036, 247143, aCurve);
	curvesConstraints[constCounter++] = aCurve;

	// body - left
	constructCurvesAsConstraints(234815, 232296, aCurve);
	curvesConstraints[constCounter++] = aCurve;

	// front_right_leg
	constructCurvesAsConstraints(75468, 7716, aCurve);
	curvesConstraints[constCounter++] = aCurve;

	// front_left_leg
	constructCurvesAsConstraints(231495, 77171, aCurve);
	curvesConstraints[constCounter++] = aCurve;

	// tail
	constructCurvesAsConstraints(230301, 113500, aCurve);
	curvesConstraints[constCounter++] = aCurve;

	/* Manual set-up for Armadillo */
	///// Head
	///constructCurvesAsConstraints(68818,6278, aCurve);
	///curvesConstraints[constCounter++] = aCurve;
	///// Stomach
	///constructCurvesAsConstraints(56965, 41616, aCurve);
	///curvesConstraints[constCounter++] = aCurve;
	///// Leg/Foot (R then L)
	///constructCurvesAsConstraints(28590, 16119, aCurve);
	///curvesConstraints[constCounter++] = aCurve;
	///constructCurvesAsConstraints(25037, 571, aCurve);
	///curvesConstraints[constCounter++] = aCurve;
	///// Arm/Hand
	///constructCurvesAsConstraints(55454, 6877, aCurve);
	///curvesConstraints[constCounter++] = aCurve;
	///constructCurvesAsConstraints(49059, 36423, aCurve);
	///curvesConstraints[constCounter++] = aCurve;
	///// Back
	///constructCurvesAsConstraints(68331, 72522, aCurve);
	///curvesConstraints[constCounter++] = aCurve;
	///// Tail
	///constructCurvesAsConstraints(24056, 1075, aCurve);
	///curvesConstraints[constCounter++] = aCurve;

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
	c.resize(2 * numConstraints);
	C.resize(2 * numConstraints, B2D.cols());

	int counter = 0;
	int elem;
	vector<Eigen::Triplet<double>> CTriplet; 
	for (int i = 0; i < curvesConstraints.size(); i++)
	{
		for (int j = 0; j < curvesConstraints[i].size() - 1; j++)
		{
			elem = curvesConstraints[i][j];
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * elem + 0, 1.0));
			c(counter++) = constraintVect2D[i][j](0);
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * elem + 1, 1.0));
			c(counter++) = constraintVect2D[i][j](1);
		}
	}

	C.setFromTriplets(CTriplet.begin(), CTriplet.end());
}

/* The path will be reversed, from init to end */
void VectorFields::constructCurvesAsConstraints(const int& init, const int& end, vector<int>& curve)
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

void VectorFields::projectCurvesToFrame()
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
		for (int j = 0; j < curvesConstraints[i].size()-1; j++)
		{
			face1 = curvesConstraints[i][j];
			face2 = curvesConstraints[i][j + 1];
			ALoc = A.block(3 * face1, 2 * face1, 3, 2);
			if (j < curvesConstraints[i].size() - 2)
			{
				face3 = curvesConstraints[i][j + 2];
				//face3 = curvesConstraints[i][curveSize-1];
				vec3D = (FC.row(face3) - FC.row(face1)).transpose();
			}
			else
			{
				vec3D = (FC.row(face2) - FC.row(face1)).transpose();
			}
			
			vec2D = ALoc.transpose() * vec3D;
			vec2D.normalize();
			constraintVect2D[i][j] = vec2D;
			//cout << "vec2D= " << vec2D << endl; 
		}
	}

	cout << "Fields are projected to 2d frame " << endl; 
}

//void VectorFields::constructSpecifiedConstraintsWithSingularities()  ==> Version 2.0
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
//	
//	int numSingConstraints = 0;
//	for (int i = 0; i < SingNeighCC.size(); i++) {
//		// Use only n-1 neighboring faces as constraints
//		for (int j = 0; j < (SingNeighCC[i].size()-1); j++) {
//			numSingConstraints++;
//		}
//	}
//
//	// Setting up matrix C and vector c
//	c.resize(2 * (globalConstraints.size()+numSingConstraints));
//
//
//	// HARD CONSTRAINTS
//	Eigen::SparseMatrix<double> CTemp;
//	vector<Eigen::Triplet<double>> CTriplet;
//	CTriplet.reserve(2 * globalConstraints.size() +  2 * 4 * 7 * SingNeighCC.size());
//	int counter = 0;
//	for (int i = 0; i < globalConstraints.size(); i++) {
//		// Matrix C
//		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * globalConstraints[i] + 0, 1.0));
//		c(counter++, 0) = sqrt(2.0);
//		//c(counter++, 1) = sqrt(2.0);
//		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * globalConstraints[i] + 1, 1.0));
//		c(counter++, 0) = sqrt(2.0);
//		//c(counter++, 1) = sqrt(2.0);
//	}
//
//	
//	// SINGULARITIES CONSTRAINTS
//	for (int id = 0; id < SingNeighCC.size(); id++) {
//		// For testing
//		sharedEdgesVect[id].resize(2*SingNeighCC[id].size()-2);
//
//		// 4. Compute rotation of among its valence
//		const double rotAngle = 2 * M_PI / (double)SingNeighCC[id].size();
//		const double cosA = cos(rotAngle);
//		const double sinA = sin(rotAngle);
//
//		// Which case? => determining which edge is the common edge
//		enum class SharedEdgeCase {Case1, Case2, Case3};
//		SharedEdgeCase edgeCase1, edgeCase2; 
//		for (int i = 0; i < (SingNeighCC[id].size() - 1); i++) {
//			// 1. Find shared edge (naively)
//			Eigen::RowVector3d es;
//			for (int f1 = 0; f1 < F.cols(); f1++) {
//				for (int f2 = 0; f2 < F.cols(); f2++) {
//					bool b1 = F(SingNeighCC[id][i], (f1+1)%F.cols()) == F(SingNeighCC[id][i + 1], f2);
//					bool b2 = F(SingNeighCC[id][i], f1) == F(SingNeighCC[id][i + 1], (f2 + 1) % F.cols());
//					if (b1 && b2) {
//						sharedEdgesVect[id][2 * i + 0] = F(SingNeighCC[id][i], f1);
//						sharedEdgesVect[id][2 * i + 1] = F(SingNeighCC[id][i], (f1 + 1) % F.cols());
//						es = V.row(F(SingNeighCC[id][i], (f1 + 1) % F.cols())) - V.row(F(SingNeighCC[id][i], f1));
//						printf("Shared edge=%d->%d\n", F(SingNeighCC[id][i], f1), F(SingNeighCC[id][i], (f1 + 1) % F.cols()));
//						
//						if (f1 == 0)		edgeCase1 = SharedEdgeCase::Case1;	// => edge V0->V1 is the shared edge => it takes 0 step to reach v0
//						else if(f1==1)		edgeCase1 = SharedEdgeCase::Case3;	// => edge V1->V2 is the shared edge => it takes 2 step to reach v0
//						else if(f1==2)		edgeCase1 = SharedEdgeCase::Case2;	// => edge V2->V0 is the shared edge => it takes 1 step to reach v0
//
//						if (f2 == 0)		edgeCase2 = SharedEdgeCase::Case1;
//						else if (f2 == 1)	edgeCase2 = SharedEdgeCase::Case3;
//						else if (f2 == 2)	edgeCase2 = SharedEdgeCase::Case2;
//					}
//				}
//			}
//			// 2. Find angles between basis1 and es
//			Eigen::VectorXd eVect;
//			Eigen::RowVector3d b11, b12;
//			eVect = A.block(3 * SingNeighCC[id][i], 2 * SingNeighCC[id][i]+0, 3, 1);
//			b11 << eVect(0), eVect(1), eVect(2);
//			eVect = A.block(3 * SingNeighCC[id][i], 2 * SingNeighCC[id][i]+1, 3, 1);
//			b12 << eVect(0), eVect(1), eVect(2);
//			cout << "______B11: " << b11 << ", B12: " << b12 << endl;
//
//			// Basis 1, Frame 1
//			double cosR12 = (b11.dot(es)) / (b11.norm()*es.norm());
//			if (cosR12 > 1.0) cosR12 = 1.0; 
//			if (cosR12 <-1.0) cosR12 = -1.0;
//			const double angleR12_1 = (edgeCase1==SharedEdgeCase::Case2 ? 2*M_PI - acos(cosR12) : acos(cosR12));
//			const double cosR12_1 = cos(angleR12_1);
//			const double sinR12_1 = sin(angleR12_1);
//			printf("______[%.2f] Rotation matrix R12_1 = [%.3f,%.3f; %.3f, %.3f]\n", angleR12_1*180.0/M_PI, cosR12_1, -sinR12_1, sinR12_1, cosR12_1);
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
//			eVect = A.block(3 * SingNeighCC[id][i+1], 2 * SingNeighCC[id][i+1] + 0, 3, 1);
//			b21 << eVect(0), eVect(1), eVect(2);
//			eVect = A.block(3 * SingNeighCC[id][i+1], 2 * SingNeighCC[id][i+1] + 1, 3, 1);
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
//			printf("______[%.2f] Rotation matrix R22_1 = [%.2f,%.2f; %.2f, %.2f]\n", angleR21_1*180.0/M_PI, cosR21_1, -sinR21_1, sinR21_1, cosR21_1);
//
//			// Basis 2, Frame 2
//			cosR21 = (b22.dot(es)) / (b22.norm()*es.norm());
//			if (cosR21 > 1.0) cosR21 = 1.0;
//			if (cosR21 < -1.0) cosR21 = -1.0;
//			const double angleR21_2 = (edgeCase2 == SharedEdgeCase::Case1 ? 2 * M_PI - acos(cosR21) : acos(cosR21));
//			const double cosR21_2 = cos(angleR21_2);
//			const double sinR21_2 = sin(angleR21_2);
//			printf("______[%.2f] Rotation matrix R22_2 = [%.2f,%.2f; %.2f, %.2f]\n", angleR21_2*180.0/M_PI, cosR21_2, -sinR21_2, sinR21_2, cosR21_2);
//
//			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, cosR12_1));
//			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, -sinR12_1));
//			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, -cosR21_1));
//			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, sinR21_1));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, cosR21_1*cosR12_1+sinR21_1*sinR12_1));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, -cosR21_1*sinR12_1+sinR21_1*cosR12_1));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, -1.0));
//			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, 0.0));
//			c(counter) = 0.0;
//			counter++;
//			
//			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, sinR12_1));
//			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, cosR12_1 ));
//			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, -sinR21_1));
//			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, -cosR21_1 ));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, -sinR21_1*cosR12_1+cosR21_1*sinR12_1));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, sinR21_1*sinR12_1+cosR21_1*cosR21_1));
//			//CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, 0.0));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, -1.0));
//			c(counter) = 0.0;
//			counter++;
//		}
//	}
//
//	C.resize(2 * (globalConstraints.size()+numSingConstraints), B2D.rows());
//	C.setFromTriplets(CTriplet.begin(), CTriplet.end());
//	//printf("Cp=%dx%d\n", C.rows(), C.cols());	
//}

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
void VectorFields::setupGlobalProblem(const Eigen::Vector3d& lambda)
{	
	Eigen::VectorXd					b, g, h, vEst;
	Eigen::SparseMatrix<double>		A_LHS;
	//Eigen::VectorXd					vEst;
	//double lambda2 = 0.4; 
	//Eigen::Vector3d lambda;
	//lambda(0) = 1.0; // 100 * MF2D.coeff(0, 0) / SF2D.coeff(0, 0);		// on harmonic energy
	//lambda(1) = 0.1; // 100 * MF2D.coeff(0, 0) / B2D.coeff(0, 0);		// on bi-harmonic energy
	//lambda(2) = 0.4;											// on the constraint
	
	constructConstraints();
	//setupRHSGlobalProblemMapped(g, h, vEst, b);
	//setupLHSGlobalProblemMapped(A_LHS);
	//solveGlobalSystemMappedLDLT(vEst, A_LHS, b);
	//solveGlobalSystemMappedLU_GPU();

	setupRHSGlobalProblemSoftConstraints(lambda, b);
	setupLHSGlobalProblemSoftConstraints(lambda, A_LHS);		
	solveGlobalSystemMappedLDLTSoftConstraints(A_LHS, b);
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

void VectorFields::setupRHSGlobalProblemSoftConstraints(const Eigen::Vector3d& lambda, Eigen::VectorXd& b)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Setting up the RHS of the system... ";

	b = lambda(2) * C.transpose() * c;

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}

void VectorFields::setupLHSGlobalProblemSoftConstraints(const Eigen::Vector3d& lambda, Eigen::SparseMatrix<double>& A_LHS)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Setting up the LHS of the system... ";

	//A_LHS = SF2D + lambda*(C.transpose()*C);
	A_LHS = lambda(0)*SF2D + lambda(1)*B2D + lambda(2)*C.transpose()*C;

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


	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
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

void VectorFields::solveGlobalSystemMappedLDLTSoftConstraints(Eigen::SparseMatrix<double>& A_LHS, Eigen::VectorXd& b)
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
		if (sparseSolver.info() == Eigen::InvalidInput)
			cout << "Input is Invalid. " << endl;
		cout << sparseSolver.info() << endl;
		return;
	}

	Xf = x;

	printf("____Xf size is %dx%d\n", Xf.rows(), Xf.cols());


	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}

// RANK-2 TENSOR
void VectorFields::constructMappingMatrix_TensorR2()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Constructing Mapping matrices (Global/World-Coord to Local Frame)... ";


	AT2R.resize(3 * F.rows(), 3 * F.rows());
	vector<Eigen::Triplet<double>> ATriplet;
	ATriplet.reserve(3 * 3 * F.rows());
	Eigen::Vector3d e, f, n;
	Eigen::Vector3d eeT, efT, feT, efTfeT, ffT;

	for (int i = 0; i < F.rows(); i++) {
		/* Computing the basic elements */
		e = V.row(F(i, 1)) - V.row(F(i, 0));
		e.normalize();

		n = NF.row(i);
		n.normalize();

		f = n.cross(e);
		f.normalize();

		/* Computing the values for tensor */


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

void VectorFields::constructStiffnessMatrixSF2D_TensorR2(Eigen::SparseMatrix<double>& LapCurl3D, Eigen::SparseMatrix<double>& LapCurl2D, Eigen::SparseMatrix<double>& LapDiv3D, Eigen::SparseMatrix<double>& LapDiv2D)
{

}

void VectorFields::constructStiffnessMatrixCurlPart2D_TensorR2(Eigen::SparseMatrix<double>& LapCurl3D, Eigen::SparseMatrix<double>& LapCurl2D)
{

}

void VectorFields::constructStiffnessMatrixDivPart2D_TensorR2(Eigen::SparseMatrix<double>& LapDiv3D, Eigen::SparseMatrix<double>& LapDiv2D)
{

}

// APPLICATIONS ON GLOBAL SYSTEM
void VectorFields::computeSmoothing(const double& mu, const Eigen::VectorXd& v_in, Eigen::VectorXd& v_out)
{
	/* First flavour */
	//Eigen::SparseMatrix<double> A = MF2D + mu*B2D;
	Eigen::SparseMatrix<double> A = MF2D + mu*SF2D;

	/* Second flavour */
	//Eigen::SparseMatrix<double> A = MF2D + mu*SF2D*MF2Dinv*SF2D;
	Eigen::VectorXd b = MF2D*v_in;

	Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> sparseSolver(A);
	v_out = sparseSolver.solve(b);

	double in_length = v_in.transpose()*MF2D*v_in;
	double out_length = v_out.transpose()*MF2D*v_out;
	cout << "IN Length= " << in_length << endl;
	cout << "OUT length " << out_length << endl; 
	

	/* Computing the L2-norm of the smoothed fields */
	double diff1 = (v_out - v_in).transpose()*MF2D*(v_out - v_in);
	double diff2 = v_in.transpose()*MF2D*v_in;
	double sqrt_norm = sqrt(diff1 / diff2);
	printf("The diff of v_out and v_in is %.10f \n", sqrt_norm);

	/* Computing the energy */
	double energy1 = v_in.transpose() * ((B2D) * v_in);
	double energy2 = v_out.transpose() * ((B2D) * v_out);
	printf("The energy is=%.4f ==> %.4f.\n", energy1, energy2);
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
	//Sample[0] = rand() % F.rows();
	Sample[0] = 0;
	//Sample[0] = 70267; // Arma 43k
	//Sample[0] = 5461;	// For Armadilo of 10k vertices

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

	sampleDistance = D; 
}

void VectorFields::constructBasis()
{
	// Select the types of basis construction
	Eigen::SparseMatrix<double> BasisFunctions;

	constructBasis_LocalEigenProblem();
	//constructBasis_LocalEigenProblem10();
	//constructBasis_OptProblem();
	//constructBasis_GradOfLocalFunction(BasisFunctions);
	//constructBasis_EigenPatch(BasisFunctions);
}

void VectorFields::constructBasis_LocalEigenProblem()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Constructing Basis...\n";

	double	coef = sqrt(pow(1.5, 2) + pow(1.5, 2));
	double distRatio = coef * sqrt((double)V.rows() / (double)Sample.size());

	// Setup sizes of each element to construct basis
	try {
		BasisTemp.resize(2 * F.rows(), 2 * Sample.size());
	}
	catch (string &msg) {
		cout << "Cannot allocate memory for basis.." << endl;
	}

	Basis.resize(BasisTemp.rows(), BasisTemp.cols());
	vector<vector<Eigen::Triplet<double>>> UiTriplet(Sample.size());

	/* Set-UP Laplace Matrix */
	Eigen::SparseMatrix<double> LapForBasis = MF2Dinv * SF2DAsym;

	cout << "....Constructing and solving local systems...";
	const int NUM_PROCESS = 4;
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

	omp_set_num_threads(1);
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

		printf("num threads=%d, iproc=%d, ID=%d, start=%d, to end=%d, num els=%d\n", ntids, iproc, tid, istart, istart + ipts, ipts);

		Eigen::VectorXd				D(F.rows());
		for (int i = 0; i < F.rows(); i++) {
			D(i) = numeric_limits<double>::infinity();
		}

		//cout << "[" << tid << "] Number of processors " << iproc << ", with " << ntids << " threads." << endl;

		UiTriplet[tid].reserve(2.0 * ((double)ipts / (double)Sample.size()) * 2 * 10.0 * F.rows());

		// Computing the values of each element
		for (id = istart; id < (istart + ipts); id++) {
		//for (id = istart; id < (istart + ipts) && id < 10; id++) {
			if (id >= Sample.size()) break;

			vector<Eigen::Triplet<double>> BTriplet, C1Triplet, C2Triplet;

			LocalFields localField(id);
			t1 = chrono::high_resolution_clock::now();
			//localField.constructSubdomain(Sample[id], V, F, avgEdgeLength, AdjMF3N, distRatio);
			localField.constructSubdomain(Sample[id], V, F, avgEdgeLength, AdjMF2Ring, distRatio);
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
			localField.constructLocalEigenProblemWithSelector(SF2DAsym, MF2D, AdjMF2Ring, 2, doubleArea, UiTriplet[id]);
			//localField.constructLocalEigenProblem(SF2D, AdjMF3N, doubleArea, UiTriplet[id]);
			t2 = chrono::high_resolution_clock::now();
			durations[3] += t2 - t1;


			if (id == 0)
			{
				SubDomain = localField.SubDomain;
				Boundary = localField.Boundary;
				//patchDijkstraDist = localField.dijksFaceDistMapped;
			}

			// To get local elements for visualizing subdomain
			if (id == 0 || id == 46) {
				cout << "Getting element of ID " << id << endl;

				for (int fid : localField.SubDomain) {
					localSystem(fid) = 0.3;
				}

				for (int fid : localField.Boundary) {
					localSystem(fid) = 0.7;
				}

				//localSystem(localField.sampleID) = 1.0;


			}

			/* Localized eigenproblems */
			//if (id == 15)
			{
				//Eigen::SparseMatrix<double> MTempStiff;
				//localField.obtainLocalMatrixPatch2D(SF2D, MTempStiff);
				//visualizeSparseMatrixInMatlab(MTempStiff);
				//localField.constructLocalEigenProblem(SF2D, AdjMF2Ring, doubleArea, eigFieldsLocal);
				//localField.constructLocalEigenProblemWithSelector(SF2D, AdjMF2Ring, doubleArea, eigFieldsLocal);
			}

		}
	}
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;

	cout << "....Gathering local elements as basis matrix... ";
	t1 = chrono::high_resolution_clock::now();
	gatherBasisElements(UiTriplet,2);
	Basis = BasisTemp;
	//normalizeBasisAbs();

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
	printf("....[3] Solving local eigenvalue problems: %.8f seconds.\n", durations[3].count());

	// Information about Basis
	printf("> Basis Structure information \n");
	printf("....Size = %dx%d\n", Basis.rows(), Basis.cols());
	printf("....NNZ per row = %.2f\n", (double)Basis.nonZeros() / (double)Basis.rows());
}

void VectorFields::constructBasis_OptProblem()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Constructing Basis...\n";

	double	coef = sqrt(pow(1.3, 2) + pow(1.1, 2));
	double distRatio = coef * sqrt((double)V.rows() / (double)Sample.size());

	// Setup sizes of each element to construct basis
	try {
		BasisTemp.resize(2 * F.rows(), 2 * Sample.size());
	}
	catch (string &msg) {
		cout << "Cannot allocate memory for basis.." << endl;
	}

	Basis.resize(BasisTemp.rows(), BasisTemp.cols());
	vector<vector<Eigen::Triplet<double>>> UiTriplet(Sample.size());

	/* Set-UP Laplace Matrix */
	Eigen::SparseMatrix<double> LapForBasis = MF2Dinv * SF2DAsym;

	cout << "....Constructing and solving local systems...";
	const int NUM_PROCESS = 8;
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

	//omp_set_num_threads(1);
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

		printf("num threads=%d, iproc=%d, ID=%d, start=%d, to end=%d, num els=%d\n", ntids, iproc, tid, istart, istart + ipts, ipts);

		Eigen::VectorXd				D(F.rows());
		for (int i = 0; i < F.rows(); i++) {
			D(i) = numeric_limits<double>::infinity();
		}

		UiTriplet[tid].reserve(2.0 * ((double)ipts / (double)Sample.size()) * 2 * 10.0 * F.rows());

		// Computing the values of each element
		for (id = istart; id < (istart + ipts); id++) {
			if (id >= Sample.size()) break;

			vector<Eigen::Triplet<double>> BTriplet, C1Triplet, C2Triplet;

			LocalFields localField(id);
			t1 = chrono::high_resolution_clock::now();
			//localField.constructSubdomain(Sample[id], V, F, avgEdgeLength, AdjMF3N, distRatio);
			localField.constructSubdomain(Sample[id], V, F, avgEdgeLength, AdjMF2Ring, distRatio);
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
			//localField.constructMatrixBLocal(B2D, AdjMF2Ring, BTriplet);			
			t2 = chrono::high_resolution_clock::now();
			durations[3] += t2 - t1;

			t1 = chrono::high_resolution_clock::now();
			localField.constructLocalConstraints(C1Triplet, C2Triplet);
			//localField.constructLocalConstraintsWithLaplacian(doubleArea, AdjMF2Ring, SF2D, C1Triplet, C2Triplet);
			///localField.constructLocalConstraintsWithLaplacian(doubleArea, LapForBasis, C1Triplet, C2Triplet);
			t2 = chrono::high_resolution_clock::now();
			durations[4] += t2 - t1;

			t1 = chrono::high_resolution_clock::now();
			localField.setupRHSLocalProblemMapped();
			t2 = chrono::high_resolution_clock::now();
			durations[5] += t2 - t1;

			t1 = chrono::high_resolution_clock::now();
			//localField.setupLHSLocalProblemMapped();
			localField.setupLHSLocalProblemMapped(BTriplet, C1Triplet, C2Triplet);
			t2 = chrono::high_resolution_clock::now();
			durations[6] += t2 - t1;

			localField.computeDijkstraFaceDistance(V, F, FC, AdjMF3N);

			t1 = chrono::high_resolution_clock::now();
			localField.solveLocalSystemMappedLDLT(UiTriplet[id]);
			t2 = chrono::high_resolution_clock::now();
			durations[7] += t2 - t1;


			/* To get local elements for visualizing subdomain */
			//if (id == 0 || id == 46) {
			//	cout << "Getting element of ID " << id << endl;
			//
			//	for (int fid : localField.SubDomain) {
			//		localSystem(fid) = 0.3;
			//	}
			//
			//	for (int fid : localField.Boundary) {
			//		localSystem(fid) = 0.7;
			//	}/
			//	//localSystem(localField.sampleID) = 1.0;/
			//}
		}
	}
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;

	cout << "....Gathering local elements as basis matrix... ";
	t1 = chrono::high_resolution_clock::now();
	gatherBasisElements(UiTriplet, 2);
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;

	cout << "....Partition of unity of the basis matrix... ";
	t1 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	//Basis = BasisTemp;
	//normalizeBasis();
	normalizeBasisAbs(2);

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
	printf("....NNZ per row = %.2f\n", (double)Basis.nonZeros() / (double)Basis.rows());
}


void VectorFields::constructBasis_GradOfLocalFunction(Eigen::SparseMatrix<double>& BasisFunctions)
{
	//string filename = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/MATLAB Implementation/Data/CDragon_Basis_1000_full.mat";
	//string filename = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/MATLAB Implementation/Data/CDragon_Basis_1000_full.mat";
	string filename = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/MATLAB Implementation/Data/Cube_Sharp_Basis_1000.txt";
	Eigen::SparseMatrix<double> BasisGrad, BasisCoGrad;
	Eigen::MatrixXd BasisTemp;
	

	/* Obtaining matrix from old data */
	if (false) {			/* MATLAB DATA*/
	//ReadSparseMatrixFromMatlab(BasisFunctions, filename);
		ReadDenseMatrixFromMatlab(BasisTemp, filename);
		int basCols = min((int)BasisTemp.cols(), 1000);

		/* Storing the data as sparse matrix */
		vector<Eigen::Triplet<double>> ATriplet;
		ATriplet.reserve(40 * BasisTemp.rows());
		for (int j = 0; j < basCols; j++)
		{
			for (int i = 0; i < BasisTemp.rows(); i++)
			{
				if (BasisTemp(i, j) > 0.0000000001)
				{
					ATriplet.push_back(Eigen::Triplet<double>(i, j, -BasisTemp(i, j)));
				}
			}
		}
		BasisFunctions.resize(BasisTemp.rows(), basCols);	
		BasisFunctions.setFromTriplets(ATriplet.begin(), ATriplet.end());
	}
	else 
	{
		readEigenSparseMatrixFromBinary(filename, BasisFunctions);
	}
	printf("Size of the basis function=%dx%d; ", BasisFunctions.rows(), BasisFunctions.cols());
	printf("__with %d nonzeros (%.5f nnz per row) \n", BasisFunctions.nonZeros(), (double)BasisFunctions.nonZeros() / (double)BasisFunctions.rows());

	/* Computing the graident and co-gradient of the basis functions */
	BasisGrad = A.transpose() * GF3D * BasisFunctions;
	BasisCoGrad = A.transpose() * J3D * GF3D * BasisFunctions;
	printf("Grad =%dx%d\n", BasisGrad.rows(), BasisGrad.cols());
	printf("Co-Grad =%dx%d\n", BasisCoGrad.rows(), BasisCoGrad.cols());

	/* Storing them as new basis */
	vector<Eigen::Triplet<double>> BTriplet;
	BTriplet.reserve(BasisGrad.nonZeros() + BasisCoGrad.nonZeros());

	/* Getting the elements of the Gradient part */
	for (int i = 0; i < BasisGrad.outerSize(); i++) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(BasisGrad, i); it; ++it) {
			BTriplet.push_back(Eigen::Triplet<double>(it.row(), 2 * it.col(), it.value()));
		}
	}

	/* Getting the elements of the Co-Gradient part */
	for (int i = 0; i < BasisCoGrad.outerSize(); i++) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(BasisCoGrad, i); it; ++it) {
			BTriplet.push_back(Eigen::Triplet<double>(it.row(), 2 * it.col() + 1, it.value()));
		}
	}

	Basis.resize(BasisGrad.rows(), BasisGrad.cols() + BasisCoGrad.cols());
	Basis.setFromTriplets(BTriplet.begin(), BTriplet.end());
	printf("Size of the BASIS function=%dx%d; ", Basis.rows(), Basis.cols());
	printf("__with %d nonzeros (%.5f nnz per row\n", Basis.nonZeros(), (double)Basis.nonZeros() / (double)Basis.rows());

	//normalizeBasisAbs(2);
}

void VectorFields::constructBasis_EigenPatch(Eigen::SparseMatrix<double>& BasisFunctions)
{
	/* Check the basis function*/
	if (BasisFunctions.cols() < 1)
	{
		constructBasis_GradOfLocalFunction(BasisFunctions);
	}

	/* Check if we have eigenfunctions already */
	if (eigFieldFull2D.cols() < 1)
	{
		computeEigenMatlab(SF2DAsym, MF2D, 2, eigFieldFull2D, eigValuesFull, "hello");
		string filename = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/Cube_Round_2_Full_eigFields.mat";
		//ReadDenseMatrixFromMatlab(eigFieldFull2D, filename);
	}

	/* Setting up the matrices */
	vector<Eigen::Triplet<double>> BTriplet;
	BTriplet.reserve(2 * BasisFunctions.nonZeros());
	//Basis.resize()

	/* Get the mapping matrix */
	vector<set<int>> VFMap(V.rows());
	for (int i = 0; i < F.rows(); i++)
	{
		for (int j = 0; j < F.cols(); j++)
		{
			VFMap[F(i, j)].insert(i);
		}
	}

	/* Showing the result of mapping*/
	//for (int i = 0; i < 100; i++)
	//{
	//	printf("V[%d] :", i);
	//	for (int j : VFMap[i])
	//	{
	//		printf("%d |", j);
	//	}
	//	printf("\n");
	//}

	for (int i = 0; i < BasisFunctions.outerSize(); i++) 
	{
		/* Throwing out values at each vertex to be at every face */
		vector<set<double>> BasisWeight(F.rows());
		for (Eigen::SparseMatrix<double>::InnerIterator it(BasisFunctions, i); it; ++it) {
			for (int j : VFMap[it.row()])
			{
				BasisWeight[j].insert(it.value());
			}
		}
		
		/*Setting the weight on each face * eigenVectors on each face */
		for (int j = 0; j < F.rows(); j++)
		{
			if (BasisWeight[j].size() < 3)continue;

			double sum = 0;
			for (double k : BasisWeight[j])
			{
				sum += k;
			}
			double weight = sum / 3.0; 

			//if (i == 0)
			//{
			//	printf("Weight of face %d = %.5f, eVect=[%.3f;%.3f] => [%.3f;%.3f]\n", j, weight, eigFieldFull2D(2 * j + 0, 0), eigFieldFull2D(2 * j + 0, 1), weight*eigFieldFull2D(2 * j + 0, 0), weight*eigFieldFull2D(2 * j + 0, 1));
			//}

			BTriplet.push_back(Eigen::Triplet<double>(2 * j + 0, 2 * i + 0, weight * eigFieldFull2D(2 * j + 0, 0)));
			BTriplet.push_back(Eigen::Triplet<double>(2 * j + 1, 2 * i + 0, weight * eigFieldFull2D(2 * j + 1, 0)));
			BTriplet.push_back(Eigen::Triplet<double>(2 * j + 0, 2 * i + 1, weight * eigFieldFull2D(2 * j + 0, 1)));
			BTriplet.push_back(Eigen::Triplet<double>(2 * j + 1, 2 * i + 1, weight * eigFieldFull2D(2 * j + 1, 1)));
		}
	}

	Basis.resize(0, 0);
	Basis.resize(2 * F.rows(), 2 * BasisFunctions.cols());
	Basis.setFromTriplets(BTriplet.begin(), BTriplet.end());	
}

void VectorFields::constructBasis_LocalEigenProblem10()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Constructing Basis...\n";

	const int EIG_NUM = 10; 
	double	coef = sqrt(pow(0.5, 2) + pow(0.7, 2));
	double distRatio = coef * sqrt((double)V.rows() / (double)Sample.size());

	// Setup sizes of each element to construct basis
	try {
		BasisTemp.resize(2 * F.rows(), EIG_NUM * Sample.size());
		printf("Basis size = %dx%d\n", BasisTemp.rows(), BasisTemp.cols());
	}
	catch (string &msg) {
		cout << "Cannot allocate memory for basis.." << endl;
	}

	Basis.resize(BasisTemp.rows(), BasisTemp.cols());
	vector<vector<Eigen::Triplet<double>>> UiTriplet(Sample.size());

	/* Set-UP Laplace Matrix */
	Eigen::SparseMatrix<double> LapForBasis = MF2Dinv * SF2DAsym;

	cout << "....Constructing and solving local systems...";
	const int NUM_PROCESS = 4;
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

	omp_set_num_threads(1);
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

		printf("num threads=%d, iproc=%d, ID=%d, start=%d, to end=%d, num els=%d\n", ntids, iproc, tid, istart, istart + ipts, ipts);

		Eigen::VectorXd				D(F.rows());
		for (int i = 0; i < F.rows(); i++) {
			D(i) = numeric_limits<double>::infinity();
		}

		UiTriplet[tid].reserve(2.0 * ((double)ipts / (double)Sample.size()) * EIG_NUM * 10.0 * F.rows());

		// Computing the values of each element
		for (id = istart; id < (istart + ipts); id++) {
		//for (id = istart; id < (istart + ipts) && id < 10; id++) {
			if (id >= Sample.size()) break;

			vector<Eigen::Triplet<double>> BTriplet, C1Triplet, C2Triplet;

			LocalFields localField(id);
			t1 = chrono::high_resolution_clock::now();
			//localField.constructSubdomain(Sample[id], V, F, avgEdgeLength, AdjMF3N, distRatio);
			localField.constructSubdomain(Sample[id], V, F, avgEdgeLength, AdjMF2Ring, distRatio);
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
			localField.constructLocalEigenProblemWithSelector(SF2DAsym, MF2D, AdjMF2Ring, EIG_NUM, doubleArea, UiTriplet[id]);
			//localField.constructLocalEigenProblem(SF2D, AdjMF3N, doubleArea, UiTriplet[id]);
			t2 = chrono::high_resolution_clock::now();
			durations[3] += t2 - t1;


			if (id == 0)
			{
				SubDomain = localField.SubDomain;
				Boundary = localField.Boundary;
				//patchDijkstraDist = localField.dijksFaceDistMapped;
			}

			// To get local elements for visualizing subdomain
			if (id == 0 || id == 46) {
				cout << "Getting element of ID " << id << endl;

				for (int fid : localField.SubDomain) {
					localSystem(fid) = 0.3;
				}

				for (int fid : localField.Boundary) {
					localSystem(fid) = 0.7;
				}			

			}
		}
	}
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;

	cout << "....Gathering local elements as basis matrix... ";
	t1 = chrono::high_resolution_clock::now();
	gatherBasisElements(UiTriplet, EIG_NUM);
	Basis = BasisTemp;
	//normalizeBasis();
	//normalizeBasisAbs(EIG_NUM);

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
	printf("....[3] Solving local eigenvalue problems: %.8f seconds.\n", durations[3].count());

	// Information about Basis
	printf("> Basis Structure information \n");
	printf("....Size = %dx%d\n", Basis.rows(), Basis.cols());
	printf("....NNZ per row = %.2f\n", (double)Basis.nonZeros() / (double)Basis.rows());
}


void VectorFields::constructBasisEigenVects()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Constructing Basis...\n";

	double	coef = sqrt(pow(1.1, 2) + pow(1.3, 2));
	double distRatio = coef * sqrt((double)V.rows() / (double)Sample.size());

	// Setup sizes of each element to construct basis
	try {
		BasisTemp.resize(2 * F.rows(), 2 * Sample.size());
	}
	catch (string &msg) {
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
		tid = omp_get_thread_num();
		ntids = omp_get_num_threads();
		ipts = (int)ceil(1.00*(double)Sample.size() / (double)ntids);
		istart = tid * ipts;
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

			vector<Eigen::Triplet<double>> BTriplet, C1Triplet, C2Triplet;

			LocalFields localField(id);
			t1 = chrono::high_resolution_clock::now();
			localField.constructSubdomain(Sample[id], V, F, avgEdgeLength, AdjMF3N, distRatio);
			t2 = chrono::high_resolution_clock::now();
			durations[0] += t2 - t1;
			
			t1 = chrono::high_resolution_clock::now();
			localField.constructLocalElements(F);
			t2 = chrono::high_resolution_clock::now();
			durations[2] += t2 - t1;			

			// Get the patch for id 0
			if(id==0)
				localPatchElements = localField.InnerElements;

			t1 = chrono::high_resolution_clock::now();
			//localField.solveLocalSystemMappedLDLT(UiTriplet[id]);
			//localField.constructLocalEigenProblem(SF2D, AdjMF2Ring, doubleArea, UiTriplet[id]);
			localField.constructLocalEigenProblem(SF2D, AdjMF2Ring, doubleArea, UiTriplet[id]);
			//localField.constructLocalEigenProblem(SF2DAsym, AdjMF3N, doubleArea, UiTriplet[id]);
			t2 = chrono::high_resolution_clock::now();
			durations[7] += t2 - t1;
			//cout << "System " << id << " ( " << XfLoc.rows() << ") is solved." << endl; 
			//printf("System %d (%d) is solved.\n", id, XfLoc.rows());


			//localField.measureXF(doubleArea, J);

		}
	}
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;

	cout << "....Gathering local elements as basis matrix... ";
	t1 = chrono::high_resolution_clock::now();
	gatherBasisElements(UiTriplet, 2);
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;

	cout << "....Partition of unity of the basis matrix... ";
	t1 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	//Basis = BasisTemp; 
	//normalizeBasis();
	normalizeBasisAbs(2);

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
	printf("....[7] Solving Eigenvalue problem: %.8f seconds.\n", durations[7].count());

	// Information about Basis
	printf("> Basis Structure information \n");
	printf("....Size = %dx%d\n", Basis.rows(), Basis.cols());
	printf("....NNZ per row = %.2f\n", (double)Basis.nonZeros() / (double)Basis.rows());
}

void VectorFields::gatherBasisElements(const vector<vector<Eigen::Triplet<double>>> &UiTriplet, const int& NUM_EIGEN)
{
	/* Set the basis sum to be zero on each column-group */
	vector<Eigen::Triplet<double>> BTriplet;
	BasisSum.resize(2 * F.rows(), NUM_EIGEN);
	for (int i = 0; i < BasisSum.rows(); i++) 
	{
		for (int j = 0; j < NUM_EIGEN; j++)
		{
			BasisSum(i, j) = 0.0;
		}
	}

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
		std::copy(UiTriplet[j].begin(), UiTriplet[j].end(), BTriplet.begin() + tripSize);
	}	
	BasisTemp.setFromTriplets(BTriplet.begin(), BTriplet.end());

	/* Computing the basis sum -> cor normalization */
	for (int k = 0; k < BasisTemp.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(BasisTemp, k); it; ++it) {
			BasisSum(it.row(), it.col() % NUM_EIGEN) += it.value();
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

void VectorFields::normalizeBasisAbs(const int& stride)
{
	Eigen::MatrixXd normSum(F.rows(), stride), normSumN(F.rows(), stride);
	BasisSumN.resize(BasisTemp.rows(), stride);
	vector<Eigen::Triplet<double>> BNTriplet;
	BNTriplet.reserve(BasisTemp.nonZeros());

	Eigen::MatrixXd BasisNorm(F.rows(), stride);

	for (int i = 0; i < normSum.rows(); i++) {
		for (int j = 0; j < normSum.cols(); j++) {
			normSum(i, j) = 0.0;
			//normSumN(i, j) = 0.0;
			//BasisSumN(2 * i + 0, j) = 0.0;
			//BasisSumN(2 * i + 1, j) = 0.0;
		}
	}

	// Getting the sum of norm on each pair on each frame 
	for (int k = 0; k < BasisTemp.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(BasisTemp, k); it; ++it) {
			if (it.row() % 2 == 1) continue;

			double a = it.value();
			double b = BasisTemp.coeff(it.row() + 1, it.col());
			double norm = sqrt(a*a + b*b);
			normSum(it.row() / 2, it.col() % stride) += norm;
		}
	}
	
	// Normalize the system
	for (int k = 0; k < BasisTemp.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(BasisTemp, k); it; ++it) {
			double newValue = it.value() / normSum(it.row()/2, it.col()%stride);
			// newValue *= (2.0);
			// To have the basis with norm 2.0
			BNTriplet.push_back(Eigen::Triplet<double>(it.row(), it.col(), newValue));
			BasisSumN(it.row(), it.col() % stride) += newValue;
		}
	}	

	Basis.setFromTriplets(BNTriplet.begin(), BNTriplet.end());	
}

void VectorFields::storeBasis(const string& filename)
{
	writeEigenSparseMatrixToBinary(Basis, filename);
}

void VectorFields::retrieveBasis(const string& filename)
{
	readEigenSparseMatrixFromBinary(filename, Basis);

	//BasisTemp = Basis; 
	//normalizeBasisAbs();
	printf("Basis size=%dx%d (nnz=%.4f)\n", Basis.rows(), Basis.cols(), (double)Basis.nonZeros()/(double)Basis.rows());
}


void VectorFields::setAndSolveUserSystem(const Eigen::Vector3d& lambda)
{
	// Declare function-scoped variables
	Eigen::VectorXd					bBar, gBar, hBar, vEstBar;
	Eigen::SparseMatrix<double>		A_LHSBar;
	//const double lambda2 = 0.4; 
	//Eigen::Vector3d lambda;
	//lambda(0) = 1.0;  // 100 * MF2D.coeff(0, 0) / SF2D.coeff(0, 0);		// on harmonic energy
	//lambda(1) = 0.1; // 100 * MF2D.coeff(0, 0) / B2D.coeff(0, 0);		// on bi-harmonic energy
	//lambda(2) = 0.4;											// on the constraint

	//setupReducedBiLaplacian();
	getUserConstraints();
	//setupRHSUserProblemMapped(gBar, hBar, vEstBar, bBar);
	//setupLHSUserProblemMapped(A_LHSBar);
	//solveUserSystemMappedLDLT(vEstBar, A_LHSBar, bBar);
	setupRHSUserProblemMappedSoftConstraints(lambda, bBar);
	setupLHSUserProblemMappedSoftConstraints(lambda, A_LHSBar);
	solveUserSystemMappedLDLTSoftConstraints(A_LHSBar, bBar);

	mapSolutionToFullRes();
}

void VectorFields::setupReducedBiLaplacian()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Computign Reduced Bi-Laplacian...";

	B2DBar = Basis.transpose() * B2D * Basis; 

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;
	printf(".... Local Basis = %dx%d\n", B2DBar.rows(), B2DBar.cols());

	/* Getting the information about nonzeros */
	double nnz_num = (double)B2DBar.nonZeros() / (double)B2DBar.rows();
	double nnz_perc = nnz_num / (double)B2DBar.cols();
	printf(".... NNZ per row = %.2f\n", nnz_num);
	printf(".... Percentage of NNZ = %.20f \n", nnz_perc*100);

}
void VectorFields::getUserConstraints()
{	
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Obtaining user constraints ";

	constructConstraints();

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

void VectorFields::setupRHSUserProblemMappedSoftConstraints(const Eigen::Vector3d& lambda, Eigen::VectorXd& bBar)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Constructing RHS (mapped)...";

	bBar = lambda(2)*CBar.transpose() * cBar; 

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;
}

void VectorFields::setupLHSUserProblemMappedSoftConstraints(const Eigen::Vector3d& lambda, Eigen::SparseMatrix<double>& A_LHSBar)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Constructing LHS (mapped)...";

	Eigen::SparseMatrix<double> SF2DBar = Basis.transpose() * SF2D * Basis; 
	Eigen::SparseMatrix<double> B2DBar = Basis.transpose() * B2D * Basis;
	//A_LHSBar = SF2DBar + lambda*CBar.transpose()*CBar; 
	A_LHSBar = lambda(0)*SF2DBar +  lambda(1)*B2DBar + lambda(2)*CBar.transpose()*CBar;

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;
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

void VectorFields::solveUserSystemMappedLDLTSoftConstraints(Eigen::SparseMatrix<double>& A_LHSBar, Eigen::VectorXd& bBar)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Solving reduced system...\n";

	//XLowDim.resize(B2DBar.rows());
	//XFullDim.resize(Basis.rows());
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

	XLowDim = x; 

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "..in Total of " << duration.count() << " seconds." << endl;
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

void VectorFields::computeEigenFields(const int &numEigs, const string& filename)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Computing reference eigenproblem (in Matlab)... ";

	//computeEigenMatlab(SF2DAsym, MF2D, eigFieldFull2D, eigValuesFull);
	computeEigenMatlab(SF2DAsym, MF2D, numEigs, eigFieldFull2D, eigValuesFull, filename);
	//computeEigenMatlab(SF2D, MF2D, numEigs, eigFieldFull2D, eigValuesFull, "hello");
	//cout << "::::: Eigen Values (Full Res) \n" << eigValuesFull << endl;
	//WriteSparseMatrixToMatlab(MF2D, "hello");

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}

void VectorFields::retrieveEigenFields(const string& filename)
{
	ReadDenseMatrixFromMatlab(eigFieldFull2D, filename);
	ReadVectorFromMatlab(eigValuesFull, filename);
}

void VectorFields::computeApproxEigenFields(const int &numEigs, const string& filename)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Computing restricted eigenproblem (in Matlab)...\n ";

	cout << "____Computing reduced system... ";
	Eigen::SparseMatrix<double> Mbar = Basis.transpose() * MF2D * Basis;
	Eigen::SparseMatrix<double> SFAsymbar = Basis.transpose() * SF2DAsym * Basis;




	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;


	t1 = chrono::high_resolution_clock::now();
	//computeEigenGPU(Sbar, Mbar, eigFieldReduced2D, eigValuesReduced);
	//computeEigenMatlab(Sbar, Mbar, eigFieldReduced2D, eigValuesReduced);
	computeEigenMatlab(SFAsymbar, Mbar, numEigs, eigFieldReduced2D, eigValuesReduced, filename);
	//cout << "::::: Eigen Values (Reduced) \n" << eigValuesReduced << endl;

	//WriteSparseMatrixToMatlab(Basis, "hello");

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "____Computing the eigenproblem in " << duration.count() << " seconds" << endl;
}

void VectorFields::retrieveApproxEigenFields() 
{
	//ReadSparseMatrixFromMatlab(Basis, "Hello");
	ReadDenseMatrixFromMatlab(eigFieldReduced2D, "Hello");
	ReadVectorFromMatlab(eigValuesReduced, "hello");
}

void VectorFields::measureApproxAccuracyL2Norm()
{
	/* Normalizing the fields on each patch */
	//Eigen::Vector2d xfRef2, xfApp2;
	//for (int i = 0; i < F.rows(); i++)
	//{
	//	xfRef2 = Xf.block(2 * i, 0, 2, 1); xfRef2.normalize();
	//	Xf.block(2 * i, 0, 2, 1) = xfRef2;
	//	xfApp2 = XFullDim.block(2 * i, 0, 2, 1); xfApp2.normalize();
	//	XFullDim.block(2 * i, 0, 2, 1) = xfApp2;
	//}

	Eigen::VectorXd diffV = Xf - XFullDim;
	double xf = Xf.transpose() * MF2D * Xf;
	double diff = diffV.transpose() * MF2D * diffV;
	const double L2norm = sqrt(diff / xf);

	printf("L2norm 0 = %.10f (%.10f / %.10f) \n", L2norm, diff, xf); 	
	printf("Max Error = %.3f \n", diffV.maxCoeff() / xf);

	/* Computing the energy */
	double refHarmEnergy = Xf.transpose()       * SF2D * Xf;
	double appHarmEnergy = XFullDim.transpose() * SF2D * XFullDim;
	double refBiHarmEnergy = Xf.transpose()       * B2D * Xf;
	double appBiHarmEnergy = XFullDim.transpose() * B2D * XFullDim;

	cout << ">> [REF] Energy: Harm=" << refHarmEnergy << ", Biharmonic=" << refBiHarmEnergy << endl;
	cout << ">> [APP] Energy: Harm=" << appHarmEnergy << ", Biharmonic=" << appBiHarmEnergy << endl;
	cout << "         Relative Harm-energy =" << abs(refHarmEnergy - appHarmEnergy) / refHarmEnergy << endl;
	cout << "         Relative Biharm-energy =" << abs(refBiHarmEnergy - appBiHarmEnergy) / refBiHarmEnergy << endl;


	//printf("MSE = %.3f \n", (diffV.sum()/diffV.size()) / xf);
}

void VectorFields::measureDirichletEnergy()
{
	double dirichlet = Xf.transpose() * ((B2D * MF2D) * Xf);
	cout << "__Dirichlet Energy\n \t__FullRes: " << dirichlet; 
	dirichlet = XFullDim.transpose() * ((B2D * MF2D) * XFullDim); 
	cout << ": Reduced: " << dirichlet << endl; 

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

void VectorFields::measureL2NormEigVectors()
{
	/* Getting overview of first few values*/
	printf("FULL \t REDUCED\n");
	for (int i = 0; i < 20; i++)
	{
		printf("%.4f \t %.4f \n", eigFieldFull2D(0,i), eigFieldReduced2D(0,i));
	}

	/* Computing the norm*/
	for (int i = 0; i < 100; i++)
	{
		Eigen::VectorXd diff = eigFieldFull2D.col(i) - eigFieldReduced2D.col(i);
		double norm1 = diff.transpose()*MF2D*diff; 
		double norm2 = eigFieldFull2D.col(i).transpose()*MF2D*eigFieldFull2D.col(i); 
		double l2norm = sqrt(norm1 / norm2);
		printf("Norm of %d diff is %.10f \n", i, l2norm);
	}
}

void VectorFields::vectorFieldsDesignTest()
{
	const int NUM_TESTS = 1;
	Eigen::VectorXd l2norm_error(NUM_TESTS);

	for (int i = 0; i < NUM_TESTS; i++)
	{	
		/* Setting up the constraints */
		constructConstraints();

		/* Solving the full resolution (incl the constraints) */
		Eigen::Vector3d lambda;
		lambda(0) = 1.0; 	// on harmonic energy
		lambda(1) = 0.1; 	// on bi-harmonic energy
		lambda(2) = 0.4;
		setupGlobalProblem(lambda);

		/* Solving the reduced system*/
		setupReducedBiLaplacian();
		setAndSolveUserSystem(lambda);			

		/* Computing L2-norm error*/
		Eigen::VectorXd diffV = Xf - XFullDim;
		double xf = Xf.transpose() * MF2D * Xf;
		double diff = diffV.transpose() * MF2D * diffV;
		const double L2norm = sqrt(diff / xf);
		l2norm_error(i) = L2norm;
		cout << "Error " << i << " = " << L2norm << endl;

		/* Computing the energy */
		double refHarmEnergy   = Xf.transpose()       * SF2D * Xf; 
		double appHarmEnergy   = XFullDim.transpose() * SF2D * XFullDim;
		double refBiHarmEnergy = Xf.transpose()       * B2D * Xf;
		double appBiHarmEnergy = XFullDim.transpose() * B2D * XFullDim;

		cout << ">> [REF] Energy: Harm=" << refHarmEnergy << ", Biharmonic=" << refBiHarmEnergy << endl;		
		cout << ">> [APP] Energy: Harm=" << appHarmEnergy << ", Biharmonic=" << appBiHarmEnergy << endl;
		cout << "         Relative Harm-energy =" << abs(refHarmEnergy - appHarmEnergy) / refHarmEnergy << endl;
		cout << "         Relative Biharm-energy =" << abs(refBiHarmEnergy - appBiHarmEnergy) / refBiHarmEnergy << endl;
	}
	cout << "ERRORS: \n" << l2norm_error << endl; 
}

void VectorFields::vectorFieldsDesignTest_Normalized()
{
	const int NUM_TESTS = 1;
	Eigen::VectorXd l2norm_error(NUM_TESTS);

	for (int i = 0; i < NUM_TESTS; i++)
	{
		/* Setting up the constraints */
		constructConstraints();

		/* Solving the full resolution (incl the constraints) */
		Eigen::Vector3d lambda;
		lambda(0) = 1.0; 	// on harmonic energy
		lambda(1) = 0.1; 	// on bi-harmonic energy
		lambda(2) = 0.4;
		setupGlobalProblem(lambda);

		/* Solving the reduced system*/
		setupReducedBiLaplacian();
		setAndSolveUserSystem(lambda);

		/* 'Normalize' each component of the vector fields */
		Eigen::VectorXd VFieldsRef(2*F.rows()), VFieldsApp(2*F.rows());
		Eigen::Vector2d vf2Ref, vf2App;

		for (int i = 0; i < F.rows(); i++)
		{
			/* Normalize the full res vector fields */
			vf2Ref = Xf.block(2 * i, 0, 2, 1);
			VFieldsRef.block(2 * i, 0, 2, 1) = vf2Ref.normalized();

			/* Normalize the approximated vectorfields */
			vf2App = XFullDim.block(2 * i, 0, 2, 1);
			VFieldsApp.block(2 * i, 0, 2, 1) = vf2App.normalized();
		}

		/* old values = new values */
		Xf = VFieldsRef;
		XFullDim = VFieldsApp;

		/* Computing L2-norm error*/
		Eigen::VectorXd diffV = Xf - XFullDim;
		double xf = Xf.transpose() * /*MF2D **/ Xf;
		double diff = diffV.transpose() * /*MF2D **/ diffV;
		const double L2norm = sqrt(diff / xf);
		l2norm_error(i) = L2norm;
		cout << "Error " << i << " = " << L2norm << endl;
	}
	cout << "ERRORS: \n" << l2norm_error << endl;
}

/* ====================== APPLICATIONS ON REDUCED SYSTEM ============================*/
void VectorFields::computeSmoothingApprox(const double& mu, const Eigen::VectorXd& v_in, Eigen::VectorXd& v_out)
{
	/* Reduced Matrices */
	Eigen::SparseMatrix<double> MF2DBar = (Basis.transpose()*MF2D)*Basis; 
	Eigen::SparseMatrix<double> SF2DBar = (Basis.transpose()*SF2D)*Basis;
	Eigen::VectorXd v_inBar = Basis.transpose()*v_in;
	Eigen::VectorXd v_outBar;
	cout << "The reduced matrices are set up\n"; 
	/* First flavour */
	Eigen::SparseMatrix<double> AL = MF2DBar + mu*SF2DBar ;

	/* Second flavour */
	//Eigen::SparseMatrix<double> A = MF2DBar + mu*SF2DBar*(Basis.transpose()*MF2Dinv*Basis)*SF2DBar;
	Eigen::VectorXd b = MF2DBar*v_inBar;

	Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> sparseSolver(AL);
	v_outBar = sparseSolver.solve(b);
	v_out = Basis*v_outBar;
	cout << "VOUT \n" << v_out.block(0, 0, 100, 1) << endl; 

	/* Computing the L2-norm of the smoothed fields w.r.t input*/
	double diff1 = (v_out - v_in).transpose()*MF2D*(v_out - v_in);
	double diff2 = v_in.transpose()*MF2D*v_in;
	double sqrt_norm = sqrt(diff1 / diff2);
	printf("The diff of v_out and v_in is %.10f \n", sqrt_norm);

	/* Computing the energy */
	double energy1 = v_in.transpose() * ((B2D * MF2D) * v_in);
	double energy2 = v_out.transpose() * ((B2D * MF2D) * v_out);
	printf("The energy is=%.4f ==> %.4f.\n", energy1, energy2);
}

void VectorFields::ConstructCurvatureTensor(igl::opengl::glfw::Viewer &viewer)
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
	CurvatureTensor2D.resize(2 * F.rows(), 2 * F.rows());
	CurvatureTensor2D.reserve(2 * 2 * F.rows());			// 2*F rows, each with 2 non-zero entries
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
		mT = ( m1 + m2 + m3) / (2.0*doubleArea(i)*doubleArea(i));

		/* Inserting the 2x2 matrix*/
		Eigen::MatrixXd ALoc = A.block(3 * i, 2 * i, 3, 2);
		mT2D = ALoc.transpose() * mT * ALoc;
		CTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, mT2D(0, 0)));
		CTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 0, mT2D(1, 0)));
		CTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 1, mT2D(0, 1)));
		CTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, mT2D(1, 1)));

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
	CurvatureTensor2D.setFromTriplets(CTriplet.begin(), CTriplet.end());
}

// [OLD] Wrong implementation
//void VectorFields::ConstructCurvatureTensor()
//{
//	/* Obtain the principal curvature using LibIGL (vertex-based) */
//	Eigen::MatrixXd PD1, PD2;
//	Eigen::VectorXd PV1, PV2;
//	igl::principal_curvature(V, F, PD1, PD2, PV1, PV2);
//
//	/* Test on curvature and mean curvature*/
//	double meanCurve1 = 0.5*(PV1(0) + PV2(0));
//	double curveDir1 = PD1.row(0).dot(PD2.row(0));
//	printf("Mean=%.4f | dot=%.4f \n", meanCurve1, curveDir1);
//	meanCurve1 = 0.5*(PV1(1) + PV2(1));
//	curveDir1 = PD1.row(1).dot(PD2.row(1));
//	printf("Mean=%.4f | dot=%.4f \n", meanCurve1, curveDir1);
//
//	/* Covert the vertex-based to face-based principal curvatures */
//	Eigen::MatrixXd CurvatureTensor3D;
//	CurvatureTensor3D.setZero(3 * F.rows(), 2);
//	for (int i = 0; i < F.rows(); i++)
//	{
//		for (int j = 0; j < F.cols(); j++)
//		{
//			/* Maximum curvature direction */
//			CurvatureTensor3D.block(3 * i, 0, 3, 1) += (PD1.row(F(i, j))).transpose() / 3.0;
//			/* Minimum curvature direction */
//			CurvatureTensor3D.block(3 * i, 1, 3, 1) += (PD2.row(F(i, j))).transpose() / 3.0;
//		}
//		//CurvatureTensor3D.row(i) /= double(F.cols());		
//	}
//
//	CurvatureTensor = A.transpose() * CurvatureTensor3D;
//
//}


/* Temporary functions*/


void sortEigenIndex(double eig1, double eig2, double eig3, int& smallest, int& middle, int& largest)
{
	if (eig1 > eig2)
	{
		if (eig1 > eig3)
		{
			largest = 0;
			if (eig2 > eig3)
			{
				middle = 1;
				smallest = 2; 
			} 
			else
			{
				middle = 2; 
				smallest = 1; 
			}
		}
		else
		{
			largest = 2; 
			middle = 0;
			smallest = 1; 
		}
	}
	else if (eig2 > eig3)
	{
		largest = 1;
		if (eig1 > eig3)
		{
			middle = 0;
			smallest = 2; 
		} 
		else
		{
			middle = 2;
			smallest = 0;
		}
	} 
	else
	{
		largest = 2; 
		middle = 1;
		smallest = 0;
	}

	printf("Eig1=%.4f, Eig2=%.4f, Eig3=%.4f  | smallest: %d, middle: %d, largest: %d\n", eig1, eig2, eig3, smallest, middle, largest);
}

void VectorFields::ComputeCurvatureFields()
{
	/* Resizing the principal curvatures */
	CurvatureTensorField2D.resize(2 * F.rows(), 2);



	
	/* Computing the eigenvectors => principal curvatures */

	/* Making the construction parallel*/
	int id, tid, ntids, ipts, istart, iproc;
#pragma omp parallel private(tid,ntids,ipts,istart,id)	
	{
		iproc = omp_get_num_procs();
		//iproc = 1; 
		tid = omp_get_thread_num();
		ntids = omp_get_num_threads();
		//ipts = (int)ceil(1.00*(double)F.rows() / (double)ntids);
		ipts = F.rows() / ntids;
		istart = tid * ipts;
		if (tid == ntids - 1) ipts = F.rows() - istart;
		if (ipts <= 0) ipts = 0;

		//cout << "[" << tid << "] Number of processors " << iproc << ", with " << ntids << " threads." << endl;
		// Computing the values of each element
		for (id = istart; id < (istart + ipts); id++) {
			printf("Row of %d \n", id);
			if (id >= F.rows()) "This is an ERROR!!!\n";

			/* Local variables */
			Eigen::MatrixXd MLoc2D(2, 2), TensorLoc2D(2, 2), EigFields2D;
			Eigen::SparseMatrix<double> MLoc(2, 2), TensorLoc(2, 2);
			Eigen::VectorXd eigVals2D;

			/* Dummy mass matrix*/
			MLoc2D << 1.0, 0.0, 0.0, 1.0;
			MLoc.coeffRef(0, 0) = 1.0;
			MLoc.coeffRef(1, 1) = 1.0;
			int smallest, middle, largest;


			TensorLoc2D = CurvatureTensor2D.block(2 * id, 2 * id, 2, 2);
			vector<Eigen::Triplet<double>> TTriplet;
			TTriplet.reserve(4);
			TTriplet.push_back(Eigen::Triplet<double>(0, 0, CurvatureTensor2D.coeff(2 * id + 0, 2 * id + 0)));
			TTriplet.push_back(Eigen::Triplet<double>(1, 0, CurvatureTensor2D.coeff(2 * id + 1, 2 * id + 0)));
			TTriplet.push_back(Eigen::Triplet<double>(0, 1, CurvatureTensor2D.coeff(2 * id + 0, 2 * id + 1)));
			TTriplet.push_back(Eigen::Triplet<double>(1, 1, CurvatureTensor2D.coeff(2 * id + 1, 2 * id + 1)));
			TensorLoc.setFromTriplets(TTriplet.begin(), TTriplet.end());

			//if (i == 0) cout << "HEYYYYY\n" << TensorLoc << endl << endl;
			//if (i == 0) cout << "HEYYYYY\n" << TensorLoc2D << endl << endl;

			TensorLoc2D = CurvatureTensor2D.block(2 * id, 2 * id, 2, 2);
			//computeEigenGPU(TensorLoc2D, MLoc2D, EigFields2D, eigVals2D);
			//computeEigenMatlab(TensorLoc, MLoc, EigFields2D, eigVals2D);
			computeEigenMatlab(TensorLoc, MLoc, (int)MLoc.rows(), EigFields2D, eigVals2D, "hello there");
			//sortEigenIndex(abs(eigVals(0)), abs(eigVals(1)), abs(eigVals(2)), smallest, middle, largest);
			if ((eigVals2D(0)) > (eigVals2D(1))) { largest = 0; smallest = 1; }
			else { largest = 1; smallest = 0; }
			printf("__[%d] eigVal1=%.4f, eigVec[%.4f;%.4f]  \t eigVal2=%.4f, eigVec[%.4f;%.4f]\n",
					id, eigVals2D(smallest), EigFields2D(0, smallest), EigFields2D(1, smallest),
					    eigVals2D(largest),  EigFields2D(0, largest),  EigFields2D(1, largest));

			CurvatureTensorField2D.block(2 * id, 0, 2, 1) = EigFields2D.col(largest);
			CurvatureTensorField2D.block(2 * id, 1, 2, 1) = EigFields2D.col(smallest);
		}
	}
}

//
//<<<<<<< HEAD
//		TensorLoc2D = CurvatureTensor2D.block(2 * i, 2 * i, 2, 2);
//		vector<Eigen::Triplet<double>> TTriplet;
//		TTriplet.reserve(4);
//		TTriplet.push_back(Eigen::Triplet<double>(0, 0, CurvatureTensor2D.coeff(2 * i + 0, 2 * i + 0)));
//		TTriplet.push_back(Eigen::Triplet<double>(1, 0, CurvatureTensor2D.coeff(2 * i + 1, 2 * i + 0)));
//		TTriplet.push_back(Eigen::Triplet<double>(0, 1, CurvatureTensor2D.coeff(2 * i + 0, 2 * i + 1)));
//		TTriplet.push_back(Eigen::Triplet<double>(1, 1, CurvatureTensor2D.coeff(2 * i + 1, 2 * i + 1)));
//		TensorLoc.setFromTriplets(TTriplet.begin(), TTriplet.end());
//		
//		if (i == 0) cout << "HEYYYYY\n" << TensorLoc << endl << endl; 
//		if (i == 0) cout << "HEYYYYY\n" << TensorLoc2D << endl << endl;
//
//		//computeEigenGPU(TensorLoc2D, MLoc2D, EigFields2D, eigVals2D);
//		//computeEigenMatlab(TensorLoc, MLoc, EigFields2D, eigVals2D);
//		computeEigenMatlab(TensorLoc, MLoc, (int) MLoc.rows(), EigFields2D, eigVals2D, "hello there");
//		//sortEigenIndex(abs(eigVals(0)), abs(eigVals(1)), abs(eigVals(2)), smallest, middle, largest);
//		if ((eigVals2D(0)) > (eigVals2D(1))) { largest = 0; smallest = 1; }
//		else { largest = 1; smallest = 0; }
//		printf("__[%d] eigVal1=%.4f, eigVec[%.4f;%.4f]  \t eigVal2=%.4f, eigVec[%.4f;%.4f]\n",
//				i, eigVals2D(0), EigFields2D(0, 0), EigFields2D(1, 0),
//				   eigVals2D(1), EigFields2D(0, 1), EigFields2D(1, 1));
//
//		CurvatureTensorField2D.block(2 * i, 0, 2, 1) = EigFields2D.col(largest);
//		CurvatureTensorField2D.block(2 * i, 1, 2, 1) = EigFields2D.col(smallest);
//=======		

	
	//for (int i = 0; i < F.rows(); i++)
	//{
	//	TensorLoc2D = CurvatureTensor2D.block(2 * i, 2 * i, 2, 2);
	//	computeEigenGPU(TensorLoc2D, MLoc2D, EigFields2D, eigVals2D);
	//	//sortEigenIndex(abs(eigVals(0)), abs(eigVals(1)), abs(eigVals(2)), smallest, middle, largest);
	//	if ((eigVals2D(0)) > (eigVals2D(1))) { largest = 0; smallest = 1; }
	//	else { largest = 1; smallest = 0; }
	//	//printf("__[%d] eigVal1=%.4f, eigVec[%.4f;%.4f]  \t eigVal2=%.4f, eigVec[%.4f;%.4f]\n",
	//	//		i, eigVals2D(0), EigFields2D(0, 0), EigFields2D(1, 0),
	//	//		   eigVals2D(1), EigFields2D(0, 1), EigFields2D(1, 1));
	//
	//	CurvatureTensorField2D.block(2 * i, 0, 2, 1) = EigFields2D.col(largest);
	//	CurvatureTensorField2D.block(2 * i, 1, 2, 1) = EigFields2D.col(smallest);
	//}


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

	igl::edges(F, E);
	const int genus = (2 - V.rows() + E.rows() - F.rows()) / 2;
	printf("This model is of genus %d\n", genus);

}

void VectorFields::scaleMesh()
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
	V = V * (1.0/scaleFactor);


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

void VectorFields::readArrowMesh(const string &meshFile)
{
	// For actual work of reading mesh object
	VArrow.resize(0, 0);
	FArrow.resize(0, 0);

	if (meshFile.substr(meshFile.find_last_of(".") + 1) == "off") {
		igl::readOFF(meshFile, VArrow, FArrow);
	}
	else if (meshFile.substr(meshFile.find_last_of(".") + 1) == "obj") {
		igl::readOBJ(meshFile, VArrow, FArrow);
	}
	else {
		cout << "Error! File type can be either .OFF or .OBJ only." << endl;
		cout << "Program will exit in 2 seconds." << endl;
		Sleep(2000);
		exit(10);
	}

	printf("....V=%dx%d\n", VArrow.rows(), VArrow.cols());
	printf("....F=%dx%d\n", FArrow.rows(), FArrow.cols());
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

