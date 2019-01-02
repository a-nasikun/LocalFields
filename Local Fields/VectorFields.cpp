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
	//constructSpecifiedConstraints();
	
	constructSingularities();
	constructSpecifiedConstraintsWithSingularities();

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
	const int numConstraints = 2;
	set<int> constraints;
	//vector<int> globalConstraints(numConstraints);
	globalConstraints.resize(numConstraints);
	Eigen::VectorXd D;
	D.resize(F.rows());

	// Initialize the value of D
	for (int i = 0; i < F.rows(); i++) {
		D(i) = numeric_limits<double>::infinity();
	}

	constraints.insert(0);
	int curPoint = 0;

	do {
		Eigen::VectorXi::Index maxIndex;
		computeDijkstraDistanceFaceForSampling(curPoint, D);
		D.maxCoeff(&maxIndex);
		constraints.insert(maxIndex);
		curPoint = maxIndex;
	} while (constraints.size() <= numConstraints);

	int counter1 = 0;
	for (int i : constraints) {
		globalConstraints[counter1++] = i;
	}
	printf("Constraints = %d\n", globalConstraints.size());

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
	printf("Cp=%dx%d\n", C.rows(), C.cols());


	// Setting up vector c (There are 2 vector c)
	srand(time(NULL));
	c.resize(2 * globalConstraints.size(), 2);
	for (int i = 0; i < globalConstraints.size(); i++) {
		c(2 * i + 0, 0) = sqrt(2.0);
		c(2 * i + 1, 0) = sqrt(2.0);
		c(2 * i + 0, 1) = sqrt(2.0);
		c(2 * i + 1, 1) = sqrt(2.0);
	}
	printf("cBar=%dx%d\n", c.rows(), c.cols());
}

void VectorFields::constructSingularities()
{
	const int NUM_SINGS = 2;
	singularities.resize(NUM_SINGS);
	SingNeighCC.resize(NUM_SINGS);

	time_t t;
	srand((unsigned)time(&t));
	//srand(time(NULL));
	for (int id = 0; id < NUM_SINGS; id++) {
		const int SingLocation = rand() % V.rows();
		//const int SingLocation = id * (int)(V.rows() / NUM_SINGS) + 5;
		singularities[id] = SingLocation;
		const int SingNeighNum = VFNeighFull[SingLocation].size();
		const int firstNeigh = VFNeighFull[SingLocation].begin()->fId;

		SingNeighCC[id].resize(SingNeighNum);
		SingNeighCC[id][0] = firstNeigh;
		int curNeigh = firstNeigh;
		int vertex1 = SingLocation;

		//printf("Vertex %d has %d neighbors.\n", SingLocation, SingNeighNum);
		
		for (int i2 = 1; i2<SingNeighNum; i2++) {
			//printf("number of neighbors of vertex %d = %d (%d)\n", SingLocation, i2 /*SingNeighCC[id].size()*/, VFNeighFull[SingLocation].size());
			int vertex2;
			for (int i = 0; i < F.cols(); i++) {
				if (F(curNeigh, i%F.cols()) == vertex1) {
					vertex2 = F(curNeigh, (i + F.cols() - 1) % F.cols());
				}
			}
			for (std::set<VtoFPair>::iterator it1 = next(VFNeighFull[SingLocation].begin(), 1); it1 != VFNeighFull[SingLocation].end(); ++it1) {
				for (int i = 0; i < F.cols(); i++) {
					if (F(it1->fId, i) == vertex1 && F(it1->fId, (i + 1) % F.cols()) == vertex2) {
						SingNeighCC[id][i2] = it1->fId;
						//printf("     Inserting %d as neighbor\n", it1->fId);
						curNeigh = it1->fId;
					}
				}
			}
		}
	}
}

void VectorFields::constructSpecifiedConstraintsWithSingularities()
{
	// Define the constraints
	const int numConstraints = 30;
	set<int> constraints;

	globalConstraints.resize(numConstraints);
	Eigen::VectorXd D;
	D.resize(F.rows());

	// Initialize the value of D
	for (int i = 0; i < F.rows(); i++) {
		D(i) = numeric_limits<double>::infinity();
	}

	constraints.insert(0);
	int curPoint = 0;

	// Creating random constraints
	do {
		Eigen::VectorXi::Index maxIndex;
		computeDijkstraDistanceFaceForSampling(curPoint, D);
		D.maxCoeff(&maxIndex);
		constraints.insert(maxIndex);
		curPoint = maxIndex;
	} while (constraints.size() <= numConstraints);

	int counter1 = 0;
	for (int i : constraints) {
		globalConstraints[counter1++] = i;
	}
	//printf("Constraints = %d\n", globalConstraints.size());

	//constructSingularities();

	int numSingConstraints = 0;
	for (int i = 0; i < SingNeighCC.size(); i++) {
		// Use only n-1 neighboring faces as constraints
		for (int j = 0; j < (SingNeighCC[i].size()-1); j++) {
			numSingConstraints++;
		}
	}

	// Setting up matrix C and vector c
	c.resize(2 * (globalConstraints.size()+numSingConstraints), 2);
	//printf("cBar=%dx%d\n", c.rows(), c.cols());


	// HARD CONSTRAINTS
	Eigen::SparseMatrix<double> CTemp;
	vector<Eigen::Triplet<double>> CTriplet;
	CTriplet.reserve(2 * globalConstraints.size() + 2 * 3 * 7 * SingNeighCC.size());
	int counter = 0;
	for (int i = 0; i < globalConstraints.size(); i++) {
		// Matrix C
		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * globalConstraints[i] + 0, 1.0));
		c(counter, 0) = sqrt(2.0);
		c(counter++, 1) = sqrt(2.0);
		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * globalConstraints[i] + 1, 1.0));
		c(counter, 0) = sqrt(2.0);
		c(counter++, 1) = sqrt(2.0);
	}

	
	// SINGULARITIES CONSTRAINTS
	for (int id = 0; id < SingNeighCC.size(); id++) {
		const double rotAngle = 2 * M_PI / (double)SingNeighCC[id].size();
		const double cosA = cos(rotAngle);
		const double sinA = sin(rotAngle);
		for (int i = 0; i < (SingNeighCC[id].size() - 1); i++) {
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, cosA));
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, -sinA));
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, -1.0));
			c(counter, 0) = 0.0;
			c(counter, 1) = 0.0;
			counter++;

			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, sinA));
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, cosA));
			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, -1.0));
			c(counter, 0) = 0.0;
			c(counter, 1) = 0.0;
			counter++;
		}
	}

	C.resize(2 * (globalConstraints.size()+numSingConstraints), B2D.rows());
	C.setFromTriplets(CTriplet.begin(), CTriplet.end());
	//printf("Cp=%dx%d\n", C.rows(), C.cols());	
}

void VectorFields::setupGlobalProblem()
{	
	constructConstraints();	
	setupRHSGlobalProblemMapped();
	setupLHSGlobalProblemMapped();
	//setupRHSGlobalProblem();
	//setupLHSGlobalProblem();	

	//solveGlobalSystem();
	solveGlobalSystemMappedLDLT();
}

void VectorFields::setupRHSGlobalProblem()
{
	Eigen::VectorXd		zeroElements(B2D.rows());
	for (int i = 0; i < B2D.rows(); i++) zeroElements(i) = 0.0;

	b.resize(B2D.rows() + c.rows(),2);
	b.col(0) << zeroElements, c.col(0); 
	b.col(1) << zeroElements, c.col(1);

	//cout << "b: " << endl << b << endl << endl; 
}

void VectorFields::setupRHSGlobalProblemMapped()
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
	b.resize(B2D.rows() + c.rows(), 2);

	// First column of b
	h = C * vEst - c.col(0);
	b.col(0) << g, h;

	// Second Column of b
	h = C * vEst - c.col(1);
	b.col(1) << g, h;

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}

void VectorFields::setupLHSGlobalProblem()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Constructing LHS... ";

	A_LHS.resize(B2D.rows() + C.rows(), B2D.cols() + C.rows());

	vector<Eigen::Triplet<double>>	ATriplet;
	ATriplet.reserve(10 * B2D.rows());

	for (int k = 0; k < B2D.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(B2D, k); it; ++it) {
			ATriplet.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
		}
	}

	for (int k = 0; k < C.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(C, k); it; ++it) {
			ATriplet.push_back(Eigen::Triplet<double>(B2D.rows() + it.row(), it.col(), it.value()));
			ATriplet.push_back(Eigen::Triplet<double>(it.col(), B2D.cols() + it.row(), -it.value()));
		}
	}
	A_LHS.setFromTriplets(ATriplet.begin(), ATriplet.end());

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}

void VectorFields::setupLHSGlobalProblemMapped()
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

void VectorFields::solveGlobalSystem()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Solving the global system... \n";

	Xf.resize(B2D.rows(), 2);

	// Setting up the solver
	Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
	Eigen::VectorXd x; 
	//solver.compute(A);
	
	solver.analyzePattern(A_LHS);

	cout << "....Factorizing LHS..." << endl;
	solver.factorize(A_LHS);
	if (solver.info() != Eigen::Success) {
		cout << "Factorization Failed." << endl;
		return;
	}

	//printf("Size of B=%d\n", b.rows());
	// FIRST BASIS
	cout << "....Solving first problem (first frame)." << endl;
	x = solver.solve(b.col(0));
	if (solver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system." << endl;
		return;
	}

	Xf.col(0) = x.block(0, 0, B2D.rows(), 1);	

	// SECOND BASIS
	cout << "....Solving second problem. (second frame)" << endl;
	x = solver.solve(b.col(1));
	if (solver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system." << endl;
		return;
	}
	Xf.col(1) = x.block(0, 0, B2D.rows(), 1);
	
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "..in " << duration.count() << " seconds" << endl;
}

void VectorFields::solveGlobalSystemMappedLDLT()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Solving the global system (Pardiso LDLT)... \n";

	//cout << "Starting to solve problem." << endl;
	Xf.resize(B2D.rows(), 2);
	
	// Setting up the solver
	//Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> sparseSolver(A_LHS);
	Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> sparseSolver(A_LHS);	
	//Eigen::PastixLDLT<Eigen::SparseMatrix<double>,1> sparseSolver(A_LHS);

	// FIRST BASIS
	cout << "....Solving first problem (first frame)..." << endl;
	Eigen::VectorXd x = sparseSolver.solve(b.col(0));
	
	if (sparseSolver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system. " << endl;
		if (sparseSolver.info() == Eigen::NumericalIssue)
			cout << "NUMERICAL ISSUE. " << endl;
		cout << sparseSolver.info() << endl;
		return;
	}
	
	Xf.col(0) = -x.block(0, 0, B2D.rows(), 1) + vEst;
	
	// SECOND BASIS
	cout << "....Solving second problem (second frame)." << endl;
	x = sparseSolver.solve(b.col(1));
	if (sparseSolver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system." << endl;
		return;
	}
	
	Xf.col(1) = -x.block(0, 0, B2D.rows(), 1) + vEst;
	
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

	double	coef = sqrt(pow(1.1, 2) + pow(0.7, 2));
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
	cout << "Average non-zeros is " << (double)numNonZeroes / (double)numElements << endl; 
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
	setupUserBasis();
	getUserConstraints();
	setupRHSUserProblemMapped();
	setupLHSUserProblemMapped();
	solveUserSystemMappedLDLT();
	mapSolutionToFullRes();
	//obtainUserVectorFields();
}

void VectorFields::setupUserBasis()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Constructing Local Basis...";

	B2Dbar = Basis.transpose() * B2D * Basis; 

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;
	printf(".... Local Basis = %dx%d\n", B2Dbar.rows(), B2Dbar.cols());
}
void VectorFields::getUserConstraints()
{	
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Obtaining user constraints ";

	//getUserConstraintsRandom();

	//const vector<vector<int>> selectedFaces;
	//getUserConstraintsGUI(selectedFaces);
	//getUserConstraintsSpecified();
	userConstraints = globalConstraints; 
	Cbar = C * Basis;
	cBar = c;

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;

	printf(".... C_Lobal = %dx%d\n", Cbar.rows(), Cbar.cols());
	printf(".... c_Lobal = %dx%d\n", cBar.rows(), cBar.cols());
}

void VectorFields::getUserConstraintsRandom()
{
	// Define the constraints
	const int numConstraints = 20;

	srand(time(NULL));
	set<int> constraints;
	userConstraints.resize(numConstraints);
	do {
		int c = rand() % F.rows();
		constraints.insert(c);
	} while (constraints.size() <= numConstraints);

	int counter1 = 0;
	for (int i : constraints) {
		userConstraints[counter1++] = i;
	}
	printf("Constraints = %d\n", userConstraints.size());

	// Setting up matrix C
	Eigen::SparseMatrix<double> CTemp;
	vector<Eigen::Triplet<double>> CTriplet;
	CTriplet.reserve(2 * userConstraints.size());

	int counter = 0;
	for (int i = 0; i < userConstraints.size(); i++) {
		CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * userConstraints[i] + 0, 1.0));
		CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * userConstraints[i] + 1, 1.0));
	}
	CTemp.resize(2 * userConstraints.size(), Basis.rows());
	CTemp.setFromTriplets(CTriplet.begin(), CTriplet.end());
	printf("CTemp=%dx%d\n", CTemp.rows(), CTemp.cols());

	Cbar = CTemp * Basis;
	printf("Cbar=%dx%d\n", Cbar.rows(), Cbar.cols());
	// Setting up vector c (There are 2 vector c)
	srand(time(NULL));
	cBar.resize(2 * userConstraints.size(), 2);
	for (int i = 0; i < userConstraints.size(); i++) {
		cBar(2 * i + 0, 0) = ((double)(rand() % 1000) - 500) / 1000.0;
		cBar(2 * i + 1, 0) = ((double)(rand() % 1000) - 500) / 1000.0;

		cBar(2 * i + 0, 1) = ((double)(rand() % 1000) - 500) / 1000.0;
		cBar(2 * i + 1, 1) = ((double)(rand() % 1000) - 500) / 1000.0;
	}
	printf("cBar=%dx%d\n", cBar.rows(), cBar.cols());
}

void VectorFields::getUserConstraintsGUI(const vector<vector<int>> &selectedFaces)
{

}

void VectorFields::getUserConstraintsSpecified()
{
	// Define the constraints
	const int numConstraints = 20;
	set<int> constraints;
	userConstraints.resize(numConstraints);
	Eigen::VectorXd D;
	D.resize(F.rows());

	// Initialize the value of D
	for (int i = 0; i < F.rows(); i++) {
		D(i) = numeric_limits<double>::infinity();
	}

	constraints.insert(0);	
	int curPoint = 0;
		
	do {
		Eigen::VectorXi::Index maxIndex;
		computeDijkstraDistanceFaceForSampling(curPoint, D);
		D.maxCoeff(&maxIndex);
		constraints.insert(maxIndex);
		curPoint = maxIndex;
	} while (constraints.size() <= numConstraints);

	int counter1 = 0;
	for (int i : constraints) {
		userConstraints[counter1++] = i;
	}
	printf("Constraints = %d\n", userConstraints.size());

	// Setting up matrix C
	Eigen::SparseMatrix<double> CTemp;
	vector<Eigen::Triplet<double>> CTriplet;
	CTriplet.reserve(2 * userConstraints.size());
	int counter = 0;
	for (int i = 0; i < userConstraints.size(); i++) {
		CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * userConstraints[i] + 0, 1.0));
		CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * userConstraints[i] + 1, 1.0));
	}
	CTemp.resize(2 * userConstraints.size(), Basis.rows());
	CTemp.setFromTriplets(CTriplet.begin(), CTriplet.end());
	printf("CTemp=%dx%d\n", CTemp.rows(), CTemp.cols());

	Cbar = CTemp * Basis;
	printf("Cbar=%dx%d\n", Cbar.rows(), Cbar.cols());

	// Setting up vector c (There are 2 vector c)
	srand(time(NULL));
	cBar.resize(2 * userConstraints.size(), 2);
	for (int i = 0; i < userConstraints.size(); i++) {
		cBar(2 * i + 0, 0) = sqrt(2.0);
		cBar(2 * i + 1, 0) = sqrt(2.0);
		cBar(2 * i + 0, 1) = sqrt(2.0);
		cBar(2 * i + 1, 1) = sqrt(2.0);
	}
	printf("cBar=%dx%d\n", cBar.rows(), cBar.cols());
}

void VectorFields::setupRHSUserProblemMapped()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Constructing RHS (mapped)...";
	
	vEstUser.resize(B2Dbar.rows());
	for (int i = 0; i < vEstUser.rows(); i++) {
		vEstUser(i) = 0.5;
	}

	gbar = B2Dbar * vEstUser;
	//printf("gbar=%dx%d\n", gbar.rows(), gbar.cols());
	bBar.resize(B2Dbar.rows() + cBar.rows(), 2);
	//printf("bbar=%dx%d\n",bBar.rows(), bBar.cols());

	// First column of b
	hbar = Cbar * vEstUser - cBar.col(0);
	//printf("hbar=%dx%d\n", hbar.rows(), hbar.cols());

	bBar.col(0) << gbar, hbar;
	//printf("bbar=%dx%d\n", bBar.rows(), bBar.cols());
	// Second Column of b
	hbar = Cbar * vEstUser - cBar.col(1);
	//printf("hbar=%dx%d\n", hbar.rows(), hbar.cols());

	bBar.col(1) << gbar, hbar;


	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;
}

void VectorFields::setupLHSUserProblemMapped()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Constructing LHS (mapped)...";


	A_LHSbar.resize(B2Dbar.rows() + Cbar.rows(), B2Dbar.cols() + Cbar.rows());
	vector<Eigen::Triplet<double>>	ATriplet;
	ATriplet.reserve(B2Dbar.nonZeros() + 2 * Cbar.nonZeros());

	for (int k = 0; k < B2Dbar.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(B2Dbar, k); it; ++it) {
			ATriplet.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
		}
	}

	for (int k = 0; k < Cbar.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(Cbar, k); it; ++it) {
			ATriplet.push_back(Eigen::Triplet<double>(B2Dbar.rows() + it.row(), it.col(), it.value()));
			ATriplet.push_back(Eigen::Triplet<double>(it.col(), B2Dbar.cols() + it.row(), it.value()));
		}
	}
	A_LHSbar.setFromTriplets(ATriplet.begin(), ATriplet.end());

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;
	printf("....Local LHS = %dx%d\n", A_LHSbar.rows(), A_LHSbar.cols());
}

void VectorFields::solveUserSystemMappedLDLT()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Solving reduced system...\n";

	XLowDim.resize(B2Dbar.rows(), 2);
	XFullDim.resize(Basis.rows(), 2);
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> sparseSolver(A_LHSbar);

	cout << "....Solving for the first frame.\n";
	Eigen::VectorXd x = sparseSolver.solve(bBar.col(0));

	if (sparseSolver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system. " << endl;
		if (sparseSolver.info() == Eigen::NumericalIssue)
			cout << "NUMERICAL ISSUE. " << endl;
		cout << sparseSolver.info() << endl;
		return;
	}

	XLowDim.col(0) = -x.block(0, 0, B2Dbar.rows(), 1) + vEstUser;

	// SECOND BASIS
	cout << "....Solving for the second frame.\n";
	x = sparseSolver.solve(bBar.col(1));
	if (sparseSolver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system." << endl;
		return;
	}
	XLowDim.col(1) = -x.block(0, 0, B2Dbar.rows(), 1) + vEstUser;
	
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
	Eigen::VectorXd diff = Xf.col(0) - XFullDim.col(0);
	double xf = Xf.col(0).transpose() * MF2D * Xf.col(0);
	double L2norm = diff.transpose() * MF2D * diff; 
	printf("Diff 0 = %.10f\n", sqrt(L2norm / xf)); 

	diff = Xf.col(1) - XFullDim.col(1);
	xf = Xf.col(1).transpose() * MF2D * Xf.col(1);
	L2norm = diff.transpose() * MF2D * diff;
	printf("Diff 1 = %.10f\n", sqrt(L2norm / xf));
}


void VectorFields::constructLocalElements()
{
	LocalElements.resize(SubDomain.size() + Boundary.size());
	GlobToLocMap.resize(F.rows());
	for (int i = 0; i < F.rows(); i++) GlobToLocMap[i] = -1; 

	int counter = 0;
	for (int face : SubDomain) {
		LocalElements[counter] = face; 
		GlobToLocMap[face] = counter;
		counter++;
	}

	for (int face : Boundary) {
		LocalElements[counter] = face; 
		GlobToLocMap[face] = counter;
		counter++;
	}
}

void VectorFields::constructMatrixBLocal()
{
	BLoc.resize(2 * LocalElements.size(), 2 * LocalElements.size());
	vector<Eigen::Triplet<double>> BTriplet;
	BTriplet.reserve(10 * BLoc.rows());

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
}

void VectorFields::constructLocalConstraints()
{
	// Setting up matrix C
	CLoc.resize(2 * (1+LocalElements.size()), BLoc.cols());
	vector<Eigen::Triplet<double>> CTriplet;
	CTriplet.reserve(2 * (1 + LocalElements.size()));

	int counter = 0; 
	CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * GlobToLocMap[sample] + 0, 1.0));
	CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * GlobToLocMap[sample] + 1, 1.0));
	for(int bound : Boundary){
		CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * GlobToLocMap[bound] + 0, 1.0));
		CTriplet.push_back(Eigen::Triplet<double>(counter++, 2 * GlobToLocMap[bound] + 1, 1.0));
	}
	CLoc.setFromTriplets(CTriplet.begin(), CTriplet.end());
	cout << "Selection of faces is done." << endl;

	// Setting up vector c (There are 2 vector c)
	cLoc.resize(2 * (1 + LocalElements.size()), 2);
	Eigen::VectorXd zeroElements(2 * LocalElements.size());
	for (int i = 0; i < zeroElements.size(); i++) zeroElements(i) = 0.0;
	cLoc.col(0) << 1.0, 0.0, zeroElements;
	cLoc.col(1) << 0.0, 1.0, zeroElements;
	cout << "Definition of constraints is done." << endl;
}

void VectorFields::setupRHSLocalProblem()
{
	Eigen::VectorXd		zeroElements(BLoc.rows());
	for (int i = 0; i < zeroElements.size(); i++) zeroElements(i) = 0.0;

	bLoc.resize(BLoc.rows() + cLoc.rows(), 2);
	bLoc.col(0) << zeroElements, cLoc.col(0);
	bLoc.col(1) << zeroElements, cLoc.col(1);
	cout << "RHS is constructed." << endl;
}

void VectorFields::setupLHSLocalProblem()
{
	ALoc.resize(BLoc.rows() + CLoc.rows(), BLoc.cols() + CLoc.rows());	
	vector<Eigen::Triplet<double>>	ATriplet;
	ATriplet.reserve(BLoc.nonZeros() + 2 * CLoc.nonZeros());

	for (int k = 0; k < BLoc.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(BLoc, k); it; ++it) {
			ATriplet.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
		}
	}

	for (int k = 0; k < CLoc.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(CLoc, k); it; ++it) {
			ATriplet.push_back(Eigen::Triplet<double>(BLoc.rows() + it.row(), it.col(), it.value()));
			ATriplet.push_back(Eigen::Triplet<double>(it.col(), BLoc.cols() + it.row(), -it.value()));
		}
	}
	ALoc.setFromTriplets(ATriplet.begin(), ATriplet.end());
	cout << "LHS is constructed." << endl;
}

void VectorFields::solveLocalSystem()
{
	cout << "Starting to solve problem." << endl;
	XfLoc.resize(BLoc.rows(), 2);

	// Setting up the solver
	Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
	Eigen::VectorXd x;
	//solver.compute(A);

	solver.analyzePattern(ALoc);

	cout << "Factorizing LHS." << endl;
	solver.factorize(ALoc);
	if (solver.info() != Eigen::Success) {
		cout << "Factorization Failed." << endl;
		return;
	}

	printf("Size of B=%d\n", bLoc.rows());
	// FIRST BASIS
	cout << "Solving first problem." << endl;
	x = solver.solve(bLoc.col(0));
	if (solver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system." << endl;
		return;
	}
	

	XfLoc.col(0) = x.block(0, 0, BLoc.rows(), 1);

	// SECOND BASIS
	cout << "Solving second problem." << endl;
	x = solver.solve(bLoc.col(1));
	if (solver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system." << endl;
		return;
	}

	XfLoc.col(1) = x.block(0, 0, BLoc.rows(), 1);
	//for (int i = 0; i < B.rows(); i++) {
	//	//Xf.block(columnOrderMap(i), 0,1,1) = xx.block(i, 0, 1, 1);
	//}
	cout << "Problems solved." << endl;
	//cout << XfLoc << endl;
	printf("XfLoc=%dx%d\n", XfLoc.rows(), XfLoc.cols());
}

/* ====================== SETTING UP UTILITY MATRICES ============================*/ 
void VectorFields::constructSubdomainSingle(const int &source)
{
	sample = source; 

	priority_queue<VertexPair, std::vector<VertexPair>, std::greater<VertexPair>> DistPQueue;
	Eigen::VectorXd D(F.rows());
	const double maxDist = 4.0 * avgEdgeLength; 

	// Computing distance for initial sample points S
	for (int i = 0; i < F.rows(); i++) {
		D(i) = numeric_limits<double>::infinity();
	}

	D(source) = 0.0f;
	VertexPair vp{ source,D(source) };
	DistPQueue.push(vp);
	SubDomain.insert(source);

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
	} while (distFromCenter < maxDist );
	//} while (!DistPQueue.empty());

	cout << "There are " << SubDomain.size() << " elements in this sub-domain." << endl; 
}

void VectorFields::constructBoundary()
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


	cout << "Such subdomain has " << Boundary.size() << " elements in its boundary." << endl;
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