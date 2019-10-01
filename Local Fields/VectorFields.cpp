#include "VectorFields.h"

#include <igl/per_vertex_normals.h>
#include <igl/gaussian_curvature.h>
#include <igl/invert_diag.h>
#include <Eigen/Eigenvalues>
#include <random>
#include <Eigen/OrderingMethods>
#include <Eigen/CholmodSupport>
#include <suitesparse/cholmod.h>



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
	constructRandomHardConstraints();
	///constructSoftConstraints();
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
	C.resize(2 * constNum, B2DAsym.cols());

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
	C.resize(2 * constNum, B2DAsym.cols());
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
	const bool readFromFile = false;
	string filename = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Constraints/Constraints_Arma_Farthest_20.txt";;
	//string filename = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Constraints/Constraints_Cube_Rand_25.txt";

	if (readFromFile)
	{
		cout << "____Loading constraints from file \n";
		LoadSTDVectorFromTxtFile(filename, globalConstraints);
	}
	else {

		// Define the constraints
		const int numConstraints = 4;
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
		int curPoint = 0; 
		//int curPoint = dis(gen);
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

		WriteSTDVectorToTxtFile(globalConstraints, filename);
	}
	//////////////////////
		
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
		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * globalConstraints[i] + 0, 1.0));
		c(counter, 0) = cRand(0);
		counter++;

		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * globalConstraints[i] + 1, 1.0));
		c(counter, 0) = cRand(1);
		counter++;
	}
	C.resize(2 * globalConstraints.size(), B2DAsym.rows());
	C.setFromTriplets(CTriplet.begin(), CTriplet.end());
	//printf("Cp=%dx%d\n", C.rows(), C.cols());
}

void VectorFields::constructRandomHardConstraints()
{
	// Define the constraints
	const bool readFromFile = true;			/// IMPORTANT!!!!!!!!!!!!!!!!!!!!
	bool lineNotFound = true;
	//string filename = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Constraints/Constraints_CDragon_Rand_20.txt";;
	//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Armadillo_randConstraints.txt";
	string resultFile = "D:/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Fertility_randConstraints.txt";
	//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/CDragon_randConstraints.txt";
	//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Kitten_randConstraints.txt";
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
						//cout << oneWord << "|";
						getline(iStream, oneWord, ',');
					}
					//cout << endl;

					globalConstraints.resize(constraint.size());

					for (int i = 0; i < constraint.size();i++) {
						globalConstraints[i] = constraint[i];
					}

					lineNotFound = false; 
					//return;
				}

				//cout << "line =" << lineNow << endl; 
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
		std::uniform_int_distribution<> disConst(10, 50); // From 0 to F.rows()-1
		int numConstraints = disConst(gen);
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
	c.resize(2 * globalConstraints.size());
	Eigen::Vector2d cRand;

	int counter = 0;
	for (int i = 0; i < globalConstraints.size(); i++) {
		//cRand(0) = (double)(rand() % F.rows()) / (double) F.rows();
		//cRand(1) = (double)(rand() % F.rows()) / (double)F.rows();
		cRand << 1.0, 0.0;
		cRand.normalize();

		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * globalConstraints[i] + 0, 1.0));
		c(counter, 0) = cRand(0);
		counter++;

		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * globalConstraints[i] + 1, 1.0));
		c(counter, 0) = cRand(1);
		counter++;
	}
	C.resize(2 * globalConstraints.size(), B2DAsym.rows());
	C.setFromTriplets(CTriplet.begin(), CTriplet.end());
}

void VectorFields::pushNewUserConstraints(const int& fInit, const int& fEnd)
{
	userVisualConstraints.push_back(fInit);
	userVisualConstraints.push_back(fEnd);
}
void VectorFields::constructInteractiveConstraints()
{
	
}

void VectorFields::addHardConstraints()
{
	int oldNumConstr, newNumConstr;
	oldNumConstr = C.rows();
	
	/* Define the constraints */
	const int CRows = c.rows();
	const int curConstSize = userVisualConstraints.size();
	int constFid;
	Eigen::Vector2d normDir;

	printf("CRows: %d \n", CRows);

	/* Global constraints from user input */
	for (int i = curConstSize-2; i < curConstSize; i += 2)
	{
		/* Location of constraints */
		constFid = userVisualConstraints[i];
		globalConstraints.push_back(constFid);
		/* Getting the constraints + making them into local coordinates */
		Eigen::RowVector3d dir = FC.row(userVisualConstraints[i + 1]) - FC.row(constFid);
		cout << "Dir 3d: " << dir << endl;
		Eigen::MatrixXd ALoc(3, 2);
		ALoc = A.block(3 * constFid, 2 * constFid, 3, 2);
		normDir = ALoc.transpose() * dir.transpose();
		normDir.normalize();
		cout << "Dir 2d: " << normDir.transpose() << endl;
	}

	/* Setting up matrix C and column vector c */
	c.conservativeResize(CRows + 2);

	/* Putting the constraints into action */
	ConstrTriplet.push_back(Eigen::Triplet<double>(CRows + 0, 2 * constFid + 0, 1.0));
	c(CRows + 0) = normDir(0);
	ConstrTriplet.push_back(Eigen::Triplet<double>(CRows + 1, 2 * constFid + 1, 1.0));
	c(CRows + 1) = normDir(1);

	C.resize(0, 0);
	C.resize(CRows+2, B2DAsym.rows());
	C.setFromTriplets(ConstrTriplet.begin(), ConstrTriplet.end());

	printf("C:%dx%d | constraints: %d->%d | size: %d | entries: %d \n ", C.rows(),  C.cols(),constFid, userVisualConstraints[curConstSize - 1], curConstSize, ConstrTriplet.size());
	cout << "c: " << c.transpose() << endl << endl;

	newNumConstr = C.rows();
	deltaConstraints = newNumConstr - oldNumConstr;
}

void VectorFields::constructInteractiveConstraintsWithSingularities(igl::opengl::glfw::Viewer &viewer)
{
	/* Define the constraints */
	//cout << "The hard constraints part \n";

	int oldNumConstr, newNumConstr;
	oldNumConstr = C.rows();

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

	/* Counting the singularities*/
	int numSingConstraints = 0;
	for (int i = 0; i < SingNeighCC.size(); i++) {		
		for (int j = 0; j < (SingNeighCC[i].size()); j++) {		// Use only n-1 neighboring faces as constraints
			numSingConstraints++;
		}
	}
	
	/* Setting up matrix C and column vector c */
	vector<Eigen::Triplet<double>> CTriplet;
	CTriplet.reserve(2 * globalConstraints.size() + 2 * 4 * numSingConstraints);
	c.resize(2 * (globalConstraints.size() + numSingConstraints)); 

	/* Putting the constraints into action */
	for (int i = 0; i < globalConstraints.size(); i++) {
		CTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * globalConstraints[i] + 0, 1.0));
		c(2 * i + 0) = constraintValues[i](0);

		CTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * globalConstraints[i] + 1, 1.0));
		c(2 * i + 1) = constraintValues[i](1);
	}


	/* The singularity parts */
	//cout << "Constraints: The singularity parts \n";
	int counter = 2*globalConstraints.size();
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
			int i2 = (i < (SingNeighCC[id].size()) ? i + 1 : 0);

			// [a] obtain shared edges
			for (int f = 0; f < F.cols(); f++) {
				if (F(SingNeighCC[id][i], f) == singularities[id])
				{
					// [b] get the two edges
					edge1 = V.row(F(SingNeighCC[id][i], (f == 0 ? 2 : f - 1))) - V.row(F(SingNeighCC[id][i], f));
					edge2 = V.row(F(SingNeighCC[id][i], (f == 2 ? 0 : f + 1))) - V.row(F(SingNeighCC[id][i], f));
					angle = edge2.dot(edge1) / (edge1.norm()*edge2.norm());
					angle = acos(angle);

					// [c] get the angle
					internalAngle[id][i] = angle;
					gaussAngle[id] += angle;
				}
			}
		}
	}

	// Show the angles
	viewer.data().points.resize(0, 6);
	viewer.data().points.resize(SingNeighCC.size(), 6);
	for (int id = 0; id < SingNeighCC.size(); id++)
	{			
		//printf("> Singularity is in point %d \n", userSingularConstraints[id]);
		viewer.data().points.row(id) << V.row(userSingularConstraints[id]), Eigen::RowVector3d(0.0, 0.1, 0.9);
		//viewer.data().add_points(V.row(singularities[id]), Eigen::RowVector3d(0.0, 0.1, 0.9));
		
		////printf("__Gauss angle = %.4f\n", gaussAngle[id] * 180.0 / M_PI);
		//for (int i = 0; i < (SingNeighCC[id].size()); i++) {
		//	//printf("______ angle %d (F %d) = %.3f \n", i, SingNeighCC[id][i], internalAngle[id][i] * 180.0 / M_PI);
		//	Eigen::Vector3d firstBasis_ = A.block(3 * SingNeighCC[id][i], 2 * SingNeighCC[id][i], 3, 1);
		//	//viewer.data().add_edges(FC.row(SingNeighCC[id][i]), FC.row(SingNeighCC[id][i]) + firstBasis_.transpose().normalized()*avgEdgeLength*1.0, Eigen::RowVector3d(0.0, 1.0, 0.7));
		//}
	}

	// Show the angles
	for (int id = 0; id < SingNeighCC.size(); id++)					// For each singularity
	{
		double rotAngle = gaussAngle[id] / (double)SingNeighCC[id].size();
		Eigen::VectorXd inpCol(SingNeighCC.size()); inpCol.setLinSpaced(0.0, 1.0);
		Eigen::MatrixXd edgeCol; igl::jet(inpCol, true, edgeCol);

		for (int i = 0; i < (SingNeighCC[id].size()); i++) {		// iterate over all faces in one singularity
			int curFace = SingNeighCC[id][i];
			int nextFace;
			if (i == (SingNeighCC[id].size() - 1)) nextFace = SingNeighCC[id][0];
			else nextFace = SingNeighCC[id][i + 1];
			//if (i == SingNeighCC[id].size() - 2) SingNeighCC[id][SingNeighCC[id].size() - 1] = SingNeighCC[id][0];

			/* Iterate to find shared edge between two neighboring faces */
			int sharedEdgeID;
			for (int e1 = 0; e1 < 3; e1++) {						// iterate over edges sharing that face
				for (int e2 = 0; e2 < 3; e2++) {
					if (FE(curFace, e1) == FE(nextFace, e2))
					{
						sharedEdgeID = FE(curFace, e1);
					}
				}
			}
			///printf("|%d:%d => %d", curFace, nextFace, sharedEdgeID);

			/* Obtaining the transport angle (bring the T{i+1} to T{i}) */
			double angTarget, angSource;
			if (EF(sharedEdgeID, 0) == curFace)
			{
				angTarget = FrameRot(sharedEdgeID, 0);
				angSource = FrameRot(sharedEdgeID, 1);
				//printf("[0] Angle T1: %.3f | T2: %.3f \n", 180.0/M_PI*FrameRot(sharedEdgeID, 0), 180.0/M_PI*FrameRot(sharedEdgeID, 1));			
			}
			else if (EF(sharedEdgeID, 1) == curFace)
			{
				angTarget = FrameRot(sharedEdgeID, 1);
				angSource = FrameRot(sharedEdgeID, 0);
				//printf("[1] Angle T1: %.3f | T2: %.3f \n", 180.0 / M_PI*FrameRot(sharedEdgeID, 1), 180.0 / M_PI*FrameRot(sharedEdgeID, 0));
			}
			///printf("Angle T1: %.3f | T2: %.3f \n", 180.0 / M_PI*angTarget, 180.0 / M_PI*angSource);


			double totalRot = -rotAngle + angTarget - angSource + M_PI;

			///printf("RotAngle: %.4f | transportAngle: %.4f | target: %.4f | source: %.4f \n", -rotAngle*180.0 / M_PI, (angTarget - angSource + M_PI)*180.0 / M_PI, angTarget*180.0 / M_PI, angSource*180.0 / M_PI);

			Eigen::Matrix2d transfRotMat; transfRotMat << cos(totalRot), -sin(totalRot), sin(totalRot), cos(totalRot);

			// vector for visualization
			// -- vector in my (next) neighbor
			Eigen::Vector2d vn; vn << 1.0, 0.0;
			Eigen::Vector2d vm = transfRotMat*vn;

			Eigen::MatrixXd An; An = A.block(3 * nextFace, 2 * nextFace, 3, 2);
			Eigen::MatrixXd Am; Am = A.block(3 * curFace, 2 * curFace, 3, 2);

			Eigen::Vector3d edgeN = An*vn;
			Eigen::Vector3d edgeM = Am*vm;

			CTriplet.push_back(Eigen::Triplet<double>(counter + 0, 2 * nextFace + 0, transfRotMat(0, 0)));	// the reference triangle
			CTriplet.push_back(Eigen::Triplet<double>(counter + 1, 2 * nextFace + 0, transfRotMat(1, 0)));
			CTriplet.push_back(Eigen::Triplet<double>(counter + 0, 2 * nextFace + 1, transfRotMat(0, 1)));
			CTriplet.push_back(Eigen::Triplet<double>(counter + 1, 2 * nextFace + 1, transfRotMat(1, 1)));

			CTriplet.push_back(Eigen::Triplet<double>(counter + 0, 2 * curFace + 0, -1.0));					// the neighbor (next, CCW)
			CTriplet.push_back(Eigen::Triplet<double>(counter + 1, 2 * curFace + 1, -1.0));

			c(counter + 0) = 0.0;
			c(counter + 1) = 0.0;
			counter += 2;

			viewer.data().add_edges(FC.row(curFace), FC.row(curFace) + edgeM.transpose().normalized()*avgEdgeLength, edgeCol.row(i));
			viewer.data().add_edges(FC.row(nextFace), FC.row(nextFace) + edgeN.transpose().normalized()*avgEdgeLength, edgeCol.row(i));
		}
	}
	cout << "Settin gup C\n";
	C.resize(2 * (globalConstraints.size() + numSingConstraints), B2DAsym.rows());
	///printf("counter: %d | C: %dx%d \n", counter, C.rows(), C.cols());

	C.setFromTriplets(CTriplet.begin(), CTriplet.end());

	newNumConstr = C.rows();
	deltaConstraints = newNumConstr - oldNumConstr;
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

	C.resize(4 * globalConstraints.size(), B2DAsym.rows());
	C.setFromTriplets(CTriplet.begin(), CTriplet.end());
	//visualizeSparseMatrixInMatlab(C);
}

void VectorFields::resetInteractiveConstraints()
{
	userVisualConstraints.clear();
	userVisualConstraints.shrink_to_fit();

	userSingularConstraints.clear();
	userSingularConstraints.shrink_to_fit();

	C.resize(0, 0);
	c.resize(0);

	globalConstraints.clear(); 
	globalConstraints.shrink_to_fit();

	ConstrTriplet.clear();
	ConstrTriplet.shrink_to_fit();
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

void VectorFields::constructInteractiveSingularities()
{
	singularities.clear();
	SingNeighCC.clear();
	mappedBasis.clear();
	mappedBasis2.clear();
	sharedEdgesVect.clear();

	//cout << "Getting vertex singularities from user input \n";
	const int NUM_SINGS = userSingularConstraints.size();

	//cout << "Singular points: ";
	//for (int i : userSingularConstraints) cout << " " << i << "|";
	//cout << endl; 

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
		const int SingLocation = userSingularConstraints[id];
		singularities[id] = SingLocation;
		const int SingNeighNum = VFAdjacency.col(SingLocation).nonZeros();
		Eigen::SparseMatrix<bool>::InnerIterator it0(VFAdjacency, SingLocation);
		const int firstNeigh = it0.row();


		// Inserting the first neighbor (the one with lowest index/row number)
		SingNeighCC[id].resize(SingNeighNum);
		SingNeighCC[id][0] = firstNeigh;
		int curNeigh = firstNeigh;
		int vertex1 = SingLocation;
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

}

void VectorFields::addSingularityConstraints()
{
	/* Defining LOCATION for the singularities */
	//time_t t;
	//srand((unsigned)time(&t));
	srand(time(NULL));

	int oldNumConstr, newNumConstr;
	oldNumConstr = C.rows();

	const int id = userSingularConstraints.size() - 1;
	int CRows = C.rows();
	
	// Defining varaibles for singularities
	const int SingLocation = userSingularConstraints[id];
	singularities.push_back(SingLocation);
	const int SingNeighNum = VFAdjacency.col(SingLocation).nonZeros();
	Eigen::SparseMatrix<bool>::InnerIterator it0(VFAdjacency, SingLocation);
	const int firstNeigh = it0.row();


	// Inserting the first neighbor (the one with lowest index/row number)
	vector<int> emptyVect;
	SingNeighCC.push_back(emptyVect);
	SingNeighCC[id].resize(SingNeighNum);
	SingNeighCC[id][0] = firstNeigh;
	int curNeigh = firstNeigh;
	int vertex1 = SingLocation;

	// Getting the neighboring valence triangles in order
	for (int i2 = 1; i2 < SingNeighNum; i2++) {
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

	/* Defining the VALUES for singularities */
	int numSingConstraints = 0;
	
	//for (int j = 0; j < (SingNeighCC[i].size()-1); j++) {		// no consntraints on the last triangle
	for (int j = 0; j < (SingNeighCC[id].size()); j++) {
		numSingConstraints++;
	}
	//printf("Sing vert: %d | %d faces \n", SingLocation, numSingConstraints);

	// Setting up matrix C and vector c
	c.conservativeResize(C.rows() + 2*numSingConstraints);

	
	// Setting up hard constraints for neighboring faces
	Eigen::MatrixXd		ALoc(3, 2);
	Eigen::RowVector3d	c1, c2, field3D;
	Eigen::Vector2d		field2D;
	vector<vector<double>>		internalAngle(SingNeighCC.size());
	vector<double>				gaussAngle(SingNeighCC.size());

	// Local variables to compute angle on each triangle
	//cout << "Local variables to compute angle on each triangle \n";
	Eigen::Vector3d		edge1, edge2;
	double				angle;
	int					singLoc;
	
	/* Computing the gauss angles */
	internalAngle[id].resize(SingNeighCC[id].size());
	gaussAngle[id] = 0.0;

	for (int i = 0; i < (SingNeighCC[id].size()); i++)
	{
		int i2 = (i < (SingNeighCC[id].size()) ? i + 1 : 0);

		// [a] obtain shared edges
		for (int f = 0; f < F.cols(); f++) {
			if (F(SingNeighCC[id][i], f) == singularities[id])
			{
				// [b] get the two edges
				edge1 = V.row(F(SingNeighCC[id][i], (f == 0 ? 2 : f - 1))) - V.row(F(SingNeighCC[id][i], f));
				edge2 = V.row(F(SingNeighCC[id][i], (f == 2 ? 0 : f + 1))) - V.row(F(SingNeighCC[id][i], f));
				angle = edge2.dot(edge1) / (edge1.norm()*edge2.norm());
				angle = acos(angle);

				// [c] get the angle
				internalAngle[id][i] = angle;
				gaussAngle[id] += angle;
			}
		}
	}

	// Show the angles
	//for (int id = 0; id < SingNeighCC.size(); id++)
	//{
	//	viewer.data().add_points(V.row(singularities[id]), Eigen::RowVector3d(1.0, 0.1, 0.1));
	//	printf("__Gauss angle = %.4f\n", gaussAngle[id] * 180.0 / M_PI);
	//	for (int i = 0; i < (SingNeighCC[id].size()); i++) {
	//		printf("______ angle %d (F %d) = %.3f \n", i, SingNeighCC[id][i], internalAngle[id][i] * 180.0 / M_PI);
	//		Eigen::Vector3d firstBasis_ = A.block(3 * SingNeighCC[id][i], 2 * SingNeighCC[id][i], 3, 1);
	//		//viewer.data().add_edges(FC.row(SingNeighCC[id][i]), FC.row(SingNeighCC[id][i]) + firstBasis_.transpose().normalized()*avgEdgeLength*1.0, Eigen::RowVector3d(0.0, 1.0, 0.7));
	//	}
	//}


	///cout << "Start assigning values \n";
	int counter = C.rows();
	double rotAngle = gaussAngle[id] / (double)SingNeighCC[id].size();
	Eigen::VectorXd inpCol(SingNeighCC.size()); inpCol.setLinSpaced(0.0, 1.0);
	Eigen::MatrixXd edgeCol; igl::jet(inpCol, true, edgeCol);

	for (int i = 0; i < (SingNeighCC[id].size()); i++) {		// iterate over all faces in one singularity
		int curFace = SingNeighCC[id][i];
		int nextFace;
		if (i == (SingNeighCC[id].size() - 1)) nextFace = SingNeighCC[id][0];
		else nextFace = SingNeighCC[id][i + 1];
		//if (i == SingNeighCC[id].size() - 2) SingNeighCC[id][SingNeighCC[id].size() - 1] = SingNeighCC[id][0];

		/* Iterate to find shared edge between two neighboring faces */
		int sharedEdgeID;
		for (int e1 = 0; e1 < 3; e1++) {						// iterate over edges sharing that face
			for (int e2 = 0; e2 < 3; e2++) {
				if (FE(curFace, e1) == FE(nextFace, e2))
				{
					sharedEdgeID = FE(curFace, e1);
				}
			}
		}
		//printf("|%d:%d => %d \n", curFace, nextFace, sharedEdgeID);

		/* Obtaining the transport angle (bring the T{i+1} to T{i}) */
		double angTarget, angSource;
		if (EF(sharedEdgeID, 0) == curFace)
		{
			angTarget = FrameRot(sharedEdgeID, 0);
			angSource = FrameRot(sharedEdgeID, 1);
			//printf("[0] Angle T1: %.3f | T2: %.3f \n", 180.0/M_PI*FrameRot(sharedEdgeID, 0), 180.0/M_PI*FrameRot(sharedEdgeID, 1));			
		}
		else if (EF(sharedEdgeID, 1) == curFace)
		{
			angTarget = FrameRot(sharedEdgeID, 1);
			angSource = FrameRot(sharedEdgeID, 0);
			//printf("[1] Angle T1: %.3f | T2: %.3f \n", 180.0 / M_PI*FrameRot(sharedEdgeID, 1), 180.0 / M_PI*FrameRot(sharedEdgeID, 0));
		}

		double totalRot = -rotAngle + angTarget - angSource + M_PI;

		Eigen::Matrix2d transfRotMat; transfRotMat << cos(totalRot), -sin(totalRot), sin(totalRot), cos(totalRot);

		// vector for visualization
		// -- vector in my (next) neighbor
		Eigen::Vector2d vn; vn << 1.0, 0.0;
		Eigen::Vector2d vm = transfRotMat*vn;

		Eigen::MatrixXd An; An = A.block(3 * nextFace, 2 * nextFace, 3, 2);
		Eigen::MatrixXd Am; Am = A.block(3 * curFace, 2 * curFace, 3, 2);

		Eigen::Vector3d edgeN = An*vn;
		Eigen::Vector3d edgeM = Am*vm;

		ConstrTriplet.push_back(Eigen::Triplet<double>(counter + 0, 2 * nextFace + 0, transfRotMat(0, 0)));	// the reference triangle
		ConstrTriplet.push_back(Eigen::Triplet<double>(counter + 1, 2 * nextFace + 0, transfRotMat(1, 0)));
		ConstrTriplet.push_back(Eigen::Triplet<double>(counter + 0, 2 * nextFace + 1, transfRotMat(0, 1)));
		ConstrTriplet.push_back(Eigen::Triplet<double>(counter + 1, 2 * nextFace + 1, transfRotMat(1, 1)));

		ConstrTriplet.push_back(Eigen::Triplet<double>(counter + 0, 2 * curFace + 0, -1.0));					// the neighbor (next, CCW)
		ConstrTriplet.push_back(Eigen::Triplet<double>(counter + 1, 2 * curFace + 1, -1.0));
		
		c(counter + 0) = 0.0;
		c(counter + 1) = 0.0;
		counter += 2;
	}

		//cout << "Setting up the matrix for constraint \n";
	C.resize(0, 0);
	C.resize(CRows + 2 * numSingConstraints, B2DAsym.rows());
	printf("C: %dx%d \n", C.rows(), C.cols());
	C.setFromTriplets(ConstrTriplet.begin(), ConstrTriplet.end());

	newNumConstr = C.rows();
	deltaConstraints = newNumConstr - oldNumConstr;
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

	C.resize(2 * (globalConstraints.size() + numSingConstraints), B2DAsym.rows());
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

	C.resize(2 * (globalConstraints.size() + numSingConstraints), B2DAsym.rows());
	C.setFromTriplets(CTriplet.begin(), CTriplet.end());
}

void VectorFields::constructHardConstraintsWithSingularitiesWithGauss(igl::opengl::glfw::Viewer &viewer)
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
		//for (int j = 0; j < (SingNeighCC[i].size()-1); j++) {		// no consntraints on the last triangle
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
		//c(counter++, 0) = sqrt(2.0);
		c(counter++, 0) = sqrt(1.0);
		//c(counter++, 1) = sqrt(2.0);
		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * globalConstraints[i] + 1, 1.0));
		//c(counter++, 0) = sqrt(2.0);
		c(counter++, 0) = sqrt(0.0);
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
			int i2 = (i < (SingNeighCC[id].size()) ? i + 1 : 0);

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
		viewer.data().add_points(V.row(singularities[id]), Eigen::RowVector3d(1.0, 0.1, 0.1));
		printf("__Gauss angle = %.4f\n", gaussAngle[id]*180.0/M_PI);
		for (int i = 0; i < (SingNeighCC[id].size()); i++) {
			printf("______ angle %d (F %d) = %.3f \n", i, SingNeighCC[id][i], internalAngle[id][i] * 180.0 / M_PI);
			Eigen::Vector3d firstBasis_ = A.block(3 * SingNeighCC[id][i], 2 * SingNeighCC[id][i], 3, 1);
			//viewer.data().add_edges(FC.row(SingNeighCC[id][i]), FC.row(SingNeighCC[id][i]) + firstBasis_.transpose().normalized()*avgEdgeLength*1.0, Eigen::RowVector3d(0.0, 1.0, 0.7));
		}
	}

	// Show the angles
	for (int id = 0; id < SingNeighCC.size(); id++)					// For each singularity
	{		
		double rotAngle = gaussAngle[id] / (double)SingNeighCC[id].size();
		Eigen::VectorXd inpCol(SingNeighCC.size()); inpCol.setLinSpaced(0.0, 1.0);
		Eigen::MatrixXd edgeCol; igl::jet(inpCol, true, edgeCol);

		for (int i = 0; i < (SingNeighCC[id].size()); i++) {		// iterate over all faces in one singularity
			int curFace = SingNeighCC[id][i];
			int nextFace;
			if (i == (SingNeighCC[id].size() - 1)) nextFace = SingNeighCC[id][0];
			else nextFace = SingNeighCC[id][i + 1];
			//if (i == SingNeighCC[id].size() - 2) SingNeighCC[id][SingNeighCC[id].size() - 1] = SingNeighCC[id][0];

			/* Iterate to find shared edge between two neighboring faces */
			int sharedEdgeID;
			for (int e1 = 0; e1 < 3; e1++) {						// iterate over edges sharing that face
				for (int e2 = 0; e2 < 3; e2++) {
					if (FE(curFace, e1) == FE(nextFace, e2))
					{
						sharedEdgeID = FE(curFace, e1);
					}
				}
			}
			printf("|%d:%d => %d", curFace, nextFace, sharedEdgeID);

			/* Obtaining the transport angle (bring the T{i+1} to T{i}) */
			double angTarget, angSource;
			if (EF(sharedEdgeID, 0) == curFace)
			{
				angTarget = FrameRot(sharedEdgeID, 0);
				angSource = FrameRot(sharedEdgeID, 1);
				//printf("[0] Angle T1: %.3f | T2: %.3f \n", 180.0/M_PI*FrameRot(sharedEdgeID, 0), 180.0/M_PI*FrameRot(sharedEdgeID, 1));			
			}
			else if (EF(sharedEdgeID, 1) == curFace)
			{
				angTarget = FrameRot(sharedEdgeID, 1);
				angSource = FrameRot(sharedEdgeID, 0);
				//printf("[1] Angle T1: %.3f | T2: %.3f \n", 180.0 / M_PI*FrameRot(sharedEdgeID, 1), 180.0 / M_PI*FrameRot(sharedEdgeID, 0));
			}	
			///printf("Angle T1: %.3f | T2: %.3f \n", 180.0 / M_PI*angTarget, 180.0 / M_PI*angSource);


			double totalRot = -rotAngle + angTarget - angSource + M_PI;

			///printf("RotAngle: %.4f | transportAngle: %.4f | target: %.4f | source: %.4f \n", -rotAngle*180.0 / M_PI, (angTarget - angSource + M_PI)*180.0 / M_PI, angTarget*180.0 / M_PI, angSource*180.0 / M_PI);

			Eigen::Matrix2d transfRotMat; transfRotMat << cos(totalRot), -sin(totalRot), sin(totalRot), cos(totalRot);

			// vector for visualization
			// -- vector in my (next) neighbor
			Eigen::Vector2d vn; vn << 1.0, 0.0;
			Eigen::Vector2d vm = transfRotMat*vn;

			Eigen::MatrixXd An; An = A.block(3 * nextFace, 2 * nextFace, 3, 2);
			Eigen::MatrixXd Am; Am = A.block(3 * curFace, 2 * curFace, 3, 2);
			
			Eigen::Vector3d edgeN = An*vn;
			Eigen::Vector3d edgeM = Am*vm;

			CTriplet.push_back(Eigen::Triplet<double>(counter + 0, 2 * nextFace + 0, transfRotMat(0, 0)));	// the reference triangle
			CTriplet.push_back(Eigen::Triplet<double>(counter + 1, 2 * nextFace + 0, transfRotMat(1, 0)));
			CTriplet.push_back(Eigen::Triplet<double>(counter + 0, 2 * nextFace + 1, transfRotMat(0, 1)));
			CTriplet.push_back(Eigen::Triplet<double>(counter + 1, 2 * nextFace + 1, transfRotMat(1, 1)));

			CTriplet.push_back(Eigen::Triplet<double>(counter + 0, 2 * curFace + 0, -1.0));					// the neighbor (next, CCW)
			CTriplet.push_back(Eigen::Triplet<double>(counter + 1, 2 * curFace + 1, -1.0));

			c(counter + 0) = 0.0;
			c(counter + 1) = 0.0;
			counter+=2;

			viewer.data().add_edges(FC.row(curFace), FC.row(curFace) + edgeM.transpose().normalized()*avgEdgeLength, edgeCol.row(i));
			viewer.data().add_edges(FC.row(nextFace), FC.row(nextFace) + edgeN.transpose().normalized()*avgEdgeLength, edgeCol.row(i));
		}

		//for (int i = 0; i < (SingNeighCC[id].size()); i++) {
		//	printf("Edges in %d : ", SingNeighCC[id][i]);
		//	for (int j = 0; j < 3; j++) {
		//		printf("|%d ", FE(SingNeighCC[id][i], j));
		//	}
		//	cout << endl; 
		//}
		
	}

	C.resize(2 * (globalConstraints.size() + numSingConstraints), B2DAsym.rows());
	C.setFromTriplets(CTriplet.begin(), CTriplet.end());

	//// SINGULARITIES CONSTRAINTS
	//for (int id = 0; id < SingNeighCC.size(); id++) 
	//{
	//	// Getting the shared-edges of two neighboring faces For testing
	//	sharedEdgesVect[id].resize(2 * SingNeighCC[id].size() - 2);
	//
	//	// 4. Compute rotation of among its valence
	//	const double rotAngle = gaussAngle[id] / (double)SingNeighCC[id].size();	// all edges with similar angle => next iter: relative angle
	//	const double cosConst = cos(/*2*M_PI - */rotAngle);
	//	const double sinConst = sin(/*2*M_PI - */rotAngle);
	//	printf("Angle=%.3f, sin=%.3f, cos=%.3f\n", rotAngle*180.0 / M_PI, sinConst, cosConst);
	//
	//	/* Give hard constraint on face 1 */
	//	Eigen::MatrixXd ALoc(3, 2);
	//	for (int f = 0; f < F.cols(); f++) {
	//		if (F(SingNeighCC[id][0], f) == singularities[id])
	//		{
	//			Eigen::Vector3d edge = V.row(F(SingNeighCC[id][0], (f==0 ? 2 : f-1))) - V.row(F(SingNeighCC[id][0], (f == 2 ? 0 : f + 1)));
	//			ALoc = A.block(3 * SingNeighCC[id][0], 2 * SingNeighCC[id][0], 3, 2);
	//			Eigen::Vector2d edge2D = ALoc.transpose() * edge; 
	//			edge2D = edge2D.normalized() / 4.0; 
	//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][0]+0, 1.0));
	//			c(counter) = edge2D(0);
	//			counter++;
	//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][0]+1, 1.0));
	//			c(counter) = edge2D(1);
	//			counter++;
	//		}
	//	}
	//
	//	// Which case? => determining which edge is the common edge
	//	enum class SharedEdgeCase { Case1, Case2, Case3 };
	//	SharedEdgeCase edgeCase1, edgeCase2;
	//	for (int i = 0; i < (SingNeighCC[id].size() - 1); i++) 
	//	{
	//		// 1. Find shared edge (naively)
	//		Eigen::RowVector3d es;
	//		for (int f1 = 0; f1 < F.cols(); f1++) {
	//			for (int f2 = 0; f2 < F.cols(); f2++) {
	//				bool b1 = F(SingNeighCC[id][i], (f1 + 1) % F.cols()) == F(SingNeighCC[id][i + 1], f2);
	//				bool b2 = F(SingNeighCC[id][i], f1) == F(SingNeighCC[id][i + 1], (f2 + 1) % F.cols());
	//				if (b1 && b2) {
	//					sharedEdgesVect[id][2 * i + 0] = F(SingNeighCC[id][i], f1);
	//					sharedEdgesVect[id][2 * i + 1] = F(SingNeighCC[id][i], (f1 + 1) % F.cols());
	//					es = V.row(F(SingNeighCC[id][i], (f1 + 1) % F.cols())) - V.row(F(SingNeighCC[id][i], f1));
	//					printf("Shared edge=%d->%d\n", F(SingNeighCC[id][i], f1), F(SingNeighCC[id][i], (f1 + 1) % F.cols()));
	//
	//					if (f1 == 0)		edgeCase1 = SharedEdgeCase::Case1;	// => edge V0->V1 is the shared edge => it takes 0 step to reach v0
	//					else if (f1 == 1)	edgeCase1 = SharedEdgeCase::Case3;	// => edge V1->V2 is the shared edge => it takes 2 step to reach v0
	//					else if (f1 == 2)	edgeCase1 = SharedEdgeCase::Case2;	// => edge V2->V0 is the shared edge => it takes 1 step to reach v0
	//
	//					if (f2 == 0)		edgeCase2 = SharedEdgeCase::Case1;
	//					else if (f2 == 1)	edgeCase2 = SharedEdgeCase::Case3;
	//					else if (f2 == 2)	edgeCase2 = SharedEdgeCase::Case2;
	//				}
	//			}
	//		}
	//
	//		// 2. Find angles between basis1 and shared_edges es
	//		Eigen::VectorXd eVect;
	//		Eigen::RowVector3d b11, b12;
	//		eVect = A.block(3 * SingNeighCC[id][i], 2 * SingNeighCC[id][i] + 0, 3, 1);
	//		b11 << eVect(0), eVect(1), eVect(2);
	//		eVect = A.block(3 * SingNeighCC[id][i], 2 * SingNeighCC[id][i] + 1, 3, 1);
	//		b12 << eVect(0), eVect(1), eVect(2);
	//		//cout << "______B11: " << b11 << ", B12: " << b12 << endl;
	//
	//		// Basis 1, Frame 1
	//		double cosR12 = (b11.dot(es)) / (b11.norm()*es.norm());
	//		if (cosR12 > 1.0) cosR12 = 1.0;
	//		if (cosR12 <-1.0) cosR12 = -1.0;
	//		const double angleR12_1 = (edgeCase1 == SharedEdgeCase::Case2 ? 2 * M_PI - acos(cosR12) : acos(cosR12));
	//		printf("______[%.2f] Rotation matrix R12_1\n", angleR12_1*180.0 / M_PI);
	//
	//		// 3. Find angles between basis2 and es
	//		Eigen::RowVector3d b21, b22;
	//		eVect = A.block(3 * SingNeighCC[id][i + 1], 2 * SingNeighCC[id][i + 1] + 0, 3, 1);
	//		b21 << eVect(0), eVect(1), eVect(2);
	//		eVect = A.block(3 * SingNeighCC[id][i + 1], 2 * SingNeighCC[id][i + 1] + 1, 3, 1);
	//		b22 << eVect(0), eVect(1), eVect(2);
	//		//cout << "______B21: " << b21 << ", B22: " << b22 << endl;
	//
	//		// Basis 2, Frame 1
	//		double cosR21 = (b21.dot(es)) / (b21.norm()*es.norm());
	//		if (cosR21 > 1.0) cosR21 = 1.0;
	//		if (cosR21 < -1.0) cosR21 = -1.0;
	//		double angleR21_1 = (edgeCase2 == SharedEdgeCase::Case3 ? 2 * M_PI - acos(cosR21) : acos(cosR21));
	//		angleR21_1 = 2 * M_PI - angleR21_1;
	//		printf("______[%.2f] Rotation matrix R22_1 = [%.2f]\n", angleR21_1*180.0 / M_PI);
	//		
	//		const double RotAngle = (angleR12_1 + angleR21_1 > 2 * M_PI ? (angleR12_1 + angleR21_1) - 2 * M_PI : angleR12_1 + angleR21_1);
	//		const double cosBasis = cos(RotAngle);
	//		const double sinBasis = sin(RotAngle);
	//		printf("____ To map basis1 -> basis2: rotate by %.2f degree (cos=%.2f, cin=%.2f))\n", (RotAngle)*180.0 / M_PI, cosBasis, sinBasis);
	//
	//
	//		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, cosBasis));
	//		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, -sinBasis));
	//		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, -cosConst));
	//		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, sinConst));
	//		c(counter) = 0.0;
	//		counter++;
	//
	//		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, sinBasis));
	//		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, cosBasis));
	//		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, -sinConst));
	//		CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, -cosConst));
	//		c(counter) = 0.0;
	//		counter++;
	//	}
	//}

	
}

//void VectorFields::constructHardConstraintsWithSingularitiesWithGauss()
//{
//	// Define the constraints
//	const int numConstraints = 20;
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
//		//for (int j = 0; j < (SingNeighCC[i].size()-1); j++) {
//		for (int j = 0; j < (SingNeighCC[i].size()); j++) {
//			numSingConstraints++;
//		}
//	}
//
//	// Setting up matrix C and vector c
//	c.resize(2 * (globalConstraints.size() + numSingConstraints));
//
//
//	// HARD CONSTRAINTS
//	Eigen::SparseMatrix<double> CTemp;
//	vector<Eigen::Triplet<double>> CTriplet;
//	CTriplet.reserve(2 * globalConstraints.size() + 2 * 4 * 7 * SingNeighCC.size());
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
//	// Setting up hard constraints for neighboring faces
//	Eigen::MatrixXd		ALoc(3, 2);
//	Eigen::RowVector3d	c1, c2, field3D;
//	Eigen::Vector2d		field2D;
//	vector<vector<double>>		internalAngle(SingNeighCC.size());
//	vector<double>				gaussAngle(SingNeighCC.size());
//
//	// Local variables to compute angle on each triangle
//	Eigen::Vector3d		edge1, edge2;
//	double				angle;
//	int					singLoc;
//	for (int id = 0; id < SingNeighCC.size(); id++)
//	{
//		internalAngle[id].resize(SingNeighCC[id].size());
//		gaussAngle[id] = 0.0;
//
//		for (int i = 0; i < (SingNeighCC[id].size()); i++)
//		{
//			int i2 = (i < (SingNeighCC[id].size() - 1) ? i + 1 : 0);
//
//			// [a] obtain shared edges
//			for (int f = 0; f < F.cols(); f++) {
//				if (F(SingNeighCC[id][i], f) == singularities[id])
//				{
//					// [b] get the two edges
//					edge1 = V.row(F(SingNeighCC[id][i], (f == 0 ? 2 : f - 1))) - V.row(F(SingNeighCC[id][i], f));
//					edge2 = V.row(F(SingNeighCC[id][i], (f == 2 ? 0 : f + 1))) - V.row(F(SingNeighCC[id][i], f));
//					angle = edge2.dot(edge1) / (edge1.norm()*edge2.norm());
//					angle = acos(angle);
//
//					// [c] get the angle
//					internalAngle[id][i] = angle;
//					gaussAngle[id] += angle;
//				}
//			}
//		}
//	}
//
//	// Show the angles
//	for (int id = 0; id < SingNeighCC.size(); id++)
//	{
//		printf("__Gauss angle = %.4f\n", gaussAngle[id] * 180.0 / M_PI);
//		for (int i = 0; i < (SingNeighCC[id].size()); i++) {
//			printf("______ angle %d = %.3f \n", i, internalAngle[id][i] * 180.0 / M_PI);
//		}
//	}
//
//	// SINGULARITIES CONSTRAINTS
//	for (int id = 0; id < SingNeighCC.size(); id++)
//	{
//		// Getting the shared-edges of two neighboring faces For testing
//		sharedEdgesVect[id].resize(2 * SingNeighCC[id].size() - 2);
//
//		// 4. Compute rotation of among its valence
//		const double rotAngle = gaussAngle[id] / (double)SingNeighCC[id].size();	// all edges with similar angle => next iter: relative angle
//		const double cosConst = cos(/*2*M_PI - */rotAngle);
//		const double sinConst = sin(/*2*M_PI - */rotAngle);
//		printf("Angle=%.3f, sin=%.3f, cos=%.3f\n", rotAngle*180.0 / M_PI, sinConst, cosConst);
//
//		/* Give hard constraint on face 1 */
//		Eigen::MatrixXd ALoc(3, 2);
//		for (int f = 0; f < F.cols(); f++) {
//			if (F(SingNeighCC[id][0], f) == singularities[id])
//			{
//				Eigen::Vector3d edge = V.row(F(SingNeighCC[id][0], (f == 0 ? 2 : f - 1))) - V.row(F(SingNeighCC[id][0], (f == 2 ? 0 : f + 1)));
//				ALoc = A.block(3 * SingNeighCC[id][0], 2 * SingNeighCC[id][0], 3, 2);
//				Eigen::Vector2d edge2D = ALoc.transpose() * edge;
//				edge2D = edge2D.normalized() / 4.0;
//				CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][0] + 0, 1.0));
//				c(counter) = edge2D(0);
//				counter++;
//				CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][0] + 1, 1.0));
//				c(counter) = edge2D(1);
//				counter++;
//			}
//		}
//
//		// Which case? => determining which edge is the common edge
//		enum class SharedEdgeCase { Case1, Case2, Case3 };
//		SharedEdgeCase edgeCase1, edgeCase2;
//		for (int i = 0; i < (SingNeighCC[id].size() - 1); i++)
//		{
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
//						if (f1 == 0)		edgeCase1 = SharedEdgeCase::Case1;	// => edge V0->V1 is the shared edge => it takes 0 step to reach v0
//						else if (f1 == 1)	edgeCase1 = SharedEdgeCase::Case3;	// => edge V1->V2 is the shared edge => it takes 2 step to reach v0
//						else if (f1 == 2)	edgeCase1 = SharedEdgeCase::Case2;	// => edge V2->V0 is the shared edge => it takes 1 step to reach v0
//
//						if (f2 == 0)		edgeCase2 = SharedEdgeCase::Case1;
//						else if (f2 == 1)	edgeCase2 = SharedEdgeCase::Case3;
//						else if (f2 == 2)	edgeCase2 = SharedEdgeCase::Case2;
//					}
//				}
//			}
//
//			// 2. Find angles between basis1 and shared_edges es
//			Eigen::VectorXd eVect;
//			Eigen::RowVector3d b11, b12;
//			eVect = A.block(3 * SingNeighCC[id][i], 2 * SingNeighCC[id][i] + 0, 3, 1);
//			b11 << eVect(0), eVect(1), eVect(2);
//			eVect = A.block(3 * SingNeighCC[id][i], 2 * SingNeighCC[id][i] + 1, 3, 1);
//			b12 << eVect(0), eVect(1), eVect(2);
//			//cout << "______B11: " << b11 << ", B12: " << b12 << endl;
//
//			// Basis 1, Frame 1
//			double cosR12 = (b11.dot(es)) / (b11.norm()*es.norm());
//			if (cosR12 > 1.0) cosR12 = 1.0;
//			if (cosR12 <-1.0) cosR12 = -1.0;
//			const double angleR12_1 = (edgeCase1 == SharedEdgeCase::Case2 ? 2 * M_PI - acos(cosR12) : acos(cosR12));
//			printf("______[%.2f] Rotation matrix R12_1\n", angleR12_1*180.0 / M_PI);
//
//			// 3. Find angles between basis2 and es
//			Eigen::RowVector3d b21, b22;
//			eVect = A.block(3 * SingNeighCC[id][i + 1], 2 * SingNeighCC[id][i + 1] + 0, 3, 1);
//			b21 << eVect(0), eVect(1), eVect(2);
//			eVect = A.block(3 * SingNeighCC[id][i + 1], 2 * SingNeighCC[id][i + 1] + 1, 3, 1);
//			b22 << eVect(0), eVect(1), eVect(2);
//			//cout << "______B21: " << b21 << ", B22: " << b22 << endl;
//
//			// Basis 2, Frame 1
//			double cosR21 = (b21.dot(es)) / (b21.norm()*es.norm());
//			if (cosR21 > 1.0) cosR21 = 1.0;
//			if (cosR21 < -1.0) cosR21 = -1.0;
//			double angleR21_1 = (edgeCase2 == SharedEdgeCase::Case3 ? 2 * M_PI - acos(cosR21) : acos(cosR21));
//			angleR21_1 = 2 * M_PI - angleR21_1;
//			printf("______[%.2f] Rotation matrix R22_1 = [%.2f]\n", angleR21_1*180.0 / M_PI);
//
//			const double RotAngle = (angleR12_1 + angleR21_1 > 2 * M_PI ? (angleR12_1 + angleR21_1) - 2 * M_PI : angleR12_1 + angleR21_1);
//			const double cosBasis = cos(RotAngle);
//			const double sinBasis = sin(RotAngle);
//			printf("____ To map basis1 -> basis2: rotate by %.2f degree (cos=%.2f, cin=%.2f))\n", (RotAngle)*180.0 / M_PI, cosBasis, sinBasis);
//
//
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, cosBasis));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, -sinBasis));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, -cosConst));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, sinConst));
//			c(counter) = 0.0;
//			counter++;
//
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 0, sinBasis));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i] + 1, cosBasis));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 0, -sinConst));
//			CTriplet.push_back(Eigen::Triplet<double>(counter, 2 * SingNeighCC[id][i + 1] + 1, -cosConst));
//			c(counter) = 0.0;
//			counter++;
//		}
//	}
//
//	C.resize(2 * (globalConstraints.size() + numSingConstraints), B2D.rows());
//	C.setFromTriplets(CTriplet.begin(), CTriplet.end());
//}

void VectorFields::constructCurves_antiqueHead()
{
	const int NUM_CURVES = 36;
	curvesConstraints.resize(NUM_CURVES);

	srand(time(NULL));
	int init_, end_;
	vector<int> aCurve;

	int constCounter = 0;
	// top strokes
	constructCurvesAsConstraints(882082, 915190, aCurve);	// right
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(883547, 263397, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(884451, 910057, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(883937, 909749, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(885304, 911206, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(894021, 916842, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(240426, 256750, aCurve);	// left
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(889899, 295945, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(238188, 295255, aCurve);
	curvesConstraints[constCounter++] = aCurve;

	// nose + lips + cheeks 
	constructCurvesAsConstraints(130107, 58014, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(717267, 71571, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(81667, 741126, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(705473, 72860, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(60333, 94137, aCurve);
	curvesConstraints[constCounter++] = aCurve;

	// top (exc curl)
	constructCurvesAsConstraints(208690, 1002301, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(351026, 1003722, aCurve);
	curvesConstraints[constCounter++] = aCurve;

	// curl right
	constructCurvesAsConstraints(846933, 858428, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(207165, 858454, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(214019, 865303, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(214038, 865379, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(214115, 214467, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(865740, 865804, aCurve); ////
	curvesConstraints[constCounter++] = aCurve;
	//constructCurvesAsConstraints(865804, 214540, aCurve);
	//curvesConstraints[constCounter++] = aCurve;
	//constructCurvesAsConstraints(214540, 214913, aCurve);
	//curvesConstraints[constCounter++] = aCurve;
	//constructCurvesAsConstraints(866176, 866240, aCurve);
	//curvesConstraints[constCounter++] = aCurve;
	//constructCurvesAsConstraints(214976, 207336, aCurve);
	//curvesConstraints[constCounter++] = aCurve;

	// curl left
	constructCurvesAsConstraints(845207, 193244, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(844490, 203236, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(854500, 218977, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(870241, 218393, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	//constructCurvesAsConstraints(869657, 217933, aCurve);
	//curvesConstraints[constCounter++] = aCurve;
	//constructCurvesAsConstraints(869196, 204084, aCurve);
	//curvesConstraints[constCounter++] = aCurve;
	//constructCurvesAsConstraints(855348, 855628, aCurve);
	//curvesConstraints[constCounter++] = aCurve;

	// stroke above the eyes
	constructCurvesAsConstraints(793784, 145642, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(795980, 127924, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(796342, 796670, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(787811, 786977, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(788966, 974301, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(322926, 1173159, aCurve);
	curvesConstraints[constCounter++] = aCurve;

	// eye brows
	constructCurvesAsConstraints(676102, 47368, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(700788, 33018, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(710473, 12457, aCurve);
	curvesConstraints[constCounter++] = aCurve;
	constructCurvesAsConstraints(665321, 40900, aCurve);
	curvesConstraints[constCounter++] = aCurve;

}

void VectorFields::constructSoftConstraints()
{
	constructCurves_antiqueHead();

	if (false)
	{

		const int NUM_CURVES = 8;
		curvesConstraints.resize(NUM_CURVES);

		srand(time(NULL));
		int init_, end_;
		vector<int> aCurve;

		int constCounter = 0;


		////* Manual set-up for Chinese Dragon */
		///// Face
		///constructCurvesAsConstraints(152474, 51474, aCurve);
		///curvesConstraints[constCounter++] = aCurve;
		///
		///// Back
		///constructCurvesAsConstraints(44109, 68907, aCurve);
		///curvesConstraints[constCounter++] = aCurve;
		///
		///// body - bottom
		///constructCurvesAsConstraints(13471, 195817, aCurve);
		///curvesConstraints[constCounter++] = aCurve;
		///
		///// body - right
		///constructCurvesAsConstraints(123036, 247143, aCurve);
		///curvesConstraints[constCounter++] = aCurve;
		///
		///// body - left
		///constructCurvesAsConstraints(234815, 232296, aCurve);
		///curvesConstraints[constCounter++] = aCurve;
		///
		///// front_right_leg
		///constructCurvesAsConstraints(75468, 7716, aCurve);
		///curvesConstraints[constCounter++] = aCurve;
		///
		///// front_left_leg
		///constructCurvesAsConstraints(231495, 77171, aCurve);
		///curvesConstraints[constCounter++] = aCurve;
		///
		///// tail
		///constructCurvesAsConstraints(230301, 113500, aCurve);
		///curvesConstraints[constCounter++] = aCurve;

		/* Manual set-up for Armadillo */
		// Head
		constructCurvesAsConstraints(68818, 6278, aCurve);
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
	}

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
	C.resize(2 * numConstraints, B2DAsym.cols());

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

void VectorFields::measureSoftConstraintError(const Eigen::Vector3d& lambda)
{
	setupGlobalProblem(lambda);
	setAndSolveUserSystem(lambda);

	Eigen::VectorXd diff = (Xf - XFullDim);
	double error = diff.transpose()*MF2D*diff;
	double ref = Xf.transpose()*MF2D*Xf;
	double relError = sqrt(error / ref);
	cout << "The l2-norm error is " << relError << endl; 
	cout << "__The rel error is: " << error / ref << endl; 
	cout << "__The apprxo length is: " << XFullDim.transpose()*MF2D*XFullDim << endl;
	cout << "__The difference length is: " << error << endl;
	cout << "__The reference length is: " << ref << endl;

	error = diff.transpose()*SF2DAsym*diff;
	ref = Xf.transpose()*SF2DAsym*Xf;
	relError = sqrt(error / ref);
	cout << "The rel. energy error is " << relError << endl;
	cout << "__The rel energy is: " << error / ref << endl;
	cout << "__The apprxo energy is: " << XFullDim.transpose()*SF2DAsym*XFullDim << endl;
	cout << "__The difference energy is: " << error << endl;
	cout << "__The reference energy is: " << ref << endl;
	cout << "__The apprxo energy is: " << XFullDim.transpose()*SF2DAsym*XFullDim << endl;
}

void VectorFields::precomputeForSoftConstraints(const Eigen::Vector3d& lambda)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2, t3;
	chrono::duration<double>					duration;
	cout << "> Precompute for Soft-constraints... \n";
	t1 = chrono::high_resolution_clock::now();

	t2 = chrono::high_resolution_clock::now();
	SBbar = Basis.transpose()* (lambda(0)*SF2DAsym + lambda(1)*B2DAsym) * Basis;
	t3 = chrono::high_resolution_clock::now();
	duration = t3 - t2;
	cout << "....UT * (a*S + b*B) * U in " << duration.count() << " seconds" << endl;

	t2 = chrono::high_resolution_clock::now();
	CTCbar = Basis.transpose()* (C.transpose()*C) * Basis;
	t3 = chrono::high_resolution_clock::now();
	duration = t3 - t2;
	cout << "....UT * CT * C * U in " << duration.count() << " seconds" << endl;

	t2 = chrono::high_resolution_clock::now();
	CTcbar = Basis.transpose()* C.transpose()*c;
	t3 = chrono::high_resolution_clock::now();
	duration = t3 - t2;
	cout << "....UT * CT * c in " << duration.count() << " seconds" << endl;

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}

void VectorFields::projectCurvesToFrame()
{
	cout << "Projecting the constraints...";
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

	cout << "... is done " << endl; 
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
	//Eigen::SparseMatrix<double>		tempB2D = B2D;
	//B2D = B2DAsym;

	// lambda 0: dirichlet
	// lambda 1: bilaplacian
	// lambda 2: (soft-) constraint	

	//constructConstraints();
	//setupRHSGlobalProblemMapped(g, h, vEst, b);
	//setupLHSGlobalProblemMapped(A_LHS);
	//solveGlobalSystemMappedLDLT(vEst, A_LHS, b);

	arbField2D = Xf; 
	//solveGlobalSystemMappedLU_GPU();

	setupRHSGlobalProblemSoftConstraints(lambda, b);
	setupLHSGlobalProblemSoftConstraints(lambda, A_LHS);		
	solveGlobalSystemMappedLDLTSoftConstraints(A_LHS, b);

	//B2D = tempB2D;
}
void VectorFields::setupGlobalProblem(const Eigen::Vector3d& lambda, Eigen::MatrixXd& M)
{
	setupGlobalProblem(lambda);
	M.col(testID) = Xf;
}

void VectorFields::setupRHSGlobalProblemMapped(Eigen::VectorXd& g, Eigen::VectorXd& h, Eigen::VectorXd& vEst, Eigen::VectorXd& b)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Constructing RHS... ";

	vEst.resize(B2DAsym.cols());
	for (int i = 0; i < vEst.rows(); i++) {
		vEst(i) = 0.5;
	}

	g = B2DAsym * vEst;
	b.resize(B2DAsym.rows() + c.rows(), c.cols());

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

	A_LHS.resize(B2DAsym.rows() + C.rows(), B2DAsym.cols() + C.rows());

	vector<Eigen::Triplet<double>>	ATriplet;
	ATriplet.reserve(10 * B2DAsym.rows());		// It should be #rows x 4 blocks @ 2 elements (8) + #constraints,
											// but made it 10 for safety + simplicity

	for (int k = 0; k < B2DAsym.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(B2DAsym, k); it; ++it) {
			ATriplet.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
		}
	}

	for (int k = 0; k < C.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(C, k); it; ++it) {
			ATriplet.push_back(Eigen::Triplet<double>(B2DAsym.rows() + it.row(), it.col(), it.value()));
			ATriplet.push_back(Eigen::Triplet<double>(it.col(), B2DAsym.cols() + it.row(), it.value()));
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

	//Eigen::SparseMatrix<double> Mconst = C*MF2D*C.transpose();
	//b = lambda(2) * C.transpose() * Mconst * c;
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

	//const double lambda_1 = 10000 / B2D.coeff(0, 0);

	//Eigen::SparseMatrix<double> Mconst = C*MF2D*C.transpose();	
	//A_LHS = lambda(0)*SF2DAsym + lambda(1)*B2D + lambda(2)*C.transpose()*C;
	//A_LHS = lambda(0)*SF2DAsym + lambda(2)*C.transpose()*Mconst*C;
	//A_LHS = lambda(0)*SF2DAsym + lambda(2)*C.transpose()*C;
	A_LHS = lambda(0)*SF2DAsym + lambda(1)*B2DAsym + lambda(2)*C.transpose()*C;

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}


void VectorFields::solveGlobalSystemMappedLDLT(Eigen::VectorXd& vEst, Eigen::SparseMatrix<double>& A_LHS, Eigen::VectorXd& b)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2, t3;
	chrono::duration<double>					duration;
	cout << "> Solving the global system (Pardiso LDLT)... \n";
	t1 = chrono::high_resolution_clock::now();

	//cout << "Starting to solve problem." << endl;
	Xf.resize(B2DAsym.rows());
	
	// Setting up the solver
	//Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> sparseSolver(A_LHS);
	Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> sparseSolver;			
	//Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>> sparseSolver(A_LHS);
	//Eigen::PastixLDLT<Eigen::SparseMatrix<double>,1> sparseSolver(A_LHS);

	sparseSolver.analyzePattern(A_LHS);
	sparseSolver.factorize(A_LHS);
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "....Factorization in " << duration.count() << " seconds" << endl;


	// FIRST BASIS
	cout << "....Solvingthe problem..." << endl;
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
	Xf = -x.block(0, 0, B2DAsym.rows(), 1) + vEst;

	t2 = chrono::high_resolution_clock::now();
	printf("____Xf size is %dx%d\n", Xf.rows(), Xf.cols());	

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
	Xf.resize(B2DAsym.rows());

	Eigen::MatrixXd X;
	solveLUinCUDA(A_LHS, b, X);


	Xf.col(0) = -X.block(0, 0, B2DAsym.rows(), 1) + vEst;
	Xf.col(1) = -X.block(0, 1, B2DAsym.rows(), 1) + vEst;
	//cout << Xf.block(0, 0, 100, 2) << endl; 

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "..in " << duration.count() << " seconds" << endl;
}

void VectorFields::solveGlobalSystemMappedLDLTSoftConstraints(Eigen::SparseMatrix<double>& A_LHS, Eigen::VectorXd& b)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2, t3;
	chrono::duration<double>					duration;
	
	cout << "> Solving the global system (Pardiso LDLT)... \n";
	t1 = chrono::high_resolution_clock::now();
	//cout << "Starting to solve problem." << endl;
	Xf.resize(B2DAsym.rows());


	// Setting up the solver
	//Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> sparseSolver(A_LHS);
	Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> sparseSolver;
	//Eigen::PastixLDLT<Eigen::SparseMatrix<double>,1> sparseSolver(A_LHS);
	sparseSolver.analyzePattern(A_LHS);
	sparseSolver.factorize(A_LHS);
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "....Factorizing the LHS in " << duration.count() << " seconds" << endl;


	// FIRST BASIS
	//cout << "....Solving first problem (first frame)..." << endl;
	t2 = chrono::high_resolution_clock::now();
	Xf = sparseSolver.solve(b);
	t3 = chrono::high_resolution_clock::now();

	duration = t3 - t2;
	cout << "....Solving the LHS in " << duration.count()*1000 << " m.seconds" << endl;


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

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;

	printf("____Xf size is %dx%d\n", Xf.rows(), Xf.cols());
}

// Alignment fields (Maximal curvature direction) 
void VectorFields::computeMaximalPrincipalCurvature(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::VectorXd &PD, Eigen::VectorXd& PV)
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
	double energy1 = v_in.transpose() * ((B2DAsym) * v_in);
	double energy2 = v_out.transpose() * ((B2DAsym) * v_out);
	printf("The energy is=%.4f ==> %.4f.\n", energy1, energy2);
}

void VectorFields::selectAdaptiveRegions(igl::opengl::glfw::Viewer &viewer)
{
	cout << "Dealing with adaptivity\n";
	
	//const int face_id1 = 5213; 	const int face_id2 = 44893;	// arma43k
	//const int face_id1 = 26806; 	const int face_id2 = 29748;	// arma43k
	//const int face_id1 = 421007; 	const int face_id2 = 397117;	// arma43k
	//const int face_id1 = 0; const int face_id2 = 5000;
	const int face_id1 = 520374; const int face_id2 = 137423;			//bimba 1m
	Eigen::VectorXd dist;
	faceScale.resize(F.rows());	
	fieldScale.resize(F.rows());
	dist.resize(F.rows());

	Eigen::VectorXd faceColor(F.rows());
	Eigen::MatrixXd fCol(F.rows(),3);

	cout << "COmputing the dijkstra distance \n";
	computeDijkstraDistanceFace(face_id1, dist);

	double upBound = 1.25*(FC.row(face_id1) - FC.row(face_id2)).norm();
	double scaleBound = 1.5*upBound;
	double coef = (-1.0) / (scaleBound - upBound);

	for(int i=0; i<F.rows(); i++)
	{
		//printf("ID=%d distance=%.4f (c.t. %.4f) \n", i, dist(i), upBound);
		if (dist(i) < upBound) {
			faceScale(i) = 6.0;
			faceColor(i) = 0.7;
			fieldScale(i) = 1.0;
			fCol.row(i) = Eigen::RowVector3d(232.0 / 255.0, 232.0 / 255.0, 232.0 / 255.0);
		}
		else {
			faceScale(i) = 1.0;
			faceColor(i) = 0.3;
			fieldScale(i) = max(0.0, coef*(dist(i) - upBound)+ 1.0);
			//fCol.row(i) = Eigen::RowVector3d(152.0 / 255.0, 152.0 / 255.0, 152.0 / 255.0);
			fCol.row(i) = Eigen::RowVector3d(112.0 / 255.0, 161.0 / 255.0, 215.0 / 255.0);
		}
	
	}

	//Eigen::MatrixXd fCol;
	//igl::jet(faceColor, false, fCol);
	viewer.data().set_colors(fCol);
}

void VectorFields::selectAdaptiveRegions_Curvature(igl::opengl::glfw::Viewer &viewer)
{
	cout << "Dealing with adaptivity\n";

	//const int face_id1 = 5213; 	const int face_id2 = 44893;	// arma43k
	const int face_id1 = 26806; 	const int face_id2 = 29748;	// arma43k
	Eigen::VectorXd dist;
	faceScale.resize(F.rows());
	dist.resize(F.rows());

	Eigen::VectorXd faceColor(F.rows());

	/* Computing the distance */
	cout << "COmputing the dijkstra distance \n";
	computeDijkstraDistanceFace(face_id1, dist);
	double upBound = 1.25*(FC.row(face_id1) - FC.row(face_id2)).norm();

	/* Computing curvature */
	Eigen::VectorXd kV, kF;
	Eigen::SparseMatrix<double> M_, Minv_;
	igl::gaussian_curvature(V, F, kV);
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M_);
	igl::invert_diag(M_, Minv_);
	kV = (Minv_*kV).eval();
	kF.resize(F.rows());

	for (int i = 0; i < F.rows(); i++)
	{
		kF(i) = (kV(F(i, 0)) + kV(F(i, 1)) + kV(F(i, 2))) / 3.0;
	}

	int counterN = 0;
	for (int i = 0; i<F.rows(); i++)
	{
		
		//if (dist(i) < upBound) {
			double val = exp(abs(kF(i))/25.0);
			// substract by 0.5 to make it from 0.0 - 0.5
			// mult by 2 to make it from 0.0 to 1.0
			// mult by 2.5 to scale toward desired result in this non-uniform sampling
			// add by 1 to scale everything from 1.0 to 3.5
			val = ((val / (val + 1.0))-0.5)*2*2.5 + 1.0;
			//printf("ID=%d | curv=%.3f | sigm: %.4f \n", i, kF(i), val);
			faceScale(i) = val;
			faceColor(i) = 0.7;
			counterN++;
		//}
		//else {
		//	faceScale(i) = 1.0;
		//	faceColor(i) = 0.3;
		//}

	}

	cout << "There are " << counterN << " entries in the domain\n";

	cout << "(curv Vert) Min: " << kV.minCoeff() << " | max: " << kV.maxCoeff() << endl;
	cout << "(curv Face) Min: " << kF.minCoeff() << " | max: " << kF.maxCoeff() << endl;
	cout << "(scale Face) Min: " << faceScale.minCoeff() << " | max: " << faceScale.maxCoeff() << endl;

	Eigen::MatrixXd fCol;
	igl::jet(faceColor, false, fCol);
	//igl::jet(kV, true, fCol);
	viewer.data().set_colors(fCol);
}

void VectorFields::constructSamples(const int &n)
{
	numSample = n; 
	cout << "> Constructing " << n << " samples in ";
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;

	t1 = chrono::high_resolution_clock::now();
	farthestPointSampling();
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;

	cout << duration.count() << "seconds" << endl;	

	///testViennaCL2(SF2DAsym, MF2Dinv, eigFieldFull2D, eigValuesFull);

	/* Counting samples inside and outside the regions */
	int counterB = 0, counterNB = 0;
	for (int sample : Sample)
	{
		if (faceScale(sample) > 1.01)
		{
			counterB++;
		}
		else {
			counterNB++;
		}
	}
	printf("From %d samples, %d are inside the region and %d are outside \n", Sample.size(), counterB, counterNB);
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
	//Sample[0] = 69298; // Arma 43k
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

void VectorFields::constructMultiBasis()
{
	cout << "\n========================= REDUCED/LOCAL-PROBLEM =============================\n";
	vector<int> sampleSizeVect{5000, 1000, 500};
	//vector<int> sampleSizeVect{250, 500, 1000, 2500, 5000, 10000, 25000};
	numSupport = 40.0;
	for (int sample : sampleSizeVect)
	{	
		constructSamples(sample);		
		constructBasis();
		///loadAndConstructBasis();		
		string filename_basis = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Neptune_" + to_string(2 * sample) + "_Eigfields_" + to_string((int)numSupport) + "sup_spectra";
		//storeBasis(filename_basis);		
	}
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

	//double	coef = sqrt(pow(1.7, 2) + pow(1.9, 2));			// regular
	//double	coef = 1.5*sqrt(pow(1.7, 2) + pow(1.9, 2));			// adaptive
	//double	coef = sqrt(pow(1.1, 2) + pow(1.3, 2));
	//double distRatio = coef * sqrt((double)V.rows() / (double)Sample.size());

	double avg_area = 0;
	for (int i = 0; i < F.rows(); i++) avg_area += doubleArea(i) / 2.0;
	avg_area /= (double)F.rows();

	this->numSupport = 40.0;
	bool adaptiveBasis = false;			// IMPORTANT FLAG!!!!!
	if (adaptiveBasis) this->numSupport *= 4.0; 

	double distRatio = sqrt((this->numSupport*(double)F.rows()*avg_area) / (M_PI*2.0*(double)Sample.size()));

	// Setup sizes of each element to construct basis
	try {
		Basis.resize(1, 1);
		Basis.data().clear();
		Basis.data().squeeze();
		Basis.resize(2 * F.rows(), 2 * Sample.size());
	}
	catch (string &msg) {
		cout << "Cannot allocate memory for basis.." << endl;
	}

	//Basis.resize(BasisTemp.rows(), BasisTemp.cols());
	vector<vector<Eigen::Triplet<double>>> UiTriplet(Sample.size());
	cout << "Transforming the dirichlet enregy \n";
	Eigen::SparseMatrix<double> Sh = MF2DhNeg*SF2DAsym*MF2DhNeg;
	printf("Sh=%dx%d (%d) | Mh=%dx%d (%d nnzs) \n", Sh.rows(), Sh.cols(), Sh.nonZeros(), MF2DhNeg.rows(), MF2DhNeg.cols(), MF2DhNeg.nonZeros());

	/* Set-UP Laplace Matrix */
	//Eigen::SparseMatrix<double> LapForBasis = MF2Dinv * SF2DAsym;

	///cout << "....Constructing and solving local systems...";
	const int NUM_PROCESS = 4;
	durations.resize(NUM_PROCESS);
	const int Num_fields = 2;
	const int NUM_EIG = 2;

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

	///omp_set_num_threads(1);

	cout << "Setup the timing parameters: \n";
	//const int NUM_THREADS = omp_get_num_procs();
	const int NUM_THREADS = omp_get_num_procs()/2;
	//const int NUM_THREADS = 8;
	//const int NUM_THREADS = 1;
	vector<chrono::high_resolution_clock::time_point> t_end(NUM_THREADS);
	vector<chrono::duration<double>> eigen_dur(NUM_THREADS);
	vector<chrono::duration<double>> subdom_dur(NUM_THREADS);
	vector<chrono::duration<double>> boundary_dur(NUM_THREADS);
	vector<chrono::duration<double>> localEls_dur(NUM_THREADS);
	cout << "Setup the timing parameters: DONE! \n";

	duration = t0 - t0;
	for (int i = 0; i < NUM_THREADS; i++)
	{
		t_end[i] = chrono::high_resolution_clock::now();
		eigen_dur[i]    = duration;
		subdom_dur[i]   = duration;
		boundary_dur[i] = duration;
		localEls_dur[i] = duration;
	}
	
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
		
		//std::vector<Engine*> ep;
		//ep.resize(ntids);
		for (int i = 0; i < ntids; i++) 
		{
			//ep[i] = engOpen(NULL);
			//void *vpDcom = NULL;
			//int iret;
			//ep[tid] = engOpenSingleUse(NULL, vpDcom, &iret);
		}

		chrono::duration<double> dur_;

		printf("num threads=%d, iproc=%d, ID=%d, start=%d, to end=%d, num els=%d\n", ntids, iproc, tid, istart, istart + ipts, ipts);

		Eigen::VectorXd				D(F.rows());
		vector<bool>				visitedFaces(F.rows());
		for (int i = 0; i < F.rows(); i++) {
			D(i) = numeric_limits<double>::infinity();
			visitedFaces[i] = false;
		}

		//cout << "[" << tid << "] Number of processors " << iproc << ", with " << ntids << " threads." << endl;

		//UiTriplet[tid].reserve(2.0 * ((double)ipts / (double)Sample.size()) * 2 * 10.0 * F.rows());

		// Computing the values of each element
		for (id = istart; id < (istart + ipts); id++) {
		//for (id = istart; id < (istart + ipts) && id < 10; id++) {
			if (id >= Sample.size()) break;

			//vector<Eigen::Triplet<double>> BTriplet, C1Triplet, C2Triplet;
			//UiTriplet[id].reserve(2.0 * ((double)ipts / (double)Sample.size()) * 2 * 10.0 * F.rows());


			LocalFields localField(id);
			t1 = chrono::high_resolution_clock::now();
			//localField.constructSubdomain(Sample[id], V, F, avgEdgeLength, AdjMF3N, distRatio);
			localField.constructSubdomain(Sample[id], V, F, avgEdgeLength, faceScale, AdjMF2Ring, distRatio);
			//localField.constructSubdomain(Sample[id], V, F, D, AdjMF2Ring, Sample.size(), this->numSupport, NUM_EIG);
			t2 = chrono::high_resolution_clock::now();
			dur_ = t2 - t1;
			durations[0] += t2 - t1;
			subdom_dur[tid] += dur_;

			t1 = chrono::high_resolution_clock::now();
			localField.constructBoundary(F, visitedFaces, AdjMF3N, AdjMF2Ring);
			//localField.constructSelectorMatrix(F, doubleArea);
			t2 = chrono::high_resolution_clock::now();
			dur_ = t2 - t1;
			durations[1] += t2 - t1;
			boundary_dur[tid] += dur_;

			t1 = chrono::high_resolution_clock::now();
			localField.constructLocalElements(Num_fields, F);
			t2 = chrono::high_resolution_clock::now();
			dur_ = t2 - t1;
			durations[2] += t2 - t1;	
			localEls_dur[tid] += dur_;

			UiTriplet[id].reserve(2 * localField.InnerElements.size());

			t1 = chrono::high_resolution_clock::now();
			//ep[tid] = engOpen(NULL);
			if (id % 50 == 0)
			printf("Running element %d on thread %d. \n", id, tid);
			///if(id%((int)(Sample.size()/4))==0)
			///	cout << "[" << id << "] Constructing local eigen problem\n ";

			localField.constructLocalEigenProblemWithSelector(Num_fields, Sh, MF2DhNeg, AdjMF2Ring, 2, doubleArea, UiTriplet[id]);			
			//localField.constructLocalEigenProblemWithSelectorMatrix(Num_fields, SF2DAsym, MF2D, AdjMF2Ring, 2, doubleArea, UiTriplet[id]);
			//localField.constructLocalEigenProblemWithSelectorRotEig(ep[tid], tid, SF2DAsym, MF2D, AdjMF2Ring, 2, doubleArea, UiTriplet[id]);		// 2nd basis: 90 rotation of the first basis
			//engClose(ep[tid]);
			//localField.constructLocalEigenProblem(SF2D, AdjMF3N, doubleArea, UiTriplet[id]);
			t2 = chrono::high_resolution_clock::now();
			dur_ = t2 - t1;
			durations[3] += t2 - t1;
			eigen_dur[tid] += dur_;


			//if (id == 0)
			//{
			//	SubDomain = localField.SubDomain;
			//	Boundary = localField.Boundary;
			//	//patchDijkstraDist = localField.dijksFaceDistMapped;
			//}

			// To get local elements for visualizing subdomain
			if (id == 0 || id == 46) {
				cout << "Getting element of ID " << id << endl;
			
				for (int fid : localField.SubDomain) {
					localSystem(fid) = 0.3;
				}
			
				//for (int fid : localField.Boundary) {
				//	localSystem(fid) = 0.7;
				//}			
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

			if (id == istart + ipts - 1)
			{
				t_end[tid] = chrono::high_resolution_clock::now();
				//cout << "Thread " << tid << " is ALL DONE!!!\n";
			}

		}

		//for (int i = 0; i < ntids; i++)
		//{
		//	engClose(ep[i]);
		//}
	}
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;

	cout << "....Gathering local elements as basis matrix... ";
	t1 = chrono::high_resolution_clock::now();

	///checkBasisSupport(UiTriplet);

	bool writeBasisCompsToFile = false;
	if(writeBasisCompsToFile)
		writeBasisElementsToFile(UiTriplet, 2);
	else 		
		gatherBasisElements(UiTriplet, 2);

	//Basis = BasisTemp;
	//normalizeBasisAbs(2);

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
	printf("....NNZ=%d, per row = %.4f\n", Basis.nonZeros(),  (double)Basis.nonZeros() / (double)Basis.rows());

	//cout << "TIMING on each process for BASIS. \n";
	//chrono::duration<double> tot1 = t0-t0;
	//for (int i = 0; i < NUM_THREADS; i++)
	//{
	//	chrono::duration<double> d = t_end[i] - t0;
	//	cout << "[" << i << "] process. Subdomain: " << subdom_dur[i].count() << " \t | Boundary: " << boundary_dur[i].count()
	//		 << " \t | Local entries: " << localEls_dur[i].count() << " | Eigen: " << eigen_dur[i].count() << " | tend: " << d.count() << endl; 
	//
	//	tot1 += subdom_dur[i];
	//}
	//cout << "Total duration: " << tot1.count() << endl; 
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
			localField.constructSubdomain(Sample[id], V, F, avgEdgeLength, faceScale, AdjMF2Ring, distRatio);
			t2 = chrono::high_resolution_clock::now();
			durations[0] += t2 - t1;

			t1 = chrono::high_resolution_clock::now();
			localField.constructBoundary(F, AdjMF3N, AdjMF2Ring);
			t2 = chrono::high_resolution_clock::now();
			durations[1] += t2 - t1;

			t1 = chrono::high_resolution_clock::now();
			localField.constructLocalElements(2, F);
			t2 = chrono::high_resolution_clock::now();
			durations[2] += t2 - t1;

			t1 = chrono::high_resolution_clock::now();
			//localField.constructMatrixBLocal(B2D);
			localField.constructMatrixBLocal(B2DAsym, AdjMF2Ring);
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

			//localField.computeDijkstraFaceDistance(V, F, FC, AdjMF3N);

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
	string filename = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/MATLAB Implementation/Data/Genus2_Basis_1000_full.txt";
	//string filename = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/MATLAB Implementation/Data/CDragon_Basis_1000_full.mat";
	//string filename = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/MATLAB Implementation/Data/Cube_Sharp_Basis_1000.txt";
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
		string filename = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/Arma_2_Ref_eigfields";
		computeEigenMatlab(SF2DAsym, MF2D, 2, eigFieldFull2D, eigValuesFull, "hello");
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

void VectorFields::constructBasis_Coarsening(igl::opengl::glfw::Viewer &viewer)
{
	string meshFile = "D:/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/AIM_Kitten-watertight/366_kitten_500.obj";
	Eigen::MatrixXd V2;
	Eigen::MatrixXi F2;

	/* Loading small scale mesh */
	cout << "Reading the 2nd mesh "; 
	if (meshFile.substr(meshFile.find_last_of(".") + 1) == "off") {
		igl::readOFF(meshFile, V2, F2);
	}
	else if (meshFile.substr(meshFile.find_last_of(".") + 1) == "obj") {
		igl::readOBJ(meshFile, V2, F2);
	}
	else {
		cout << "Error! File type can be either .OFF or .OBJ only." << endl;
		Sleep(2000); exit(10);
	}
	printf(" with %d faces and %d vertices \n", F2.rows(), V2.rows());

	/* Loading the mesh */
	scaleMesh(V2, F2);
	//viewer.append_mesh();
	//for (int i = 0; i < V2.rows(); i++)
	//{
	//	V2.row(i) = V2.row(i) + Eigen::RowVector3d(0.25, 0.0, 0.0);
	//}
	//viewer.data().set_mesh(V2, F2);
	//viewer.selected_data_index = 1;

	cout << "Creating the basis for each local frame \n";
	/* Create the basis matrix for the 2nd mesh */
	Eigen::SparseMatrix<double> A2;
	A2.resize(3 * F2.rows(), 2 * F2.rows());
	vector<Eigen::Triplet<double>> ATriplet;
	ATriplet.reserve(3 * 2 * F2.rows());
	Eigen::Vector3d e, f, n;
	Eigen::MatrixXd NF2;
	igl::per_face_normals(V2, F2, NF2);

	for (int i = 0; i < F2.rows(); i++) {
		e = V2.row(F2(i, 1)) - V2.row(F2(i, 0));
		e.normalize();

		n = NF2.row(i);
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
	A2.setFromTriplets(ATriplet.begin(), ATriplet.end());

	/* Computing the centroid of every triangle faces */
	Eigen::MatrixXd FC2(F2.rows(), F2.cols());
	for (int i = 0; i < F2.rows(); i++)
	{
		FC2.row(i) = (V2.row(F2(i, 0)) + V2.row(F2(i, 1)) + V2.row(F2(i, 2))) / 3.0; 
	}

	/* Creating the hash table */
	cout << "Preparing the hash table\n";
	const int hashSize = 20;
	vector<vector<vector<vector<int>>>> hashTable(hashSize);	// z-size is 10
	for (int i = 0; i < hashSize; i++)	{
		hashTable[i].resize(hashSize);							// y-size is 10
		for (int j = 0; j < hashSize; j++) {
			hashTable[i][j].resize(hashSize);					// x-size is 10
			for (int k = 0; k < hashSize; k++) {
				hashTable[i][j][k].reserve(100);
			}
		}
	}
	
	cout << "Max, min, length and grid \n";
	Eigen::RowVectorXd minV(V2.rows(), 3);
	Eigen::RowVectorXd maxV(V2.rows(), 3);
	Eigen::RowVector3d length, grid;

	/* Get the min and max coefficients */
	for (int i = 0; i < V2.cols(); i++)
	{
		minV(i) = V2.col(i).minCoeff();
		maxV(i) = V2.col(i).maxCoeff();
		length(i) = maxV(i) - minV(i);
	}
	grid = length / (double)hashSize;

	// populating the hash table
	cout << "Populating hash table \n";
	vector<int> ids(3);
	for (int i = 0; i < F2.rows(); i++)	{
		for(int j=0; j<3; j++)		ids[j] = (int)floor((FC2(i, j) - minV(j)) / grid(j));
		hashTable[ids[2]][ids[1]][ids[0]].push_back(i);
	}

	// Traversing the hash table
	//for (int z_ = 0; z_ < hashSize; z_++) {
	//	for (int y_ = 0; y_ < hashSize; y_++) {
	//		for (int x_ = 0; x_ < hashSize; x_++) {
	//			if (hashTable[z_][y_][x_].size() > 0)	{
	//				printf("The idx[%d][%d][%d] has %d entries: ", z_, y_, x_, hashTable[z_][y_][x_].size());
	//				for (int k : hashTable[z_][y_][x_]) {
	//					cout << k << " | ";
	//				}
	//				cout << endl; 
	//			}
	//		}
	//	}
	//}

	/* Selecting the closest triangle on the 2nd mesh from the 1st mesh*/
	cout << "Selecting closest entries \n"; 
	vector<int> largeToSmallCorr(F.rows()); 
	vector<set<int>> smallToLargeCorr(F2.rows());
	Eigen::VectorXd corrColor(F.rows());
	for (int i = 0; i < F.rows(); i++) {
		for (int j = 0; j < 3; j++)
		{
			ids[j] = (int)floor((FC(i, j) - minV(j)) / grid(j));
			if (ids[j] < 0) ids[j] = 0;
			if (ids[j] >= hashSize) ids[j] = hashSize-1;
		}

		// if the table has a content, find the triangle (in mesh 2) minimum distance to mesh 1
		vector<int> neigh;
		if (hashTable[ids[2]][ids[1]][ids[0]].size() > 0) {
			neigh.reserve(100);
			for (int fID : hashTable[ids[2]][ids[1]][ids[0]])
			{
				neigh.push_back(fID);
			}
		}
		else {
			int m = 1;
			do {
				for (int k_ = -m; k_ <=m; k_++) {
					if (ids[2] + k_ < 0 || ids[2] + k_ >= hashSize) continue;
					for (int j_ = -m; j_ <= m; j_++) {
						if (ids[1] + j_ < 0 || ids[1] + j_ >= hashSize) continue;
						for (int i_ = -m; i_ <=m ; i_++) {
							if (ids[0] + i_ < 0 || ids[0] + i_ >= hashSize) continue;

							for (int fID : hashTable[ids[2] + k_][ids[1] + j_][ids[0] + i_])
							{
								neigh.push_back(fID);
							}
						}
					}
				}
				m++;
			} while (neigh.size() < 1);
		}

		double minDist = std::numeric_limits<double>::max();
		int minID;
		for (int fID : neigh)
		{
			double d_ = (FC.row(i) - FC2.row(fID)).norm();
			if (d_ < minDist) {
				minID = fID;
				minDist = d_;
			}
		}

		largeToSmallCorr[i] = minID;
		smallToLargeCorr[minID].insert(i);
		//corrColor(i) = (double) (50*minID % (int) F.rows());
		corrColor(i) = FC(minID,0)*FC(minID, 1)*FC(minID, 2);

		if (i == 147312 || i == 139270) {
			printf("Face %d (in mesh 1) is closest to face %d (in mesh 2) \n", i, minID);
		}

		//printf("Face %d (in mesh 1) is closest to face %d (in mesh 2) \n", i, minID);
		//printf("Face %d is located at hash[%d][%d][%d] with %d entries \n", i, ids[2], ids[1], ids[0], hashTable[ids[2]][ids[1]][ids[0]].size());		
	}

	Eigen::MatrixXd fCol;	
	igl::jet(corrColor, true, fCol);
	viewer.data().set_colors(fCol);
	viewer.selected_data_index = 0;


	/* Mapping from coarse to full resolution mesh */
	Basis.resize(2 * F.rows(), 2 * F2.rows());
	vector<Eigen::Triplet<double>> BTriplet;
	BTriplet.reserve(2 * F.rows());

	for (int i = 0; i < F2.rows(); i++)	
	{
		Eigen::MatrixXd base(3,2);
		Eigen::MatrixXd entry(2, 2);
		base = A2.block(3 * i, 2 * i, 3, 2);
		for (int face : smallToLargeCorr[i])
		{
			Eigen::MatrixXd ALoc(3, 2);
			
			ALoc = A.block(3 * face, 2 * face, 3, 2);
			entry = ALoc.transpose()*base;

			BTriplet.push_back(Eigen::Triplet<double>(2 * face + 0, 2 * i + 0, entry(0, 0)));
			BTriplet.push_back(Eigen::Triplet<double>(2 * face + 1, 2 * i + 0, entry(1, 0)));
			BTriplet.push_back(Eigen::Triplet<double>(2 * face + 0, 2 * i + 1, entry(0, 1)));
			BTriplet.push_back(Eigen::Triplet<double>(2 * face + 1, 2 * i + 1, entry(1, 1)));
		}
	}
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
	double	coef = sqrt(pow(1.5, 2) + pow(1.7, 2));
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
	const int Num_fields = 2;

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

	Eigen::SparseMatrix<double> Sh = MF2DhNeg*SF2DAsym*MF2DhNeg;

	omp_set_num_threads(16);
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

		//UiTriplet[tid].reserve(2.0 * ((double)ipts / (double)Sample.size()) * EIG_NUM * 10.0 * F.rows());

		// Computing the values of each element
		for (id = istart; id < (istart + ipts); id++) {
		//for (id = istart; id < (istart + ipts) && id < 10; id++) {
			if (id >= Sample.size()) break;

			vector<Eigen::Triplet<double>> BTriplet, C1Triplet, C2Triplet;

			LocalFields localField(id);
			t1 = chrono::high_resolution_clock::now();
			//localField.constructSubdomain(Sample[id], V, F, avgEdgeLength, AdjMF3N, distRatio);
			//localField.constructSubdomain(Sample[id], V, F, avgEdgeLength, faceScale, AdjMF2Ring, distRatio);
			localField.constructSubdomain(Sample[id], V, F, D, AdjMF2Ring, Sample.size(), this->numSupport, EIG_NUM);
			t2 = chrono::high_resolution_clock::now();
			durations[0] += t2 - t1;

			t1 = chrono::high_resolution_clock::now();
			//localField.constructBoundary(F, AdjMF3N, AdjMF2Ring);
			localField.constructSelectorMatrix(F, doubleArea);
			t2 = chrono::high_resolution_clock::now();
			durations[1] += t2 - t1;

			t1 = chrono::high_resolution_clock::now();
			localField.constructLocalElements(Num_fields, F);
			t2 = chrono::high_resolution_clock::now();
			durations[2] += t2 - t1;

			printf("ID=%d has %d entries \n", id, localField.InnerElements.size());
			UiTriplet[id].reserve(2 * localField.InnerElements.size());

			t1 = chrono::high_resolution_clock::now();
			//localField.constructLocalEigenProblemWithSelector(Num_fields, Sh, MF2DhNeg, AdjMF2Ring, EIG_NUM, doubleArea, UiTriplet[id]);
			localField.constructLocalEigenProblemWithSelectorMatrix(Num_fields, SF2DAsym, MF2D, AdjMF2Ring, EIG_NUM, doubleArea, UiTriplet[id]);
			//localField.constructLocalEigenProblemWithSelector(SF2DAsym, MF2D, AdjMF2Ring, EIG_NUM, doubleArea, UiTriplet[id]);
			//localField.constructLocalEigenProblem(SF2D, AdjMF3N, doubleArea, UiTriplet[id]);
			t2 = chrono::high_resolution_clock::now();
			durations[3] += t2 - t1;


			//if (id == 0)
			//{
			//	SubDomain = localField.SubDomain;
			//	Boundary = localField.Boundary;
			//	//patchDijkstraDist = localField.dijksFaceDistMapped;
			//}
			//
			//// To get local elements for visualizing subdomain
			//if (id == 0 || id == 46) {
			//	cout << "Getting element of ID " << id << endl;
			//
			//	for (int fid : localField.SubDomain) {
			//		localSystem(fid) = 0.3;
			//	}
			//
			//	for (int fid : localField.Boundary) {
			//		localSystem(fid) = 0.7;
			//	}			
			//
			//}
		}
	}
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;

	cout << "....Gathering local elements as basis matrix... ";
	t1 = chrono::high_resolution_clock::now();
	gatherBasisElements(UiTriplet, EIG_NUM);
	printf("Basis 1st row nnz: %d \n", Basis.col(0).nonZeros());
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
	printf("> Basis Structure information (%d nnz)\n", Basis.nonZeros());
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
			localField.constructLocalElements(2, F);
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
		//std::copy(UiTriplet[j].begin(), UiTriplet[j].end(), BTriplet.begin() + tripSize);
		std::move(UiTriplet[j].begin(), UiTriplet[j].end(), BTriplet.begin() + tripSize);
	}	
	//BTriplet.reserve(totalElements);
	//for (int i = 0; i < Sample.size(); i++)
	//{
	//	for (int j = 0; j < UiTriplet[i].size(); j++)
	//	{
	//		BTriplet.push_back(UiTriplet[i][j]);
	//	}
	//	printf("Triplet %d has %d nnz \n", i, UiTriplet[i].size());
	//}
	Basis.setFromTriplets(BTriplet.begin(), BTriplet.end());
	printf("Basis space: %dx%d | BTriplet size: %d \n", Basis.rows(), Basis.cols(), BTriplet.size());
	vector<int> checkPoint = { 1, 190000, 400000, 600000 };
	//for(int i : checkPoint)
	//printf("%d data: [%d,  %d]=%.4f \n", i,  BTriplet[i].row(), BTriplet[i].col(), BTriplet[i].value());

	BTriplet.clear(); BTriplet.shrink_to_fit();
	//for (int k = 0; k < Basis.outerSize(); ++k) {
	//int k = 0;
	//	for (Eigen::SparseMatrix<double>::InnerIterator it(Basis, k); it; ++it) {
	//		//if(it.row()<1000) 
	//		printf("[%d,%d]=%.3f \n", it.row(), it.col(), it.value());
	//		//ATriplet.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
	//	}
	//}


	// empty the data
	//BTriplet.clear();
	//BTriplet.shrink_to_fit();

	/* Computing the basis sum -> cor normalization */
	//for (int k = 0; k < BasisTemp.outerSize(); ++k) {
	//	for (Eigen::SparseMatrix<double>::InnerIterator it(BasisTemp, k); it; ++it) {
	//		BasisSum(it.row(), it.col() % NUM_EIGEN) += it.value();
	//	}
	//}
}
void VectorFields::checkBasisSupport(const vector<vector<Eigen::Triplet<double>>> &UiTriplet)
{

	cout << "Checking the support \n";
	Eigen::VectorXd sumPerRow; sumPerRow.resize(2 * F.rows()); sumPerRow.setZero();
	Eigen::VectorXi countPerRow; countPerRow.resize(2 * F.rows()); countPerRow.setZero();

	for (int j = 0; j < (Sample.size()-50); j++) {
		//printf("Triplet %d has %d entries \n", j, UiTriplet[j].size());
		{
			for (int i = 0; i < UiTriplet[j].size(); i++)
			{
				sumPerRow[UiTriplet[j][i].row()] += abs(UiTriplet[j][i].value());
				countPerRow[UiTriplet[j][i].row()] += 1;
			}
		}
	}

	printf("Size of coutn per row: %d with last entry: %d \n", countPerRow.size(), countPerRow(2 * F.rows() - 1));

	for (int i = 0; i < countPerRow.size(); i++)
	{
		//if (i % 10 == 0)
		//{
		//	printf("The %d-th row has %d entries (total val: %.4f) \n", i, countPerRow(i), sumPerRow(i));
		//}
		
		if (countPerRow(i) < 25)
		{
			//printf("PROBLEM! The %d-th row has only %d entries \n", i, countPerRow(i));
			lowSupportFaces.push_back(i);
		}
	}

	/* Getting the IDs of the faces with low support */
	for (int i = 0; i < lowSupportFaces.size(); i++)
	{
		printf("[%d] ID:%d has only %d support. Problem! \n ", i, lowSupportFaces[i], countPerRow(lowSupportFaces[i]));
	}

	/* Finding the distance between them */
	Eigen::VectorXd faceDistance(F.rows());
	faceDistance.setConstant(std::numeric_limits<double>::max());
	computeDijkstraDistanceFaceForSampling(lowSupportFaces[0], faceDistance);

	/* Showing distance among faces */
	for (int i = 0; i < lowSupportFaces.size(); i++)
	{
		printf("Distance from Face %d to Face %d is %.32f \n", lowSupportFaces[0], lowSupportFaces[i], faceDistance(lowSupportFaces[i]));
	}

}

void VectorFields::writeBasisElementsToFile(const vector<vector<Eigen::Triplet<double>>> &UiTriplet, const int& NUM_EIGEN)
{
	string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/BasisComponents_Fertility_"+ to_string(2* Sample.size()) + "_EigFields_160sup.txt"; 
	
	std::ofstream ofs;
	ofs.open(resultFile, std::ofstream::out | std::ofstream::app);

	for (int j = 0; j < Sample.size(); j++) {
		for (Eigen::Triplet<double> trip : UiTriplet[j])
		{
			ofs << trip.row() << "," << trip.col() << "," << trip.value() << "\n";
		}
	}	
	
	ofs.close();
	cout << "Basis written to file \n";
}

void VectorFields::loadAndConstructBasis()
{
	string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/BasisComponents_Fertility_" + to_string(2 * Sample.size()) + "_EigFields_160sup.txt";
	ifstream file(resultFile);
	string oneLine, oneWord;
	int i = 0, j;
	double v;
	int counter;

	vector<Eigen::Triplet<double>> MTriplet;
	MTriplet.reserve(2*Sample.size() * 200);

	if (file.is_open())
	{
		Basis.resize(2*F.rows(), 2*Sample.size());

		/* Obtain each member elements */
		while (getline(file, oneLine))
		{			
			istringstream iStream(oneLine);

			getline(iStream, oneWord, ',');
			int i = stoi(oneWord);
			getline(iStream, oneWord, ',');
			int j = stoi(oneWord);
			getline(iStream, oneWord, ',');
			double v = stod(oneWord);
			MTriplet.push_back(Eigen::Triplet<double>(i, j, v));			
		}
	}
	file.close();
	printf(", with %d elements\n", MTriplet.size());
	Basis.setFromTriplets(MTriplet.begin(), MTriplet.end());
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
	/* For testing only, will be cleared later */
	BasisTemp = Basis;
	Basis.resize(0, 0);
	Basis.resize(BasisTemp.rows(), BasisTemp.cols());


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
			Eigen::Vector2d ab; 
			ab << a, b; 
			//double norm = sqrt(a*a + b*b);
			double norm = ab.norm();
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
	///getUserConstraints();
	//getUserConstraintsEfficient();

	///setupRHSUserProblemMapped(gBar, hBar, vEstBar, bBar);
	///setupLHSUserProblemMapped(A_LHSBar);
	///solveUserSystemMappedLDLT(vEstBar, A_LHSBar, bBar);
	
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

	B2DBar = Basis.transpose() * B2DAsym * Basis; 

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

	//constructConstraints();

	printf("Basis: %dx%d | C: %dx%d | c:%d \n", Basis.rows(), Basis.cols(), C.rows(), C.cols(), c.rows());

	//userConstraints = globalConstraints; 
	CBar			= C * Basis;
	cBar			= c;

	printf("Basis: %dx%d | C: %dx%d | CBar: %dx%d | c:%d | cBar:%d \n", Basis.rows(), Basis.cols(), C.rows(), C.cols(), CBar.rows(), CBar.cols(), c.rows(), cBar.rows());


	/* MUCH FASTER: Alternative of CBar construction */
	//vector<Eigen::Triplet<double>> CTriplet;
	//CTriplet.reserve(40 * 2*globalConstraints.size());
	//vector<double> constraints_(2 * globalConstraints.size());
	//for (int i = 0; i < globalConstraints.size(); i++) { 
	//	constraints_[2 * i]   = 2*globalConstraints[i]; 
	//	constraints_[2 * i+1] = 2*globalConstraints[i]+1;
	//}
	////for(int k=0; k<Basis.transpose().outerSize(); ++k)
	//for(int k=0; k<constraints_.size(); k++)
	//{
	//	for (Eigen::SparseMatrix<double>::InnerIterator it(BasisT, constraints_[k]); it; ++it)
	//	{
	//		CTriplet.push_back(Eigen::Triplet<double>(k, it.row(), it.value()));
	//	}
	//}
	//CBar.resize(0, 0);
	//CBar.resize(2 * globalConstraints.size(), Basis.cols());
	//CBar.setFromTriplets(CTriplet.begin(), CTriplet.end());
	CBarT = CBar.transpose();

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	//cout << "in " << duration.count() << " seconds." << endl;

	//printf(".... C_LoCal = %dx%d\n", CBar.rows(), CBar.cols());
	//printf(".... c_LoCal = %dx%d\n", cBar.rows(), cBar.cols());
}

void VectorFields::getUserConstraintsEfficient()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Obtaining user constraints ";

	cBar = c;
	Eigen::SparseMatrix<double> CT = C.transpose();

	//printf("Basis: %dx%d | C: %dx%d | CBar: %dx%d | c:%d | cBar:%d \n", Basis.rows(), Basis.cols(), C.rows(), C.cols(), CBar.rows(), CBar.cols(), c.rows(), cBar.rows());
	
	/* MUCH FASTER: Alternative of CBar construction */
	/* For the stroke hard constraints */
	vector<Eigen::Triplet<double>> CTriplet;
	CTriplet.reserve(40 * 2*globalConstraints.size());
	
	{
		//for (int i = 2 * globalConstraints.size(); i < C.rows(); i++)
		for (int i = 0; i < C.rows(); i++)
		{
			vector<int> filledIdx;
			//printf("i:%d :", i - 2 * globalConstraints.size());
			for (Eigen::SparseMatrix<double>::InnerIterator it(CT, i); it; ++it)
			{
				filledIdx.push_back(it.row());
				//printf("%d |", it.row());
			}
			//cout << endl;

			for (int k = 0; k < Basis.cols(); k++)
			{
				double e_ik = 0;
				//e_ik = C.coeff(i, filledIdx[0])*Basis.coeff(filledIdx[0], k) + C.coeff(i, filledIdx[1])*Basis.coeff(filledIdx[1], k) + C.coeff(i, filledIdx[2])*Basis.coeff(filledIdx[2], k);
				for (int val : filledIdx) {
					e_ik += C.coeff(i, val)*Basis.coeff(val, k);
				}
				if (abs(e_ik) > 100.0*std::numeric_limits<double>::epsilon()) {
					CTriplet.push_back(Eigen::Triplet<double>(i, k, e_ik));
					//printf("data [%d,%d]=%.4f \n", i, k, e_ik);
				}
			}
		}
	}

	CBar.resize(0, 0);
	CBar.resize(C.rows(), Basis.cols());
	CBar.setFromTriplets(CTriplet.begin(), CTriplet.end());
	CBarT = CBar.transpose();

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	//cout << "in " << duration.count() << " seconds." << endl;

	printf(".... C_LoCal = %dx%d\n", CBar.rows(), CBar.cols());
	printf(".... c_LoCal = %dx%d\n", cBar.rows(), cBar.cols());

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

	printf("Rhs: %d \n", bBar.rows());
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

	/* Mass matrix of the selected faces (on constraints) */
	//Eigen::SparseMatrix<double> Mconst = C*MF2D*C.transpose();

	//printf("Siz of Mconst: %dx%d\n", Mconst.rows(), Mconst.cols());
	//printf("Siz of Cbar: %dx%d\n", CBar.rows(), CBar.cols());

	//bBar = lambda(2)*CBar.transpose() * cBar; 
	//bBar = lambda(2)*Basis.transpose()*(C.transpose()*c);
	//bBar = lambda(2)*CBar.transpose() * Mconst * cBar;
	bBar = lambda(2)*CTcbar;

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

	///Eigen::SparseMatrix<double> SF2DBar = Basis.transpose() * SF2DAsym * Basis;
	//Eigen::SparseMatrix<double> B2DBar = Basis.transpose() * B2D * Basis;
	//A_LHSBar = SF2DBar + lambda*CBar.transpose()*CBar; 


	//const double lambda_1 = 10000 / B2D.coeff(0, 0);
	//cout << "lambda_2 " << lambda(2) << endl;

	/* Local matrix */
	//Eigen::SparseMatrix<double> Mconst = C*MF2D*C.transpose();

	//printf("Siz of Mconst: %dx%d\n", Mconst.rows(), Mconst.cols());
	//printf("Siz of Cbar: %dx%d\n", CBar.rows(), CBar.cols());

	//A_LHSBar = lambda(0)*SF2DBar +  lambda(1)*B2DBar + lambda(2)*CBar.transpose()*CBar;
	//A_LHSBar = lambda(0)*SF2DBar + lambda(2)*CBar.transpose()*Mconst*CBar;
	//A_LHSBar = lambda(0)*SF2DBar + lambda(2)*CBar.transpose()*CBar;
	///Eigen::SparseMatrix<double> CTCbar;
	///CTCbar = Basis.transpose()* (C.transpose()*C)*Basis;
	//A_LHSBar = lambda(0)*SF2DBar + lambda(1)*B2DBar + lambda(2)*CBar.transpose()*CBar;
	///A_LHSBar = lambda(0)*SF2DBar + lambda(1)*B2DBar + lambda(2)*CTCbar;
	A_LHSBar = SBbar + lambda(2)*CTCbar;

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
	//Eigen::SparseLU<Eigen::SparseMatrix<double>> sparseSolver;

	//sparseSolver.analyzePattern(A_LHSBar);
	//sparseSolver.compute(A_LHSBar);

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
	cout << "> Solving reduced system...\n";
	t0 = chrono::high_resolution_clock::now();

	//XLowDim.resize(B2DBar.rows());
	//XFullDim.resize(Basis.rows());
	//Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> sparseSolver(A_LHSBar);

	Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> sparseSolver;
	sparseSolver.analyzePattern(A_LHSBar);
	sparseSolver.factorize(A_LHSBar);
	t1 = chrono::high_resolution_clock::now(); 
	duration = t1 - t0;
	cout << "....Re-factorization in " << duration.count()*1000 << " mili.secs." << endl;

	//cout << "....Solving for the first frame.\n";
	t1 = chrono::high_resolution_clock::now();
	Eigen::VectorXd x = sparseSolver.solve(bBar);
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "....Solving in red. space in " << duration.count() * 1000 << " mili.secs." << endl;

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
	cout << "___Mapping to full-resolution...";
	t0 = chrono::high_resolution_clock::now();

	//XFullDim = Basis * XLowDim;
	//SparseMatrix_Vector_Multiplication(Basis, XLowDim, XFullDim);
	//Eigen::SparseMatrix<double, Eigen::RowMajor> BasisRow = Basis; 
	//SparseMatrix_Vector_Multiplication_CSR(BasisRow, XLowDim, XFullDim);
	performLifting();

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << " in " << duration.count() << " seconds." << endl;

	//cout << "Fields \n";
	//cout << XFullDim.block(0, 0, 100, 1) << endl; 


	//printf("....XFull (%dx%d) =  Basis (%dx%d) * XLowDim (%dx%d) \n", XFullDim.rows(), XFullDim.cols(), Basis.rows(), Basis.cols(), XLowDim.rows(), XLowDim.cols());
}

void VectorFields::initializeParametersForLifting()
{
	cudaError_t			cudaStat1 = cudaSuccess;

	// initialize the system
	cusparseCreateMatDescr(&descrA);
	cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO);
	cusparseCreate(&handle);

	// Matrix variables
	BasisRow = Basis; 
	BasisRow.makeCompressed();
	const int nnz = BasisRow.nonZeros();
	const int m = BasisRow.rows();
	const int n = BasisRow.cols();

	// Populating the matrix in CPU
	double* h_csrVal	= (double*)malloc(nnz * sizeof(double));
	int* h_csrRowPtr	= (int*)malloc((m + 1) * sizeof(int));
	int* h_csrColInd	= (int*)malloc(nnz * sizeof(int));
	h_csrVal			= BasisRow.valuePtr();
	h_csrRowPtr			= BasisRow.outerIndexPtr();
	h_csrColInd			= BasisRow.innerIndexPtr();	

	// Allocating memory in device/GPU
	cudaStat1 = cudaMalloc(&d_csrColInd, nnz * sizeof(int));     //cout << "__col:alloc_status:" << cudaStat1 << endl;
	cudaStat1 = cudaMalloc(&d_csrRowPtr, (m + 1) * sizeof(int)); //cout << "__row:alloc_status:" << cudaStat1 << endl;
	cudaStat1 = cudaMalloc(&d_csrVal, nnz * sizeof(double));     //cout << "__val:alloc_status:" << cudaStat1 << endl;

	cudaStat1 = cudaMemcpy(d_csrRowPtr, h_csrRowPtr, (m + 1) * sizeof(int), cudaMemcpyHostToDevice);  //cout << "__rows: status:" << cudaStat1 << endl;
	cudaStat1 = cudaMemcpy(d_csrColInd, h_csrColInd, nnz * sizeof(int), cudaMemcpyHostToDevice);      //cout << "__col: status:" << cudaStat1 << endl;
	cudaStat1 = cudaMemcpy(d_csrVal, h_csrVal, nnz * sizeof(double), cudaMemcpyHostToDevice);         //cout << "__val: status:" << cudaStat1 << endl;

}

void VectorFields::performLifting()
{
	// Setting up some variable
	cudaError_t			cudaStat1 = cudaSuccess;
	const int nnz = BasisRow.nonZeros();
	const int m = BasisRow.rows();
	const int n = BasisRow.cols();

	// Populating data in CPU
	//double* h_a = (double*)malloc(n * sizeof(double));
	//h_a = XLowDim.data();
	double* h_b = (double*)malloc(m * sizeof(double));
	for (int i = 0; i < m; i++) h_b[i] = 0.5;

	// Allocating memory in device/GPU
	double *d_a;  cudaStat1 = cudaMalloc(&d_a, n * sizeof(double));				 //cout << "__alloc_status:" << cudaStat1 << endl;
	double *d_b;  cudaStat1 = cudaMalloc(&d_b, m * sizeof(double));				 //cout << "__alloc_status:" << cudaStat1 << endl;

	// Copying data to CUDA/GPU
	cudaStat1 = cudaMemcpy(d_a, XLowDim.data(), n * sizeof(double), cudaMemcpyHostToDevice);					//cout << "__alloc_status:" << cudaStat1 << endl;
	cudaStat1 = cudaMemcpy(d_b, h_b, m * sizeof(double), cudaMemcpyHostToDevice);					//cout << "__alloc_status:" << cudaStat1 << endl;

	// The multiciplication
	double alpha = 1.0;
	double beta = 0.0;
 	cusparseStatus_t cusparseStat1 = cusparseDcsrmv(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, m, n, nnz, &alpha, descrA, d_csrVal, d_csrRowPtr, d_csrColInd, d_a, &beta, d_b);
	//cout << "__status:" << cusparseStat1 << endl;

	// Copying to CPU
	cudaMemcpy(h_b, d_b, m * sizeof(double), cudaMemcpyDeviceToHost);

	//XFullDim.resize(2 * F.rows());
	XFullDim = Eigen::Map<Eigen::VectorXd>(h_b, m);
}


void VectorFields::obtainUserVectorFields()
{

	cout << "Hello there \n" << endl; 
}

// INTERACTIVE/REAL-TIME SYSTEM VIA SCHUR COMPLEMENT
void VectorFields::setAndSolveInteractiveSystem(const Eigen::Vector3d& lambda)
{
	// Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	cout << "Set constraint, set system, solve system, and lift it up in :";
	t0 = chrono::high_resolution_clock::now();

	obtainConstraints();
	solveInteractiveSystem();

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << " in " << duration.count() << " seconds." << endl;
}

void VectorFields::obtainConstraints()
{
	//getUserConstraints();
	getUserConstraintsEfficient();
}

void VectorFields::preComputeReducedElements()
{
	// Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "Precompute B2D, vEst, and B2D*vEst ...";


	// Factorization of B2DBar;
	B2DBarFactor.analyzePattern(B2DBar);
	B2DBarFactor.factorize(B2DBar);

	// Set up the estimation variable vAdd
	vAdd.resize(B2DBar.rows());
	for (int i = 0; i < vAdd.rows(); i++) {
		vAdd(i) = 0.5;
	}

	// setup the Bv
	BvBar = B2DBar * vAdd;

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << " in " << duration.count() << " seconds." << endl;
}

void VectorFields::solveInteractiveSystem()
{
	// Timing
	///chrono::high_resolution_clock::time_point	t0, t1, t2;
	///chrono::duration<double>					duration;
	///cout << "Solve interactive system ...\n";
	///
	////* ================== 1. Setting up LHS ================== */
	///cout << "__Create LHS: ";
	///t0 = chrono::high_resolution_clock::now();
	///t1 = chrono::high_resolution_clock::now();
	//Eigen::MatrixXd BC(CBar.cols(), CBar.rows());
	if (CBar.rows() == 2)
		BC.resize(CBar.cols(), CBar.rows());
	else
		BC.conservativeResize(Eigen::NoChange, BC.cols() + deltaConstraints);
		//BC.conservativeResize(Eigen::NoChange, BC.cols() + 2);
	Eigen::MatrixXd LHS;
	Eigen::VectorXd bc;
	//for (int i = 2; i > 0; i--)
	for (int i = deltaConstraints; i > 0; i--)
	{
		bc = CBarT.col(BC.cols()-i);
		//bc.transposeInPlace();
		BC.col(BC.cols()-i) = B2DBarFactor.solve(bc);
	}
	LHS = CBar*BC;
	///t2 = chrono::high_resolution_clock::now();
	///duration = t2 - t1;
	///cout << " in " << duration.count() << " seconds." << endl;
	///
	////* ================== 1. Setting up RHS ================== */
	///cout << "__Create rhs: ";
	///t1 = chrono::high_resolution_clock::now();
	Eigen::VectorXd bbv = B2DBarFactor.solve(BvBar);
	Eigen::VectorXd cbbv = CBar*bbv;
	Eigen::VectorXd cvc = CBar*vAdd - cBar; 
	Eigen::VectorXd rhs = cbbv - cvc;
	///t2 = chrono::high_resolution_clock::now();
	///duration = t2 - t1;
	///cout << " in " << duration.count() << " seconds." << endl;
	///
	////* ================== 1. Solve the 1st System  ================== */
	///cout << "__Solve Schur-complement system: ";
	///t1 = chrono::high_resolution_clock::now();
	Eigen::LDLT<Eigen::MatrixXd> LHS_Fact;
	LHS_Fact.compute(LHS);
	Eigen::VectorXd lambda = LHS_Fact.solve(rhs);
	///t2 = chrono::high_resolution_clock::now();
	///duration = t2 - t1;
	///cout << " in " << duration.count() << " seconds." << endl;
	////* ================== 2. Setting up LHS ================== */
	///
	///
	////* ================== 2. Setting up RHS ================== */
	///cout << "__dense rhs: ";
	///t1 = chrono::high_resolution_clock::now();
	rhs = CBar.transpose()*lambda - BvBar;
	///t2 = chrono::high_resolution_clock::now();
	///duration = t2 - t1;
	///cout << " in " << duration.count() << " seconds." << endl;

	/* ================== 2. Solve the 2nd System  ================== */
	///cout << "__solve 2nd system: ";
	///t1 = chrono::high_resolution_clock::now();
	Eigen::VectorXd x = B2DBarFactor.solve(rhs);
	///t2 = chrono::high_resolution_clock::now();
	///duration = t2 - t1;
	///cout << " in " << duration.count() << " seconds." << endl;

	/* ================== 3. Map x to xStar  ================== */
	///cout << "__map x to x*: ";
	///t1 = chrono::high_resolution_clock::now();
	XLowDim = x + vAdd;
	///t2 = chrono::high_resolution_clock::now();
	///duration = t2 - t1;
	///cout << " in " << duration.count() << " seconds." << endl;
	///
	///t2 = chrono::high_resolution_clock::now();
	///duration = t2 - t0;
	///cout << " in " << duration.count() << " seconds." << endl;

	/* ================== 4. Map to full resolution  ================== */
	mapSolutionToFullRes();
}

void VectorFields::computeEigenFields(const int &numEigs, const string& filename)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Computing reference eigenproblem... ";

	//computeEigenMatlab(SF2DAsym, MF2D, eigFieldFull2D, eigValuesFull);
	//computeEigenMatlab(SF2DAsym, MF2D, numEigs, eigFieldFull2D, eigValuesFull, filename);
	//computeEigenSpectra_GenSym(SF2DAsym, MF2D, numEigs, eigFieldFull2D, eigValuesFull, filename);
	//computeEigenSpectra_RegNSym(SF2DAsym, MF2Dinv, numEigs, eigFieldFull2D, eigValuesFull, filename);
	Eigen::SparseMatrix<double> Sh = MF2DhNeg*SF2DAsym*MF2DhNeg;
	//computeEigenSpectra_RegSym_Transf(Sh, MF2DhNeg, numEigs, eigFieldFull2D, eigValuesFull, filename);
	//computeEigenSpectra_RegSym_Custom(Sh, MF2DhNeg, numEigs, eigFieldFull2D, eigValuesFull, filename);
	computeEigenMatlab(SF2DAsym, MF2D, numEigs, eigFieldFull2D, eigValuesFull, "hello");
	//cout << "::::: Eigen Values (Full Res) \n" << eigValuesFull << endl;
	//WriteSparseMatrixToMatlab(MF2D, "hello");

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;
}

void VectorFields::retrieveEigenFields(const string& filename)
{
	ReadDenseMatrixFromMatlab(eigFieldFull2D, filename);
	ReadVectorFromMatlab(eigValuesFull, filename, 2*F.rows());
}

void VectorFields::computeApproxEigenFields(const int &numEigs, const string& filename)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;

	bool obtainEigenBasis = false;
	if (obtainEigenBasis)
	{
		Eigen::MatrixXd EigenBasis;
		string filenameBasis = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/Armadillo_1000_eigenfields_Ref_2";
		ReadDenseMatrixFromMatlab(EigenBasis, filenameBasis, 172964, 1000);
		//EigenBasis = EigenBasis.block(0, 0, 172964, 500);
		Eigen::MatrixXd MbarD = EigenBasis.transpose() * MF2D * EigenBasis;
		Eigen::MatrixXd SFAsymbarD = EigenBasis.transpose() * SF2DAsym * EigenBasis;

		vector<Eigen::Triplet<double>> MTriplet;
		vector<Eigen::Triplet<double>> STriplet;
		MTriplet.reserve(MbarD.rows()*MbarD.cols());
		STriplet.reserve(SFAsymbarD.rows()*SFAsymbarD.cols());
		for (int i = 0; i < MbarD.cols(); i++) {
			for (int j = 0; j < MbarD.rows(); j++) {
				MTriplet.push_back(Eigen::Triplet<double>(j, i, MbarD(j, i)));
				STriplet.push_back(Eigen::Triplet<double>(j, i, SFAsymbarD(j, i)));
			}
		}
		Eigen::SparseMatrix<double> Mbar;
		Eigen::SparseMatrix<double> SFAsymbar;
		Mbar.resize(MbarD.rows(), MbarD.cols());
		SFAsymbar.resize(SFAsymbarD.rows(), SFAsymbarD.cols());
		Mbar.setFromTriplets(MTriplet.begin(), MTriplet.end());
		SFAsymbar.setFromTriplets(STriplet.begin(), STriplet.end());
		computeEigenMatlab(SFAsymbar, Mbar, numEigs, eigFieldReduced2D, eigValuesReduced, filename);
	}
	else {		
		cout << "> Computing restricted eigenproblem (in Matlab)...\n ";

		cout << "____Computing reduced system... ";
		t1 = chrono::high_resolution_clock::now();
		Eigen::SparseMatrix<double> Mbar = Basis.transpose() * MF2D * Basis;
		Eigen::SparseMatrix<double> SFAsymbar = Basis.transpose() * SF2DAsym * Basis;
		t2 = chrono::high_resolution_clock::now();
		duration = t2 - t1;
		cout << "in " << duration.count() << " seconds" << endl;

		cout << "____Computing restricted eigenproblem ... ";
		t1 = chrono::high_resolution_clock::now();
		computeEigenGPU(SFAsymbar, Mbar, eigFieldReduced2D, eigValuesReduced);
		//computeEigenMatlab(Sbar, Mbar, eigFieldReduced2D, eigValuesReduced);
		//computeEigenMatlab(SFAsymbar, Mbar, numEigs, eigFieldReduced2D, eigValuesReduced, filename);
		//cout << "::::: Eigen Values (Reduced) \n" << eigValuesReduced << endl;

		t2 = chrono::high_resolution_clock::now();
		duration = t2 - t1;
		cout << "in " << duration.count() << " seconds" << endl;
		//WriteSparseMatrixToMatlab(Basis, "hello");
	}

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "____Computing the eigenproblem in " << duration.count() << " seconds" << endl;
}

void VectorFields::computeApproxEigenFields_Mult()
{
	vector<string> modelFile;
	vector<string> basisFile;

	// populating the model File
	modelFile.push_back("../LocalFields/Models/Armadillo/Armadillo_43243.obj");
	modelFile.push_back("../LocalFields/Models/AIM894_Chinese Dragon/894_dragon_tris.obj");
	modelFile.push_back("../LocalFields/Models/AIM_fertility_watertight/fertility.obj");
	modelFile.push_back("../LocalFields/Models/AIM_Ramesses_clean_watertight/814_Ramesses_1.5Mtriangles_clean.off");
	modelFile.push_back("D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/models/AIM_Neptune_clean__watertight_4M triangles/803_neptune_4Mtriangles_manifold.off");

	// populating the model File
	basisFile.push_back("D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Arma_2000_EigFields_35sup");
	basisFile.push_back("D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_CDragon_2000_Eigfields_40sup");
	basisFile.push_back("D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Fertility_2000_Eigfields_40sup");
	basisFile.push_back("D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Ramses_2000_Eigfields_40sup_spectra");
	basisFile.push_back("D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Neptune_2000_Eigfields_40sup_spectra");
	for (int i = 0; i < modelFile.size(); i++)
	{
		readMesh(modelFile[i]);
		scaleMesh(V, F);

		computeEdges();
		computeAverageEdgeLength();
		computeFaceCenter();
		computeFaceNormal();
		constructVFNeighbors();

		constructVertexAdjacencyMatrix();
		constructFaceAdjacency3NMatrix();
		constructFaceAdjacency2RingMatrix();
		constructEVList();
		constructEFList();
		selectFaceToDraw(2000);

		/* MATRIX CONSTRUCTIONS */
		constructMassMatrices();
		constructRotationMatrix();
		constructMappingMatrix();

		constructGradient3D();
		constructGradientStar3D();
		constructStiffnessMatrices_Implicit();

		/* Construct sample */
		numSample = 1000;
		numSupport = 40.0;
		faceScale.resize(F.rows()); 
		faceScale.setConstant(1.0);
		constructSamples(numSample);

		/* Retrive basis*/
		if(i<3)
			constructBasis();
		else 
			retrieveBasis(basisFile[i]);

		const int eigsToCompute = 500;
		computeApproxEigenFields(eigsToCompute, "hello");

		/* Lifting to full res */
		chrono::high_resolution_clock::time_point	t1, t2;
		chrono::duration<double>					duration;
		cout << "Lifting to the full resolution...";
		t1 = chrono::high_resolution_clock::now();

		Eigen::MatrixXd Eig_Full = Basis*eigFieldReduced2D.block(0,0, Basis.cols(), eigsToCompute);

		t2 = chrono::high_resolution_clock::now();
		duration = t2 - t1;
		cout << "in " << duration.count() << " seconds" << endl;
	}
}

void VectorFields::retrieveApproxEigenFields() 
{
	//ReadSparseMatrixFromMatlab(Basis, "Hello");
	ReadDenseMatrixFromMatlab(eigFieldReduced2D, "Hello");
	ReadVectorFromMatlab(eigValuesReduced, "hello", 2*F.rows());
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
	double refHarmEnergy = Xf.transpose()       * SF2DAsym * Xf;
	double appHarmEnergy = XFullDim.transpose() * SF2DAsym * XFullDim;
	double refBiHarmEnergy = Xf.transpose()       * B2DAsym * Xf;
	double appBiHarmEnergy = XFullDim.transpose() * B2DAsym * XFullDim;

	cout << ">> [REF] Energy: Harm=" << refHarmEnergy << ", Biharmonic=" << refBiHarmEnergy << endl;
	cout << ">> [APP] Energy: Harm=" << appHarmEnergy << ", Biharmonic=" << appBiHarmEnergy << endl;
	cout << "         Relative Harm-energy =" << abs(refHarmEnergy - appHarmEnergy) / refHarmEnergy << endl;
	cout << "         Relative Biharm-energy =" << abs(refBiHarmEnergy - appBiHarmEnergy) / refBiHarmEnergy << endl;


	//printf("MSE = %.3f \n", (diffV.sum()/diffV.size()) / xf);
}

void VectorFields::measureDirichletEnergy()
{
	double dirichlet = Xf.transpose() * ((SF2DAsym) * Xf);
	cout << "__Dirichlet Energy\n \t__FullRes: " << dirichlet; 
	dirichlet = XFullDim.transpose() * SF2DAsym * XFullDim; 
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
		double refHarmEnergy   = Xf.transpose()       * SF2DAsym * Xf; 
		double appHarmEnergy   = XFullDim.transpose() * SF2DAsym * XFullDim;
		double refBiHarmEnergy = Xf.transpose()       * B2DAsym * Xf;
		double appBiHarmEnergy = XFullDim.transpose() * B2DAsym * XFullDim;

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
	double energy1 = v_in.transpose() * ((B2DAsym * MF2D) * v_in);
	double energy2 = v_out.transpose() * ((B2DAsym * MF2D) * v_out);
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

void VectorFields::scaleMesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F)
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

