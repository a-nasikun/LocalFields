#include "VectorFields.h"
#include <fstream>
#include <random>

/* ====================== ITEMS FOR TESTING ONLY ============================*/
void VectorFields::constructArbitraryField()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Constructing arbitrary field...";

	// Construct the Field
	arbField.resize(V.rows());
		
	// Dijstra-based Arbitrary Field
	//computeDijkstraDistanceVertex(0, arbField);

	// Position based arbitrary scalar field
	for (int i = 0; i < V.rows(); i++) {
		arbField(i) = V(i, 0) *  V(i, 1) *  V(i, 2);
		//arbField(i) = V(i, 1);
	}


	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;
}

void VectorFields::constructArbitraryField2D()
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Constructing gradient 2D of arbitrary field...";


	Eigen::SparseMatrix<double> Grad3D, Grad2D;
	igl::grad(V, F, Grad3D);
	printf("Size of grad3d: %dx%d\n", Grad3D.rows(), Grad3D.cols());
	rearrangeGradient3D(Grad3D);
	printf("Size of grad3d [2]: %dx%d\n", Grad3D.rows(), Grad3D.cols());
	Grad2D = A.transpose()*Grad3D; 
	printf("Grad3d [2]: %dx%d || arbFields: %d. \n", Grad2D.rows(), Grad2D.cols(), arbField.rows());
	arbField2D = 10.0*Grad2D * arbField;

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;
}
void VectorFields::testProjection_MyBasis_NoRegularizer(const Eigen::SparseMatrix<double>& U, const Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> &sparseSolver, const Eigen::SparseMatrix<double>& B, const Eigen::VectorXd& a, const Eigen::VectorXd& v, double &error)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> [Testing basis 2D of arbitrary field without regularizer...]\n";

	// Construct matrices for Test
	//cout << "____Assigning variables\n";
	//Eigen::SparseMatrix<double> U = Basis;// BasisTemp;
	//Eigen::VectorXd				v = inputFields; 
	//Eigen::VectorXd				a = (U.transpose()*(MF2D*v));
	//Eigen::SparseMatrix<double> B = U.transpose() * MF2D * U;
	//
	//cout << "____Solving linear system variables\n";
	//Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> sparseSolver(B);

	Eigen::VectorXd w = sparseSolver.solve(a);

	if (sparseSolver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system. " << endl;
		if (sparseSolver.info() == Eigen::NumericalIssue)
			cout << "NUMERICAL ISSUE. " << endl;
		cout << sparseSolver.info() << endl;
		return;
	}

	//Eigen::VectorXd wb;
	//wb.resize(U.rows());

	//cout << "____Getting total SUM(wi*bi) \n";
	wb = U*w; 

	// Compare their L2-Norm
	cout << "____Computing L2-norm \n";
	Eigen::VectorXd diff = v - wb; 
	double length1 = wb.transpose()*MF2D*wb;
	double norm1 = diff.transpose()*MF2D*diff;
	double norm2 = v.transpose()*MF2D*v;
	double normL2 = sqrt(norm1 / norm2); 
	error = normL2; 

	//cout << "____The L2 Norm is << " << normL2 << ": sqrt(" << norm1 << "/" << norm2 << ")" << endl;

	/* Measuring the 'length' of each vector */
	cout << "ENERGY ASYM \n";
	double harm_energy1 = v.transpose()*SF2DAsym*v;
	double harm_energy2 = wb.transpose()*SF2DAsym*wb;
	double harm_relEnergy = abs(harm_energy1 - harm_energy2) / harm_energy1;

	//cout << "____Harmonic Energy => Ref=" << harm_energy1 << ", Approx:" << harm_energy2 << endl;
	//cout << "____Relative energy: " << harm_relEnergy << endl;

	/* Measuring the 'length' of each vector */
	//cout << "ENERGY SYM \n";
	double harm_energy1_sym = v.transpose()*SF2D*v;
	double harm_energy2_sym = wb.transpose()*SF2D*wb;
	double harm_relEnergy_sym = abs(harm_energy1_sym - harm_energy2_sym) / harm_energy1_sym;
	//cout << "____Harmonic Energy => Ref=" << harm_energy1_sym << ", Approx:" << harm_energy2_sym << endl;
	//cout << "____Relative energy: " << harm_relEnergy_sym << endl;
	
	/* Measuring the energy */	
	double biharm_energy1 = v.transpose()*B2DAsym*v;
	double biharm_energy2 = wb.transpose()*B2DAsym*wb;
	double biharm_relEnergy = abs(biharm_energy1 - biharm_energy2) / biharm_energy1;
	//cout << "____Bi-Harmonic ENERGY => Ref=" << biharm_energy1 << ", Approx:" << biharm_energy2 << endl;
	//cout << "____aRelative energy: " << biharm_relEnergy << endl;

	/* Measuring the 'length' of each vector */
	//cout << "ENERGY SYM \n";
	double biharm_energy1_sym = v.transpose()*B2D*v;
	double biharm_energy2_sym = wb.transpose()*B2D*wb;
	double biharm_relEnergy_sym = abs(biharm_energy1_sym - biharm_energy2_sym) / biharm_energy1_sym;
	//cout << "____Harmonic Energy => Ref=" << biharm_energy1_sym << ", Approx:" << biharm_energy2_sym << endl;
	//cout << "____Relative energy: " << biharm_relEnergy_sym << endl;


	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	//cout << "in " << duration.count() << " seconds." << endl;

	bool writeToFile = true;
	if (writeToFile) {
		std::ofstream ofs;
		//string resultFile = "D:/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Kitten_L2projection_eigenFields_" + to_string(Basis.cols()) + "_coarsening.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Kitten_L2projection_eigenFields_" + to_string(Basis.cols()) + ".txt";
		
		string resultFile = "D:/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Fertility_L2projection_variousBases_combined.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/CDragon_L2projection_eigenFields_40sup.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/CDragon_L2projection_optAlg.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/CDragon_L2projection_eigFields10.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/CDragon_L2projection_eigPatch.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/CDragon_L2projection_grad.txt";

		ofs.open(resultFile, std::ofstream::out | std::ofstream::app);

		ofs << globalConstraints.size() << "\t"
			<< harm_energy1 << "\t" << harm_energy1_sym << "\t" << biharm_energy1 << "\t" << biharm_energy1_sym << "\t" << norm2 << "\t"	// reference (harmAsym, harmSym, biharmAsym, biharmSym)
			<< harm_energy2 << "\t" << harm_relEnergy << "\t" << harm_energy2_sym << "\t" << harm_relEnergy_sym << "\t"						// approx (harmAsym, harmRelEnergyAsym,harmSym, harmRelEnergySym)
			<< biharm_energy2 << "\t" << biharm_relEnergy << "\t" << biharm_energy2_sym << "\t" << biharm_relEnergy_sym << "\t"				// approx (biharmAsym, biharmRelEnergyAsym, biharmSym, biharmRelEnergySym)
			<< length1 << "\t" << normL2 << "\n";																												// l2-norm

		ofs.close();
	}
}

void VectorFields::testProjection_EigenBasis_NoRegularizer(const Eigen::MatrixXd& U, const Eigen::LDLT<Eigen::MatrixXd>& denseSolver, const Eigen::VectorXd& a,  const Eigen::VectorXd& v, double &error)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> [Testing EIGENBASIS without regularizer...]\n";

	// Construct matrices for Test
	cout << "____Assigning variables\n";
	//Eigen::MatrixXd U = Basis;// BasisTemp;
	//Eigen::VectorXd				v = inputFields;
	//Eigen::VectorXd				a = (U.transpose()*(MF2D*v));
	//Eigen::MatrixXd				B = U.transpose() * MF2D * U;

	cout << "____Solving linear system variables\n";
	//Eigen::VectorXd w = B.ldlt().solve(a);
	Eigen::VectorXd w = denseSolver.solve(a);
		
	cout << "____Getting total SUM(wi*bi) \n";
	wbEigen = U*w;

	// Compare their L2-Norm
	cout << "____Computing L2-norm \n";
	Eigen::VectorXd diff = v - wbEigen;
	double length1 = wbEigen.transpose()*MF2D*wbEigen;
	double norm1 = diff.transpose()*MF2D*diff;
	double norm2 = v.transpose()*MF2D*v;
	double normL2 = sqrt(norm1 / norm2);
	error = normL2;

	cout << "____The L2 Norm is << " << normL2 << ": sqrt(" << norm1 << "/" << norm2 << ")" << endl;

	/* Measuring the 'length' of each vector */
	cout << "ENERGY ASYM \n";
	double harm_energy1 = v.transpose()*SF2DAsym*v;
	double harm_energy2 = wbEigen.transpose()*SF2DAsym*wbEigen;
	double harm_relEnergy = abs(harm_energy1 - harm_energy2) / harm_energy1;

	cout << "____Harmonic Energy => Ref=" << harm_energy1 << ", Approx:" << harm_energy2 << endl;
	cout << "____Relative energy: " << harm_relEnergy << endl;

	/* Measuring the 'length' of each vector */
	cout << "ENERGY SYM \n";
	double harm_energy1_sym = v.transpose()*SF2D*v;
	double harm_energy2_sym = wbEigen.transpose()*SF2D*wbEigen;
	double harm_relEnergy_sym = abs(harm_energy1_sym - harm_energy2_sym) / harm_energy1_sym;
	cout << "____Harmonic Energy => Ref=" << harm_energy1_sym << ", Approx:" << harm_energy2_sym << endl;
	cout << "____Relative energy: " << harm_relEnergy_sym << endl;

	/* Measuring the energy */
	double biharm_energy1 = v.transpose()*B2DAsym*v;
	double biharm_energy2 = wbEigen.transpose()*B2DAsym*wbEigen;
	double biharm_relEnergy = abs(biharm_energy1 - biharm_energy2) / biharm_energy1;
	cout << "____Bi-Harmonic ENERGY => Ref=" << biharm_energy1 << ", Approx:" << biharm_energy2 << endl;
	cout << "____aRelative energy: " << biharm_relEnergy << endl;

	/* Measuring the 'length' of each vector */
	cout << "ENERGY SYM \n";
	double biharm_energy1_sym = v.transpose()*B2D*v;
	double biharm_energy2_sym = wbEigen.transpose()*B2D*wbEigen;
	double biharm_relEnergy_sym = abs(biharm_energy1_sym - biharm_energy2_sym) / biharm_energy1_sym;
	cout << "____Harmonic Energy => Ref=" << biharm_energy1_sym << ", Approx:" << biharm_energy2_sym << endl;
	cout << "____Relative energy: " << biharm_relEnergy_sym << endl;


	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;

	bool writeToFile = true;
	if (writeToFile) {
		std::ofstream ofs;
		string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Armadillo_L2projection_modalBasis.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Fertility_L2projection_modalBasis_sameStorage.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/CDragon_L2projection_modalBasis_500.txt";

		ofs.open(resultFile, std::ofstream::out | std::ofstream::app);

		ofs << globalConstraints.size() << "\t"
			<< harm_energy1 << "\t" << harm_energy1_sym << "\t" << biharm_energy1 << "\t" << biharm_energy1_sym << "\t" << norm2 << "\t"	// reference (harmAsym, harmSym, biharmAsym, biharmSym)
			<< harm_energy2 << "\t" << harm_relEnergy << "\t" << harm_energy2_sym << "\t" << harm_relEnergy_sym << "\t"						// approx (harmAsym, harmRelEnergyAsym,harmSym, harmRelEnergySym)
			<< biharm_energy2 << "\t" << biharm_relEnergy << "\t" << biharm_energy2_sym << "\t" << biharm_relEnergy_sym << "\t"				// approx (biharmAsym, biharmRelEnergyAsym, biharmSym, biharmRelEnergySym)
			<< length1 << "\t" << normL2 << "\n";																												// l2-norm

		ofs.close();
	}
}

void VectorFields::testProjection_MyBasis_WithRegularizer(const Eigen::SparseMatrix<double>& Basis, const Eigen::VectorXd& inputFields, const Eigen::SparseMatrix<double>& MReg, double &error)
{

	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Testing basis 2D of arbitrary field WITH regularizer...\n";

	// Construct matrices for Test	
	cout << "__[APPROXIMATION]....\n";
	//const double				lambda = 10000 * MF2D.coeff(0,0) / MReg.coeff(0,0);
	
	Eigen::VectorXd id(MF2D.rows());
	id.setConstant(1.0);
	double en1 = id.transpose()*MF2D*id;
	double en2 = id.transpose()*MReg*id;

	Eigen::SparseMatrix<double> U = Basis;// BasisTemp;
	Eigen::VectorXd				v = inputFields;
	/* Perturbed the input fields*/
	//perturbVectorFields(v);
	//perturbVectorFieldsRegular(v);
	//pertFields = v;
	
	srand(time(NULL));	
	int l1_ = testID % 5;
	//double l2_ = (double)(l1_ + 1) / 10.0;
	double l2_ = 25 * (l1_+1);
	//int l1_ = testID;
	//double l2_ = 100+(double)(l1_+1)*10.0;

	double				lambda = l2_;
	cout << "__lamba: " << lambda << endl;
	printf("U=%dx%d | MF=%dx%d | MREg=%dx%d | v=%d \n", U.rows(), U.cols(), MF2D.rows(), MF2D.cols(), MReg.rows(), MReg.cols(), v.rows());
	lambda = lambda * en1 / en2; 
	//const double				lambda = 0.5;
	Eigen::VectorXd				a = U.transpose()*MF2D*v;
	Eigen::SparseMatrix<double> B = U.transpose() * (MF2D + lambda*MReg) * U;

	cout << "____Solving linear system variables (with lambda=" << lambda << ")\n";
	//Eigen::PardisoLLT<Eigen::SparseMatrix<double>> sparseSolver(B);
	Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> sparseSolver(B);
	
	Eigen::VectorXd w = sparseSolver.solve(a);

	if (sparseSolver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system. " << endl;
		if (sparseSolver.info() == Eigen::NumericalIssue)
			cout << "NUMERICAL ISSUE. " << endl;
		cout << sparseSolver.info() << endl;
		return;
	}
	wb = U*w;
	
	cout << "__[REFERENCE]....\n";
	a = MF2D*v;
	B = MF2D + lambda*MReg; 
	Eigen::VectorXd  wRef;
	//Eigen::PardisoLLT<Eigen::SparseMatrix<double>> sparseSolver2(B);
	Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> sparseSolver2(B);
	wRef = sparseSolver2.solve(a);
	projRef = wRef;

	// Compare their L2-Norm
	cout << "____Computing L2-norm \n";
	Eigen::VectorXd diff = wRef - wb;
	double norm1 = diff.transpose()*MF2D*diff;
	double norm2 = wRef.transpose()*MF2D*wRef;
	double normL2 = sqrt(norm1 / norm2);	
	cout << "The L-2 Norm is << " << normL2 << endl;	
	cout << "____The L2 Norm is << " << normL2 << ": sqrt(" << norm1 << "/" << norm2 << ")" << endl;

	/* Measuring the 'length' of each vector */
	cout << "HARMONIC \n";
	double harm_energy1 = wRef.transpose()*SF2DAsym*wRef;
	double harm_energy2 = wb.transpose()*SF2DAsym*wb;
	double harm_relEnergy = abs(harm_energy1 - harm_energy2) / harm_energy1;
	cout << "____Harmonic Energy => Ref=" << harm_energy1 << ", Approx:" << harm_energy2 << endl;
	cout << "____Relative energy: " << harm_relEnergy << endl;

	/* Measuring the 'length' of each vector */
	double harm_energy1_sym = wRef.transpose()*SF2D*wRef;
	double harm_energy2_sym = wb.transpose()*SF2D*wb;
	double harm_relEnergy_sym = abs(harm_energy1_sym - harm_energy2_sym) / harm_energy1_sym;
	cout << "____Harmonic Energy => Ref=" << harm_energy1_sym << ", Approx:" << harm_energy2_sym << endl;
	cout << "____Relative energy: " << harm_relEnergy_sym << endl;

	/* Measuring the energy */
	cout << "BIHARMONIC \n";
	double biharm_energy1 = wRef.transpose()*B2DAsym*wRef;
	double biharm_energy2 = wb.transpose()*B2DAsym*wb;
	double biharm_relEnergy = abs(biharm_energy1 - biharm_energy2) / biharm_energy1;
	cout << "____Bi-Harmonic ENERGY => Ref=" << biharm_energy1 << ", Approx:" << biharm_energy2 << endl;
	cout << "____aRelative energy: " << biharm_relEnergy << endl;

	/* Measuring the 'length' of each vector */
	cout << "ENERGY SYM \n";
	double biharm_energy1_sym = wRef.transpose()*B2D*wRef;
	double biharm_energy2_sym = wb.transpose()*B2D*wb;
	double biharm_relEnergy_sym = abs(biharm_energy1_sym - biharm_energy2_sym) / biharm_energy1_sym;
	cout << "____Harmonic Energy => Ref=" << biharm_energy1_sym << ", Approx:" << biharm_energy2_sym << endl;
	cout << "____Relative energy: " << biharm_relEnergy_sym << endl;	


	/* Measuring the 'Initial Energy' of input vector */
	cout << "Input vector: ENERGY SYM \n";
	double harm_energy_input = v.transpose()*SF2DAsym*v;
	double harm_energy_input_sym = v.transpose()*SF2D*v;
	double biharm_energy_input = v.transpose()*B2DAsym*v;
	double biharm_energy_input_sym = v.transpose()*B2D*v;
	cout << "____Input Harm Energy => aSym" << harm_energy_input << ", Sym:" << harm_energy_input_sym << endl;
	cout << "____Input BiHarm Energy => aSym" << biharm_energy_input << ", Sym:" << biharm_energy_input_sym << endl;

	/* Measuring the 'length' of each vector */
	double length_init = v.transpose()*MF2D*v;
	double length1 = wRef.transpose()*MF2D*wRef;
	double length2 = wb.transpose()*MF2D*wb;
	cout << "Length => Ref=" << length1 << ", Approx:" << length2 << ", initial: " << v.transpose()*MF2D*v << endl;
	cout << "Relative length: " << length2 / length1 * 100.0 << "%" << endl;

	/* Measuring the performance of target/goal */
	cout << "Goal of the minimization \n";
	double goal_energy_input = v.transpose() * B * v;
	double goal_energy_ref = wRef.transpose()*B*wRef;
	double goal_energy_app = wb.transpose()*B*wb;
	double goal_energy_rel = abs(goal_energy_app - goal_energy_ref) / goal_energy_ref;
	cout << "____Input goal" << goal_energy_input << endl;
	cout << "____Goal: Ref" << goal_energy_ref << ", App:" << goal_energy_app << endl;
	error = goal_energy_rel;

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;

	bool writeToFile = true;
	if (writeToFile) {
		std::ofstream ofs;
		string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/CDragon_L2projectionWithReg_eigenFields" + to_string(Basis.cols()) + ".txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Armadillo_L2projectionWithReg_eigenFields" + to_string(Basis.cols()) + ".txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Kitten_L2projectionWithReg_eigenFields_" + to_string(Basis.cols()) + "_40sup.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Kitten_L2projectionWithReg_eigenFieldsNRot_" + to_string(Basis.cols()) + "_40sup.txt";

		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Armadillo_L2projectionWithReg_eigFields.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Armadillo_L2projectionWithReg_optAlg.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Armadillo_L2projectionWithReg_eigFields10.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Armadillo_L2projectionWithReg_eigPatch.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Armadillo_L2projectionWithReg_grad.txt";

		ofs.open(resultFile, std::ofstream::out | std::ofstream::app);

		ofs << globalConstraints.size() << "\t" << lambda * en2 / en1 << "\t" << lambda << "\t"
			<< harm_energy_input << "\t" << harm_energy_input_sym << "\t" << biharm_energy_input << "\t" << biharm_energy_input_sym << "\t" << goal_energy_input << "\t" << length_init << "\t"	// initial input
			<< harm_energy1 << "\t" << harm_energy1_sym << "\t" << biharm_energy1 << "\t" << biharm_energy1_sym << "\t" << goal_energy_ref << "\t" << norm2 << "\t"	// reference (harmAsym, harmSym, biharmAsym, biharmSym)
			<< harm_energy2 << "\t" << harm_relEnergy << "\t" << harm_energy2_sym << "\t" << harm_relEnergy_sym << "\t"	<< goal_energy_app << "\t"					// approx (harmAsym, harmRelEnergyAsym,harmSym, harmRelEnergySym)
			<< biharm_energy2 << "\t" << biharm_relEnergy << "\t" << biharm_energy2_sym << "\t" << biharm_relEnergy_sym << "\t"				// approx (biharmAsym, biharmRelEnergyAsym, biharmSym, biharmRelEnergySym)
			<< length_init << "\t" << length1 << "\t" << length2 << "\t" << length2/length1 * 100.0 << "\t"
			<< goal_energy_rel << "\t" << normL2 << "\n";																												// l2-norm

		ofs.close();
	}
}

void VectorFields::testProjection_EigenBasis_WithRegularizer(const Eigen::MatrixXd& U, const Eigen::VectorXd& v, const Eigen::SparseMatrix<double>& MReg, double &error)
{

	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Testing basis 2D of arbitrary field WITH regularizer...\n";

	// Construct matrices for Test	
	cout << "__[APPROXIMATION]....\n";
	//const double				lambda = 10000 * MF2D.coeff(0,0) / MReg.coeff(0,0);
	
	Eigen::VectorXd id(MF2D.rows());
	id.setConstant(1.0);
	double en1 = id.transpose()*MF2D*id;
	double en2 = id.transpose()*MReg*id;

	//Eigen::MatrixXd				U = Basis;
	//Eigen::VectorXd				v = inputFields;
	
	/* Perturbed the input fields*/
	//perturbVectorFields(v);
	//pertFields = v;

	srand(time(NULL));
	int l1_ = testID % 5;
	double l2_ = 25 * (l1_ + 1);

	double				lambda = l2_;
	cout << "__lamba: " << lambda << endl;
	lambda = lambda * en1 / en2;

	Eigen::VectorXd				a = U.transpose()*MF2D*v;
	Eigen::MatrixXd				B = U.transpose() * (MF2D + lambda*MReg) * U;

	cout << "____Solving linear system variables (with lambda=" << lambda << ")\n";
	Eigen::VectorXd w = B.ldlt().solve(a);	
	wbEigen = U*w;

	cout << "__[REFERENCE]....\n";
	a = MF2D*v;
	Eigen::SparseMatrix<double> BRef = MF2D + lambda*MReg;
	Eigen::VectorXd  wRef;
	Eigen::PardisoLLT<Eigen::SparseMatrix<double>> sparseSolver2(BRef);
	wRef = sparseSolver2.solve(a);
	projRef = wRef;
	
	// Compare their L2-Norm
	cout << "____Computing L2-norm \n";
	Eigen::VectorXd diff = wRef - wbEigen;
	double norm1 = diff.transpose()*MF2D*diff;
	double norm2 = wRef.transpose()*MF2D*wRef;
	double normL2 = sqrt(norm1 / norm2);
	//error = normL2;
	cout << "The L-2 Norm is << " << normL2 << endl;
	cout << "____The L2 Norm is << " << normL2 << ": sqrt(" << norm1 << "/" << norm2 << ")" << endl;

	/* Measuring the 'length' of each vector */
	cout << "HARMONIC \n";
	double harm_energy1 = wRef.transpose()*SF2DAsym*wRef;
	double harm_energy2 = wbEigen.transpose()*SF2DAsym*wbEigen;
	double harm_relEnergy = abs(harm_energy1 - harm_energy2) / harm_energy1;
	cout << "____Harmonic Energy => Ref=" << harm_energy1 << ", Approx:" << harm_energy2 << endl;
	cout << "____Relative energy: " << harm_relEnergy << endl;

	/* Measuring the 'length' of each vector */
	double harm_energy1_sym = wRef.transpose()*SF2D*wRef;
	double harm_energy2_sym = wbEigen.transpose()*SF2D*wbEigen;
	double harm_relEnergy_sym = abs(harm_energy1_sym - harm_energy2_sym) / harm_energy1_sym;
	cout << "____Harmonic Energy => Ref=" << harm_energy1_sym << ", Approx:" << harm_energy2_sym << endl;
	cout << "____Relative energy: " << harm_relEnergy_sym << endl;

	/* Measuring the energy */
	cout << "BIHARMONIC \n";
	double biharm_energy1 = wRef.transpose()*B2DAsym*wRef;
	double biharm_energy2 = wbEigen.transpose()*B2DAsym*wbEigen;
	double biharm_relEnergy = abs(biharm_energy1 - biharm_energy2) / biharm_energy1;
	cout << "____Bi-Harmonic ENERGY => Ref=" << biharm_energy1 << ", Approx:" << biharm_energy2 << endl;
	cout << "____aRelative energy: " << biharm_relEnergy << endl;

	/* Measuring the 'length' of each vector */
	cout << "ENERGY SYM \n";
	double biharm_energy1_sym = wRef.transpose()*B2D*wRef;
	double biharm_energy2_sym = wbEigen.transpose()*B2D*wbEigen;
	double biharm_relEnergy_sym = abs(biharm_energy1_sym - biharm_energy2_sym) / biharm_energy1_sym;
	cout << "____Harmonic Energy => Ref=" << biharm_energy1_sym << ", Approx:" << biharm_energy2_sym << endl;
	cout << "____Relative energy: " << biharm_relEnergy_sym << endl;


	/* Measuring the 'Initial Energy' of input vector */
	cout << "Input vector: ENERGY SYM \n";
	double harm_energy_input = v.transpose()*SF2DAsym*v;
	double harm_energy_input_sym = v.transpose()*SF2D*v;
	double biharm_energy_input = v.transpose()*B2DAsym*v;
	double biharm_energy_input_sym = v.transpose()*B2D*v;
	cout << "____Input Harm Energy => aSym" << harm_energy_input << ", Sym:" << harm_energy_input_sym << endl;
	cout << "____Input BiHarm Energy => aSym" << biharm_energy_input << ", Sym:" << biharm_energy_input_sym << endl;

	/* Measuring the 'length' of each vector */
	double length_init = v.transpose()*MF2D*v;
	double length1 = wRef.transpose()*MF2D*wRef;
	double length2 = wbEigen.transpose()*MF2D*wbEigen;
	cout << "Length => Ref=" << length1 << ", Approx:" << length2 << ", initial: " << v.transpose()*MF2D*v << endl;
	cout << "Relative length: " << length2 / length1 * 100.0 << "%" << endl;

	/* Measuring the performance of target/goal */
	cout << "Goal of the minimization \n";
	double goal_energy_input = v.transpose() * BRef * v;
	double goal_energy_ref = wRef.transpose()*BRef*wRef;
	double goal_energy_app = wbEigen.transpose()*BRef*wbEigen;
	double goal_energy_rel = abs(goal_energy_app - goal_energy_ref) / goal_energy_ref;
	cout << "____Input goal" << goal_energy_input << endl;
	cout << "____Goal: Ref" << goal_energy_ref << ", App:" << goal_energy_app << endl;
	error = goal_energy_rel;

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;

	bool writeToFile = true;
	if (writeToFile) {
		std::ofstream ofs;
		string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Bimba_L2projectionWithReg_modalBasis_sameStorage.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Kitten_L2projectionWithReg_modalBasis_667.txt";

		ofs.open(resultFile, std::ofstream::out | std::ofstream::app);

		ofs << globalConstraints.size() << "\t" << lambda * en2 / en1 << "\t" << lambda << "\t"
			<< harm_energy_input << "\t" << harm_energy_input_sym << "\t" << biharm_energy_input << "\t" << biharm_energy_input_sym << "\t" << goal_energy_input << "\t" << length_init << "\t"	// initial input
			<< harm_energy1 << "\t" << harm_energy1_sym << "\t" << biharm_energy1 << "\t" << biharm_energy1_sym << "\t" << goal_energy_ref << "\t" << norm2 << "\t"	// reference (harmAsym, harmSym, biharmAsym, biharmSym)
			<< harm_energy2 << "\t" << harm_relEnergy << "\t" << harm_energy2_sym << "\t" << harm_relEnergy_sym << "\t" << goal_energy_app << "\t"					// approx (harmAsym, harmRelEnergyAsym,harmSym, harmRelEnergySym)
			<< biharm_energy2 << "\t" << biharm_relEnergy << "\t" << biharm_energy2_sym << "\t" << biharm_relEnergy_sym << "\t"				// approx (biharmAsym, biharmRelEnergyAsym, biharmSym, biharmRelEnergySym)
			<< length_init << "\t" << length1 << "\t" << length2 << "\t" << length2 / length1 * 100.0 << "\t"
			<< goal_energy_rel << "\t" << normL2 << "\n";

		//ofs << globalConstraints.size() << "\t" << lambda * en2 / en1 << "\t" << lambda << "\t"
		//	<< harm_energy_input << "\t" << harm_energy_input_sym << "\t" << biharm_energy_input << "\t" << biharm_energy_input_sym << "\t" << length_init << "\t"	// initial input
		//	<< harm_energy1 << "\t" << harm_energy1_sym << "\t" << biharm_energy1 << "\t" << biharm_energy1_sym << "\t" << norm2 << "\t"	// reference (harmAsym, harmSym, biharmAsym, biharmSym)
		//	<< harm_energy2 << "\t" << harm_relEnergy << "\t" << harm_energy2_sym << "\t" << harm_relEnergy_sym << "\t"						// approx (harmAsym, harmRelEnergyAsym,harmSym, harmRelEnergySym)
		//	<< biharm_energy2 << "\t" << biharm_relEnergy << "\t" << biharm_energy2_sym << "\t" << biharm_relEnergy_sym << "\t"				// approx (biharmAsym, biharmRelEnergyAsym, biharmSym, biharmRelEnergySym)
		//	<< length_init << "\t" << length1 << "\t" << length2 << "\t" << length2 / length1 * 100.0 << "\t"
		//	<< normL2 << "\n";																												// l2-norm

		ofs.close();
	}
}

void VectorFields::testProjection_EigenBasis_WithRegularizer(const Eigen::MatrixXd& Basis, const Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> &sparseSolver_Ref, const Eigen::SparseMatrix<double>& B_Ref, const Eigen::VectorXd& a_Ref, const Eigen::LDLT<Eigen::MatrixXd>& denseSolver_Red, const Eigen::MatrixXd& B_Red, const Eigen::VectorXd& a_Red, const Eigen::VectorXd& inputFields, const Eigen::SparseMatrix<double>& MReg, double &error)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Testing basis 2D of arbitrary field WITH regularizer...\n";

	// Construct matrices for Test	
	cout << "__[APPROXIMATION]....\n";
	//const double				lambda = 10000 * MF2D.coeff(0,0) / MReg.coeff(0,0);

	Eigen::VectorXd id(MF2D.rows());
	id.setConstant(1.0);
	double en1 = id.transpose()*MF2D*id;
	double en2 = id.transpose()*MReg*id;

	//Eigen::MatrixXd				U = Basis;
	Eigen::VectorXd				v = inputFields;

	/* Perturbed the input fields*/
	//perturbVectorFields(v);
	//pertFields = v;

	srand(time(NULL));
	int l1_ = testID % 5;
	double l2_ = 25 * (l1_ + 1);

	double				lambda = l2_;
	cout << "__lamba: " << lambda << endl;
	lambda = lambda * en1 / en2;

	//Eigen::VectorXd				a = U.transpose()*MF2D*v;
	//Eigen::MatrixXd				B = U.transpose() * (MF2D + lambda*MReg) * U;

	cout << "____Solving linear system variables (with lambda=" << lambda << ")\n";
	Eigen::VectorXd w = denseSolver_Red.solve(a_Red);
	wbEigen = Basis*w;

	cout << "__[REFERENCE]....\n";
	//a = MF2D*v;
	//Eigen::SparseMatrix<double> BRef = MF2D + lambda*MReg;
	//Eigen::PardisoLLT<Eigen::SparseMatrix<double>> sparseSolver2(BRef);
	Eigen::VectorXd  wRef;
	wRef = sparseSolver_Ref.solve(a_Ref);
	projRef = wRef;

	// Compare their L2-Norm
	cout << "____Computing L2-norm \n";
	Eigen::VectorXd diff = wRef - wbEigen;
	double norm1 = diff.transpose()*MF2D*diff;
	double norm2 = wRef.transpose()*MF2D*wRef;
	double normL2 = sqrt(norm1 / norm2);
	//error = normL2;
	cout << "The L-2 Norm is << " << normL2 << endl;
	cout << "____The L2 Norm is << " << normL2 << ": sqrt(" << norm1 << "/" << norm2 << ")" << endl;

	/* Measuring the 'length' of each vector */
	cout << "HARMONIC \n";
	double harm_energy1 = wRef.transpose()*SF2DAsym*wRef;
	double harm_energy2 = wbEigen.transpose()*SF2DAsym*wbEigen;
	double harm_relEnergy = abs(harm_energy1 - harm_energy2) / harm_energy1;
	cout << "____Harmonic Energy => Ref=" << harm_energy1 << ", Approx:" << harm_energy2 << endl;
	cout << "____Relative energy: " << harm_relEnergy << endl;

	/* Measuring the 'length' of each vector */
	double harm_energy1_sym = wRef.transpose()*SF2D*wRef;
	double harm_energy2_sym = wbEigen.transpose()*SF2D*wbEigen;
	double harm_relEnergy_sym = abs(harm_energy1_sym - harm_energy2_sym) / harm_energy1_sym;
	cout << "____Harmonic Energy => Ref=" << harm_energy1_sym << ", Approx:" << harm_energy2_sym << endl;
	cout << "____Relative energy: " << harm_relEnergy_sym << endl;

	/* Measuring the energy */
	cout << "BIHARMONIC \n";
	double biharm_energy1 = wRef.transpose()*B2DAsym*wRef;
	double biharm_energy2 = wbEigen.transpose()*B2DAsym*wbEigen;
	double biharm_relEnergy = abs(biharm_energy1 - biharm_energy2) / biharm_energy1;
	cout << "____Bi-Harmonic ENERGY => Ref=" << biharm_energy1 << ", Approx:" << biharm_energy2 << endl;
	cout << "____aRelative energy: " << biharm_relEnergy << endl;

	/* Measuring the 'length' of each vector */
	cout << "ENERGY SYM \n";
	double biharm_energy1_sym = wRef.transpose()*B2D*wRef;
	double biharm_energy2_sym = wbEigen.transpose()*B2D*wbEigen;
	double biharm_relEnergy_sym = abs(biharm_energy1_sym - biharm_energy2_sym) / biharm_energy1_sym;
	cout << "____Harmonic Energy => Ref=" << biharm_energy1_sym << ", Approx:" << biharm_energy2_sym << endl;
	cout << "____Relative energy: " << biharm_relEnergy_sym << endl;


	/* Measuring the 'Initial Energy' of input vector */
	cout << "Input vector: ENERGY SYM \n";
	double harm_energy_input = v.transpose()*SF2DAsym*v;
	double harm_energy_input_sym = v.transpose()*SF2D*v;
	double biharm_energy_input = v.transpose()*B2DAsym*v;
	double biharm_energy_input_sym = v.transpose()*B2D*v;
	cout << "____Input Harm Energy => aSym" << harm_energy_input << ", Sym:" << harm_energy_input_sym << endl;
	cout << "____Input BiHarm Energy => aSym" << biharm_energy_input << ", Sym:" << biharm_energy_input_sym << endl;

	/* Measuring the 'length' of each vector */
	double length_init = v.transpose()*MF2D*v;
	double length1 = wRef.transpose()*MF2D*wRef;
	double length2 = wbEigen.transpose()*MF2D*wbEigen;
	cout << "Length => Ref=" << length1 << ", Approx:" << length2 << ", initial: " << v.transpose()*MF2D*v << endl;
	cout << "Relative length: " << length2 / length1 * 100.0 << "%" << endl;

	/* Measuring the performance of target/goal */
	cout << "Goal of the minimization \n";
	double goal_energy_input = v.transpose() * B_Ref * v;
	double goal_energy_ref = wRef.transpose()*B_Ref*wRef;
	double goal_energy_app = wbEigen.transpose()*B_Ref*wbEigen;
	double goal_energy_rel = abs(goal_energy_app - goal_energy_ref) / goal_energy_ref;
	cout << "____Input goal" << goal_energy_input << endl;
	cout << "____Goal: Ref" << goal_energy_ref << ", App:" << goal_energy_app << endl;
	error = goal_energy_rel;

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;

	bool writeToFile = true;
	if (writeToFile) {
		std::ofstream ofs;
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Bimba_L2projectionWithReg_modalBasis_sameStorage.txt";
		string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Fertility_L2projectionWithReg_modalBasis_250.txt";

		ofs.open(resultFile, std::ofstream::out | std::ofstream::app);

		ofs << globalConstraints.size() << "\t" << lambda * en2 / en1 << "\t" << lambda << "\t"
			<< harm_energy_input << "\t" << harm_energy_input_sym << "\t" << biharm_energy_input << "\t" << biharm_energy_input_sym << "\t" << goal_energy_input << "\t" << length_init << "\t"	// initial input
			<< harm_energy1 << "\t" << harm_energy1_sym << "\t" << biharm_energy1 << "\t" << biharm_energy1_sym << "\t" << goal_energy_ref << "\t" << norm2 << "\t"	// reference (harmAsym, harmSym, biharmAsym, biharmSym)
			<< harm_energy2 << "\t" << harm_relEnergy << "\t" << harm_energy2_sym << "\t" << harm_relEnergy_sym << "\t" << goal_energy_app << "\t"					// approx (harmAsym, harmRelEnergyAsym,harmSym, harmRelEnergySym)
			<< biharm_energy2 << "\t" << biharm_relEnergy << "\t" << biharm_energy2_sym << "\t" << biharm_relEnergy_sym << "\t"				// approx (biharmAsym, biharmRelEnergyAsym, biharmSym, biharmRelEnergySym)
			<< length_init << "\t" << length1 << "\t" << length2 << "\t" << length2 / length1 * 100.0 << "\t"
			<< goal_energy_rel << "\t" << normL2 << "\n";

		//ofs << globalConstraints.size() << "\t" << lambda * en2 / en1 << "\t" << lambda << "\t"
		//	<< harm_energy_input << "\t" << harm_energy_input_sym << "\t" << biharm_energy_input << "\t" << biharm_energy_input_sym << "\t" << length_init << "\t"	// initial input
		//	<< harm_energy1 << "\t" << harm_energy1_sym << "\t" << biharm_energy1 << "\t" << biharm_energy1_sym << "\t" << norm2 << "\t"	// reference (harmAsym, harmSym, biharmAsym, biharmSym)
		//	<< harm_energy2 << "\t" << harm_relEnergy << "\t" << harm_energy2_sym << "\t" << harm_relEnergy_sym << "\t"						// approx (harmAsym, harmRelEnergyAsym,harmSym, harmRelEnergySym)
		//	<< biharm_energy2 << "\t" << biharm_relEnergy << "\t" << biharm_energy2_sym << "\t" << biharm_relEnergy_sym << "\t"				// approx (biharmAsym, biharmRelEnergyAsym, biharmSym, biharmRelEnergySym)
		//	<< length_init << "\t" << length1 << "\t" << length2 << "\t" << length2 / length1 * 100.0 << "\t"
		//	<< normL2 << "\n";																												// l2-norm

		ofs.close();
	}
}

void VectorFields::testProjection_MyBasis_WithRegularizer(const Eigen::SparseMatrix<double>& Basis,  const Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> &sparseSolver_Ref, const Eigen::SparseMatrix<double>& B_Ref, const Eigen::VectorXd& a_Ref, const Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> &sparseSolver_Red, const Eigen::SparseMatrix<double>& B_Red, const Eigen::VectorXd& a_Red, const Eigen::VectorXd& inputFields, const Eigen::SparseMatrix<double>& MReg, double &error)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> Testing basis 2D of arbitrary field WITH regularizer...\n";

	// Construct matrices for Test	
	cout << "__[APPROXIMATION]....\n";
	
	Eigen::VectorXd id(MF2D.rows());
	id.setConstant(1.0);
	double en1 = id.transpose()*MF2D*id;
	double en2 = id.transpose()*MReg*id;	
	int l1_ = testID % 5;
	double l2_ = 25 * (l1_ + 1);
	double				lambda = l2_;
	lambda = lambda * en1 / en2;
		
	Eigen::VectorXd				v = inputFields;	

	
	Eigen::VectorXd w = sparseSolver_Red.solve(a_Red);

	if (sparseSolver_Red.info() != Eigen::Success) {
		cout << "Cannot solve the linear system. " << endl;
		if (sparseSolver_Red.info() == Eigen::NumericalIssue)
			cout << "NUMERICAL ISSUE. " << endl;
		cout << sparseSolver_Red.info() << endl;
		return;
	}
	wb = Basis*w;

	cout << "__[REFERENCE]....\n";	
	Eigen::VectorXd  wRef;
	wRef = sparseSolver_Ref.solve(a_Ref);
	projRef = wRef;

	// Compare their L2-Norm
	cout << "____Computing L2-norm \n";
	Eigen::VectorXd diff = wRef - wb;
	double norm1 = diff.transpose()*MF2D*diff;
	double norm2 = wRef.transpose()*MF2D*wRef;
	double normL2 = sqrt(norm1 / norm2);
	cout << "The L-2 Norm is << " << normL2 << endl;
	cout << "____The L2 Norm is << " << normL2 << ": sqrt(" << norm1 << "/" << norm2 << ")" << endl;

	/* Measuring the 'length' of each vector */
	cout << "HARMONIC \n";
	double harm_energy1 = wRef.transpose()*SF2DAsym*wRef;
	double harm_energy2 = wb.transpose()*SF2DAsym*wb;
	double harm_relEnergy = abs(harm_energy1 - harm_energy2) / harm_energy1;
	cout << "____Harmonic Energy => Ref=" << harm_energy1 << ", Approx:" << harm_energy2 << endl;
	cout << "____Relative energy: " << harm_relEnergy << endl;

	/* Measuring the 'length' of each vector */
	double harm_energy1_sym = wRef.transpose()*SF2D*wRef;
	double harm_energy2_sym = wb.transpose()*SF2D*wb;
	double harm_relEnergy_sym = abs(harm_energy1_sym - harm_energy2_sym) / harm_energy1_sym;
	cout << "____Harmonic Energy => Ref=" << harm_energy1_sym << ", Approx:" << harm_energy2_sym << endl;
	cout << "____Relative energy: " << harm_relEnergy_sym << endl;

	/* Measuring the energy */
	cout << "BIHARMONIC \n";
	double biharm_energy1 = wRef.transpose()*B2DAsym*wRef;
	double biharm_energy2 = wb.transpose()*B2DAsym*wb;
	double biharm_relEnergy = abs(biharm_energy1 - biharm_energy2) / biharm_energy1;
	cout << "____Bi-Harmonic ENERGY => Ref=" << biharm_energy1 << ", Approx:" << biharm_energy2 << endl;
	cout << "____aRelative energy: " << biharm_relEnergy << endl;

	/* Measuring the 'length' of each vector */
	cout << "ENERGY SYM \n";
	double biharm_energy1_sym = wRef.transpose()*B2D*wRef;
	double biharm_energy2_sym = wb.transpose()*B2D*wb;
	double biharm_relEnergy_sym = abs(biharm_energy1_sym - biharm_energy2_sym) / biharm_energy1_sym;
	cout << "____Harmonic Energy => Ref=" << biharm_energy1_sym << ", Approx:" << biharm_energy2_sym << endl;
	cout << "____Relative energy: " << biharm_relEnergy_sym << endl;


	/* Measuring the 'Initial Energy' of input vector */
	cout << "Input vector: ENERGY SYM \n";
	double harm_energy_input = v.transpose()*SF2DAsym*v;
	double harm_energy_input_sym = v.transpose()*SF2D*v;
	double biharm_energy_input = v.transpose()*B2DAsym*v;
	double biharm_energy_input_sym = v.transpose()*B2D*v;
	cout << "____Input Harm Energy => aSym" << harm_energy_input << ", Sym:" << harm_energy_input_sym << endl;
	cout << "____Input BiHarm Energy => aSym" << biharm_energy_input << ", Sym:" << biharm_energy_input_sym << endl;

	/* Measuring the 'length' of each vector */
	double length_init = v.transpose()*MF2D*v;
	double length1 = wRef.transpose()*MF2D*wRef;
	double length2 = wb.transpose()*MF2D*wb;
	cout << "Length => Ref=" << length1 << ", Approx:" << length2 << ", initial: " << v.transpose()*MF2D*v << endl;
	cout << "Relative length: " << length2 / length1 * 100.0 << "%" << endl;

	/* Measuring the performance of target/goal */
	cout << "Goal of the minimization \n";
	double goal_energy_input = v.transpose() * B_Ref * v;
	double goal_energy_ref = wRef.transpose()*B_Ref*wRef;
	double goal_energy_app = wb.transpose()*B_Ref*wb;
	double goal_energy_rel = sqrt(abs(goal_energy_app - goal_energy_ref) / goal_energy_ref);
	cout << "____Input goal" << goal_energy_input << endl;
	cout << "____Goal: Ref" << goal_energy_ref << ", App:" << goal_energy_app << endl;
	error = goal_energy_rel;

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;

	bool writeToFile = true;
	if (writeToFile) {
		std::ofstream ofs;
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Fertility_L2projectionWithReg_eigenFields" + to_string(Basis.cols()) + "_2.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Armadillo_L2projectionWithReg_eigenFields" + to_string(Basis.cols()) + ".txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Kitten_L2projectionWithReg_eigenFields_" + to_string(Basis.cols()) + "_40sup.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Kitten_L2projectionWithReg_eigenFieldsNRot_" + to_string(Basis.cols()) + "_40sup.txt";

		string resultFile = "D:/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Fertility_L2projectionWithReg_variousBases_combined.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Armadillo_L2projectionWithReg_eigFields.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Armadillo_L2projectionWithReg_optAlg.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Armadillo_L2projectionWithReg_eigFields10.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Armadillo_L2projectionWithReg_eigPatch.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Armadillo_L2projectionWithReg_grad.txt";

		ofs.open(resultFile, std::ofstream::out | std::ofstream::app);

		ofs << globalConstraints.size() << "\t" << lambda * en2 / en1 << "\t" << lambda << "\t"
			<< harm_energy_input << "\t" << harm_energy_input_sym << "\t" << biharm_energy_input << "\t" << biharm_energy_input_sym << "\t" << goal_energy_input << "\t" << length_init << "\t"	// initial input
			<< harm_energy1 << "\t" << harm_energy1_sym << "\t" << biharm_energy1 << "\t" << biharm_energy1_sym << "\t" << goal_energy_ref << "\t" << norm2 << "\t"	// reference (harmAsym, harmSym, biharmAsym, biharmSym)
			<< harm_energy2 << "\t" << harm_relEnergy << "\t" << harm_energy2_sym << "\t" << harm_relEnergy_sym << "\t" << goal_energy_app << "\t"					// approx (harmAsym, harmRelEnergyAsym,harmSym, harmRelEnergySym)
			<< biharm_energy2 << "\t" << biharm_relEnergy << "\t" << biharm_energy2_sym << "\t" << biharm_relEnergy_sym << "\t"				// approx (biharmAsym, biharmRelEnergyAsym, biharmSym, biharmRelEnergySym)
			<< length_init << "\t" << length1 << "\t" << length2 << "\t" << length2 / length1 * 100.0 << "\t"
			<< goal_energy_rel << "\t" << normL2 << "\n";																												// l2-norm

		ofs.close();
	}
}

void VectorFields::projectionTest()
{
	cout << "PROJECTION TEST! \n";
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;	


	const int NUM_TEST = 100;
	Eigen::VectorXd errors1(NUM_TEST), errors2(NUM_TEST);
	Eigen::MatrixXd DesignedFields(2 * F.rows(), NUM_TEST);
	Eigen::MatrixXd PerturbedFields(2 * F.rows(), NUM_TEST);

	bool useEigenBasis = false; 
	//bool readFieldsFromFile = true;
	bool readDesFieldsFromFile = true;
	bool readPertFieldsFromFile = true;

	/* Loading eigen basis for Armadillo */
	Eigen::MatrixXd EigenBasis;
	string eigBasisFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/Armadillo_1000_eigenfields_Ref_2";
	if (useEigenBasis)
	{
		ReadDenseMatrixFromMatlab(EigenBasis, eigBasisFile, 172964, 1000);
	}

	/* Fields and perturbed Fields */
	string desFieldsFile  = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/CDragon_DesignedFields";
	string pertFieldsFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/CDragon_PerturbedFields";
	if (readDesFieldsFromFile) {
		ReadDenseMatrixFromMatlab(DesignedFields, desFieldsFile, 2*F.rows(), NUM_TEST);
	}
	if (readPertFieldsFromFile) {
		ReadDenseMatrixFromMatlab(PerturbedFields, pertFieldsFile, 2 * F.rows(), NUM_TEST);
	}
	//Xf = arbField2D; 

	/* Factorization of the mass matrix approx in reduced system */
	Eigen::SparseMatrix<double>							B_NR = Basis.transpose() * MF2D * Basis;
	Eigen::PardisoLDLT<Eigen::SparseMatrix<double>>		sparseSolver_NR(B_NR);

	const int iStart = 0;
	for (int i = iStart; i < iStart+NUM_TEST; i++)
	//for (int i = 37; i < 39; i++)
	{
		t1 = chrono::high_resolution_clock::now();
		testID = i; 

		
		Eigen::Vector3d lambda;
		lambda(0) = 1.0; // 100 * MF2D.coeff(0, 0) / SF2D.coeff(0, 0);		// on harmonic energy
		lambda(1) = 1e-4; // 100 * MF2D.coeff(0, 0) / B2D.coeff(0, 0);		// on bi-harmonic energy
		lambda(2) = 0.4;

		// to get the number of constraints (for output purpose only)
		if (readDesFieldsFromFile || readPertFieldsFromFile) {
			constructRandomHardConstraints();
		}
		if (readDesFieldsFromFile) {
			Xf = DesignedFields.col(i);
		}
		if(readPertFieldsFromFile){			
			pertFields = PerturbedFields.col(i);
		}

		/* In case none are read from files, constructs one! */
		if(!readPertFieldsFromFile && !readDesFieldsFromFile)
		{
			setupGlobalProblem(lambda);
			DesignedFields.col(i) = Xf;
			Eigen::VectorXd vp_ = Xf;
			perturbVectorFields(vp_);
			pertFields = vp_;
			PerturbedFields.col(i) = pertFields;
		}
		cout << "Computing projeciton \n"; 
		printf("Size of Xf: %d | pertFields=%d | MF=%dx%d | SF=%dx%d | BNR: %dx%d \n", Xf.rows(), pertFields.rows(), MF2D.rows(), MF2D.cols(), SF2D.rows(), SF2D.cols(), B_NR.rows(), B_NR.cols());
		/* Projection to the subspace */
		/* Reference results */
		//setupGlobalProblem(Eigen::Vector3d(1,1,1));
		Eigen::VectorXd						a_NR = (Basis.transpose()*(MF2D*Xf));
		///testProjection_MyBasis_NoRegularizer(Basis, sparseSolver_NR, B_NR, a_NR, Xf, errors1(i - iStart));
		//testProjection_EigenBasis_NoRegularizer(EigenBasis, Xf, errors2(i));

		testProjection_MyBasis_WithRegularizer(Basis, pertFields, SF2DAsym, errors2(i));
		
		//testProjection_EigenBasis_WithRegularizer(EigenBasis, pertFields, SF2DAsym, errors2(i));

		t2 = chrono::high_resolution_clock::now();
		duration = t2 - t1;
		//printf("[%d] run => Error=%.10f (in %.3f seconds) \n", i, errors1(i), duration.count());		
		printf("[%d] run => [My Basis] Error=%.10f\n", i, errors1(i));
		printf("            [EigenBasis] Error=%.10f (in %.3f seconds) \n", errors2(i), duration.count());
	}

	//if (!readPertFieldsFromFile && !readDesFieldsFromFile)
	//{
	//	WriteDenseMatrixToMatlab(DesignedFields, desFieldsFile);
	//	WriteDenseMatrixToMatlab(PerturbedFields, pertFieldsFile);
	//}

	cout << "ERRORS: \n" <<  errors1 << endl << "ERRORS2 \n" << errors2 << endl; 
}

void VectorFields::projectionTest(bool &readDesFieldsFromFile, bool &readPertFieldsFromFile, bool &useEigenBasis, int start, int nTests)
{
	cout << "PROJECTION TEST! \n";
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;

	bool projTestNoReg = true;
	bool projTestWithReg = true;

	const int NUM_TEST = nTests;
	Eigen::VectorXd errors1(NUM_TEST), errors2(NUM_TEST);
	Eigen::MatrixXd DesignedFields(2 * F.rows(), NUM_TEST);
	Eigen::MatrixXd PerturbedFields(2 * F.rows(), NUM_TEST);

	//bool useEigenBasis = false;

	/* Loading eigen basis for Armadillo */
	Eigen::MatrixXd EigenBasis;
	Eigen::MatrixXd EigenBasisFull, EigenBasis2;
	string eigBasisFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/CDragon_500_eigenfields_Ref";
	//string eigBasisFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/Fertility_250_eigenfields_Ref";
	string eigBasisFile2 = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/Kitten_401-667_eigenfields_Ref";

	/* Factorization of the mass matrix approx in reduced system */
	printf("___Factorizing the reduced massmatrix Md %dx%d\n", MF2D.rows(), MF2D.cols());
	Eigen::SparseMatrix<double>							B_NR = Basis.transpose() * MF2D * Basis;
	Eigen::PardisoLDLT<Eigen::SparseMatrix<double>>		sparseSolver_NR(B_NR);
	Eigen::MatrixXd										B_MB;
	Eigen::LDLT<Eigen::MatrixXd>						denseSolver_NR;

	/* Factorization for regularizer */
	cout << "Setting up factorization of the proj with regularizer\n";
	const int NUM_SOLVER = 5; 
	Eigen::VectorXd id(MF2D.rows());
	id.setConstant(1.0);
	double en1 = id.transpose()*MF2D*id;
	double en2 = id.transpose()*SF2DAsym*id;
	vector<Eigen::PardisoLDLT<Eigen::SparseMatrix<double>>> sparseSolver_WR_Ref(NUM_SOLVER);
	vector<Eigen::PardisoLDLT<Eigen::SparseMatrix<double>>> sparseSolver_WR_Red(NUM_SOLVER);
	vector<Eigen::LDLT<Eigen::MatrixXd>>					denseSolver_WR_Red(NUM_SOLVER);
	vector<Eigen::SparseMatrix<double>>						BRef(NUM_SOLVER);
	vector<Eigen::SparseMatrix<double>>						BRed_Sp (NUM_SOLVER);
	vector<Eigen::MatrixXd>									BRed_Dn(NUM_SOLVER);

	/* Test non-zeros in the reduced system */
	printf("___Factorizing the reduced massmatrix Md %dx%d\n", SF2DAsym.rows(), SF2DAsym.cols());
	t1 = chrono::high_resolution_clock::now();
	Eigen::SparseMatrix<double> SFRed = Basis.transpose()*SF2DAsym*Basis;
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;

	//printf("Time for computing reduced laplacian: %.10f\n", duration.count());
	//printf("__non-zero entries: %d \n", SFRed.nonZeros());
	//printf("__non-zero per-row: %.10f \n", (double) SFRed.nonZeros() /  (double) Basis.cols());
	//printf("__non-zero percentrage: %.10f \n", 100.0 * (double) SFRed.nonZeros() / (double) (SFRed.rows()*SFRed.cols()));
	//printf("__Basis size=%dx%d \n", Basis.rows(), Basis.cols());

	
	if (projTestWithReg)
	{
		for (int i = 0; i < NUM_SOLVER; i++)
		{
			double l2_ = 25 * (i + 1);
			double	lambda = l2_ * en1 / en2;
			//Eigen::SparseMatrix<double> BRef = MF2D + lambda*SF2DAsym;
			//Eigen::SparseMatrix<double> BRed = Basis.transpose()*BRef*Basis;
	
			cout << "____[" << i << "] Factorizing the sparse matrix (M+a*S) (Ref and Red) \n";
			BRef[i] = MF2D + lambda*SF2DAsym;
			BRed_Sp[i] = Basis.transpose()*BRef[i] * Basis;
			sparseSolver_WR_Ref[i].analyzePattern(BRef[i]);
			sparseSolver_WR_Ref[i].factorize(BRef[i]);
			sparseSolver_WR_Ref[i].pardisoParameterArray()[1] = 0;
			sparseSolver_WR_Ref[i].pardisoParameterArray()[59] = 2;
			sparseSolver_WR_Red[i].analyzePattern(BRed_Sp[i]);
			sparseSolver_WR_Red[i].factorize(BRed_Sp[i]);
			sparseSolver_WR_Red[i].pardisoParameterArray()[1] = 0;
			sparseSolver_WR_Red[i].pardisoParameterArray()[59] = 2;
		}
	}
	
	if (useEigenBasis)
	{
		cout << "Reading eigenbasis (if possible...) \n";
		//ReadDenseMatrixFromMatlab(EigenBasis2, eigBasisFile2, 2*F.rows(), 267);
		//ReadDenseMatrixFromMatlab(EigenBasisFull, eigBasisFile, 2*F.rows(), 40);
		ReadDenseMatrixFromMatlab(EigenBasis, eigBasisFile, 2 * F.rows(), 500);
	
		// same storage
		//EigenBasis.resize(EigenBasisFull.rows(), 40);
		//EigenBasis = EigenBasisFull.block(0, 0, EigenBasisFull.rows(), 40);
		//EigenBasis = EigenBasisFull;
		// same performance
		//cout << "Alloc. memory for eigen basis\n";
		//EigenBasis.resize(EigenBasisFull.rows(), EigenBasisFull.cols() + EigenBasis2.cols());
		//cout << "Combining to eigen basis\n";
		//EigenBasis.block(0, 0, EigenBasisFull.rows(), EigenBasisFull.cols()) = EigenBasisFull.block(0, 0, EigenBasisFull.rows(), EigenBasisFull.cols());
		//EigenBasis.block(0, EigenBasisFull.cols(), EigenBasis2.rows(), EigenBasis2.cols()) = EigenBasis2.block(0, 0, EigenBasis2.rows(), EigenBasis2.cols());
	
		// Free memory manually
		EigenBasisFull.resize(0, 0);
		EigenBasis2.resize(0, 0);
		
	
		B_MB = EigenBasis.transpose()*MF2D*EigenBasis;
		cout << "Factorizing using eigen basis\n";
		denseSolver_NR.compute(B_MB);
	
		/* The one with regularizer*/
		if (projTestWithReg)
		{
			for (int i = 0; i < NUM_SOLVER; i++)
			{
				double l2_ = 25 * (i + 1);
				double	lambda = l2_ * en1 / en2;
				//BRef[i] = MF2D + lambda*SF2DAsym;
				BRed_Dn[i] = EigenBasis.transpose()*BRef[i] * EigenBasis;
	
				cout << "____[" << i << "] Factorizing the dense matrix (M+a*S) (Red) \n";
				denseSolver_WR_Red[i].compute(BRed_Dn[i]);
			}
		}
	}
	
	/* Fields and perturbed Fields */
	string desFieldsFile  = "D:/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/Fertility_DesignedFields";
	string pertFieldsFile = "D:/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/Fertility_PerturbedFields";
	if (readDesFieldsFromFile) {
		ReadDenseMatrixFromMatlab(DesignedFields, desFieldsFile, 2 * F.rows(), NUM_TEST);
	}
	if (readPertFieldsFromFile) {
		ReadDenseMatrixFromMatlab(PerturbedFields, pertFieldsFile, 2 * F.rows(), NUM_TEST);
	}
	//Xf = arbField2D; 
	
	for (int i = start; i < start + NUM_TEST; i++)
		//for (int i = 37; i < 39; i++)
	{
		t1 = chrono::high_resolution_clock::now();
		testID = i;
	
		Eigen::Vector3d lambda;
		lambda(0) = 1.0; // 100 * MF2D.coeff(0, 0) / SF2D.coeff(0, 0);		// on harmonic energy
		lambda(1) = 1e-4; // 100 * MF2D.coeff(0, 0) / B2D.coeff(0, 0);		// on bi-harmonic energy
		lambda(2) = 0.4;
	
		// to get the number of constraints (for output purpose only)
		if (readDesFieldsFromFile || readPertFieldsFromFile) {
			constructRandomHardConstraints();
		}
		if (readDesFieldsFromFile) {
			pertFields = PerturbedFields.col(i-start);
		}
		if (readPertFieldsFromFile) {
			Xf = DesignedFields.col(i-start);
		}
	
		/* In case none are read from files, constructs one! */
		if (!readPertFieldsFromFile && !readDesFieldsFromFile)
		{
			setupGlobalProblem(lambda);
			DesignedFields.col(i-start) = Xf;
			Eigen::VectorXd vp_ = Xf;
			perturbVectorFields(vp_);
			pertFields = vp_;
			PerturbedFields.col(i-start) = pertFields;
		}
		/* Projection to the subspace */
		/* Reference results */
		Eigen::VectorXd										a_Ref = MF2D*Xf;
		Eigen::VectorXd										a_NR = (Basis.transpose()*a_Ref);
		testProjection_MyBasis_NoRegularizer(Basis, sparseSolver_NR, B_NR, a_NR, Xf, errors1(i-start));
		//testProjection_MyBasis_WithRegularizer(Basis, pertFields, SF2DAsym, errors2(i-start));
		testProjection_MyBasis_WithRegularizer(Basis, sparseSolver_WR_Ref[i%NUM_SOLVER], BRef[i%NUM_SOLVER], a_Ref, sparseSolver_WR_Red[i%NUM_SOLVER], BRed_Sp[i%NUM_SOLVER], a_NR, pertFields, SF2DAsym, errors1(i - start));
	
		
		//a_NR = (EigenBasis.transpose()*(MF2D*Xf));
		//testProjection_EigenBasis_NoRegularizer(EigenBasis, denseSolver_NR, a_NR, Xf, errors2(i-start));
		//testProjection_EigenBasis_WithRegularizer(EigenBasis, pertFields, SF2DAsym, errors2(i-start));
		//testProjection_EigenBasis_WithRegularizer(EigenBasis, sparseSolver_WR_Ref[i%NUM_SOLVER], BRef[i%NUM_SOLVER], a_Ref, denseSolver_WR_Red[i%NUM_SOLVER], BRed_Dn[i%NUM_SOLVER], a_NR, pertFields, SF2DAsym, errors2(i - start));
	
		t2 = chrono::high_resolution_clock::now();
		duration = t2 - t1;
		//printf("[%d] run => Error=%.10f (in %.3f seconds) \n", i, errors1(i), duration.count());		
		printf("[%d] run => [No. Reg.]  Error=%.10f\n", i, errors1(i));
		//printf("            [with Reg.] Error=%.10f (in %.3f seconds) \n", errors2(i), duration.count());
	}
	
	if (!readPertFieldsFromFile && !readDesFieldsFromFile)
	{
		/* Write the matrices to matlab:
		** this will be in a new file, so I need to combine them with previous files (ensure that we still save the old file, by renaming it) */
		WriteDenseMatrixToMatlab(DesignedFields, desFieldsFile);
		WriteDenseMatrixToMatlab(PerturbedFields, pertFieldsFile);
	
		// set to be true for the next test
		readPertFieldsFromFile = true;
		readDesFieldsFromFile = true;
	}
	
	cout << "ERRORS: \n" << errors1 << endl << "ERRORS2 \n" << errors2 << endl;
}

void VectorFields::projectionSimpleL2Test()
{
	Eigen::VectorXd										a_NR = (Basis.transpose()*(MF2D*Xf));
	Eigen::SparseMatrix<double>							B_NR = Basis.transpose() * MF2D * Basis;
	Eigen::PardisoLDLT<Eigen::SparseMatrix<double>>		sparseSolver_NR(B_NR);
	double error; 

	Eigen::VectorXd w = sparseSolver_NR.solve(a_NR);
	wb = Basis*w;

	// Compare their L2-Norm
	//cout << "____Computing L2-norm \n";
	Eigen::VectorXd diff = Xf - wb;
	double length1 = wb.transpose()*MF2D*wb;
	double norm1 = diff.transpose()*MF2D*diff;
	double norm2 = Xf.transpose()*MF2D*Xf;
	double normL2 = sqrt(norm1 / norm2);
	error = normL2;
	printf("The L-2 norm is: %.5f(sqrt(%.5f/%.5f)). \n", normL2, norm1, norm2);

	/* Construct a scaled mass matrix */
	vector<Eigen::Triplet<double>> MTriplet;
	Eigen::SparseMatrix<double> M2(2 * F.rows(), 2 * F.rows());
	for (int i = 0; i < F.rows(); i++)
	{
		MTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, faceScale(i)*doubleArea(i) / 2.0));
		MTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, faceScale(i)*doubleArea(i) / 2.0));
	}
	M2.setFromTriplets(MTriplet.begin(), MTriplet.end());
	length1 = wb.transpose()*M2*wb;
	norm1 = diff.transpose()*M2*diff;
	norm2 = Xf.transpose()*M2*Xf;
	normL2 = sqrt(norm1 / norm2);
	printf("The scaled L-2 norm is: %.5f(sqrt(%.5f/%.5f)). \n", normL2, norm1, norm2);

	/* Measuring the 'length' of each vector */
	double harm_energy1 = Xf.transpose()*SF2DAsym*Xf;
	double harm_energy2 = wb.transpose()*SF2DAsym*wb;
	double harm_relEnergy = sqrt(abs(harm_energy1 - harm_energy2) / harm_energy1);
	printf("The relative energy is: %.5f((%.5f/%.5f)). \n", harm_relEnergy, abs(harm_energy1 - harm_energy2), harm_energy1);
}

void VectorFields::convergenceTest()
{
	const int NUM_DIFF_BASIS = 6;
	vector<std::string> basisFile;
	basisFile.reserve(NUM_DIFF_BASIS);	


	basisFile.push_back("D:/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Fertility_2000_Eigfields_40sup");
	basisFile.push_back("D:/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Fertility_2000_OptAlg_40sup");
	basisFile.push_back("D:/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Fertility_400_EigFields10_40sup");
	basisFile.push_back("D:/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Fertility_2000_EigPatch_40sup");
	basisFile.push_back("D:/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Fertility_2000_Grad_40sup");

	//basisFile.push_back("D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Kitten_500_EigFields_35sup");
	//basisFile.push_back("D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Kitten_1000_EigFields_35sup");
	//basisFile.push_back("D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Kitten_2000_EigFields_35sup");
	//basisFile.push_back("D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Kitten_5000_EigFields_35sup");
	//basisFile.push_back("D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Kitten_10000_EigFields_35sup");
	//basisFile.push_back("D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Kitten_20000_EigFields_35sup");
	//basisFile.push_back("D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Kitten_50000_EigFields_35sup");

	//vector<int> subspdim{500, 1000, 2000, 5000, 10000, 20000, 50000};
	//vector<int> subspdim{ 2000 };
	//for (int i : subspdim) {
	//	basisFile.push_back("D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Fertility_" + to_string(i) + "_Eigfields_40sup");
	//}

	/* For the projection tests */
	bool readDesFields  = true;
	bool readPertFields = true;
	bool useEigenBasis  = false;
	int idStart = 0;
	int NUM_TEST = 100;

	//int counter = 0; 
	for (string file_ : basisFile)
	{
		retrieveBasis(file_);
		//if (counter > 0) { readDesFields = true; readPertFields = true; }
		projectionTest(readDesFields, readPertFields, useEigenBasis, idStart, NUM_TEST);
		//counter++;
	}
}

void VectorFields::compareModalBasis_SamePerformance()
{
	
}

void VectorFields::compareModalBasis_SameStorage()
{
	string ourBasis = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Arma_2000_EigFields_35sup";

	/* For the projection tests */
	bool readDesFields = true;
	bool readPertFields = true;
	bool useEigenBasis = true;
	int idStart = 0;
	int NUM_TEST = 50;

	retrieveBasis(ourBasis);
	projectionTest(readDesFields, readPertFields, useEigenBasis, idStart, NUM_TEST);
}

void VectorFields::testMappingMatrix()
{
	// Should be identity
	Eigen::SparseMatrix<double> ATA = A.transpose()*A;
	//visualizeSparseMatrixInMatlab(ATA);
	cout << "Block of ATA " << endl << ATA.block(5, 5, 5, 5) << endl;

	// Orthogonality of the matrix
	Eigen::Vector3d e, f;
	srand(time(NULL));
	const int idx = rand() % F.rows();
	e << A.coeffRef(3 * idx, 2 * idx), A.coeffRef(3 * idx + 1, 2 * idx), A.coeffRef(3 * idx + 2, 2 * idx);
	f << A.coeffRef(3 * idx, 2 * idx + 1), A.coeffRef(3 * idx + 1, 2 * idx + 1), A.coeffRef(3 * idx + 2, 2 * idx + 1);
	cout << idx << " => e*f = " << e.dot(f) << endl;
}

void VectorFields::testAdjMV()
{
	for (int i = 0; i < V.rows(); i++) {
		cout << i << ":";
		for (std::set<int, double>::iterator it = AdjMV[i].begin(); it != AdjMV[i].end(); ++it) {
			cout << *it << ", ";
		}
		cout << endl;
	}
}

void VectorFields::testAdjacencyAndEdges()
{
	for (int i = 0; i < F.rows(); i++) {
		//for (int i = 0; i < min(100, (int)F.rows()); i++) {
		printf("F(%d) [%d, %d, %d] :=>"
			"(0) F(%d)[%d, %d, %d] on edges(%d, %d),"
			"(1) F(%d)[%d, %d, %d] on edges(%d, %d),"
			"(2) F(%d)[%d, %d, %d] on edges(%d, %d)\n",
			i, F(i, 0), F(i, 1), F(i, 2),
			AdjMF3N(i, 0), F(AdjMF3N(i, 0), 0), F(AdjMF3N(i, 0), 1), F(AdjMF3N(i, 0), 2), EdgePairMatrix(i, 0), EdgePairMatrix(i, 1),
			AdjMF3N(i, 1), F(AdjMF3N(i, 1), 0), F(AdjMF3N(i, 1), 1), F(AdjMF3N(i, 1), 2), EdgePairMatrix(i, 2), EdgePairMatrix(i, 3),
			AdjMF3N(i, 2), F(AdjMF3N(i, 2), 0), F(AdjMF3N(i, 2), 1), F(AdjMF3N(i, 2), 2), EdgePairMatrix(i, 4), EdgePairMatrix(i, 5));
	}
}

void VectorFields::testDijkstraFace()
{
	dijkstraFace.resize(F.rows());

	// For single-sourced Dijkstra
	//const int source = *(NeighRing[0].begin());
	//computeDijkstraDistanceFace(source, dijkstraFace);

	// For multiple-sourced Dijkstra
	const int numSource = 5;
	Eigen::VectorXi source(numSource);
	srand(time(NULL));
	for (int i = 0; i < numSource; i++) {
		source(i) = rand() % F.rows();
	}
	computeDijkstraDistanceFaceMultSource(source, dijkstraFace);
}

void VectorFields::testCurlEnergy()
{
	//double gradEnergy = gradArbField3D.transpose() * LapCurl3D * gradArbField3D;
	//cout << "The energy is " << gradEnergy << endl;
}

int VectorFields::selectRandomFace()
{
	srand(time(NULL));
	//std::rand() / ((RAND_MAX) / 6);
	//int randID = rand() % F.rows();
	int randID = rand() / (RAND_MAX/F.rows());
	cout << "Selected face: " << randID << endl;
	return randID;
}

void VectorFields::checkB2DStructure()
{
	for (int i = 0; i < B2D.outerSize(); i++) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(B2D, i); it; ++it) {
			if (it.row() == 200 && it.col() % 2 == 0) {
				cout << "VALUES IN B2D :: N(0)=";
				cout << it.col() / 2 << ", " << endl;
			}
		}
	}

	cout << "NEIGHBORS :: N(0)=";
	for (int i = 0; i < AdjMF3N.cols(); i++) {
		cout << (AdjMF3N(100, i)) << ", ";
	}
	cout << endl;

	cout << "2RING:: N(0)=";
	for (int i : AdjMF2Ring[100]) {
		cout << i << ", ";
	}
	cout << endl;

}

// ITEMS FOR TESTING ONLY
void VectorFields::TEST_VECTOR(igl::opengl::glfw::Viewer &viewer, const string& meshFile)
{
	/* ========================= PRE-PROCESS ==============================*/
	cout << "========================= PRE-PROCESS ==============================\n";
	readMesh(meshFile);
	scaleMesh(V, F);

	viewer.data().set_mesh(V, F);
	//viewer.append_mesh();
	//viewer.selected_data_index = 1;
	//viewer.data().set_mesh(V, F);
	//viewer.data().show_lines = false;
	//viewer.selected_data_index = 0;

	//Eigen::SparseMatrix<double> ChrisSparseMat;
	//ReadChristopherStiffnessMatrix("D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Local Fields/Models/hodgeLaplace.txt", ChrisSparseMat);
	//WriteSparseMatrixToMatlab(ChrisSparseMat, "hello");
	//readArrowMesh("../LocalFields/Models/arrow.obj");
	computeEdges();
	computeAverageEdgeLength();
	computeFaceCenter();
	computeFaceNormal();
	constructVFNeighbors();
	//constructVFNeighborsFull();
	//constructVFAdjacency();
	//testAdjacency();
	constructVertexAdjacencyMatrix();
	constructFaceAdjacency3NMatrix();
	constructFaceAdjacency2RingMatrix();
	constructEVList();
	constructEFList(); 
	selectFaceToDraw(75000); 
	//selectFaceToDraw(max((int) round(F.rows()/20.0), 5000));
	//selectFaceToDraw(F.rows());

	/* MATRIX CONSTRUCTIONS */
	constructMassMatrices();
	constructRotationMatrix();
	constructMappingMatrix();

	/* =========== Test on PROBLEM SOLVING-related functionalities ================*/
	//constructStiffnessMatrices();
	//loadStiffnessMatrices();
	constructGradient3D();
	constructGradientStar3D();
	constructStiffnessMatrices_Implicit();	
	constructMatrixB();
	
	//constructConstraints();
	//checkB2DStructure();

	//////* ====================== GLOBAL PROBLEM ====================*/
	////////cout << "\n========================= GLOBAL PROBLEM =============================\n";
	Eigen::Vector3d lambda;
	lambda(0) = 0.25; // 100 * MF2D.coeff(0, 0) / SF2D.coeff(0, 0);		// on harmonic energy
	lambda(1) = 1e-4; // 100 * MF2D.coeff(0, 0) / B2D.coeff(0, 0);		// on bi-harmonic energy
	lambda(2) = 1;
	//setupGlobalProblem(lambda);

	/* ====================== LOCAL ELEMENTS ====================*/
	//string filename_basis = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Kitten_2000_Eigfields_40sup_Spectra";
	//string filename_basis = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Arma43k_2000_Eigfields_40sup_Matlab";
	//string filename_basis = "D:/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_CDragon_2000_OptAlgAsym_30sup";

	//string filename_basis = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Kitten_5000_Eigfields_40sup";

	//string filename_basis = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Arma_2000_EigFields_35sup";
	//string filename_basis = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Arma_2000_OptAlg_30sup";
	//string filename_basis = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Arma_2000_eigFields10_30sup";
	//string filename_basis = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Arma_2000_EigPatch_30sup";
	//string filename_basis = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Arma_2000_Grad_30sup";

	//string filename_basis = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Kitten_2000_OptAlg_30sup";
	//string filename_basis = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Kitten_2000_EigFields_35sup";
	//string filename_basis = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Kitten_2000_eigFields10_30sup";
	//string filename_basis = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Kitten_2000_EigPatch_30sup";
	//string filename_basis = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Kitten_2000_Grad_30sup";

	/* For convergence */
	//string filename_basis = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Kitten_500_EigFields_35sup";
	//string filename_basis = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Kitten_1000_EigFields_35sup";
	//string filename_basis = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Kitten_2000_EigFields_35sup";
	//string filename_basis = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Kitten_5000_EigFields_35sup";
	//string filename_basis = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Kitten_10000_EigFields_35sup";
	//string filename_basis = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Kitten_50000_EigFields_35sup";


	//string filename_basis = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_CDragon_2000_eigFields10_30sup";
	//string filename_basis = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_CDragon_2000_EigPatch_30sup";
	//string filename_basis = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_CDragon_2000_Grad_30sup";

	//string filename_vfields = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/Arma_4_Ref_eigFields_2.txt";	
	//string filename_vfields = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/Arma_constraintFields_randConst_Asym_1.txt";	//random constraint
	//string filename_vfields = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/Arma_constraintFields_randConst_Sym_1.txt";	//random constraint
	//string filename_vfields = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/Arma_constraintFields_farConst_Asym_1.txt";	//fartheset point
	//string filename_vfields = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/Arma_arbFields_xyz-axis.txt";
	//string filename_vfields = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/Arma_arbFields_y-axis.txt";

	//string filename_vfields = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/Arma_constraintFields_user_Asym_0 (from center).txt";
	//string filename_vfields = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/Arma_constraintFields_user_Asym_1 (going left).txt";
	//string filename_vfields = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/Arma_constraintFields_user_Asym_2 (going down).txt";
	//string filename_vfields = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/Arma_constraintFields_user_Asym_3 (center and down).txt";
	//string filename_vfields = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/Arma_constraintFields_user_Asym_4 (around_52).txt";
	string filename_vfields = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/Arma_constraintFields_user_Asym_5 (from right arm_35).txt";

	//testSparseMatrix();

	//string filename_vfields = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/CDragon_constraintFields_1.txt"; //farthest point constraint
	cout << "\n========================= REDUCED/LOCAL-PROBLEM =============================\n";
	numSample = 200; 
	numSupport = 40.0;
	string model = "Fertility_";
	//string filename_basis = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_" + model + to_string(numSample * 2) + "_Eigfields_" + to_string((int)numSupport) + "sup_adaptiveScale_7.5";
	//string filename_basis = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_" + model + to_string(numSample * 2) + "_Eigfields_" + to_string((int)numSupport) + "sup";
	string filename_basis = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_" + model + to_string(numSample * 2) + "_EigFields10_" + to_string((int)numSupport) + "sup";
	//selectAdaptiveRegions(viewer);
	//selectAdaptiveRegions_Curvature(viewer);
	faceScale.resize(F.rows()); faceScale.setConstant(1.0);
	///constructSamples(numSample);
	///constructBasis();	
	///storeBasis(filename_basis);			// Binary, Eigen-base
	//constructMultiBasis();
	//retrieveBasis(filename_basis);	
	///BasisT = Basis.transpose();
	//normalizeBasisAbs(2);
	//visualizeSamples(viewer);
	//visualizeSubdomain(viewer);


	//setupReducedBiLaplacian();
	//preComputeReducedElements();
	//initializeParametersForLifting();
	//setAndSolveUserSystem(lambda);
	//WriteEigenVectorToTxtFile(arbField2D, filename_vfields);
	//LoadEigenVectorFromTxtFile(filename_vfields, arbField2D);

	//cout << "Basis \n" << Basis.col(0).block(0, 0, 100, 1) << endl; 

	/* Test Spectra */
	//testSpectra();
	//testViennaCL2();

	int eigsToCompute = 10;
	string    filename_refField = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/CDragon_2k_50_eigVectorFields_Spectra_Ref";
	string filename_approxField = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/Kitten_25_Approx_EigenBasis_2000dim_30sup";
	//computeEigenFields(eigsToCompute, filename_refField);	
	//retrieveEigenFields(filename_refField);
	//computeApproxEigenFields(eigsToCompute, filename_approxField);
	//retrieveApproxEigenFields();

	// Store the eigenfields as vector fields
	//string filename_vfields = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/Arma_4_Ref_eigFields";
	//WriteEigenVectorToTxtFile(eigFieldFull2D.col(0), filename_vfields+"_0.txt");
	//WriteEigenVectorToTxtFile(eigFieldFull2D.col(1), filename_vfields+"_1.txt");
	//WriteEigenVectorToTxtFile(eigFieldFull2D.col(2), filename_vfields+"_2.txt");
	//WriteEigenVectorToTxtFile(eigFieldFull2D.col(3), filename_vfields+"_3.txt");

	//testEnergyOfLocalPatch(viewer);

	//printDataForVTK();
	//writeEigenFieldsForVTK();

	/* ====================== TESTING BASIS ====================*/
	/* _____ Projection Test ___________________________________*/
	//constructArbitraryField();
	//constructArbitraryField2D();
	//WriteEigenVectorToTxtFile(arbField2D, filename_vfields);
	//LoadEigenVectorFromTxtFile(filename_vfields, arbField2D);
	//double error; 
	//projectionTest();
	convergenceTest();
	//compareModalBasis_SameStorage();

	/* Alignment using (maximal/) principal curvature */
	//Eigen::MatrixXd CCol;
	//Eigen::VectorXd PD, PV;
	//computeMaximalPrincipalCurvature(V, F, PD, PV);
	//visualize2Dfields(viewer, PD, Eigen::RowVector3d(0.0, 0.1, 0.9), 1.0);
	//igl::jet(PV, true, CCol);
	//viewer.data().set_colors(CCol);

	/* _____ Vector fields design test __________________________*/
	//vectorFieldsDesignTest();
	//vectorFieldsDesignTest_Normalized();


	//visualizeGlobalConstraints(viewer);
	//measureDirichletEnergy();

	/* ====================== PARALLEL TRANSPORT ====================*/
	//computeDijkstraForParallelTransport(200, 5000);
	//constructParallelTransport();
	//visualizeParallelTransportPath(viewer);
	//visualizeParallelTransport(viewer);


	/* ==================== VISUALIZATION ======================== */
	/* GLOBAL  */
	//visualizeApproximatedFields(viewer);
	//visualizeGlobalConstraints(viewer);
	//visualizeSingularitiesConstraints(viewer);
	//visualizeSharedEdges(viewer);

	/* LOCAL  */
	//XFullDim = eigFieldFull2D.col(0);
	//visualizeApproxResult(viewer);	
	//visualizeUserConstraints(viewer);
	//visualizeSamples(viewer);
	//visualizeSingularitiesConstraints(viewer);

	/* VISUALIZATION FOR TESTING PURPOSE */
	//visualizeNeighboringRings(viewer);
	//visualizeDijkstraFace(viewer);
	//visualizeArbField(viewer);
	//visualizeVertexFacesNeighbors(viewer, 0);
	//testEdgesAddition(viewer);
	//visualizePatchDijkstra(viewer);

	/* SOFT CONSTRAINTS */
	//visualizeCurveConstraints(viewer);
	//visualizeSoftConstraints(viewer);
	//measureSoftConstraintError(lambda);


	/* MEASURE ACCURACY */
	//measureApproxAccuracyL2Norm();

	/* PROJECTION ON ADAPTIVE SAMPLING */

	filename_vfields = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/Fertility_constraintFields_user_7.txt";
	//loadVectorFieldsFromFile(filename_vfields, Xf);
	//visualizeApproximatedFields(viewer);
	//filename_vfields = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/Fertility_constraintFields_user_7_constraints.txt";
	//loadConstraintsFromFile(filename_vfields);
	//visualizeGlobalConstraints(viewer);

	/* PRojection on adaptive basis */
	//cout << "[1] Projection using adaptive basis \n";
	//projectionSimpleL2Test();
	//wbEigen = wb; 
	//
	///* PRojection on uniform basis */
	//filename_basis = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_" + model + to_string(numSample * 2) + "_Eigfields_" + to_string((int)numSupport) + "sup_Spectra";
	//retrieveBasis(filename_basis);
	//
	//cout << "[2] Projection using regular/isotropic basis \n";
	//projectionSimpleL2Test();

	/* Basis via Coarsening */
	//constructBasis_Coarsening(viewer);
	//printf("Basis =%dx%d \n", Basis.rows(), Basis.cols());
	//projectionTest();

	//computeEigenFields(eigsToCompute, filename_refField);		
	//computeApproxEigenFields(eigsToCompute, filename_approxField);
}

void VectorFields::constructParallelTransport()
{
	Eigen::Vector2d PTvalue;
	PTvalue << sqrt(2.0)/2.0, sqrt(2.0)/2.0; 
	//PTvalue << 1.0, 0.0;

	parallelTransport.resize(PTpath.size());
	parallelTransport[0] = PTvalue;
	PTsharedEdges.resize(2 * PTpath.size());

	for (int i = 0; i < PTpath.size() - 1; i++) {
		int face1 = PTpath[i];
		int face2 = PTpath[i + 1];
		
		// 1. Find shared edge (naively)
		enum class SharedEdgeCase { Case1, Case2, Case3 };
		SharedEdgeCase edgeCase1, edgeCase2;
		Eigen::RowVector3d es;
		for (int f1 = 0; f1 < F.cols(); f1++) {
			for (int f2 = 0; f2 < F.cols(); f2++) {
				bool b1 = F(face1, (f1 + 1) % F.cols()) == F(face2, f2);
				bool b2 = F(face1, f1) == F(face2, (f2 + 1) % F.cols());
				if (b1 && b2) {						
					es = V.row(F(face1, (f1 + 1) % F.cols())) - V.row(F(face1, f1));
					PTsharedEdges[2 * i + 0] = F(face1, f1);
					PTsharedEdges[2 * i + 1] = F(face1, (f1 + 1) % F.cols());
					
					if (f1 == 0)		edgeCase1 = SharedEdgeCase::Case1;	// => edge V0->V1 is the shared edge => it takes 0 step to reach v0
					else if (f1 == 1)	edgeCase1 = SharedEdgeCase::Case3;	// => edge V1->V2 is the shared edge => it takes 2 step to reach v0
					else if (f1 == 2)	edgeCase1 = SharedEdgeCase::Case2;	// => edge V2->V0 is the shared edge => it takes 1 step to reach v0

					if (f2 == 0)		edgeCase2 = SharedEdgeCase::Case1;
					else if (f2 == 1)	edgeCase2 = SharedEdgeCase::Case3;
					else if (f2 == 2)	edgeCase2 = SharedEdgeCase::Case2;
				}
			}
		}

		// 2. Find angles between basis1 and es
		Eigen::VectorXd eVect;
		Eigen::RowVector3d b11, b12;
		eVect = A.block(3 * face1, 2 * face1 + 0, 3, 1);
		b11 << eVect(0), eVect(1), eVect(2);

		double cosR12 = (b11.dot(es)) / (b11.norm()*es.norm());
		if (cosR12 > 1.0) cosR12 = 1.0;
		if (cosR12 <-1.0) cosR12 = -1.0;
		const double angleR12_1 = (edgeCase1 == SharedEdgeCase::Case2 ? 2 * M_PI - acos(cosR12) : acos(cosR12));
		const double cosR12_1 = cos(angleR12_1);
		const double sinR12_1 = sin(angleR12_1);
		//printf("______[%.2f] Rotation matrix R12_1 = [%.3f,%.3f; %.3f, %.3f]\n", angleR12_1*180.0 / M_PI, cosR12_1, -sinR12_1, sinR12_1, cosR12_1);

		Eigen::RowVector3d b21, b22;
		eVect = A.block(3 * face2, 2 * face2, 3, 1);
		b21 << eVect(0), eVect(1), eVect(2);

		double cosR21 = (b21.dot(es)) / (b21.norm()*es.norm());
		if (cosR21 > 1.0) cosR21 = 1.0;
		if (cosR21 < -1.0) cosR21 = -1.0;
		double angleR21_1 = (edgeCase2 == SharedEdgeCase::Case3 ? 2 * M_PI - acos(cosR21) : acos(cosR21));
		const double cosR21_1 = cos(angleR21_1);
		const double sinR21_1 = sin(angleR21_1);
		//printf("______[%.2f] Rotation matrix R22_1 = [%.2f,%.2f; %.2f, %.2f]\n", angleR21_1*180.0 / M_PI, cosR21_1, -sinR21_1, sinR21_1, cosR21_1);
		
		double anglePhi1;
		Eigen::Vector2d bas(1.0, 0.0);
		if (parallelTransport[i](1) < 0)
			anglePhi1 = 2*M_PI - acos(bas.dot(parallelTransport[i])/parallelTransport[i].norm());
		else
			anglePhi1 = acos(bas.dot(parallelTransport[i]) / parallelTransport[i].norm());

		//printf("______basis to field : [%.2f]\n", anglePhi1*180.0/M_PI);
		const double anglePhi2 = anglePhi1 - angleR12_1 + angleR21_1;
		//printf("______NEW angle : [%.2f]\n", anglePhi2*180.0 / M_PI);
		const double cosBasis = cos(anglePhi2);
		const double sinBasis = sin(anglePhi2);

		// Rotate basis
		// Obtain the mapped basis
		Eigen::Matrix2d RotMat2D;
		RotMat2D << cosBasis, -sinBasis, sinBasis, cosBasis;
		Eigen::Vector2d transported = RotMat2D * bas;// parallelTransport[i];
		parallelTransport[i + 1] = transported;		
	}
}

void VectorFields::writeBasisToFile()
{
	printf("Basis=%dx%d || #F=%d\n", Basis.rows(), Basis.cols(), F.rows());
	Eigen::SparseMatrix<double> Basis3D = A*Basis;
	printf("Basis=%dx%d || Basis3D=%dx%d || #F=%d\n", Basis.rows(), Basis.cols(), Basis3D.rows(), Basis3D.cols(), F.rows());

	WriteSparseMatrixToMatlab(Basis3D, "hello");
}

void VectorFields::writeField3DToFile()
{
	Eigen::VectorXd Field3D = A*XFullDim;
	WriteDenseMatrixToMatlab(Field3D, "hello");
}

void VectorFields::printDataForVTK()
{
	/* PRint vertices*/
	cout << "POINTS " << V.rows() << " float\n";
	cout << V << endl;

	/* Print faces */
	cout << "POLYGONS " << F.rows() << " " << F.rows() * (F.cols()+1) << endl; 
	for (int i = 0; i < F.rows(); i++)
	{
		cout << F.cols() << " " << F(i, 0) << " " << F(i, 1) << " " << F(i, 2) << endl;
	}

	/* Print eigenfields*/
	Eigen::VectorXd eig3D;
	eig3D = A * eigFieldFull2D.col(0);
	cout << "CELL_DATA " << F.rows() << endl; 
	cout << "VECTORS EigenField float\n";
	for (int i = 0; i < F.rows(); i++)
	{
		cout << F.cols() << " " << eig3D(3 * i + 0) << " " << eig3D(3 * i + 1) << " " << eig3D(3 * i + 2) << endl;
	}
}

void VectorFields::writeEigenFieldsForVTK()
{
	for (int id = 0; id < 1; id++)
	{
		string filename = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/VTK and ParaView/Test Data/Torus_4k_EigFields_face_ref_"+ to_string(id+1) +".vtk";
		
		ofstream file(filename);
		if (file.is_open())
		{
			file << "# vtk DataFile Version 2.0\n";
			file << "Torus Eigenfields\n";
			file << "ASCII\n";
			file << "DATASET POLYDATA\n";

			/* PRint vertices*/
			file << "POINTS " << V.rows() << " double\n";
			for (int i = 0; i < V.rows(); i++)
			{
				file << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << "\n";
			}

			/* Print faces */
			file << "POLYGONS " << F.rows() << " " << F.rows() * (F.cols() + 1) << "\n";
			for (int i = 0; i < F.rows(); i++)
			{
				file << F.cols() << " " << F(i, 0) << " " << F(i, 1) << " " << F(i, 2) << "\n";
			}

			/* Print eigenfields*/
			Eigen::VectorXd eig3D;
			Eigen::MatrixXd EigFields;

			//EigFields = Basis.block(0, 0, Basis.rows(), eigFieldReduced2D.cols())*eigFieldReduced2D;
			EigFields = eigFieldFull2D;
			double const scale = 1.0;
			eig3D = A * EigFields.col(id);

			/* FACE-base DATA */
			file << "CELL_DATA " << F.rows() << "\n";
			file << "VECTORS EigenField_" << id+1 << " double\n";
			for (int i = 0; i < F.rows(); i++)
			{
				file << scale*eig3D(3 * i + 0) << " " << scale*eig3D(3 * i + 1) << " " << scale*eig3D(3 * i + 2) << "\n";
				cout << "face=" << i << ": " << scale*eig3D(3 * i + 0) << " " << scale*eig3D(3 * i + 1) << " " << scale*eig3D(3 * i + 2) << "\n";
			}

			/* POINT-base DATA */
			//Eigen::VectorXd VEigFields;
			//VEigFields.setZero(3 * V.rows());
			//Eigen::VectorXi VNumNeighFaces;
			//VNumNeighFaces.setZero(V.rows());
			//
			//for (int i = 0; i < F.rows(); i++)
			//{
			//	VEigFields.block(3 * F(i, 0), 0, 3, 1) += eig3D.block(3 + i, 0, 3, 1);
			//	VEigFields.block(3 * F(i, 1), 0, 3, 1) += eig3D.block(3 + i, 0, 3, 1);
			//	VEigFields.block(3 * F(i, 2), 0, 3, 1) += eig3D.block(3 + i, 0, 3, 1);
			//	//VEigFields.row(F(i, 0)) += EigFields.row(i);
			//	//VEigFields.row(F(i, 1)) += EigFields.row(i);
			//	//VEigFields.row(F(i, 2)) += EigFields.row(i);
			//	VNumNeighFaces(F(i, 0)) += 1;
			//	VNumNeighFaces(F(i, 1)) += 1;
			//	VNumNeighFaces(F(i, 2)) += 1;
			//}
			//
			//for (int i = 0; i < V.rows(); i++)
			//{
			//	VEigFields.block(3 * i, 0, 3, 1) /= (double)VNumNeighFaces(i);
			//	VEigFields.block(3 * i, 0, 3, 1) = VEigFields.block(3 * i, 0, 3, 1).normalized();
			//}
			//
			//file << "POINT_DATA " << V.rows() << "\n";
			//file << "VECTORS EigenField_" << id+1 << " double\n";
			//for (int i = 0; i < V.rows(); i++)
			//{
			//	file << scale*VEigFields(3 * i + 0) << " " << scale*VEigFields(3 * i + 1) << " " << scale*VEigFields(3 * i + 2) << "\n";
			//}

			file.close();
		}
	}
}

void VectorFields::testEdgesAddition(igl::opengl::glfw::Viewer &viewer)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();
	cout << "> Adding edges... ";

	/* Adding edges */
	Eigen::MatrixXd transformedFields(F.rows(), F.cols());
	Eigen::VectorXd fields3D = A*arbField2D; 
	for (int i = 0; i < F.rows(); i++)
	{
		transformedFields.row(i) = (fields3D.block(3 * i, 0, 3, 1)).transpose();
	}

	viewer.data().add_edges(FC, FC + transformedFields*avgEdgeLength, Eigen::RowVector3d(0.0, 0.8, 0.1));


	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;

}

void VectorFields::testEnergyOfLocalPatch(igl::opengl::glfw::Viewer &viewer)
{
	// Identifying elements of local patch
	for (int i = 0; i < localPatchElements.size(); i++)
	{
		viewer.data().add_points(FC.row(localPatchElements[i]), Eigen::RowVector3d(0.0, 0.2, 0.9));
	}

	// Map from Global to Local
	vector<int> GlobToLocMap; 
	GlobToLocMap.resize(F.rows());
	for (int i = 0; i < F.rows(); i++)
	{
		GlobToLocMap[i] = -1;
	}

	int counter = 0;
	for (int face : localPatchElements) {
		GlobToLocMap[face] = counter;
		counter++;
	}

	// Construct Mass Matrix
	vector<Eigen::Triplet<double>> MTriplet;
	MTriplet.reserve(2 * localPatchElements.size());
	Eigen::SparseMatrix<double> MLocal(2 * localPatchElements.size(), 2 * localPatchElements.size());
	for (int i = 0; i < localPatchElements.size(); i++)
	{
		MTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, doubleArea(localPatchElements[i]) / 2.0));
		MTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, doubleArea(localPatchElements[i]) / 2.0));
	}
	MLocal.setFromTriplets(MTriplet.begin(), MTriplet.end());

	// Obtain the Stiffness matrix on the patch
	vector<Eigen::Triplet<double>> SFTriplet;
	SFTriplet.reserve(2 * localPatchElements.size());
	Eigen::SparseMatrix<double> SFLocal(2 * localPatchElements.size(), 2 * localPatchElements.size());
	for (int i = 0; i < localPatchElements.size(); i++)
	{
		// Diagonal elements
		int li = localPatchElements[i];
		SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 0, SF2D.coeff(2 * li + 0, 2 * li + 0)));
		SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 0, 2 * i + 1, SF2D.coeff(2 * li + 0, 2 * li + 1)));
		SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 0, SF2D.coeff(2 * li + 1, 2 * li + 0)));
		SFTriplet.push_back(Eigen::Triplet<double>(2 * i + 1, 2 * i + 1, SF2D.coeff(2 * li + 1, 2 * li + 1)));

		// Non-diagonal elements
		for (int j : AdjMF2Ring[li]) {
			const int neigh = j;
			if (GlobToLocMap[neigh] >= 0) {
				int neighLoc = GlobToLocMap[neigh];
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

	SFLocal.setFromTriplets(SFTriplet.begin(), SFTriplet.end());
	visualizeSparseMatrixInMatlab(SFLocal);

	// Get the components of the harmonic energy on the patch
	Eigen::VectorXd localHarmEnergy(2 * localPatchElements.size());
	for (int i = 0; i < localPatchElements.size(); i++)
	{
		Eigen::Vector2d vectLoc;
		vectLoc(0) = eigFieldFull2D(2 * localPatchElements[i] + 0, 0);
		vectLoc(1) = eigFieldFull2D(2 * localPatchElements[i] + 1, 0);
		//vectLoc.normalize();
		localHarmEnergy(2 * i + 0) = vectLoc(0);
		localHarmEnergy(2 * i + 1) = vectLoc(1);
	}

	// Compute the dirichlet energy
	double energy = localHarmEnergy.transpose()*SFLocal*localHarmEnergy;
	printf("Local Energy = %.5f\n", energy);

	double totalEnergy = eigFieldFull2D.col(0).transpose() * SF2D * eigFieldFull2D.col(0);
	printf("Total Energy = %.5f\n", totalEnergy);
}

void VectorFields::testGradients()
{
	//arbField = A*arbField2D;

	arbFieldE3D.resize(E.rows());
	for (int i = 0; i < E.rows(); i++) {
		double v1 = arbField(E(i, 0));
		double v2 = arbField(E(i, 1));
		double v = 0.5*(v1 + v2);
		arbFieldE3D(i) = v; 
	}

	Eigen::VectorXd gV3D = GF3D*arbField;
	Eigen::VectorXd gE3D = GFStar3D*arbFieldE3D;

	Eigen::VectorXd diff = gV3D - gE3D;
	double l2norm = diff.transpose()*MF3D*diff;
	double ref = gV3D.transpose()*MF3D*gV3D;
	l2norm = sqrt(l2norm / ref);
	printf(">>>The l2norm diff is %.10f\n", l2norm);

	//cout << "grad Vertex \n " << gV3D.block(0, 0, 30, 1) << endl; 
	//cout << "grad Edge \n " << gE3D.block(0, 0, 30, 1) << endl;
}

void VectorFields::testRotation() {
	cout << "[ Testing the rotation matrix ]\n";
	cout << "[ 1 TEST ]\n";
	Eigen::VectorXd v = A * arbField2D;
	Eigen::VectorXd ajv = A.transpose()*J3D * v;
	Eigen::VectorXd jav = J*A.transpose()*v;

	Eigen::VectorXd diff = ajv - jav;
	double diff1 = diff.transpose() * MF2D * diff;
	double diff2 = jav.transpose() * MF2D * jav;
	double l2norm = sqrt(diff1 / diff2);
	printf("The difference of AT*J3D*v-J*AT is %.20f\n", l2norm);
	cout << "The difference of AT*J3D*v-J*AT is " << l2norm << endl;

	cout << "[ 2 TEST ]\n";
	Eigen::VectorXd v1 = arbField2D;
	jav = J3D * A * v1;
	ajv = A * J * v1; 
	diff = ajv - jav;
	diff1 = diff.transpose() * MF3D * diff;
	diff2 = jav.transpose() * MF3D * jav;
	l2norm = sqrt(diff1 / diff2);
	printf("The difference of AT*J3D*v-J*AT is %.20f\n", l2norm);
	cout << "The difference of AT*J3D*v-J*AT is " << l2norm << endl;


	/*
	//printf("J3D=%dx%d, arbField=%d\n", J3D.rows(), J3D.cols(), arbField.size());
	//Eigen::VectorXd v = GF3D * arbField;
	Eigen::VectorXd v = A * arbField2D;

	// Dot prodduct
	Eigen::VectorXd Jv = J3D * v;
	double dotProd = Jv.dot(v);
	printf("The dot product of Jv with v is %.8f\n", dotProd);

	Jv = J3D*Jv;
	double vTemp = (Jv + v).transpose()*MF3D*(Jv + v);
	printf("J(Jv)+v is %.8f\n", vTemp);

	Jv = J3D * J3D * Jv;
	Eigen::VectorXd diff = Jv - v;
	vTemp = diff.transpose()*MF3D*diff;
	printf("J*J*J*J*v - v is %.8f\n", vTemp);

	cout << "Local-coordinate Case \n";
	v = arbField2D;

	// Dot prodduct
	Jv = J * v;
	dotProd = Jv.dot(v);
	printf("The dot product of Jv with v is %.8f\n", dotProd);

	Jv = J*Jv;
	vTemp = (Jv + v).transpose()*MF2D*(Jv + v);
	printf("J(Jv)+v is %.8f\n", vTemp);

	Jv = J * J * Jv;
	diff = Jv - v;
	vTemp = diff.transpose()*MF2D*diff;
	printf("J*J*J*J*v - v is %.8f\n", vTemp);
	*/


}

void VectorFields::testMassMatrix()
{
	cout << "Test Mass matrix \n";
	double refArea = 0.0;
	for (int i = 0; i < F.rows(); i++)
	{
		refArea += doubleArea(i) / 2.0; 
	}
	printf("The area of the surface is %.10f\n", refArea);

	// Identity matrix
	Eigen::VectorXd IdV(V.rows()); IdV.setConstant(1.0);
	Eigen::VectorXd IdE(E.rows()); IdE.setConstant(1.0);
	Eigen::VectorXd IdF(3*F.rows()); IdF.setConstant(1.0);

	double vArea = IdV.transpose()*MV*IdV;
	printf("Area [VERTEX] = %.10f\n", vArea);
	printf("Dim: MF3D=%dx%d, IdF=%d\n", MF3D.rows(), MF3D.cols(), IdF.size());
	double fArea = IdF.transpose()*MF3D*IdF;
	printf("Area [FACE  ] = %.10f\n", fArea);
	printf("Dim: MStar=%dx%d, IdE=%d\n", MStar.rows(), MStar.cols(), IdE.size());
	double eArea = IdE.transpose()*MStar*IdE;
	printf("Area [EDGE  ] = %.10f\n", eArea);
}

void VectorFields::projectionMatrixTest()
{
	Eigen::SparseMatrix<double> Id(A.cols(), A.cols());
	Id.setIdentity();

	Eigen::SparseMatrix<double> AAT = A * A.transpose();
	Eigen::SparseMatrix<double> ATA = A.transpose()*A;

	WriteSparseMatrixToMatlab(AAT, "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/AAT");
	WriteSparseMatrixToMatlab(ATA, "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/ATA");

}

void VectorFields::testCurlAndDiv()
{
	Eigen::SparseMatrix<double> Curl = MStarInv*GFStar3D.transpose()*J3D*MF3D;
	Eigen::SparseMatrix<double> Div = -MStarInv*GFStar3D.transpose()*MF3D;

	Curl = Curl * A;
	Div = Div * A;
	Eigen::SparseMatrix<double> DivJ = Div*J;

	//WriteSparseMatrixToMatlab(Curl, "D:/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/Curl");
	//WriteSparseMatrixToMatlab(DivJ, "D:/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/DivJ");

	Eigen::SparseMatrix<double> LapC = -MF3D * J3D * GFStar3D*MStarInv*GFStar3D.transpose()*J3D*MF3D;
	Eigen::SparseMatrix<double> LapD = MF3D * GFStar3D * MStarInv * GFStar3D.transpose()*MF3D;

	//LapC = A.transpose() * LapC * A;
	//LapD = A.transpose() * LapD * A;	

	Eigen::SparseMatrix<double> JLCJ = J * LapC * J;

	//WriteSparseMatrixToMatlab(LapD, "D:/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/LapD");
	//WriteSparseMatrixToMatlab(JLCJ, "D:/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/JLapCJ");
}

void VectorFields::perturbVectorFields(Eigen::VectorXd& inputFields)
{
	//const int NUM_PERTS = 50000;
	const int NUM_PERTS = 10000;
	set<int> perturbedFaces;

	/* Random number generator */
	std::random_device rd;								// Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd());								// Standard mersenne_twister_engine seeded with rd()

	do
	{
		std::uniform_int_distribution<> dist(0, F.rows()); // From 0 to F.rows()-1
		perturbedFaces.insert(dist(gen));
	} while (perturbedFaces.size() < NUM_PERTS);

	for (int i : perturbedFaces)
	{
		std::uniform_int_distribution<> dist(10,100);
		Eigen::Vector2d a;
		a(0) = (double)dist(gen) / 100;
		a(1) = (double)dist(gen) / 100;
		a.normalize();

		double norm_init = inputFields.block(2 * i, 0, 2, 1).norm();
		a = 0.5 * norm_init * a;
		inputFields.block(2 * i, 0, 2, 1) = inputFields.block(2 * i, 0, 2, 1) + a;
	}
}

void VectorFields::perturbVectorFieldsRegular(Eigen::VectorXd& inputFields)
{
	const int NUM_PERTS = (int)floor(0.05*F.rows());
	printf("___Num of perturbed fields: %d\n", NUM_PERTS);
	set<int> perturbedFaces;

	/* Random number generator */
	std::random_device rd;								// Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd());								// Standard mersenne_twister_engine seeded with rd()

	int counter = 0;
	do
	{		
		perturbedFaces.insert(counter+=20);
	} while (perturbedFaces.size() < NUM_PERTS);

	for (int i : perturbedFaces)
	{
		Eigen::Vector2d a;
		a(0) = (double) (i%50) / 100.0;
		a(1) = (double) ((i+25)%50) / 100.0;
		a.normalize();

		double norm_init = inputFields.block(2 * i, 0, 2, 1).norm();
		a = 0.5 * norm_init * a;
		inputFields.block(2 * i, 0, 2, 1) = inputFields.block(2 * i, 0, 2, 1) + a;
	}
}

void VectorFields::testSpectra()
{
	const int n = 100;
	Eigen::SparseMatrix<double> M(n, n);
	M.reserve(Eigen::VectorXi::Constant(n, 3));
	for (int i = 0; i < n; i++)
	{
		M.insert(i, i) = 10.0;
		if (i > 0)
			M.insert(i - 1, i) = 3.0;
		if (i < n - 1)
			M.insert(i + 1, i) = 2.0;
	}
	Eigen::VectorXd evalues;
	Eigen::MatrixXd evectors;

	///computeEigenSpectra(SF2DAsym, MF2D, 10, eigFieldFull2D, eigValuesFull, "hello");
}

void VectorFields::testSparseMatrix()
{
	//Eigen::SparseMatrix<double> MM;
	vector<int> matSize{100000, 20000000, 30000000, 300000, 400};

	cout << "Test on matrix size \n"; 
	for (int s : matSize) {
		Basis.resize(s, s);
		printf("Size=%dx%d \n", Basis.rows(), Basis.cols());
		//Basis.data().clear();
		Basis.data().squeeze();
		Basis.setIdentity();
	}
}