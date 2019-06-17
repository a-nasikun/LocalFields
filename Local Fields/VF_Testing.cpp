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
void VectorFields::testProjection_MyBasis_NoRegularizer(const Eigen::SparseMatrix<double>& Basis, const Eigen::VectorXd& inputFields, double &error)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> [Testing basis 2D of arbitrary field without regularizer...]\n";

	// Construct matrices for Test
	cout << "____Assigning variables\n";
	Eigen::SparseMatrix<double> U = Basis;// BasisTemp;
	Eigen::VectorXd				v = inputFields; 
	//Eigen::VectorXd				v =  arbField2D;
	Eigen::VectorXd				a = (U.transpose()*(MF2D*v));
	Eigen::SparseMatrix<double> B = U.transpose() * MF2D * U;

	cout << "____Solving linear system variables\n";
	Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> sparseSolver(B);
	//Eigen::PardisoLLT<Eigen::SparseMatrix<double>> sparseSolver(B);

	Eigen::VectorXd w = sparseSolver.solve(a);

	if (sparseSolver.info() != Eigen::Success) {
		cout << "Cannot solve the linear system. " << endl;
		if (sparseSolver.info() == Eigen::NumericalIssue)
			cout << "NUMERICAL ISSUE. " << endl;
		cout << sparseSolver.info() << endl;
		return;
	}

	//Eigen::VectorXd wb;
	wb.resize(U.rows());

	cout << "____Getting total SUM(wi*bi) \n";
	wb = U*w; 

	// Compare their L2-Norm
	cout << "____Computing L2-norm \n";
	Eigen::VectorXd diff = v - wb; 
	double length1 = wb.transpose()*MF2D*wb;
	double norm1 = diff.transpose()*MF2D*diff;
	double norm2 = v.transpose()*MF2D*v;
	double normL2 = sqrt(norm1 / norm2); 
	error = normL2; 

	cout << "____The L2 Norm is << " << normL2 << ": sqrt(" << norm1 << "/" << norm2 << ")" << endl;

	/* Measuring the 'length' of each vector */
	cout << "ENERGY ASYM \n";
	double harm_energy1 = v.transpose()*SF2DAsym*v;
	double harm_energy2 = wb.transpose()*SF2DAsym*wb;
	double harm_relEnergy = abs(harm_energy1 - harm_energy2) / harm_energy1;

	cout << "____Harmonic Energy => Ref=" << harm_energy1 << ", Approx:" << harm_energy2 << endl;
	cout << "____Relative energy: " << harm_relEnergy << endl;

	/* Measuring the 'length' of each vector */
	cout << "ENERGY SYM \n";
	double harm_energy1_sym = v.transpose()*SF2D*v;
	double harm_energy2_sym = wb.transpose()*SF2D*wb;
	double harm_relEnergy_sym = abs(harm_energy1_sym - harm_energy2_sym) / harm_energy1_sym;
	cout << "____Harmonic Energy => Ref=" << harm_energy1_sym << ", Approx:" << harm_energy2_sym << endl;
	cout << "____Relative energy: " << harm_relEnergy_sym << endl;
	
	/* Measuring the energy */	
	double biharm_energy1 = v.transpose()*B2DAsym*v;
	double biharm_energy2 = wb.transpose()*B2DAsym*wb;
	double biharm_relEnergy = abs(biharm_energy1 - biharm_energy2) / biharm_energy1;
	cout << "____Bi-Harmonic ENERGY => Ref=" << biharm_energy1 << ", Approx:" << biharm_energy2 << endl;
	cout << "____aRelative energy: " << biharm_relEnergy << endl;

	/* Measuring the 'length' of each vector */
	cout << "ENERGY SYM \n";
	double biharm_energy1_sym = v.transpose()*B2D*v;
	double biharm_energy2_sym = wb.transpose()*B2D*wb;
	double biharm_relEnergy_sym = abs(biharm_energy1_sym - biharm_energy2_sym) / biharm_energy1_sym;
	cout << "____Harmonic Energy => Ref=" << biharm_energy1_sym << ", Approx:" << biharm_energy2_sym << endl;
	cout << "____Relative energy: " << biharm_relEnergy_sym << endl;


	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;

	bool writeToFile = true;
	if (writeToFile) {
		std::ofstream ofs;
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Kitten_L2projection_eigenFields_500.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Kitten_L2projection_eigenFields_1000.txt";
		string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Kitten_L2projection_eigenFields_" + to_string(Basis.cols()) + ".txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Kitten_L2projection_eigenFields_5000.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Kitten_L2projection_eigenFields_10000.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Kitten_L2projection_eigenFields_20000.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Kitten_L2projection_eigenFields_50000.txt";

		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Armadillo_L2projection_eigenFields.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Armadillo_L2projection_optAlg.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Armadillo_L2projection_eigFields10.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Armadillo_L2projection_eigPatch.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Armadillo_L2projection_grad.txt";

		ofs.open(resultFile, std::ofstream::out | std::ofstream::app);

		ofs << globalConstraints.size() << "\t"
			<< harm_energy1 << "\t" << harm_energy1_sym << "\t" << biharm_energy1 << "\t" << biharm_energy1_sym << "\t" << norm2 << "\t"	// reference (harmAsym, harmSym, biharmAsym, biharmSym)
			<< harm_energy2 << "\t" << harm_relEnergy << "\t" << harm_energy2_sym << "\t" << harm_relEnergy_sym << "\t"						// approx (harmAsym, harmRelEnergyAsym,harmSym, harmRelEnergySym)
			<< biharm_energy2 << "\t" << biharm_relEnergy << "\t" << biharm_energy2_sym << "\t" << biharm_relEnergy_sym << "\t"				// approx (biharmAsym, biharmRelEnergyAsym, biharmSym, biharmRelEnergySym)
			<< length1 << "\t" << normL2 << "\n";																												// l2-norm

		ofs.close();
	}
}

void VectorFields::testProjection_EigenBasis_NoRegularizer(const Eigen::MatrixXd& Basis, const Eigen::VectorXd& inputFields, double &error)
{
	// For Timing
	chrono::high_resolution_clock::time_point	t0, t1, t2;
	chrono::duration<double>					duration;
	t0 = chrono::high_resolution_clock::now();
	cout << "> [Testing EIGENBASIS without regularizer...]\n";

	// Construct matrices for Test
	cout << "____Assigning variables\n";
	Eigen::MatrixXd U = Basis;// BasisTemp;
	Eigen::VectorXd				v = inputFields;
	//Eigen::VectorXd				v =  arbField2D;
	Eigen::VectorXd				a = (U.transpose()*(MF2D*v));
	Eigen::MatrixXd				B = U.transpose() * MF2D * U;

	cout << "____Solving linear system variables\n";
	Eigen::VectorXd w = B.ldlt().solve(a);
		
	cout << "____Getting total SUM(wi*bi) \n";
	wbEigen = U*w;

	//// Compare their L2-Norm
	//cout << "____Computing L2-norm \n";
	//Eigen::VectorXd diff = v - wbEigen;
	//double norm1 = diff.transpose()*MF2D*diff;
	//double norm2 = v.transpose()*MF2D*v;
	//double normL2 = sqrt(norm1 / norm2);
	//error = normL2;
	//
	//cout << "____The L2 Norm is << " << normL2 << ": sqrt(" << norm1 << "/" << norm2 << ")" << endl;
	//
	///* Measuring the energy */
	//double energy1 = v.transpose()*B2DAsym*v;
	//double energy2 = wbEigen.transpose()*B2DAsym*wbEigen;
	//cout << "____Bi-Harmonic ENERGY => Ref=" << energy1 << ", Approx:" << energy2 << endl;
	//cout << "____Relative energy: " << abs(energy1 - energy2) / energy1 << endl;
	//
	///* Measuring the 'length' of each vector */
	//energy1 = v.transpose()*SF2DAsym*v;
	//energy2 = wbEigen.transpose()*SF2DAsym*wbEigen;
	//cout << "____Harmonic Energy => Ref=" << energy1 << ", Approx:" << energy2 << endl;
	//cout << "____Relative energy: " << abs(energy1 - energy2) / energy1 << endl;
	//
	//t2 = chrono::high_resolution_clock::now();
	//duration = t2 - t0;
	//cout << "in " << duration.count() << " seconds." << endl;

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
	error = normL2;
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

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;

	bool writeToFile = true;
	if (writeToFile) {
		std::ofstream ofs;
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Kitten_L2projectionWithReg_eigenFields_500.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Kitten_L2projectionWithReg_eigenFields_1000.txt";
		string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Kitten_L2projectionWithReg_eigenFields" + to_string(Basis.cols()) + ".txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Kitten_L2projectionWithReg_eigenFields_5000.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Kitten_L2projectionWithReg_eigenFields_10000.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Kitten_L2projectionWithReg_eigenFields_20000.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Kitten_L2projectionWithReg_eigenFields_50000.txt";
		

		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Armadillo_L2projectionWithReg_eigFields.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Armadillo_L2projectionWithReg_optAlg.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Armadillo_L2projectionWithReg_eigFields10.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Armadillo_L2projectionWithReg_eigPatch.txt";
		//string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Armadillo_L2projectionWithReg_grad.txt";

		ofs.open(resultFile, std::ofstream::out | std::ofstream::app);

		ofs << globalConstraints.size() << "\t" << lambda * en2 / en1 << "\t" << lambda << "\t"
			<< harm_energy_input << "\t" << harm_energy_input_sym << "\t" << biharm_energy_input << "\t" << biharm_energy_input_sym << "\t" << length_init << "\t"	// initial input
			<< harm_energy1 << "\t" << harm_energy1_sym << "\t" << biharm_energy1 << "\t" << biharm_energy1_sym << "\t" << norm2 << "\t"	// reference (harmAsym, harmSym, biharmAsym, biharmSym)
			<< harm_energy2 << "\t" << harm_relEnergy << "\t" << harm_energy2_sym << "\t" << harm_relEnergy_sym << "\t"						// approx (harmAsym, harmRelEnergyAsym,harmSym, harmRelEnergySym)
			<< biharm_energy2 << "\t" << biharm_relEnergy << "\t" << biharm_energy2_sym << "\t" << biharm_relEnergy_sym << "\t"				// approx (biharmAsym, biharmRelEnergyAsym, biharmSym, biharmRelEnergySym)
			<< length_init << "\t" << length1 << "\t" << length2 << "\t" << length2/length1 * 100.0 << "\t"
			<< normL2 << "\n";																												// l2-norm

		ofs.close();
	}
}

void VectorFields::testProjection_EigenBasis_WithRegularizer(const Eigen::MatrixXd& Basis, const Eigen::VectorXd& inputFields, const Eigen::SparseMatrix<double>& MReg, double &error)
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

	Eigen::MatrixXd				U = Basis;
	Eigen::VectorXd				v = inputFields;
	
	/* Perturbed the input fields*/
	perturbVectorFields(v);
	pertFields = v;

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


	//// Compare their L2-Norm
	//cout << "____Computing L2-norm \n";
	//Eigen::VectorXd diff = wRef - wbEigen;
	//double norm1 = diff.transpose()*MF2D*diff;
	//double norm2 = wRef.transpose()*MF2D*wRef;
	//double normL2 = sqrt(norm1 / norm2);
	//error = normL2;
	//cout << "The L-2 Norm is << " << normL2 << endl;
	//
	///* Measuring the energy */
	//double energy1 = wRef.transpose()*B2DAsym*wRef;
	//double energy2 = wbEigen.transpose()*B2DAsym*wbEigen;
	//cout << "BIHARMONIC ENERGY => Ref=" << energy1 << ", Approx:" << energy2 << ",initial: " << v.transpose()*B2D*v << endl;
	//cout << "Relative energy: " << abs(energy1 - energy2) / energy1 << endl;
	//
	///* Measuring the energy */
	//energy1 = wRef.transpose()*SF2DAsym*wRef;
	//energy2 = wbEigen.transpose()*SF2DAsym*wbEigen;
	//cout << "HARMONIC ENERGY => Ref=" << energy1 << ", Approx:" << energy2 << ",initial: " << v.transpose()*SF2D*v << endl;
	//cout << "Relative energy: " << abs(energy1 - energy2) / energy1 << endl;
	//
	//
	///* Measuring the 'length' of each vector */
	//double length1 = wRef.transpose()*MF2D*wRef;
	//double length2 = wbEigen.transpose()*MF2D*wbEigen;
	//cout << "Length => Ref=" << length1 << ", Approx:" << length2 << ", initial: " << v.transpose()*MF2D*v << endl;
	//cout << "Relative length: " << length2 / length1 * 100.0 << "%" << endl;
	//
	//t2 = chrono::high_resolution_clock::now();
	//duration = t2 - t0;
	//cout << "in " << duration.count() << " seconds." << endl;


	// Compare their L2-Norm
	cout << "____Computing L2-norm \n";
	Eigen::VectorXd diff = wRef - wbEigen;
	double norm1 = diff.transpose()*MF2D*diff;
	double norm2 = wRef.transpose()*MF2D*wRef;
	double normL2 = sqrt(norm1 / norm2);
	error = normL2;
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

	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t0;
	cout << "in " << duration.count() << " seconds." << endl;

	bool writeToFile = true;
	if (writeToFile) {
		std::ofstream ofs;
		string resultFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Tests/Projections/Armadillo_L2projectionWithReg_modalBasis.txt";

		ofs.open(resultFile, std::ofstream::out | std::ofstream::app);

		ofs << globalConstraints.size() << "\t" << lambda * en2 / en1 << "\t" << lambda << "\t"
			<< harm_energy_input << "\t" << harm_energy_input_sym << "\t" << biharm_energy_input << "\t" << biharm_energy_input_sym << "\t" << length_init << "\t"	// initial input
			<< harm_energy1 << "\t" << harm_energy1_sym << "\t" << biharm_energy1 << "\t" << biharm_energy1_sym << "\t" << norm2 << "\t"	// reference (harmAsym, harmSym, biharmAsym, biharmSym)
			<< harm_energy2 << "\t" << harm_relEnergy << "\t" << harm_energy2_sym << "\t" << harm_relEnergy_sym << "\t"						// approx (harmAsym, harmRelEnergyAsym,harmSym, harmRelEnergySym)
			<< biharm_energy2 << "\t" << biharm_relEnergy << "\t" << biharm_energy2_sym << "\t" << biharm_relEnergy_sym << "\t"				// approx (biharmAsym, biharmRelEnergyAsym, biharmSym, biharmRelEnergySym)
			<< length_init << "\t" << length1 << "\t" << length2 << "\t" << length2 / length1 * 100.0 << "\t"
			<< normL2 << "\n";																												// l2-norm

		ofs.close();
	}
}

void VectorFields::projectionTest()
{
	cout << "PROJECTION TEST! \n";
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;	


	const int NUM_TEST = 1;
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
	string desFieldsFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/Kitten_DesignedFields";
	string pertFieldsFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/Kitten_PerturbedFields";
	if (readDesFieldsFromFile) {
		ReadDenseMatrixFromMatlab(DesignedFields, desFieldsFile, 2*F.rows(), NUM_TEST);
	}
	if (readPertFieldsFromFile) {
		ReadDenseMatrixFromMatlab(PerturbedFields, pertFieldsFile, 2 * F.rows(), NUM_TEST);
	}
	//Xf = arbField2D; 

	for (int i = 10; i < 10+NUM_TEST; i++)
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
			pertFields = PerturbedFields.col(i);
		}
		if(readPertFieldsFromFile){
			Xf = DesignedFields.col(i);
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
		/* Projection to the subspace */
		/* Reference results */
		//setupGlobalProblem(Eigen::Vector3d(1,1,1));
		testProjection_MyBasis_NoRegularizer(Basis, Xf, errors1(i));
		//testProjection_EigenBasis_NoRegularizer(EigenBasis, Xf, errors2(i));

		//testProjection_MyBasis_WithRegularizer(Basis, pertFields, B2DAsym, errors2(i));
		///testProjection_MyBasis_WithRegularizer(Basis, pertFields, SF2DAsym, errors2(i));
		
		//testProjection_EigenBasis_WithRegularizer(EigenBasis, pertFields, SF2DAsym, errors2(i));

		t2 = chrono::high_resolution_clock::now();
		duration = t2 - t1;
		//printf("[%d] run => Error=%.10f (in %.3f seconds) \n", i, errors1(i), duration.count());		
		printf("[%d] run => [My Basis] Error=%.10f\n", i, errors1(i));
		printf("            [EigenBasis] Error=%.10f (in %.3f seconds) \n", errors2(i), duration.count());
	}

	if (!readPertFieldsFromFile && !readDesFieldsFromFile)
	{
		WriteDenseMatrixToMatlab(DesignedFields, desFieldsFile);
		WriteDenseMatrixToMatlab(PerturbedFields, pertFieldsFile);
	}

	cout << "ERRORS: \n" <<  errors1 << endl << "ERRORS2 \n" << errors2 << endl; 
}

void VectorFields::projectionTest(bool &readDesFieldsFromFile, bool &readPertFieldsFromFile, int start, int nTests)
{
	cout << "PROJECTION TEST! \n";
	// For Timing
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;


	const int NUM_TEST = nTests;
	Eigen::VectorXd errors1(NUM_TEST), errors2(NUM_TEST);
	Eigen::MatrixXd DesignedFields(2 * F.rows(), NUM_TEST);
	Eigen::MatrixXd PerturbedFields(2 * F.rows(), NUM_TEST);

	bool useEigenBasis = false;

	/* Loading eigen basis for Armadillo */
	Eigen::MatrixXd EigenBasis;
	string eigBasisFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/Armadillo_1000_eigenfields_Ref_2";
	if (useEigenBasis)
	{
		ReadDenseMatrixFromMatlab(EigenBasis, eigBasisFile, 172964, 1000);
	}

	/* Fields and perturbed Fields */
	string desFieldsFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/Kitten_DesignedFields";
	string pertFieldsFile = "D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/VFields/Kitten_PerturbedFields";
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
		//setupGlobalProblem(Eigen::Vector3d(1,1,1));
		//testProjection_MyBasis_NoRegularizer(Basis, Xf, errors1(i));
		//testProjection_EigenBasis_NoRegularizer(EigenBasis, Xf, errors2(i));

		//testProjection_MyBasis_WithRegularizer(Basis, pertFields, B2DAsym, errors2(i));
		testProjection_MyBasis_WithRegularizer(Basis, pertFields, SF2DAsym, errors2(i));

		//testProjection_EigenBasis_WithRegularizer(EigenBasis, pertFields, SF2DAsym, errors2(i));

		t2 = chrono::high_resolution_clock::now();
		duration = t2 - t1;
		//printf("[%d] run => Error=%.10f (in %.3f seconds) \n", i, errors1(i), duration.count());		
		printf("[%d] run => [My Basis] Error=%.10f\n", i, errors1(i));
		printf("            [EigenBasis] Error=%.10f (in %.3f seconds) \n", errors2(i), duration.count());
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

void VectorFields::convergenceTest()
{
	const int NUM_DIFF_BASIS = 6;
	vector<std::string> basisFile;
	basisFile.reserve(NUM_DIFF_BASIS);

	basisFile.push_back("D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Kitten_500_EigFields_35sup");
	basisFile.push_back("D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Kitten_1000_EigFields_35sup");
	basisFile.push_back("D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Kitten_2000_EigFields_35sup");
	basisFile.push_back("D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Kitten_5000_EigFields_35sup");
	basisFile.push_back("D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Kitten_10000_EigFields_35sup");
	basisFile.push_back("D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Data/Basis/Basis_Kitten_50000_EigFields_35sup");

	/* For the projection tests */
	bool readDesFields = true;
	bool readPertFields = true;
	int idStart = 0;
	int NUM_TEST = 11;

	for (string file_ : basisFile)
	{
		retrieveBasis(file_);
		projectionTest(readDesFields, readPertFields, idStart, NUM_TEST);
	}
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
	l2norm = l2norm / ref;
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
	const int NUM_PERTS = 5000;
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