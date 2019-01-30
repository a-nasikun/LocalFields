#include "Utility.h"
#include <fstream>

double SparseMatrixMaxValue(const Eigen::SparseMatrix<double> &M)
{
	double maxVal = 0.0;
	for (int k = 0; k < M.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(M, k); it; ++it) {
			if (it.value() > maxVal) maxVal = it.value(); 
		}
	}

	return maxVal;
}

void WriteDenseMatrixToMatlab(const Eigen::MatrixXd& M, const string& filename)
{
	//Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");
	//ofstream file(filename.c_str());
	//file << M.format(CSVFormat);

	// Create a 'Fake Sparse' Matrix
	vector<Eigen::Triplet<double>> UTriplet;
	Eigen::SparseMatrix<double> MSparse;
	UTriplet.reserve(M.rows() * M.cols());
	MSparse.resize(M.rows(), M.cols());
	MSparse.reserve(M.rows() * M.cols());

	for (int i = 0; i < M.cols(); i++) {
		for (int j = 0; j < M.rows(); j++) {
			UTriplet.push_back(Eigen::Triplet<double>(j, i, M.coeff(j, i)));
		}
	}

	MSparse.setFromTriplets(UTriplet.begin(), UTriplet.end());
	printf("MSparse: %dx%d with %d nonzeros\n", MSparse.rows(), MSparse.cols(), MSparse.nonZeros());
	//cout << MSparse.block(300, 100, 100, 10) << endl; 

	// Writing the Sparse Matrix to matlab
	using namespace matlab::engine;
	Engine *ep;
	mxArray *MM = NULL, *MS = NULL, *result = NULL, *eigVecResult, *nEigs;

	const int NNZ_M = MSparse.nonZeros();
	int nnzMCounter = 0;

	double	*srm = (double*)malloc(NNZ_M * sizeof(double));
	mwIndex *irm = (mwIndex*)malloc(NNZ_M * sizeof(mwIndex));
	mwIndex *jcm = (mwIndex*)malloc((MSparse.cols() + 1) * sizeof(mwIndex));

	MM = mxCreateSparse(MSparse.rows(), MSparse.cols(), NNZ_M, mxREAL);
	srm = mxGetPr(MM);
	irm = mxGetIr(MM);
	jcm = mxGetJc(MM);

	// Getting matrix M
	jcm[0] = nnzMCounter;
	for (int i = 0; i < MSparse.outerSize(); i++) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(MSparse, i); it; ++it) {
			srm[nnzMCounter] = it.value();
			irm[nnzMCounter] = it.row();
			nnzMCounter++;
		}
		jcm[i + 1] = nnzMCounter;
	}

	// Start Matlab Engine
	ep = engOpen(NULL);
	if (!(ep = engOpen(""))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		cout << "CANNOT START MATLAB " << endl;
	}
	else {
		cout << "MATLAB STARTS. OH YEAH!!!" << endl;
	}

	engPutVariable(ep, "ApproxFields", MM);
	engEvalString(ep, "ApproxFields=full(ApproxFields);");
	engEvalString(ep, "save('F:/PROGRAMMING/Localized Vector Fields/Build/ForChristopher/ChineseDragon/ApproxFields','ApproxFields');");

}

void WriteSparseMatrixToMatlab(const Eigen::SparseMatrix<double>& M, const string& filename)
{
	printf("Size of M=%dx%d\n", M.rows(), M.cols());

	using namespace matlab::engine;
	Engine *ep;
	mxArray *MM = NULL, *MS = NULL, *result = NULL, *eigVecResult, *nEigs;

	const int NNZ_M = M.nonZeros();
	int nnzMCounter = 0;

	double	*srm = (double*)malloc(NNZ_M * sizeof(double));
	mwIndex *irm = (mwIndex*)malloc(NNZ_M * sizeof(mwIndex));
	mwIndex *jcm = (mwIndex*)malloc((M.cols() + 1) * sizeof(mwIndex));

	MM = mxCreateSparse(M.rows(), M.cols(), NNZ_M, mxREAL);
	srm = mxGetPr(MM);
	irm = mxGetIr(MM);
	jcm = mxGetJc(MM);

	// Getting matrix M
	jcm[0] = nnzMCounter;
	for (int i = 0; i < M.outerSize(); i++) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it) {
			srm[nnzMCounter] = it.value();
			irm[nnzMCounter] = it.row();
			nnzMCounter++;
		}
		jcm[i + 1] = nnzMCounter;
	}

	// Start Matlab Engine
	ep = engOpen(NULL);
	if (!(ep = engOpen(""))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		cout << "CANNOT START MATLAB " << endl;
	}
	else {
		cout << "MATLAB STARTS. OH YEAH!!!" << endl;
	}

	engPutVariable(ep, "Basis", MM);	
	// LAPTOP
	//engEvalString(ep, "save('F:/PROGRAMMING/Localized Vector Fields/Build/ForChristopher/Armadillo/Basis','Basis');");

	// WORKSTATION
	engEvalString(ep, "save('E:/Local Programming/Localized Bases for Vector Fields/LocalFields_build/ForChristopher/Armadillo/Basis','Basis');");

}

//template<typename Scalar>
void manuallyDestroySparseMatrix(Eigen::SparseMatrix<double> &M)
{
	M.~SparseMatrix();
	Eigen::SparseMatrix<double> K;
}

