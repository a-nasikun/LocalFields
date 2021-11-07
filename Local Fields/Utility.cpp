#include "Utility.h"
#include <fstream>
#include <cusparse_v2.h>

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

//<<<<<<< HEAD
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

	engPutVariable(ep, "data", MM);
	engEvalString(ep, "data=full(data);");
	string saveFile = "save('" + filename + "','data');";
	engEvalString(ep, saveFile.c_str());

	//engEvalString(ep, "save('F:/PROGRAMMING/Localized Vector Fields/Build/ForChristopher/ChineseDragon/ApproxFields','vectorFields');");

}

void WriteSparseMatrixToMatlab(const Eigen::SparseMatrix<double>& M, const string& filename)
{
	printf("Saving matrix with the size of M=%dx%d\n", M.rows(), M.cols());

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

	engPutVariable(ep, "MF", MM);	
	// LAPTOP
	//engEvalString(ep, "save('F:/PROGRAMMING/Localized Vector Fields/Build/ForChristopher/Armadillo/Basis','Basis');");
	//engEvalString(ep, "save('E:/Local Programming/Localized Bases for Vector Fields/LocalFields_build/ForChristopher/Armadillo/Basis','Basis');");

	// WORKSTATION
	//engEvalString(ep, "save('E:/Local Programming/Localized Bases for Vector Fields/LocalFields_build/ForChristopher/Armadillo/Basis','Basis');");
	string saveFile = "save('" + filename +"','MF');";
	engEvalString(ep, saveFile.c_str());

}

void ReadDenseMatrixFromMatlab(Eigen::MatrixXd& M, const string& filename, const int& nRows, const int& nCols, const int& nBlocks)
{
	using namespace matlab::engine;
	Engine *ep;
	mxArray *eigValM, *eigVecM;		// for Matlab
	double	*eigValE, *eigVecE;		// for Eigen
	const int NUM_EIGEN = nCols;
	const int NUM_ROWS = nRows;
	const int NUM_BLOCKS = nBlocks;

	M.resize(NUM_ROWS, NUM_BLOCKS*NUM_EIGEN);

	// Start Matlab Engine
	ep = engOpen(NULL);
	if (!(ep = engOpen(""))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		cout << "CANNOT START MATLAB " << endl;
	}
	else {
		cout << "MATLAB STARTS. OH YEAH!!!" << endl;
	}
	
	cout << "Loading Matrix" << endl;
	string location = "load('" + filename + ".mat');";
	engEvalString(ep, location.c_str());
	
	// First 2 blocks
	cout << "Retrieving the Matrix of " << nRows <<" rows, and " << nCols << "columns. " << endl; 
	//engEvalString(ep, "REVec = EigVec;");
	//eigVecM = engGetVariable(ep, "vectorFields");
	eigVecM = engGetVariable(ep, "data");
	//eigVecM = engGetVariable(ep, "BasisFull");
		cout << "Variable obtained from matlab \n";
	//eigVecM = engGetVariable(ep, "EigVec");
	//eigVecM = engGetVariable(ep, "BasisFull");
		cout << "Allocating memory in C++.\n ";
	eigVecE = (double*)malloc(NUM_ROWS * NUM_BLOCKS*NUM_EIGEN * sizeof(double));
		cout << "Copying data to C++ data \n";
	memcpy((void *)eigVecE, (void *)mxGetPr(eigVecM), NUM_ROWS * NUM_BLOCKS*NUM_EIGEN * sizeof(double));

	cout << "Converting the Matrix to Eigen format" << endl;
	//Eigen::MatrixXd REigCont;
	M.resize(NUM_ROWS, NUM_BLOCKS* NUM_EIGEN);
	M = Eigen::Map<Eigen::MatrixXd, 0, Eigen::OuterStride<>>(eigVecE, NUM_ROWS, NUM_BLOCKS*NUM_EIGEN, Eigen::OuterStride<>(NUM_ROWS));
	printf("Size of Matrix=%dx%d\n", M.rows(), M.cols());
	//M.block(0, 0, NUM_ROWS, NUM_BLOCKS * NUM_EIGEN) = REigCont;
	
	///engEvalString(ep, "clear;");
	///
	///mxDestroyArray(eigValM);
	///mxDestroyArray(eigVecM);
	///realloc(eigVecE, 0);
	///realloc(eigValE, 0);
	///
	///engClose(ep);
	cout << "FInish loading matrix \n";
}

void ReadDenseMatrixFromMatlab(Eigen::MatrixXd& M, const string& filename)
{
	//ReadDenseMatrixFromMatlab(M, filename, 500, 172946);
	ReadDenseMatrixFromMatlab(M, filename, 1000, 43243);
}

void ReadSparseMatrixFromMatlab(Eigen::SparseMatrix<double>& M, const string& filename)
{
	/// [TO DO]
}

void ReadVectorFromMatlab(Eigen::VectorXd& v, const string& filename, const int nRows)
{
	using namespace matlab::engine;
	Engine *ep;
	mxArray *eigValM, *eigVecM;		// for Matlab
	double	*eigValE, *eigVecE;		// for Eigen
	const int NUM_EIGEN = nRows;		// num rows
	const int NUM_BLOCKS = 1;

	v.resize(NUM_BLOCKS*NUM_EIGEN);	

	// Start Matlab Engine
	ep = engOpen(NULL);
	if (!(ep = engOpen(""))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		cout << "CANNOT START MATLAB " << endl;
	}
	else {
		cout << "MATLAB STARTS. OH YEAH!!!" << endl;
	}

	cout << "Loading Matrix" << endl;
	string location = "load('" + filename + ".mat');";
	engEvalString(ep, location.c_str());
	//engEvalString(ep, "load('D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/CDragon_500_REigVal.mat');");

	// Get the EIGENVALUES from Matlab => C++
	cout << "Retrieving the eigenvalues" << endl;
	eigValM = engGetVariable(ep, "data");
	eigValE = (double*)malloc(NUM_BLOCKS*NUM_EIGEN * sizeof(double));
	memcpy((void *)eigValE, (void *)mxGetPr(eigValM), NUM_BLOCKS*NUM_EIGEN * sizeof(double));
	v = Eigen::Map<Eigen::VectorXd>(eigValE, NUM_EIGEN);
	//cout << "Eigenvalues: \n" << v << endl; 
	cout << "Reading done\n";
		
	engEvalString(ep, "clear;");
	engClose(ep);
}

void visualizeSparseMatrixInMatlab(const Eigen::SparseMatrix<double> &M)
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

	engPutVariable(ep, "M", MM);
	engEvalString(ep, "spy(M)");
}


void ReadChristopherStiffnessMatrix(const string &filename, Eigen::SparseMatrix<double> &M)
{
	ifstream file(filename);
	string oneLine, oneWord;
	int i, j;
	double v;
	int counter = 0;

	vector<Eigen::Triplet<double>> MTriplet;
	MTriplet.reserve(6000 * 20);

	if (file.is_open())
	{
		getline(file, oneLine);
		while (getline(file, oneLine))
		{
			istringstream iStream(oneLine);
			getline(iStream, oneWord, ' ');
			i = stoi(oneWord);
			getline(iStream, oneWord, ' ');
			j = stoi(oneWord);
			getline(iStream, oneWord, ' ');
			v = stod(oneWord);

			MTriplet.push_back(Eigen::Triplet<double>(i, j, v));

		}
	}
	file.close();

	printf("Size=%dx%d, with %d elements\n", i + 1, i + 1, MTriplet.size());
	M.resize(i + 1, i + 1);
	M.setFromTriplets(MTriplet.begin(), MTriplet.end());
}

//=======
void LoadSparseMatrixFromTxtFile(const string& filename, Eigen::SparseMatrix<double> &M)
{
	ifstream file(filename);
	string oneLine, oneWord;
	int i=0, j;
	double v;
	int counter;

	vector<Eigen::Triplet<double>> MTriplet;
	MTriplet.reserve(6000 * 20);

	if (file.is_open())
	{
		/* Get the size of the matrix */
		getline(file, oneLine);		
		istringstream iStream1(oneLine);
		getline(iStream1, oneWord, ' ');
		const int m = stoi(oneWord);
		getline(iStream1, oneWord, ' ');
		const int n = stoi(oneWord);
		printf("Size=%d x %d ", m, n);
		M.resize(m,n);

		/* Obtain each member elements */
		while (getline(file, oneLine))
		{
			counter = 0; 

			istringstream iStream(oneLine);

			for (string word; iStream >> word; )
			{
				//cout << "line: " << word <<  endl; 
				if (counter == 0)
				{
					//i = stoi(word);
				}
				else if (counter % 2 == 1)
				{
					j = stoi(word);
				}
				else if (counter % 2 == 0)
				{
					v = stod(word);
					MTriplet.push_back(Eigen::Triplet<double>(i, j, v));
					//printf("______[%d, %d]=%.5f\n", i, j, v);
				}
				counter++;
			}
			i++;
		}
	}
	file.close();
	printf(", with %d elements\n", MTriplet.size());	
	M.setFromTriplets(MTriplet.begin(), MTriplet.end());
}

/* Obtain the inverse of Mass Matrix M, to get MInv*/
void ConstructInverseMassMatrix(Eigen::SparseMatrix<double> &M, Eigen::SparseMatrix<double> &MInv)
{
	MInv.resize(M.rows(), M.cols());
	vector<Eigen::Triplet<double>> MTriplet;
	MTriplet.reserve(M.rows());

	/* Getting the sum of every non-zero elements in a row */
	for (int k = 0; k < M.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(M, k); it; ++it) {
			MTriplet.push_back(Eigen::Triplet<double>(it.row(), it.col(), 1.0 / it.value()));
		}
	}
	MInv.setFromTriplets(MTriplet.begin(), MTriplet.end());
}

void WriteSTDVectorToTxtFile(const vector<int>& vector, const string& filename)
{
	ofstream myfile(filename.c_str());
	if (myfile.is_open())
	{
		cout << "Write file to text \n";
		for (int i = 0; i < vector.size(); i++)
		{
			myfile << vector[i] << "\n";
		}
		myfile.close();
	}
	else cout << "Unable to open file";
}

void LoadSTDVectorFromTxtFile(const string& filename, vector<int>& vector) 
{
	vector.reserve(100000);

	string line;
	ifstream myfile(filename.c_str());
	if (myfile.is_open())
	{
		while (getline(myfile, line))
		{
			//cout << "Inserting " << line << " to vector \n";
			vector.push_back(stoi(line));
		}
		myfile.close();
	}

	else cout << "Unable to open file";
}

void WriteEigenVectorToTxtFile(const Eigen::VectorXd& vector, const string& filename)
{
	ofstream myfile(filename.c_str());
	if (myfile.is_open())
	{
		cout << "Write file to text \n";
		for (int i = 0; i < vector.size(); i++)
		{
			myfile << vector(i) << "\n";
		}
		myfile.close();
	}
	else cout << "Unable to open file";
}

void LoadEigenVectorFromTxtFile(const string& filename, Eigen::VectorXd& vector)
{
	std::vector<double> v2;
	v2.reserve(500000);

	string line;
	ifstream myfile(filename.c_str());
	if (myfile.is_open())
	{
		while (getline(myfile, line))
		{
			//cout << "Inserting " << line << " to vector \n";
			v2.push_back(stod(line));
		}
		myfile.close();
		v2.shrink_to_fit();
		//printf("std:vector size=%d\n", (int)v2.size());
		//for (int i = 0; i < 10; i++) printf("%.4f\n", v2[i]);

		vector.resize(v2.size());
		//cout << "Getting the elements from vector\n";
		vector = Eigen::Map<Eigen::VectorXd>(v2.data(), (int)v2.size());
		//vector = 5 * vector; 
		//printf("eigen::vector size=%d\n", vector.rows());
		//cout << vector.block(0, 0, 10, 1) << endl; 
		//for (int k=0; k<v2.size(); k++)
		//{
		//	vector(k) = v2[k];
		//}
	}

	 

	else cout << "Unable to open file";
}

// Write matrix to file
void writeEigenSparseMatrixToBinary(Eigen::SparseMatrix<double> &m, const std::string &filename)
{
	printf("Tries to write matrix\n");
	int matSize = m.nonZeros();
	m.makeCompressed();

	fstream writeFile;
	writeFile.open(filename, ios::binary | ios::out);

	if (writeFile.is_open())
	{
		int rows, cols, nnzs, outS, innS;
		rows = m.rows();
		cols = m.cols();
		nnzs = m.nonZeros();
		outS = m.outerSize();
		innS = m.innerSize();

		writeFile.write((const char *)&(rows), sizeof(int));
		writeFile.write((const char *)&(cols), sizeof(int));
		writeFile.write((const char *)&(nnzs), sizeof(int));
		writeFile.write((const char *)&(outS), sizeof(int));
		writeFile.write((const char *)&(innS), sizeof(int));

		writeFile.write((const char*)(m.valuePtr()), sizeof(double)*m.nonZeros());
		writeFile.write((const char*)(m.outerIndexPtr()), sizeof(int)*m.outerSize());
		writeFile.write((const char*)(m.innerIndexPtr()), sizeof(int)*m.nonZeros());

		printf("Writing is over\n");
		writeFile.close();
	}
	printf("File is closed\n");
}

// Read matrix from file
void readEigenSparseMatrixFromBinary(const std::string &filename, Eigen::SparseMatrix<double> &m)
{
	//m.resize(0, 0);

	printf("Tries to read file \n");
	fstream readFile;
	readFile.open(filename, ios::binary | ios::in);
	if (readFile.is_open())
	{
		int rows, cols, nnzs, outS, innS;
		readFile.read((char*)&rows, sizeof(int));
		readFile.read((char*)&cols, sizeof(int));
		readFile.read((char*)&nnzs, sizeof(int));
		readFile.read((char*)&outS, sizeof(int));
		readFile.read((char*)&innS, sizeof(int));

		m.resize(rows, cols);
		m.makeCompressed();
		m.resizeNonZeros(nnzs);

		readFile.read((char*)(m.valuePtr()), sizeof(double)*nnzs);
		readFile.read((char*)(m.outerIndexPtr()), sizeof(int)*outS);
		readFile.read((char*)(m.innerIndexPtr()), sizeof(int)*nnzs);

		printf("Reading is over\n");

		m.finalize();
		readFile.close();
	}
	else {
		cout << "ERROR! Cannot read the basis from file \n";
	}
	printf("File is closed\n");
}

void writeEigenDenseMatrixToBinary(Eigen::MatrixXd M, const std::string &filename)
{
	//std::ofstream out(filename, std::ios::out | std::ios:: | std::ios::trunc);
	std::ofstream out(filename, std::ios::out | std::ios::trunc);
	int rows = M.rows(), cols = M.cols();
	out.write((char*)(&rows), sizeof(int));
	out.write((char*)(&cols), sizeof(int));
	out.write((char*)M.data(), rows*cols * sizeof(typename double));
	out.close();
}

void readEigenDenseMatrixFromBinary(const std::string &filename, Eigen::MatrixXd M)
{
	//std::ifstream in(filename, std::ios::in | std::ios::binary);
	std::ifstream in(filename, std::ios::in);
	int rows = 0, cols = 0;
	in.read((char*)(&rows), sizeof(int));
	in.read((char*)(&cols), sizeof(int));
	M.resize(rows, cols);
	in.read((char *)M.data(), rows*cols * sizeof(double));
	in.close();
}

void SparseMatrix_Vector_Multiplication_CSR(Eigen::SparseMatrix<double, Eigen::RowMajor> &S, Eigen::VectorXd& a, Eigen::VectorXd& b)
{
	//{
	//	cout << "Toying around with sparse matrix";
	//	const int nRows = 10;
	//	const int nCols = 8;
	//	Eigen::SparseMatrix<double> K(nRows, nCols);
	//	vector<Eigen::Triplet<double>> KTriplet; KTriplet.reserve(20);
	//	srand(time(NULL));
	//	for (int i = 0; i < nRows; i++)
	//	{
	//		KTriplet.push_back(Eigen::Triplet<double>(i, rand() % nCols, 1.0));
	//		KTriplet.push_back(Eigen::Triplet<double>(i, rand() % nCols, 1.0));
	//	}
	//	K.setFromTriplets(KTriplet.begin(), KTriplet.end());
	//	cout << "K \n" << K << endl << endl; 
	//	K.makeCompressed();
	//	K = K.transpose();
	//
	//	int nnzK = K.nonZeros();
	//	int m = K.rows();
	//	int n = K.cols();
	//	printf("Rows=%d | cols=%d \n", m, n);
	//
	//	double* K_csrVal = (double*)malloc(nnzK * sizeof(double));
	//	int* K_csrRowPtr = (int*)malloc((m + 1) * sizeof(int));
	//	int* K_csrColInd = (int*)malloc(nnzK * sizeof(int));
	//	K_csrVal = K.valuePtr();
	//	K_csrRowPtr = K.outerIndexPtr();
	//	K_csrColInd = K.innerIndexPtr();
	//
	//	for (int i = 0; i < m + 1; i++) printf("%d |", K_csrRowPtr[i]); cout << endl;
	//	for (int i = 0; i < nnzK; i++) printf("%d |", K_csrColInd[i]); cout << endl;
	//	for (int i = 0; i < nnzK; i++) printf("%f |", K_csrVal[i]); cout << endl;
	//
	//}


	S.makeCompressed();
	//cout << "Setting up the handle in GPU \n";
	cudaError_t			cudaStat1 = cudaSuccess;
	cusparseHandle_t	handle;
	cusparseMatDescr_t descrA;
	cusparseCreateMatDescr(&descrA);
	cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO);
	cusparseCreate(&handle);
	
	S.makeCompressed();
	const int nnz = S.nonZeros();
	const int m = S.rows();
	const int n = S.cols();
	//printf("The data: S=%dx%d | nnzA:%d, rowsA:%d |\n", m, n, nnz, m+1);

	/* Setting up the data for the vectors */
	//cout << "Allocating memory and setup in CPU \n";
	// Populating the matrix in CPU
	double* h_csrVal = (double*)malloc(nnz * sizeof(double));
	int* h_csrRowPtr = (int*)malloc((m + 1) * sizeof(int));
	int* h_csrColInd = (int*)malloc(nnz * sizeof(int));
	h_csrVal = S.valuePtr();
	h_csrRowPtr = S.outerIndexPtr();
	h_csrColInd = S.innerIndexPtr();
	double* h_a = (double*)malloc(n * sizeof(double));
	h_a = a.data();
	double* h_b = (double*)malloc(m * sizeof(double));
	for (int i = 0; i < m; i++) h_b[i] = 0.5;
	for (int i = 0; i < nnz; i++) {
		//if(i%1000==0)
		//if(h_csrColInd[i]<0)
		//	printf("[%d]=>[%.3f]\n", h_csrColInd[i], h_csrVal[i]);
	}
	for (int i = 0; i < (m+1); i++) {
		//if(i%1000==0)
		//if(h_csrRowPtr[i]<0)
		//	printf("[%d/%d]\n", h_csrRowPtr[i], m+1);
	}

	//cout << "Allocating memory and setup in GPU \n";
	// Allocating memory in device/GPU
	double *d_a;  cudaStat1 = cudaMalloc(&d_a, n * sizeof(double));				 //cout << "__alloc_status:" << cudaStat1 << endl;
	double *d_b;  cudaStat1 = cudaMalloc(&d_b, m * sizeof(double));				 //cout << "__alloc_status:" << cudaStat1 << endl;
	int *d_csrColInd; cudaStat1 = cudaMalloc(&d_csrColInd, nnz * sizeof(int));   //cout << "__col:alloc_status:" << cudaStat1 << endl;
	int *d_csrRowPtr; cudaStat1 = cudaMalloc(&d_csrRowPtr, (m+1) * sizeof(int)); //cout << "__row:alloc_status:" << cudaStat1 << endl;
	double *d_csrVal; cudaStat1 = cudaMalloc(&d_csrVal, nnz * sizeof(double));   //cout << "__val:alloc_status:" << cudaStat1 << endl;

	// Getting the buffer size
	cusparseStatus_t cusparseStat1;
	
	//cout << "Copying to GPU \n";
	// Copying data to CUDA
	cudaStat1 = cudaMemcpy(d_a, h_a, n * sizeof(double), cudaMemcpyHostToDevice);					//cout << "__alloc_status:" << cudaStat1 << endl;
	cudaStat1 = cudaMemcpy(d_b, h_b, m * sizeof(double), cudaMemcpyHostToDevice);					//cout << "__alloc_status:" << cudaStat1 << endl;
	cudaStat1 = cudaMemcpy(d_csrRowPtr, h_csrRowPtr, (m+1) * sizeof(int), cudaMemcpyHostToDevice);  //cout << "__rows: status:" << cudaStat1 << endl;
	cudaStat1 = cudaMemcpy(d_csrColInd, h_csrColInd, nnz * sizeof(int), cudaMemcpyHostToDevice);    //cout << "__col: status:" << cudaStat1 << endl;
	cudaStat1 = cudaMemcpy(d_csrVal, h_csrVal, nnz * sizeof(double), cudaMemcpyHostToDevice);       //cout << "__val: status:" << cudaStat1 << endl;
	
	//cout << "Checking the size\n";
	//printf("__Size of h_a: %d \n", sizeof(h_a) / sizeof(double));
	//printf("__Size of h_b: %d \n", sizeof(h_b) / sizeof(double));
	//printf("__Size of h_csrRowPtr: %d \n", sizeof(h_csrRowPtr) / sizeof(int));
	//printf("__Size of h_csrColInd: %d \n", sizeof(h_csrColInd) / sizeof(int));
	//printf("__Size of h_csrVal: %d \n", sizeof(h_csrVal) / sizeof(double));
	//printf("__rows: %d | cols: %d | nnz: %d \n", m, n, nnz);
	//
	//cout << "Doing the real multiplication \n";
	
	// The multiciplication
	double alpha = 1.0; 
	double beta = 0.0;
	
	//cusparseStat1 = cusparseDcsrmv(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, m, n, nnz, &alpha, descrA, d_csrVal, d_csrRowPtr, d_csrColInd, d_a, &beta, d_b);
	//cusparseStatus_t 
	cusparseStat1 = cusparseDbsrmv(handle, CUSPARSE_DIRECTION_ROW, CUSPARSE_OPERATION_NON_TRANSPOSE, m, n, nnz, &alpha, descrA, d_csrVal, d_csrRowPtr, d_csrColInd, 1, d_a, &beta, d_b);

	//cout << "__status:" << cusparseStat1 << endl;
	
	//cout << "Copying to CPU \n";
	cudaMemcpy(h_b, d_b, m * sizeof(double), cudaMemcpyDeviceToHost);
	//for (int i = 0; i < 10; i++) printf("%d: %.4f \n", i, h_b[i]);
	b = Eigen::Map<Eigen::VectorXd>(h_b, m);
	//cout << b.block(0, 0, 10, 1) << endl;

}

/* Matrix multiplication in Cuda b = S*a */
void SparseMatrix_Vector_Multiplication(const Eigen::SparseMatrix<double> &S, Eigen::VectorXd& a, Eigen::VectorXd& b)
{
	cout << "Setting up the handle in GPU \n";
	// --- Initialize cuSPARSE
	cusparseHandle_t	handle;
	cudaError_t			cudaStat1 = cudaSuccess;
	cusparseDirection_t dirA = CUSPARSE_DIRECTION_COLUMN;

	const int m = S.rows();
	const int n = S.cols();
	double alpha = 1.0;
	double beta = 0.0; 
	cusparseMatDescr_t descrA = 0, descrB = 0;
	int nnzA = S.nonZeros();
	int rowsA = m + 1;
	printf("The data: S=%dx%d | nnzA:%d, rowsA:%d |\n", m, n, nnzA, rowsA);

	/* Setting up the data for the vectors */
	cout << "Allocating memory and setup in CPU + GPU \n";
	double *h_a = (double*)malloc(n * sizeof(double));
	double *h_b = (double*)malloc(m * sizeof(double));
	h_a = a.data();
	double *d_a;  cudaStat1 = cudaMalloc(&d_a, n*sizeof(double)); cout << "__status:" << cudaStat1 << endl;
	double *d_b;  cudaStat1 = cudaMalloc(&d_b, m*sizeof(double)); cout << "__status:" << cudaStat1 << endl;
	cudaStat1 = cudaMemcpy(d_a, h_a, n * sizeof(double), cudaMemcpyHostToDevice);
	cudaStat1 = cudaMemcpy(d_b, h_b, m * sizeof(double), cudaMemcpyHostToDevice);

	cout << "Seting up the CSR data \n";
	/* Setting up the data for the Matrix (CSR) */
	double *h_csrVal = (double*)malloc(S.nonZeros() * sizeof(double));
	int *h_csrColInd = (int*)malloc(nnzA * sizeof(int));
	int *h_csrRowPtr = (int*)malloc(rowsA * sizeof(int));
	int *d_csrColInd; cudaMalloc(&d_csrColInd, nnzA * sizeof(int));
	int *d_csrRowPtr; cudaMalloc(&d_csrRowPtr, rowsA * sizeof(int));
	double *d_csrVal; cudaMalloc(&d_csrVal, nnzA * sizeof(double));
	cudaStat1 = cudaMemcpy(d_csrColInd, h_csrColInd, nnzA * sizeof(int), cudaMemcpyHostToDevice); cout << "__col: status:" << cudaStat1 << endl;
	cudaStat1 = cudaMemcpy(d_csrRowPtr, h_csrRowPtr, rowsA * sizeof(int), cudaMemcpyHostToDevice); cout << "__rows: status:" << cudaStat1 << endl;
	cudaStat1 = cudaMemcpy(d_csrVal, h_csrVal, nnzA * sizeof(double), cudaMemcpyHostToDevice); cout << "__val: status:" << cudaStat1 << endl;
	
	cout << "Finding nnz |";
	int *d_bsrRowPtr; cudaMalloc(&d_bsrRowPtr, rowsA * sizeof(int));
	int *d_bsrColInd;
	double *d_bsrVal;
	int nnzb = nnzA;	
	//cusparseXcsr2bsrNnz(handle, dirA, m, n,	descrA, d_csrRowPtr, d_csrColInd, 1, descrB, d_bsrRowPtr, &nnzb);	
	cout << " with " << nnzb << " non zero entries \n";

	cout << "Converting to BSR format \n";
	cudaStat1 = cudaMalloc(&d_bsrColInd, nnzb * sizeof(int)); cout << "__status:" << cudaStat1 << endl;
	cudaStat1 = cudaMalloc(&d_bsrVal, nnzb * sizeof(double)); cout << "__status:" << cudaStat1 << endl;

	cusparseDcsr2bsr(handle, dirA, m, n, descrA, d_csrVal, d_csrRowPtr, d_csrColInd, 1, descrB, d_bsrVal, d_bsrRowPtr, d_bsrColInd);
	
	
	cout << "Doing the real multiplication \n";
	//cusparseDbsrmv(handle, dirA, CUSPARSE_OPERATION_NON_TRANSPOSE, m, n, S.nonZeros(), &alpha, descrA, S.valuePtr(), S.outerIndexPtr(), S.innerIndexPtr(), 1, a.data(), &beta, b.data());
	cusparseCreate(&handle);
	int blockdim = 1; 
	cusparseStatus_t cusparseStat1 = cusparseDbsrmv(handle, dirA, CUSPARSE_OPERATION_NON_TRANSPOSE, m, n, nnzb, &alpha, descrB, d_bsrVal, d_bsrRowPtr, d_bsrColInd, blockdim, d_a, &beta, d_b);
	cout << "__status:" << cusparseStat1 << endl;

	cout << "Copying to CPU \n";
	cudaMemcpy(h_b, d_b, m * sizeof(double), cudaMemcpyDeviceToHost);
	b = Eigen::Map<Eigen::VectorXd>(h_b, m);
	cout << b.block(0, 0, 10, 1) << endl; 
	//cusparseDcsrmv(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, m, n, nnzS, &alpha, descrA, d_S_sp, d_S_RowIndices, d_S_ColIndices, d_a, &beta, d_b);
}
//void SparseMatrix_Vector_Multiplication(const Eigen::SparseMatrix<double> &S, Eigen::VectorXd& a, Eigen::VectorXd& b)
//{
//	// --- Initialize cuSPARSE
//	cusparseHandle_t	handle;
//	cudaError_t			cudaStat1 = cudaSuccess;
//
//	const int m = S.rows();
//	const int n = S.cols();
//	const int lda = m;
//	printf("Input S=%dx%d | a=%d | m=%d, n=%d, lda=%d \n", S.rows(), S.cols(), a.size(), m, n, lda);
//
//	cout << "allocating memory in the host/CPU \n";
//	// --- host-size matrices
//	double *h_S = (double*)malloc(m*n * sizeof(double));
//	double *h_a = (double*)malloc(n * sizeof(double));
//	double *h_b = (double*)malloc(m * sizeof(double));
//	Eigen::MatrixXd Snew(S);
//	h_S = Snew.data();
//	h_a = a.data();
//
//	cout << "allocating memory in the device/GPU \n";
//	// --- Create device arrays and copy host arrays to them
//	double *d_S;  cudaStat1 = cudaMalloc(&d_S, m*n*sizeof(double));
//	double *d_a;  cudaStat1 = cudaMalloc(&d_a, n*sizeof(double));
//	double *d_b;  cudaStat1 = cudaMalloc(&d_b, m*sizeof(double));
//	cout << cudaStat1 << ": copying data from CPU to GPU \n";
//	cudaStat1 = cudaMemcpy(d_S, h_S, m*n * sizeof(double), cudaMemcpyHostToDevice);
//	cudaStat1 = cudaMemcpy(d_a, h_a, n * sizeof(double), cudaMemcpyHostToDevice);
//	cudaStat1 = cudaMemcpy(d_b, h_b, m * sizeof(double), cudaMemcpyHostToDevice);
//
//	// --- Descriptor for sparse matrix A
//	cusparseMatDescr_t descrA;     
//	cusparseCreateMatDescr(&descrA);
//	cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL);
//	cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ONE);
//
//	cout << cudaStat1 << "Analyzing sparsity pattern in input matrix S : ";
//	int nnzS = 0;	
//	// --- Device side number of nonzero elements per row of matrix A
//	int *d_nnzPerVectorA;   
//	cudaMalloc(&d_nnzPerVectorA, m * sizeof(*d_nnzPerVectorA));
//	cusparseDnnz(handle, CUSPARSE_DIRECTION_ROW, m, n, descrA, d_S, lda, d_nnzPerVectorA, &nnzS);
//	cout << "with " << nnzS << " non-zero entries \n";
//
//	// --- Host side number of nonzero elements per row of matrix A
//	int *h_nnzPerVectorA = (int *)malloc(m * sizeof(*h_nnzPerVectorA));
//	cudaMemcpy(h_nnzPerVectorA, d_nnzPerVectorA, m * sizeof(*h_nnzPerVectorA), cudaMemcpyDeviceToHost);
//
//	cout << "Creating a sparse matrix S in GPU \n";
//	// --- Device side sparse matrix
//	double *d_S_sp;            
//	cudaMalloc(&d_S_sp, nnzS * sizeof(*d_S_sp));
//
//	int *d_S_RowIndices;    
//	cudaMalloc(&d_S_RowIndices, (m + 1) * sizeof(*d_S_RowIndices));
//	int *d_S_ColIndices;    
//	cudaMalloc(&d_S_ColIndices, nnzS * sizeof(*d_S_ColIndices));
//
//	cout << "Converting from Dense to Sparse \n";
//	cusparseDdense2csr(handle, m, n, descrA, d_S, lda, d_nnzPerVectorA, d_S_sp, d_S_RowIndices, d_S_ColIndices);
//
//	// --- Host side sparse matrices
//	cout << "allocating memory in the device/GPU \n";
//	double *h_S_sp = (double *)malloc(nnzS * sizeof(*h_S));
//	int *h_S_RowIndices = (int *)malloc((m + 1) * sizeof(*h_S_RowIndices));
//	int *h_S_ColIndices = (int *)malloc(nnzS * sizeof(*h_S_ColIndices));
//	cout << "copying data from CPU to GPU \n";
//	cudaMemcpy(h_S_sp, d_S_sp, nnzS * sizeof(*h_S_sp), cudaMemcpyDeviceToHost);
//	cudaMemcpy(h_S_RowIndices, d_S_RowIndices, (m + 1) * sizeof(*h_S_RowIndices), cudaMemcpyDeviceToHost);
//	cudaMemcpy(h_S_ColIndices, d_S_ColIndices, nnzS * sizeof(*h_S_ColIndices), cudaMemcpyDeviceToHost);
//
//	cout << "Performing multiplication \n";
//	const double alpha = 1.;
//	const double beta = 0.;
//	cusparseDcsrmv(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, m, n, nnzS, &alpha, descrA, d_S_sp, d_S_RowIndices, d_S_ColIndices, d_a ,&beta, d_b);
//	
//	cout << "Returning results to CPU \n";
//	cudaMemcpy(h_b, d_b, m * sizeof(double), cudaMemcpyDeviceToHost);
//
//
//}

//template<typename Scalar>
void manuallyDestroySparseMatrix(Eigen::SparseMatrix<double> &M)
{
	M.~SparseMatrix();
	//Eigen::SparseMatrix<double> M;
}

void evalSurfaceGraph() {
	using namespace matlab::engine;

	// Start MATLAB engine synchronously
	std::unique_ptr<MATLABEngine> matlabPtr = startMATLAB();
}