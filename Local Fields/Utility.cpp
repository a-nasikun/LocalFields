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

	engPutVariable(ep, "ApproxFields", MM);
	engEvalString(ep, "ApproxFields=full(ApproxFields);");
	engEvalString(ep, "save('F:/PROGRAMMING/Localized Vector Fields/Build/ForChristopher/ChineseDragon/ApproxFields','ApproxFields');");

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
	eigVecM = engGetVariable(ep, "EigVect");
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
	
	engEvalString(ep, "clear;");
	engClose(ep);
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

void ReadVectorFromMatlab(Eigen::VectorXd& v, const string& filename)
{
	using namespace matlab::engine;
	Engine *ep;
	mxArray *eigValM, *eigVecM;		// for Matlab
	double	*eigValE, *eigVecE;		// for Eigen
	const int NUM_EIGEN = 500;
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
	eigValM = engGetVariable(ep, "EigVal");
	eigValE = (double*)malloc(NUM_BLOCKS*NUM_EIGEN * sizeof(double));
	memcpy((void *)eigValE, (void *)mxGetPr(eigValM), NUM_BLOCKS*NUM_EIGEN * sizeof(double));
	v = Eigen::Map<Eigen::VectorXd>(eigValE, NUM_EIGEN);
	cout << "Eigenvalues: \n" << v << endl; 

		
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
double LoadSparseMatrixFromTxtFile(const string& filename, Eigen::SparseMatrix<double> &M)
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
double ConstructInverseMassMatrix(Eigen::SparseMatrix<double> &M, Eigen::SparseMatrix<double> &MInv)
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
	printf("File is closed\n");
}

void writeEigenDenseMatrixToBinary(Eigen::MatrixXd M, const std::string &filename)
{
	std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
	int rows = M.rows(), cols = M.cols();
	out.write((char*)(&rows), sizeof(int));
	out.write((char*)(&cols), sizeof(int));
	out.write((char*)M.data(), rows*cols * sizeof(typename double));
	out.close();
}

void readEigenDenseMatrixFromBinary(const std::string &filename, Eigen::MatrixXd M)
{
	std::ifstream in(filename, std::ios::in | std::ios::binary);
	int rows = 0, cols = 0;
	in.read((char*)(&rows), sizeof(int));
	in.read((char*)(&cols), sizeof(int));
	M.resize(rows, cols);
	in.read((char *)M.data(), rows*cols * sizeof(double));
	in.close();
}

//template<typename Scalar>
void manuallyDestroySparseMatrix(Eigen::SparseMatrix<double> &M)
{
	M.~SparseMatrix();
	Eigen::SparseMatrix<double> K;
}

void evalSurfaceGraph() {
	using namespace matlab::engine;

	// Start MATLAB engine synchronously
	std::unique_ptr<MATLABEngine> matlabPtr = startMATLAB();
}