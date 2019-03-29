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

//template<typename Scalar>
void manuallyDestroySparseMatrix(Eigen::SparseMatrix<double> &M)
{
	M.~SparseMatrix();
	Eigen::SparseMatrix<double> K;
}

