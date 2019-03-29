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
		getline(file, oneLine);
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
			//getline(iStream, oneWord, ' ');
			//i = stoi(oneWord);
			//getline(iStream, oneWord, ' ');
			//j = stoi(oneWord);
			//getline(iStream, oneWord, ' ');
			//v = stod(oneWord);

			//MTriplet.push_back(Eigen::Triplet<double>(i, j, v));
			i++;
		}
	}
	file.close();

	printf("Size=%dx%d, with %d elements\n", i + 1, i + 1, MTriplet.size());
	M.resize(i + 1, i + 1);
	M.setFromTriplets(MTriplet.begin(), MTriplet.end());
}

//template<typename Scalar>
void manuallyDestroySparseMatrix(Eigen::SparseMatrix<double> &M)
{
	M.~SparseMatrix();
	Eigen::SparseMatrix<double> K;
}

