#include "Utility.h"

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

//template<typename Scalar>
void manuallyDestroySparseMatrix(Eigen::SparseMatrix<double> &M)
{
	M.~SparseMatrix();
	Eigen::SparseMatrix<double> K;
}

