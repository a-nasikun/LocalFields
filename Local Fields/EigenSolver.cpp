#include "EigenSolver.h"


/* For spectra*/
#include <GenEigsSolver.h>
#include <SymEigsSolver.h>
#include <SymEigsShiftSolver.h>
#include <GenEigsRealShiftSolver.h>
#include <SymGEigsSolver.h>
#include <MatOp/SparseGenMatProd.h>
#include <MatOp/SparseSymMatProd.h>
#include <MatOp/SparseCholesky.h>
#include <MatOp/SparseSymShiftSolve.h>
#include <MatOp/SparseGenRealShiftSolve.h>

///#include "viennacl/scalar.hpp"
///#include "viennacl/vector.hpp"
///#include "viennacl/matrix.hpp"
///#include "viennacl/matrix_proxy.hpp"
///#include "viennacl/compressed_matrix.hpp"
///
///#include "viennacl/linalg/lanczos.hpp"
///#include "viennacl/io/matrix_market.hpp"


/* Computing Eigenstructure in GPU */
void computeEigenGPU(cusolverDnHandle_t& cusolverH, Eigen::SparseMatrix<double> &S_, Eigen::SparseMatrix<double> &M_, Eigen::MatrixXd &LDEigVec, Eigen::VectorXd &LDEigVal)
{
	//cusolverDnHandle_t	cusolverH = NULL;
	cusolverH = NULL; 
	cusolverStatus_t	cusolver_status = CUSOLVER_STATUS_SUCCESS;
	cudaError_t			cudaStat1 = cudaSuccess;
	cudaError_t			cudaStat2 = cudaSuccess;
	cudaError_t			cudaStat3 = cudaSuccess;
	cudaError_t			cudaStat4 = cudaSuccess;

	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					time_span;

	t1 = chrono::high_resolution_clock::now();

	const int m = S_.rows();
	const int lda = m;

	double *A = (double*)std::malloc(m*lda * sizeof(double));
	double *B = (double*)std::malloc(m*lda * sizeof(double));
	Eigen::MatrixXd Snew(S_);
	Eigen::MatrixXd Mnew(M_);
	A = Snew.data();
	B = Mnew.data();

	double	*V = (double*)std::malloc(m*lda * sizeof(double)); // [lda*m];		// eigenvectors 
	double	*W = (double*)std::malloc(m * sizeof(double)); 			// eigenvalues 
	double	*d_A = NULL;
	double	*d_B = NULL;
	double	*d_W = NULL;
	int		*devInfo = NULL;
	double	*d_work = NULL;
	int		lwork = 0;
	int		info_gpu = 0;
	
	// step 1: create cusolver/cublas handle 
	cusolver_status = cusolverDnCreate(&cusolverH);
	assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
	t2 = chrono::high_resolution_clock::now();
	time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
	printf("Time to create handle in GPU is %.5f. \n", time_span);

	// step 2: copy A and B to device 
	cudaStat1 = cudaMalloc((void**)&d_A, sizeof(double) * lda * m);
	cudaStat2 = cudaMalloc((void**)&d_B, sizeof(double) * lda * m);
	cudaStat3 = cudaMalloc((void**)&d_W, sizeof(double) * m);
	cudaStat4 = cudaMalloc((void**)&devInfo, sizeof(int));
	cout << "Status of Allocating Memory" << endl;
	cout << cudaStat1 << ", " << cudaStat2 << ", " << cudaStat3 << endl;
	
	assert(cudaSuccess == cudaStat1);
	assert(cudaSuccess == cudaStat2);
	assert(cudaSuccess == cudaStat3);
	assert(cudaSuccess == cudaStat4);
	cudaStat1 = cudaMemcpy(d_A, A, sizeof(double) * lda * m, cudaMemcpyHostToDevice);
	cudaStat2 = cudaMemcpy(d_B, B, sizeof(double) * lda * m, cudaMemcpyHostToDevice);
	assert(cudaSuccess == cudaStat1);
	assert(cudaSuccess == cudaStat2);
	cout << "Status of Copying CPU => GPU" << endl;
	cout << cudaStat1 << ", " << cudaStat2 << endl;

	// step 3: query working space of sygvd 
	cusolverEigType_t	itype = CUSOLVER_EIG_TYPE_1;		// A*x = (lambda)*B*x 
	cusolverEigMode_t	jobz = CUSOLVER_EIG_MODE_VECTOR;	// compute eigenvalues and eigenvectors. 
	cublasFillMode_t	uplo = CUBLAS_FILL_MODE_UPPER;
	cusolver_status = cusolverDnDsygvd_bufferSize(cusolverH, itype, jobz, uplo, m, d_A, lda, d_B, lda, d_W, &lwork);
	assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);
	cudaStat1 = cudaMalloc((void**)&d_work, sizeof(double)*lwork);
	assert(cudaSuccess == cudaStat1);
	cout << "Status of asking working space" << endl;
	cout << cudaStat1 << endl;
	
	// step 4: compute spectrum of (A,B) 
	cusolver_status = cusolverDnDsygvd(cusolverH, itype, jobz, uplo, m, d_A, lda, d_B, lda, d_W, d_work, lwork, devInfo);
	cudaStat1 = cudaDeviceSynchronize();
	assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
	assert(cudaSuccess == cudaStat1);

	cout << "Status of eigensolver" << endl;
	cout << cusolver_status << ", " << cudaStat1 << endl;

	cudaStat1 = cudaMemcpy(W, d_W, sizeof(double)*m, cudaMemcpyDeviceToHost);
	cudaStat2 = cudaMemcpy(V, d_A, sizeof(double)*lda*m, cudaMemcpyDeviceToHost);
	cudaStat3 = cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost);
	
	cout << "Status of Copying" << endl;
	cout << cudaStat1 << ", " << cudaStat2 << ", " << cudaStat3 << endl;

	assert(cudaSuccess == cudaStat1);
	assert(cudaSuccess == cudaStat2);
	assert(cudaSuccess == cudaStat3);
	printf("after sygvd: info_gpu = %d\n", info_gpu);
	assert(0 == info_gpu);

	t2 = chrono::high_resolution_clock::now();
	time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
	printf("Time to compute eigenstructure is %.5f. \n", time_span);

	// Copying back to EigenFormat (dense) for Eigenvalues and Eigenvectors
	LDEigVec = Eigen::Map<Eigen::MatrixXd, 0, Eigen::OuterStride<>>(V, m, lda, Eigen::OuterStride<>(m));
	LDEigVal = Eigen::Map<Eigen::VectorXd>(W, m);
	//Eigen::MatrixXd newMatA(a.rows(), a.cols());
	//newMatA = Eigen::Map<Eigen::MatrixXd, 0, Eigen::OuterStride<>>(stdFormat, a.rows(), a.cols(), Eigen::OuterStride<>(a.rows()));

	t2 = chrono::high_resolution_clock::now();
	time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
	printf("Time to copy data to convert to EigenFormat is %.5f. \n", time_span);

	// free resources 
	if (d_A) cudaFree(d_A);
	if (d_B) cudaFree(d_B);
	if (d_W) cudaFree(d_W);
	if (devInfo) cudaFree(devInfo);
	if (d_work) cudaFree(d_work);
	if (cusolverH) cusolverDnDestroy(cusolverH);
	cudaDeviceReset();

	t2 = chrono::high_resolution_clock::now();
	time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
	printf("Getting out from CUDA: %.5f. \n", time_span);
}

/* Computing Eigenstructure in GPU */
void computeEigenGPU(Eigen::SparseMatrix<double> &S_, Eigen::SparseMatrix<double> &M_, Eigen::MatrixXd &LDEigVec, Eigen::VectorXd &LDEigVal)
{
	cusolverDnHandle_t	cusolverH = NULL;
	cusolverStatus_t	cusolver_status = CUSOLVER_STATUS_SUCCESS;
	cudaError_t			cudaStat1 = cudaSuccess;
	cudaError_t			cudaStat2 = cudaSuccess;
	cudaError_t			cudaStat3 = cudaSuccess;
	cudaError_t			cudaStat4 = cudaSuccess;

	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					time_span;

	t1 = chrono::high_resolution_clock::now();

	const int m = S_.rows();
	const int lda = m;

	double *A = (double*)std::malloc(m*lda * sizeof(double));
	double *B = (double*)std::malloc(m*lda * sizeof(double));
	Eigen::MatrixXd Snew(S_);
	Eigen::MatrixXd Mnew(M_);
	A = Snew.data();
	B = Mnew.data();

	double	*V = (double*)std::malloc(m*lda * sizeof(double)); // [lda*m];		// eigenvectors 
	double	*W = (double*)std::malloc(m * sizeof(double)); 			// eigenvalues 
	double	*d_A = NULL;
	double	*d_B = NULL;
	double	*d_W = NULL;
	int		*devInfo = NULL;
	double	*d_work = NULL;
	int		lwork = 0;
	int		info_gpu = 0;

	// step 1: create cusolver/cublas handle 
	cusolver_status = cusolverDnCreate(&cusolverH);
	assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
	t2 = chrono::high_resolution_clock::now();
	time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
	printf("Time to create handle in GPU is %.5f. \n", time_span);

	// step 2: copy A and B to device 
	cudaStat1 = cudaMalloc((void**)&d_A, sizeof(double) * lda * m);
	cudaStat2 = cudaMalloc((void**)&d_B, sizeof(double) * lda * m);
	cudaStat3 = cudaMalloc((void**)&d_W, sizeof(double) * m);
	cudaStat4 = cudaMalloc((void**)&devInfo, sizeof(int));
	cout << "Status of Allocating Memory" << endl;
	cout << cudaStat1 << ", " << cudaStat2 << ", " << cudaStat3 << endl;

	assert(cudaSuccess == cudaStat1);
	assert(cudaSuccess == cudaStat2);
	assert(cudaSuccess == cudaStat3);
	assert(cudaSuccess == cudaStat4);
	cudaStat1 = cudaMemcpy(d_A, A, sizeof(double) * lda * m, cudaMemcpyHostToDevice);
	cudaStat2 = cudaMemcpy(d_B, B, sizeof(double) * lda * m, cudaMemcpyHostToDevice);
	assert(cudaSuccess == cudaStat1);
	assert(cudaSuccess == cudaStat2);
	cout << "Status of Copying CPU => GPU" << endl;
	cout << cudaStat1 << ", " << cudaStat2 << endl;

	// step 3: query working space of sygvd 
	cusolverEigType_t	itype = CUSOLVER_EIG_TYPE_1;		// A*x = (lambda)*B*x 
	cusolverEigMode_t	jobz = CUSOLVER_EIG_MODE_VECTOR;	// compute eigenvalues and eigenvectors. 
	cublasFillMode_t	uplo = CUBLAS_FILL_MODE_UPPER;
	cusolver_status = cusolverDnDsygvd_bufferSize(cusolverH, itype, jobz, uplo, m, d_A, lda, d_B, lda, d_W, &lwork);
	assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);
	cudaStat1 = cudaMalloc((void**)&d_work, sizeof(double)*lwork);
	assert(cudaSuccess == cudaStat1);
	cout << "Status of asking working space" << endl;
	cout << cudaStat1 << endl;

	// step 4: compute spectrum of (A,B) 
	cusolver_status = cusolverDnDsygvd(cusolverH, itype, jobz, uplo, m, d_A, lda, d_B, lda, d_W, d_work, lwork, devInfo);
	cudaStat1 = cudaDeviceSynchronize();
	assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
	assert(cudaSuccess == cudaStat1);

	cout << "Status of eigensolver" << endl;
	cout << cusolver_status << ", " << cudaStat1 << endl;

	cudaStat1 = cudaMemcpy(W, d_W, sizeof(double)*m, cudaMemcpyDeviceToHost);
	cudaStat2 = cudaMemcpy(V, d_A, sizeof(double)*lda*m, cudaMemcpyDeviceToHost);
	cudaStat3 = cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost);

	cout << "Status of Copying" << endl;
	cout << cudaStat1 << ", " << cudaStat2 << ", " << cudaStat3 << endl;

	assert(cudaSuccess == cudaStat1);
	assert(cudaSuccess == cudaStat2);
	assert(cudaSuccess == cudaStat3);
	printf("after sygvd: info_gpu = %d\n", info_gpu);
	assert(0 == info_gpu);

	t2 = chrono::high_resolution_clock::now();
	time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
	printf("Time to compute eigenstructure is %.5f. \n", time_span);

	// Copying back to EigenFormat (dense) for Eigenvalues and Eigenvectors
	LDEigVec = Eigen::Map<Eigen::MatrixXd, 0, Eigen::OuterStride<>>(V, m, lda, Eigen::OuterStride<>(m));
	LDEigVal = Eigen::Map<Eigen::VectorXd>(W, m);
	//Eigen::MatrixXd newMatA(a.rows(), a.cols());
	//newMatA = Eigen::Map<Eigen::MatrixXd, 0, Eigen::OuterStride<>>(stdFormat, a.rows(), a.cols(), Eigen::OuterStride<>(a.rows()));

	t2 = chrono::high_resolution_clock::now();
	time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
	printf("Time to copy data to convert to EigenFormat is %.5f. \n", time_span);

	// free resources 
	if (d_A) cudaFree(d_A);
	if (d_B) cudaFree(d_B);
	if (d_W) cudaFree(d_W);
	if (devInfo) cudaFree(devInfo);
	if (d_work) cudaFree(d_work);
	if (cusolverH) cusolverDnDestroy(cusolverH);
	cudaDeviceReset();

	t2 = chrono::high_resolution_clock::now();
	time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
	printf("Getting out from CUDA: %.5f. \n", time_span);
}

void computeEigenGPU(Eigen::MatrixXd &S_, Eigen::MatrixXd &M_, Eigen::MatrixXd &LDEigVec, Eigen::VectorXd &LDEigVal)
{
	cusolverDnHandle_t	cusolverH = NULL;
	cusolverStatus_t	cusolver_status = CUSOLVER_STATUS_SUCCESS;
	cudaError_t			cudaStat1 = cudaSuccess;
	cudaError_t			cudaStat2 = cudaSuccess;
	cudaError_t			cudaStat3 = cudaSuccess;
	cudaError_t			cudaStat4 = cudaSuccess;

	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					time_span;

	t1 = chrono::high_resolution_clock::now();

	const int m = S_.rows();
	const int lda = m;

	double *A = (double*)std::malloc(m*lda * sizeof(double));
	double *B = (double*)std::malloc(m*lda * sizeof(double));
	A = S_.data();
	B = M_.data();

	double	*V = (double*)std::malloc(m*lda * sizeof(double)); // [lda*m];		// eigenvectors 
	double	*W = (double*)std::malloc(m * sizeof(double)); 			// eigenvalues 
	double	*d_A = NULL;
	double	*d_B = NULL;
	double	*d_W = NULL;
	int		*devInfo = NULL;
	double	*d_work = NULL;
	int		lwork = 0;
	int		info_gpu = 0;

	// step 1: create cusolver/cublas handle 
	cusolver_status = cusolverDnCreate(&cusolverH);
	assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
	t2 = chrono::high_resolution_clock::now();
	time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
	printf("Time to create handle in GPU is %.5f. \n", time_span);

	// step 2: copy A and B to device 
	cudaStat1 = cudaMalloc((void**)&d_A, sizeof(double) * lda * m);
	cudaStat2 = cudaMalloc((void**)&d_B, sizeof(double) * lda * m);
	cudaStat3 = cudaMalloc((void**)&d_W, sizeof(double) * m);
	cudaStat4 = cudaMalloc((void**)&devInfo, sizeof(int));
	//cout << "Status of Allocating Memory" << endl;
	//cout << cudaStat1 << ", " << cudaStat2 << ", " << cudaStat3 << endl;

	assert(cudaSuccess == cudaStat1);
	assert(cudaSuccess == cudaStat2);
	assert(cudaSuccess == cudaStat3);
	assert(cudaSuccess == cudaStat4);
	cudaStat1 = cudaMemcpy(d_A, A, sizeof(double) * lda * m, cudaMemcpyHostToDevice);
	cudaStat2 = cudaMemcpy(d_B, B, sizeof(double) * lda * m, cudaMemcpyHostToDevice);
	assert(cudaSuccess == cudaStat1);
	assert(cudaSuccess == cudaStat2);
	//cout << "Status of Copying CPU => GPU" << endl;
	//cout << cudaStat1 << ", " << cudaStat2 << endl;

	// step 3: query working space of sygvd 
	cusolverEigType_t	itype = CUSOLVER_EIG_TYPE_1;		// A*x = (lambda)*B*x 
	cusolverEigMode_t	jobz = CUSOLVER_EIG_MODE_VECTOR;	// compute eigenvalues and eigenvectors. 
	cublasFillMode_t	uplo = CUBLAS_FILL_MODE_UPPER;
	cusolver_status = cusolverDnDsygvd_bufferSize(cusolverH, itype, jobz, uplo, m, d_A, lda, d_B, lda, d_W, &lwork);
	assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);
	cudaStat1 = cudaMalloc((void**)&d_work, sizeof(double)*lwork);
	assert(cudaSuccess == cudaStat1);
	//cout << "Status of asking working space" << endl;
	//cout << cudaStat1 << endl;

	// step 4: compute spectrum of (A,B) 
	cusolver_status = cusolverDnDsygvd(cusolverH, itype, jobz, uplo, m, d_A, lda, d_B, lda, d_W, d_work, lwork, devInfo);
	cudaStat1 = cudaDeviceSynchronize();
	assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
	assert(cudaSuccess == cudaStat1);

	//cout << "Status of eigensolver" << endl;
	//cout << cusolver_status << ", " << cudaStat1 << endl;

	cudaStat1 = cudaMemcpy(W, d_W, sizeof(double)*m, cudaMemcpyDeviceToHost);
	cudaStat2 = cudaMemcpy(V, d_A, sizeof(double)*lda*m, cudaMemcpyDeviceToHost);
	cudaStat3 = cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost);

	//cout << "Status of Copying" << endl;
	//cout << cudaStat1 << ", " << cudaStat2 << ", " << cudaStat3 << endl;

	assert(cudaSuccess == cudaStat1);
	assert(cudaSuccess == cudaStat2);
	assert(cudaSuccess == cudaStat3);
	//printf("after sygvd: info_gpu = %d\n", info_gpu);
	assert(0 == info_gpu);

	// Copying back to EigenFormat (dense) for Eigenvalues and Eigenvectors
	LDEigVec = Eigen::Map<Eigen::MatrixXd, 0, Eigen::OuterStride<>>(V, m, lda, Eigen::OuterStride<>(m));
	LDEigVal = Eigen::Map<Eigen::VectorXd>(W, m);
	//Eigen::MatrixXd newMatA(a.rows(), a.cols());
	//newMatA = Eigen::Map<Eigen::MatrixXd, 0, Eigen::OuterStride<>>(stdFormat, a.rows(), a.cols(), Eigen::OuterStride<>(a.rows()));

	t2 = chrono::high_resolution_clock::now();
	time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
	printf("Time to copy data to convert to EigenFormat is %.5f. \n", time_span);

	// free resources 
	if (d_A) cudaFree(d_A);
	if (d_B) cudaFree(d_B);
	if (d_W) cudaFree(d_W);
	if (devInfo) cudaFree(devInfo);
	if (d_work) cudaFree(d_work);
	if (cusolverH) cusolverDnDestroy(cusolverH);
	cudaDeviceReset();
}

/* Computing Eigenstructure in Matlab */
void computeEigenMatlab(Eigen::SparseMatrix<double> &S, Eigen::SparseMatrix<double> &M, Eigen::MatrixXd &EigVec, Eigen::VectorXd &EigVal)
{
	//printf("Size of S = %dx%d\n", S.rows(), S.cols());
	using namespace matlab::engine;
	Engine *ep;
	mxArray *MM = NULL, *MS = NULL, *result = NULL, *eigVecResult, *nEigs;

	chrono::high_resolution_clock::time_point	t1, t2, t3, t4;
	chrono::duration<double>					time_span, ts2;

	const int NNZ_S = S.nonZeros();
	const int NNZ_M = M.nonZeros();
	double *eigVal, *eigVec;

	// Allocate memory for S and M (sparse representation)
	double	*srs = (double*)malloc(NNZ_S * sizeof(double));
	mwIndex *irs = (mwIndex*)malloc(NNZ_S * sizeof(mwIndex));
	mwIndex *jcs = (mwIndex*)malloc((S.cols() + 1) * sizeof(mwIndex));

	double	*srm = (double*)malloc(NNZ_M * sizeof(double));
	mwIndex *irm = (mwIndex*)malloc(NNZ_M * sizeof(mwIndex));
	mwIndex *jcm = (mwIndex*)malloc((M.cols() + 1) * sizeof(mwIndex));

	// Bind MM with M, and MS with S
	MS = mxCreateSparse(S.rows(), S.cols(), NNZ_S, mxREAL);
	srs = mxGetPr(MS);
	irs = mxGetIr(MS);
	jcs = mxGetJc(MS);

	MM = mxCreateSparse(M.rows(), M.cols(), NNZ_M, mxREAL);
	srm = mxGetPr(MM);
	irm = mxGetIr(MM);
	jcm = mxGetJc(MM);

	// Setting initial variable value
	int nnzSCounter = 0;
	int nnzMCounter = 0;

	t1 = chrono::high_resolution_clock::now();

	// Getting matrix S
	jcs[0] = nnzSCounter;
	for (int i = 0; i < S.outerSize(); i++) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(S, i); it; ++it) {
			srs[nnzSCounter] = it.value();
			irs[nnzSCounter] = it.row();
			nnzSCounter++;
		}
		jcs[i + 1] = nnzSCounter;
	}

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
		//cout << "MATLAB STARTS. OH YEAH!!!" << endl;
	}	 

	// Compute Eigenvalue in Matlab
	int NUM_EIGEN = 2;

	engPutVariable(ep, "MS", MS);
	engPutVariable(ep, "MM", MM);

	t3 = chrono::high_resolution_clock::now();
	engEvalString(ep, "[EigVec, EigVal]=eigs(MS,MM,2,'smallestabs');");
	//engEvalString(ep, "[EigVec, EigVal]=eigs(MS,MM);");
	engEvalString(ep, "EigVal=diag(EigVal);");
	engEvalString(ep, "hold on; plot(1:2, EigVal(1:2),'LineWidth',1.5);"); // has to do it this way for "correct" plot
	t4 = chrono::high_resolution_clock::now();

	//engEvalString(ep, "save('D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/CDragon_500_LDEigVect_1000samples','EigVec');");
	//engEvalString(ep, "save('D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/CDragon_500_LDEigVal_1000samples','EigVal');");
	//engEvalString(ep, "save('D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/Genus2_20_REigVect','EigVec');");
	//engEvalString(ep, "save('D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/Genus2_20_REigVal','EigVal');");
	
	result = engGetVariable(ep, "EigVal");
	eigVal = (double*)malloc(NUM_EIGEN * sizeof(double));
	memcpy((void *)eigVal, (void *)mxGetPr(result), NUM_EIGEN * sizeof(double));

	eigVecResult = engGetVariable(ep, "EigVec");
	eigVec = (double*)malloc(M.rows() * NUM_EIGEN * sizeof(double));
	memcpy((void *)eigVec, (void *)mxGetPr(eigVecResult), M.rows() * NUM_EIGEN * sizeof(double));

	EigVal = Eigen::Map<Eigen::VectorXd>(eigVal, NUM_EIGEN);
	EigVec = Eigen::Map<Eigen::MatrixXd, 0, Eigen::OuterStride<>>(eigVec, M.rows(), NUM_EIGEN, Eigen::OuterStride<>(M.rows()));

	t2 = chrono::high_resolution_clock::now();
	time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
	ts2 = chrono::duration_cast<chrono::duration<double>>(t4 - t3);
	//printf("Time to compute %d first eigenstructure is %.5f. \n", NUM_EIGEN, time_span);
	//printf("Time to compute %d first eigenstructure (in MATLAB ONLY) is %.5f. \n", NUM_EIGEN, ts2);

	// Testing out
	//printf("Check orthogonality: \n Equal=%.5f, Not-Equal=%.5f\n", EigVec.col(0).transpose()*M*EigVec.col(0), EigVec.col(0).transpose()*M*EigVec.col(50));
	//cout << "Eigen Vectors" << endl << EigVec.block(0, 0, 2, 10) << endl;
	//printf("EigenVector dimension: %dx%d\n", EigVec.rows(), EigVec.cols());
}
void computeEigenMatlab(Eigen::SparseMatrix<double> &S, Eigen::SparseMatrix<double> &M, const int& numEigs, Eigen::MatrixXd &EigVec, Eigen::VectorXd &EigVal, const string& filename)
{
	//printf("Size of S = %dx%d\n", S.rows(), S.cols());
	using namespace matlab::engine;
	Engine *ep;
	mxArray *MM = NULL, *MS = NULL, *result = NULL, *eigVecResult, *nEigsBuff;
	
	/* Storing element to determine how many eigenpairs to compute*/
	if (numEigs > M.rows())
	{
		cout << "ERROR! The number of the eigenpairs cannot exceed the size of input matrix(-ces)." << endl;
		return; 
	}
	nEigsBuff = mxCreateDoubleMatrix((mwSize)1, (mwSize)1, mxREAL);
	double *numEigsPtr = mxGetPr(nEigsBuff);
	*numEigsPtr = numEigs; 

	chrono::high_resolution_clock::time_point	t1, t2, t3, t4;
	chrono::duration<double>					time_span, ts2;

	const int NNZ_S = S.nonZeros();
	const int NNZ_M = M.nonZeros();
	double *eigVal, *eigVec;

	// Allocate memory for S and M (sparse representation)
	double	*srs = (double*)malloc(NNZ_S * sizeof(double));
	mwIndex *irs = (mwIndex*)malloc(NNZ_S * sizeof(mwIndex));
	mwIndex *jcs = (mwIndex*)malloc((S.cols() + 1) * sizeof(mwIndex));

	double	*srm = (double*)malloc(NNZ_M * sizeof(double));
	mwIndex *irm = (mwIndex*)malloc(NNZ_M * sizeof(mwIndex));
	mwIndex *jcm = (mwIndex*)malloc((M.cols() + 1) * sizeof(mwIndex));

	// Bind MM with M, and MS with S
	MS = mxCreateSparse(S.rows(), S.cols(), NNZ_S, mxREAL);
	srs = mxGetPr(MS);
	irs = mxGetIr(MS);
	jcs = mxGetJc(MS);

	MM = mxCreateSparse(M.rows(), M.cols(), NNZ_M, mxREAL);
	srm = mxGetPr(MM);
	irm = mxGetIr(MM);
	jcm = mxGetJc(MM);

	// Setting initial variable value
	int nnzSCounter = 0;
	int nnzMCounter = 0;

	t1 = chrono::high_resolution_clock::now();

	// Getting matrix S
	jcs[0] = nnzSCounter;
	for (int i = 0; i < S.outerSize(); i++) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(S, i); it; ++it) {
			srs[nnzSCounter] = it.value();
			irs[nnzSCounter] = it.row();
			nnzSCounter++;
		}
		jcs[i + 1] = nnzSCounter;
	}

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
	int *matlabStatus;
	ep = engOpen(NULL);
	if (!(ep = engOpen(""))) {
	//ep = engOpenSingleUse(NULL, NULL, matlabStatus);
	//if (!(ep = engOpenSingleUse(NULL, NULL, matlabStatus))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		cout << "CANNOT START MATLAB " << endl;
	}
	else {
		//cout << "MATLAB STARTS. OH YEAH!!!" << endl;
	}

	//cout << "Status => " << matlabStatus << endl; 

	// Compute Eigenvalue in Matlab
	int NUM_EIGEN = numEigs;

	engPutVariable(ep, "MS", MS);
	engPutVariable(ep, "MM", MM);
	engPutVariable(ep, "Num", nEigsBuff);

	t3 = chrono::high_resolution_clock::now();
	engEvalString(ep, "[data, EigVal]=eigs(MS,MM, Num(1,1),'smallestabs');");
	//engEvalString(ep, "[EigVec, EigVal]=eigs(MS,MM);");
	engEvalString(ep, "EigVal=diag(EigVal);");
	if (numEigs > 2)
	{
		engEvalString(ep, "hold on; plot(1:Num(1,1), EigVal(1:Num(1,1)),'LineWidth',1.5);"); // has to do it this way for "correct" plot		
		string approxFile = "save('" + filename + "_eigFields','data','EigVal');";		
		//string approxFile = "save('" + filename + "_eigFields','EigVec');";
		//string approxFile = "save('" + filename + "_eigvalues','EigVal');";
		cout << "Saving the eigen problem\n";
		engEvalString(ep, approxFile.c_str());
	}
	t4 = chrono::high_resolution_clock::now();

	
	result = engGetVariable(ep, "EigVal");
	eigVal = (double*)malloc(NUM_EIGEN * sizeof(double));
	memcpy((void *)eigVal, (void *)mxGetPr(result), NUM_EIGEN * sizeof(double));

	eigVecResult = engGetVariable(ep, "data");
	eigVec = (double*)malloc(M.rows() * NUM_EIGEN * sizeof(double));
	memcpy((void *)eigVec, (void *)mxGetPr(eigVecResult), M.rows() * NUM_EIGEN * sizeof(double));

	EigVal = Eigen::Map<Eigen::VectorXd>(eigVal, NUM_EIGEN);
	EigVec = Eigen::Map<Eigen::MatrixXd, 0, Eigen::OuterStride<>>(eigVec, M.rows(), NUM_EIGEN, Eigen::OuterStride<>(M.rows()));

	t2 = chrono::high_resolution_clock::now();
	time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
	ts2 = chrono::duration_cast<chrono::duration<double>>(t4 - t3);

	engEvalString(ep, "clear;");

	mxDestroyArray(MM);
	mxDestroyArray(MS);
	mxDestroyArray(result);
	mxDestroyArray(eigVecResult);
	mxDestroyArray(nEigsBuff);

	engClose(ep);

/*
	eigVal = nullptr; 
	eigVal = nullptr; 
	eigVec = nullptr; 
	srs    = nullptr;
	irs    = nullptr;
	jcs    = nullptr;
	srm    = nullptr;
	irm    = nullptr;
	jcm    = nullptr;	
*/
}

void computeEigenMatlab(Eigen::SparseMatrix<double> &S, const int& numEigs, Eigen::MatrixXd &EigVec, Eigen::VectorXd &EigVal, const string& filename)
{
	//printf("Size of S = %dx%d\n", S.rows(), S.cols());
	using namespace matlab::engine;
	Engine *ep;
	mxArray *MS = NULL, *result = NULL, *eigVecResult, *nEigsBuff;

	/* Storing element to determine how many eigenpairs to compute*/	
	nEigsBuff = mxCreateDoubleMatrix((mwSize)1, (mwSize)1, mxREAL);
	double *numEigsPtr = mxGetPr(nEigsBuff);
	*numEigsPtr = numEigs;

	chrono::high_resolution_clock::time_point	t1, t2, t3, t4;
	chrono::duration<double>					time_span, ts2;

	const int NNZ_S = S.nonZeros();
	double *eigVal, *eigVec;

	// Allocate memory for S and M (sparse representation)
	double	*srs = (double*)malloc(NNZ_S * sizeof(double));
	mwIndex *irs = (mwIndex*)malloc(NNZ_S * sizeof(mwIndex));
	mwIndex *jcs = (mwIndex*)malloc((S.cols() + 1) * sizeof(mwIndex));


	// Bind MM with M, and MS with S
	MS = mxCreateSparse(S.rows(), S.cols(), NNZ_S, mxREAL);
	srs = mxGetPr(MS);
	irs = mxGetIr(MS);
	jcs = mxGetJc(MS);

	// Setting initial variable value
	int nnzSCounter = 0;
	int nnzMCounter = 0;

	t1 = chrono::high_resolution_clock::now();

	// Getting matrix S
	jcs[0] = nnzSCounter;
	for (int i = 0; i < S.outerSize(); i++) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(S, i); it; ++it) {
			srs[nnzSCounter] = it.value();
			irs[nnzSCounter] = it.row();
			nnzSCounter++;
		}
		jcs[i + 1] = nnzSCounter;
	}

	// Start Matlab Engine
	int *matlabStatus;
	ep = engOpen(NULL);
	if (!(ep = engOpen(""))) {
		//ep = engOpenSingleUse(NULL, NULL, matlabStatus);
		//if (!(ep = engOpenSingleUse(NULL, NULL, matlabStatus))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		cout << "CANNOT START MATLAB " << endl;
	}
	else {
		//cout << "MATLAB STARTS. OH YEAH!!!" << endl;
	}

	//cout << "Status => " << matlabStatus << endl; 

	// Compute Eigenvalue in Matlab
	int NUM_EIGEN = numEigs;

	engPutVariable(ep, "MS", MS);
	engPutVariable(ep, "Num", nEigsBuff);

	t3 = chrono::high_resolution_clock::now();
	engEvalString(ep, "[EigVec, EigVal]=eigs(MS, Num(1,1),'smallestabs');");
	//engEvalString(ep, "[EigVec, EigVal]=eigs(MS,MM);");
	engEvalString(ep, "EigVal=diag(EigVal);");
	if (numEigs > 2)
	{
		cout << "Should be written: to " << filename << endl;
		engEvalString(ep, "hold on; plot(1:Num(1,1), EigVal(1:Num(1,1)),'LineWidth',1.5);"); // has to do it this way for "correct" plot		
		string approxFile = "save('" + filename + "_eigFields','data','EigVal');";
		//string approxFile = "save('" + filename + "_eigFields','EigVec');";
		//string approxFile = "save('" + filename + "_eigvalues','EigVal');";
		cout << "Saving the eigen problem\n";
		///engEvalString(ep, approxFile.c_str());
	}
	t4 = chrono::high_resolution_clock::now();


	result = engGetVariable(ep, "EigVal");
	eigVal = (double*)malloc(NUM_EIGEN * sizeof(double));
	memcpy((void *)eigVal, (void *)mxGetPr(result), NUM_EIGEN * sizeof(double));

	eigVecResult = engGetVariable(ep, "EigVec");
	eigVec = (double*)malloc(S.rows() * NUM_EIGEN * sizeof(double));
	memcpy((void *)eigVec, (void *)mxGetPr(eigVecResult), S.rows() * NUM_EIGEN * sizeof(double));

	EigVal = Eigen::Map<Eigen::VectorXd>(eigVal, NUM_EIGEN);
	EigVec = Eigen::Map<Eigen::MatrixXd, 0, Eigen::OuterStride<>>(eigVec, S.rows(), NUM_EIGEN, Eigen::OuterStride<>(S.rows()));

	t2 = chrono::high_resolution_clock::now();
	time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
	ts2 = chrono::duration_cast<chrono::duration<double>>(t4 - t3);

	engEvalString(ep, "clear; clc;");
	
	mxDestroyArray(MS);
	mxDestroyArray(result);
	mxDestroyArray(eigVecResult);
	mxDestroyArray(nEigsBuff);

	engClose(ep);
}

void computeEigenMatlab(Engine*& ep, const int tid, Eigen::SparseMatrix<double> &S, Eigen::SparseMatrix<double> &M, const int& numEigs, Eigen::MatrixXd &EigVec, Eigen::VectorXd &EigVal, const string& filename)
{
	//printf("Size of S = %dx%d\n", S.rows(), S.cols());
	using namespace matlab::engine;
	//Engine *ep;
	mxArray *MM = NULL, *MS = NULL, *result = NULL, *eigVecResult, *nEigsBuff;

	/* Storing element to determine how many eigenpairs to compute*/
	if (numEigs > M.rows())
	{
		cout << "ERROR! The number of the eigenpairs cannot exceed the size of input matrix(-ces)." << endl;
		return;
	}
	nEigsBuff = mxCreateDoubleMatrix((mwSize)1, (mwSize)1, mxREAL);
	double *numEigsPtr = mxGetPr(nEigsBuff);
	*numEigsPtr = numEigs;

	chrono::high_resolution_clock::time_point	t1, t2, t3, t4;
	chrono::duration<double>					time_span, ts2;

	const int NNZ_S = S.nonZeros();
	const int NNZ_M = M.nonZeros();
	double *eigVal, *eigVec;

	// Allocate memory for S and M (sparse representation)
	double	*srs = (double*)malloc(NNZ_S * sizeof(double));
	mwIndex *irs = (mwIndex*)malloc(NNZ_S * sizeof(mwIndex));
	mwIndex *jcs = (mwIndex*)malloc((S.cols() + 1) * sizeof(mwIndex));

	double	*srm = (double*)malloc(NNZ_M * sizeof(double));
	mwIndex *irm = (mwIndex*)malloc(NNZ_M * sizeof(mwIndex));
	mwIndex *jcm = (mwIndex*)malloc((M.cols() + 1) * sizeof(mwIndex));

	// Bind MM with M, and MS with S
	MS = mxCreateSparse(S.rows(), S.cols(), NNZ_S, mxREAL);
	srs = mxGetPr(MS);
	irs = mxGetIr(MS);
	jcs = mxGetJc(MS);

	MM = mxCreateSparse(M.rows(), M.cols(), NNZ_M, mxREAL);
	srm = mxGetPr(MM);
	irm = mxGetIr(MM);
	jcm = mxGetJc(MM);

	// Setting initial variable value
	int nnzSCounter = 0;
	int nnzMCounter = 0;

	t1 = chrono::high_resolution_clock::now();

	// Getting matrix S
	jcs[0] = nnzSCounter;
	for (int i = 0; i < S.outerSize(); i++) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(S, i); it; ++it) {
			srs[nnzSCounter] = it.value();
			irs[nnzSCounter] = it.row();
			nnzSCounter++;
		}
		jcs[i + 1] = nnzSCounter;
	}

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
	int *matlabStatus;
	ep = engOpen(NULL);
	if (!(ep = engOpen(""))) {
		//ep = engOpenSingleUse(NULL, NULL, matlabStatus);
		//if (!(ep = engOpenSingleUse(NULL, NULL, matlabStatus))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		cout << "CANNOT START MATLAB " << endl;
	}
	else {
		//cout << "The " << tid<< " MATLAB engine starts. OH YEAH!!!" << endl;
	}

	//cout << "Status => " << matlabStatus << endl; 

	// Compute Eigenvalue in Matlab
	int NUM_EIGEN = numEigs;

	engPutVariable(ep, "MS", MS);
	engPutVariable(ep, "MM", MM);
	engPutVariable(ep, "Num", nEigsBuff);

	t3 = chrono::high_resolution_clock::now();
	engEvalString(ep, "[EigVec, EigVal]=eigs(MS,MM, Num(1,1),'smallestabs');");
	//engEvalString(ep, "[EigVec, EigVal]=eigs(MS,MM);");
	engEvalString(ep, "EigVal=diag(EigVal);");
	if (numEigs > 2)
	{
		engEvalString(ep, "hold on; plot(1:Num(1,1), EigVal(1:Num(1,1)),'LineWidth',1.5);"); // has to do it this way for "correct" plot		
		string approxFile = "save('" + filename + "_eigFields','EigVec','EigVal');";
		//string approxFile = "save('" + filename + "_eigFields','EigVec');";
		//string approxFile = "save('" + filename + "_eigvalues','EigVal');";
		cout << "Saving the eigen prblem\n";
		engEvalString(ep, approxFile.c_str());

		//engEvalString(ep, "save('D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/Cube_Round_500_Approx_Eigenpatch_2000dim_40sup','EigVec', 'EigVal');");

		//engEvalString(ep, "save('D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/Kitten_Small_5000vert_EigFields_Ref','EigVec');");
		//engEvalString(ep, "save('D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/Kitten_Small_5000vert_EigValues_Ref','EigVal');");
		//engEvalString(ep, "save('D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/Genus2_20_REigVect','EigVec');");
		//engEvalString(ep, "save('D:/Nasikun/4_SCHOOL/TU Delft/Research/Projects/LocalFields/Matlab Prototyping/Data/Genus2_20_REigVal','EigVal');");
	}
	t4 = chrono::high_resolution_clock::now();


	result = engGetVariable(ep, "EigVal");
	eigVal = (double*)malloc(NUM_EIGEN * sizeof(double));
	memcpy((void *)eigVal, (void *)mxGetPr(result), NUM_EIGEN * sizeof(double));

	eigVecResult = engGetVariable(ep, "EigVec");
	eigVec = (double*)malloc(M.rows() * NUM_EIGEN * sizeof(double));
	memcpy((void *)eigVec, (void *)mxGetPr(eigVecResult), M.rows() * NUM_EIGEN * sizeof(double));

	EigVal = Eigen::Map<Eigen::VectorXd>(eigVal, NUM_EIGEN);
	EigVec = Eigen::Map<Eigen::MatrixXd, 0, Eigen::OuterStride<>>(eigVec, M.rows(), NUM_EIGEN, Eigen::OuterStride<>(M.rows()));

	t2 = chrono::high_resolution_clock::now();
	time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
	ts2 = chrono::duration_cast<chrono::duration<double>>(t4 - t3);

	mxDestroyArray(MM);
	mxDestroyArray(MS);
	mxDestroyArray(result);
	mxDestroyArray(eigVecResult);
	mxDestroyArray(nEigsBuff);

	engClose(ep);
}

/* ============ Computing eigenstructure of 2x2 matrix explicitly/analytically ============= */
void computeEigenExplicit(const Eigen::Matrix2d& M, Eigen::Vector2d& EigVal, Eigen::Matrix2d& EigVect)
{
	double trace = M.trace();
	//double det = M.determinant();
	double det = M(0, 0)*M(1, 1) - M(0, 1)*M(1, 0);
	//printf("The determinant is %.5f \n", det);
	double a_ = trace / 2.0;
	double b_ = sqrt(trace*trace/4.0 - det);

	/* Computing the eigenvalues */
	EigVal(0) = a_ + b_;
	EigVal(1) = a_ - b_;

	
	/* Computing the eigenvectors for 2x2 symmetric matrix*/
	if (abs(M(1, 0)) > 1.0*std::numeric_limits<double>::epsilon())
	{
		//cout << "Reguler " << EigVal.transpose() << endl;
		Eigen::Vector2d v1(EigVal(0) - M(1, 1), M(1, 0)); v1.normalize();
		Eigen::Vector2d v2(EigVal(1) - M(1, 1), M(1, 0)); v2.normalize();
		//EigVect.col(0) = EigVal(0)*v1; 
		//EigVect.col(1) = EigVal(1)*v2;
		EigVect.col(0) = v1;
		EigVect.col(1) = v2;
	}
	else 
	{
		///cout << "Off diagonal zeros | eigVal: " << EigVal.transpose() << endl;
		///cout << "__M=" << M << "\t\t , its trace: " << trace << endl; 
		//EigVect.col(0) = EigVal(0)*Eigen::Vector2d(1.0, 0.0);
		//EigVect.col(1) = EigVal(1)*Eigen::Vector2d(0.0, 1.0);
		EigVect.col(0) = Eigen::Vector2d(1.0, 0.0);
		EigVect.col(1) = Eigen::Vector2d(0.0, 1.0);
	}
}

/* Computing Eigenstructure in Spectra */
void computeEigenSpectra(Eigen::SparseMatrix<double> &S, Eigen::SparseMatrix<double> &M, const int& numEigs, Eigen::MatrixXd &EigVec, Eigen::VectorXd &EigVal, const string& filename)
{
	/* Getting the inverse of M*/
	vector<Eigen::Triplet<double>> MTrip;
	MTrip.reserve(M.rows());
	Eigen::SparseMatrix<double> Minv_(M.rows(), M.cols());
	for (int i = 0; i < M.outerSize(); i++) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it) {
			MTrip.push_back(Eigen::Triplet<double>(it.row(), it.col(), 1.0 / it.value()));
		}
	}
	Minv_.setFromTriplets(MTrip.begin(), MTrip.end());

	/* Computing the local laplace*/
	Eigen::SparseMatrix<double> L = Minv_*S;

	/* Computing the eigenvectors */
	//Spectra::SparseGenMatProd<double> op(S);
	Spectra::SparseSymShiftSolve<double> op(L);
	//Spectra::SparseSymMatProd<double> op(L);
	//Spectra::SparseCholesky<double> op(L);
	Spectra::SymEigsShiftSolver<double, Spectra::LARGEST_MAGN, Spectra::SparseSymShiftSolve<double>> eigs(&op, numEigs, 2*numEigs, 0.0);
	eigs.init();
	int nconv = eigs.compute();
	cout << "Result: " << eigs.info() << endl; 
	if (eigs.info() == Spectra::SUCCESSFUL) {
		EigVal = eigs.eigenvalues().real();
		EigVec = eigs.eigenvectors().real();
		cout << "I can compute the eigenproblem \n";
		cout << "Eigenvalues \n " <<  EigVal << endl; 
	}
	else {
		cout << "This results in a wrong system to solve\n";
	}
}

void computeEigenSpectra_GenSym(Eigen::SparseMatrix<double> &S, Eigen::SparseMatrix<double> &M, const int& numEigs, Eigen::MatrixXd &EigVec, Eigen::VectorXd &EigVal, const string& filename)
{
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					time_span;

	t1 = chrono::high_resolution_clock::now();

	Spectra::SparseSymMatProd<double> Sop(S);
	Spectra::SparseCholesky<double> Mop(M);

	Spectra::SymGEigsSolver<double, Spectra::SMALLEST_ALGE, Spectra::SparseSymMatProd<double>, Spectra::SparseCholesky<double>, Spectra::GEIGS_CHOLESKY> eig_solver(&Sop, &Mop, numEigs, 2 * numEigs);

	eig_solver.init();
	int nconv = eig_solver.compute();

	EigVec = eig_solver.eigenvectors();
	EigVal = eig_solver.eigenvalues();

	t2 = chrono::high_resolution_clock::now();
	time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);

	WriteDenseMatrixToMatlab(EigVec, filename);
	cout << "Eigenvalues (Spectra): \n" << EigVal << endl << endl; 
}

void computeEigenSpectra_RegNSym(Eigen::SparseMatrix<double> &S, Eigen::SparseMatrix<double> &Minv, const int& numEigs, Eigen::MatrixXd &EigVec, Eigen::VectorXd &EigVal, const string& filename)
{
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					time_span;

	t1 = chrono::high_resolution_clock::now();

	Eigen::SparseMatrix<double> L = Minv*S;

	Spectra::SparseGenRealShiftSolve<double> Lop(L);
	//Spectra::SparseGenMatProd<double> Lop(L);

	//Spectra::GenEigsSolver<double, Spectra::SMALLEST_MAGN, Spectra::SparseGenMatProd<double>> eig_solver(&Lop, numEigs, 2 * numEigs);
	Spectra::GenEigsRealShiftSolver<double, Spectra::LARGEST_MAGN, Spectra::SparseGenRealShiftSolve<double>> eig_solver(&Lop, numEigs, 2 * numEigs, 0.0);

	//Spectra::SparseSymMatProd<double> Sop(S);
	//Spectra::SparseCholesky<double> Mop(M);
	//
	//Spectra::SymGEigsSolver<double, Spectra::LARGEST_MAGN, Spectra::SparseSymMatProd<double>, Spectra::SparseCholesky<double>, Spectra::GEIGS_CHOLESKY> eig_solver(&Sop, &Mop, numEigs, 2 * numEigs);

	eig_solver.init();
	eig_solver.compute();

	EigVec = eig_solver.eigenvectors().real();
	EigVal = eig_solver.eigenvalues().real();

	t2 = chrono::high_resolution_clock::now();
	time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);

	//WriteDenseMatrixToMatlab(eig_solver.eigenvectors(), filename);
	///cout << "Eigenvalues (Spectra): \n" << EigVal << endl << endl;
}

//void testViennaCL2()
//{
//	 // If you GPU does not support double precision, use `float` instead of `double`:
//     typedef double     ScalarType;
//   
//     std::vector< std::map<unsigned int, ScalarType> > host_A;
//     if (!viennacl::io::read_matrix_market_file(host_A, "D:/Nasikun/4_SCHOOL/TU Delft/Programming/Libraries/ViennaCL-1.7.1/examples/testdata/mat65k.mtx"))
//     {
//       std::cout << "Error reading Matrix file" << std::endl;
//	   return; // EXIT_FAILURE;
//     }
//   
//     viennacl::compressed_matrix<ScalarType> A;
//     viennacl::copy(host_A, A);
//   
//     viennacl::linalg::lanczos_tag ltag(0.75,    // Select a power of 0.75 as the tolerance for the machine precision.
//                                        10,      // Compute (approximations to) the 10 largest eigenvalues
//                                        viennacl::linalg::lanczos_tag::partial_reorthogonalization, // use partial reorthogonalization
//                                        30);   // Maximum size of the Krylov space
//   
//     std::cout << "Running Lanczos algorithm (eigenvalues only). This might take a while..." << std::endl;
//     std::vector<ScalarType> lanczos_eigenvalues = viennacl::linalg::eig(A, ltag);
//   
//     std::cout << "Running Lanczos algorithm (with eigenvectors). This might take a while..." << std::endl;
//     viennacl::matrix<ScalarType> approx_eigenvectors_A(A.size1(), ltag.num_eigenvalues());
//     lanczos_eigenvalues = viennacl::linalg::eig(A, approx_eigenvectors_A, ltag);
//   
//     for (std::size_t i = 0; i< lanczos_eigenvalues.size(); i++)
//     {
//       //std::cout << "Approx. eigenvalue " << std::setprecision(7) << lanczos_eigenvalues[i];
//		 std::cout << "Approx. eigenvalue " << lanczos_eigenvalues[i];
//   
//       // test approximated eigenvector by computing A*v:
//       viennacl::vector<ScalarType> approx_eigenvector = viennacl::column(approx_eigenvectors_A, static_cast<unsigned int>(i));
//       viennacl::vector<ScalarType> Aq = viennacl::linalg::prod(A, approx_eigenvector);
//       std::cout << " (" << viennacl::linalg::inner_prod(Aq, approx_eigenvector) << " for <Av,v> with approx. eigenvector v)" << std::endl;
//     }
//}
//
//void testViennaCL2(const Eigen::SparseMatrix<double> &S, const Eigen::SparseMatrix<double> &Minv, Eigen::MatrixXd &EigVects, Eigen::VectorXd &EigVals)
//{
//	// If you GPU does not support double precision, use `float` instead of `double`:
//	typedef double     ScalarType;
//	Eigen::SparseMatrix<double> L = Minv*S;
//
//	std::vector< std::map<unsigned int, ScalarType> > host_A(L.rows());
//	for (int i = 0; i < L.outerSize(); i++) {
//		for (Eigen::SparseMatrix<double>::InnerIterator it(L, i); it; ++it) {
//			host_A[it.row()][it.col()] = it.value();
//		}
//	}
//
//
//	//if (!viennacl::io::read_matrix_market_file(host_A, "D:/Nasikun/4_SCHOOL/TU Delft/Programming/Libraries/ViennaCL-1.7.1/examples/testdata/mat65k.mtx"))
//	//{
//	//	std::cout << "Error reading Matrix file" << std::endl;
//	//	return; // EXIT_FAILURE;
//	//}
//
//	viennacl::compressed_matrix<ScalarType> A;
//	viennacl::copy(host_A, A);
//	//viennacl::copy(L, A);
//
//	viennacl::linalg::lanczos_tag ltag(0.75,    // Select a power of 0.75 as the tolerance for the machine precision.
//		4,      // Compute (approximations to) the 10 largest eigenvalues
//		viennacl::linalg::lanczos_tag::partial_reorthogonalization, // use partial reorthogonalization
//		30);   // Maximum size of the Krylov space
//
//	std::cout << "Running Lanczos algorithm (eigenvalues only). This might take a while..." << std::endl;
//	std::vector<ScalarType> lanczos_eigenvalues = viennacl::linalg::eig(A, ltag);
//
//	std::cout << "Running Lanczos algorithm (with eigenvectors). This might take a while..." << std::endl;
//	viennacl::matrix<ScalarType> approx_eigenvectors_A(A.size1(), ltag.num_eigenvalues());
//	lanczos_eigenvalues = viennacl::linalg::eig(A, approx_eigenvectors_A, ltag);
//
//	for (std::size_t i = 0; i< lanczos_eigenvalues.size(); i++)
//	{
//		//std::cout << "Approx. eigenvalue " << std::setprecision(7) << lanczos_eigenvalues[i];
//		std::cout << "Approx. eigenvalue " << lanczos_eigenvalues[i];
//
//		// test approximated eigenvector by computing A*v:
//		viennacl::vector<ScalarType> approx_eigenvector = viennacl::column(approx_eigenvectors_A, static_cast<unsigned int>(i));
//		viennacl::vector<ScalarType> Aq = viennacl::linalg::prod(A, approx_eigenvector);
//		std::cout << " (" << viennacl::linalg::inner_prod(Aq, approx_eigenvector) << " for <Av,v> with approx. eigenvector v)" << std::endl;
//	}
//}