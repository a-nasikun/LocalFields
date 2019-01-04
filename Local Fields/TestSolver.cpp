#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "TestSolver.h"




/* ============================================================================
** ============================ INTEL MKL =====================================
** ============================================================================*/

/* Auxiliary routine: printing a matrix */
void print_matrix(char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda) {
	MKL_INT i, j;
	printf("\n %s\n", desc);
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) printf(" %6.2f", a[i + j*lda]);
		printf("\n");
	}
}

/* Auxiliary routine: printing a vector of integers */
void print_int_vector(char* desc, MKL_INT n, MKL_INT* a) {
	MKL_INT j;
	printf("\n %s\n", desc);
	for (j = 0; j < n; j++) printf(" %6i", a[j]);
	printf("\n");
}

/* Testing LAPACKE LDLT Solver */
void testLAPACKEdsysv()
{
	const int N = 5;
	const int NRHS = 3;
	const int LDA = N;
	const int LDB = N;

	/* Locals */
	MKL_INT n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info;
	/* Local arrays */
	MKL_INT ipiv[N];
	/*double a[LDA*N] = {
		-5.86,  3.99, -5.93, -2.82,  7.69,
		0.00,  4.46,  2.58,  4.42,  4.61,
		0.00,  0.00, -8.52,  8.57,  7.69,
		0.00,  0.00,  0.00,  3.72,  8.07,
		0.00,  0.00,  0.00,  0.00,  9.83
	};*/
	double a[LDA*N] = {
		-5.86,  3.99, -5.93, -2.82,  7.69,
		3.99,  4.46,  2.58,  4.42,  4.61,
		-5.93,  2.58, -8.52,  8.57,  7.69,
		-2.82,  4.42,  8.57,  3.72,  8.07,
		7.69,  4.61,  7.69,  8.07,  9.83
	};
	double b[LDB*NRHS] = {
		1.32,  2.22,  0.12, -6.41,  6.33,
		-6.33,  1.69, -1.56, -9.49, -3.67,
		-8.77, -8.33,  9.54,  9.56,  7.48
	};
	/* Executable statements */
	printf("LAPACKE_dsysv (column-major, high-level) Example Program Results\n");

	/* Print vector b*/
	print_matrix("Vector B", ldb, nrhs, b, ldb);

	/* Solve the equations A*X = B */
	info = LAPACKE_dsysv(LAPACK_COL_MAJOR, 'L', n, nrhs, a, lda, ipiv,
		b, ldb);
	/* Check for the exact singularity */
	if (info > 0) {
		printf("The element of the diagonal factor ");
		printf("D(%i,%i) is zero, so that D is singular;\n", info, info);
		printf("the solution could not be computed.\n");
		exit(1);
	}
	/* Print A */
	print_matrix("Matrix X", n, n, a, lda);
	
	/* Print solution */
	print_matrix("Solution", n, nrhs, b, ldb);
	/* Print details of factorization */
	print_matrix("Details of factorization", n, n, a, lda);
	/* Print pivot indices */
	print_int_vector("Pivot indices", n, ipiv);
}

/* ============================================================================
** ================================== MKL PARDISO =============================
** ============================================================================*/
void testMKL_Pardiso()
{
	/* Matrix data. */
	MKL_INT n = 8;
	MKL_INT ia[9] = { 1, 5, 8, 10, 12, 15, 17, 18, 19 };
	MKL_INT ja[18] =
	{ 1,    3,       6, 7,
		2, 3,    5,
		3,             8,
		4,       7,
		5, 6, 7,
		6,    8,
		7,
		8
	};
	double a[18] =
	{ 7.0,      1.0,           2.0, 7.0,
		-4.0, 8.0,      2.0,
		1.0,                     5.0,
		7.0,           9.0,
		5.0, 1.0, 5.0,
		-1.0,      5.0,
		11.0,
		5.0
	};
	MKL_INT mtype = -2;       /* Real symmetric matrix */
							  /* RHS and solution vectors. */
	double b[8], x[8];
	MKL_INT nrhs = 1;     /* Number of right hand sides. */
						  /* Internal solver memory pointer pt, */
						  /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
						  /* or void *pt[64] should be OK on both architectures */
	void *pt[64];
	/* Pardiso control parameters. */
	MKL_INT iparm[64];
	MKL_INT maxfct, mnum, phase, error, msglvl;
	/* Auxiliary variables. */
	MKL_INT i;
	double ddum;          /* Double dummy */
	MKL_INT idum;         /* Integer dummy. */
						  /* -------------------------------------------------------------------- */
						  /* .. Setup Pardiso control parameters. */
						  /* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++)
	{
		iparm[i] = 0;
	}
	iparm[0] = 1;         /* No solver default */
	iparm[1] = 2;         /* Fill-in reordering from METIS */
	iparm[3] = 0;         /* No iterative-direct algorithm */
	iparm[4] = 0;         /* No user fill-in reducing permutation */
	iparm[5] = 0;         /* Write solution into x */
	iparm[6] = 0;         /* Not in use */
	iparm[7] = 2;         /* Max numbers of iterative refinement steps */
	iparm[8] = 0;         /* Not in use */
	iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
	iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
	iparm[11] = 0;        /* Not in use */
	iparm[12] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
	iparm[13] = 0;        /* Output: Number of perturbed pivots */
	iparm[14] = 0;        /* Not in use */
	iparm[15] = 0;        /* Not in use */
	iparm[16] = 0;        /* Not in use */
	iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
	iparm[18] = -1;       /* Output: Mflops for LU factorization */
	iparm[19] = 0;        /* Output: Numbers of CG Iterations */
	maxfct = 1;           /* Maximum number of numerical factorizations. */
	mnum = 1;         /* Which factorization to use. */
	msglvl = 1;           /* Print statistical information in file */
	error = 0;            /* Initialize error flag */
						  /* -------------------------------------------------------------------- */
						  /* .. Initialize the internal solver memory pointer. This is only */
						  /* necessary for the FIRST call of the PARDISO solver. */
						  /* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++)
	{
		pt[i] = 0;
	}
	/* -------------------------------------------------------------------- */
	/* .. Reordering and Symbolic Factorization. This step also allocates */
	/* all memory that is necessary for the factorization. */
	/* -------------------------------------------------------------------- */
	phase = 11;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		printf("\nERROR during symbolic factorization: %d", error);
		exit(1);
	}
	printf("\nReordering completed ... ");
	printf("\nNumber of nonzeros in factors = %d", iparm[17]);
	printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
	/* -------------------------------------------------------------------- */
	/* .. Numerical factorization. */
	/* -------------------------------------------------------------------- */
	phase = 22;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		printf("\nERROR during numerical factorization: %d", error);
		exit(2);
	}
	printf("\nFactorization completed ... ");
	/* -------------------------------------------------------------------- */
	/* .. Back substitution and iterative refinement. */
	/* -------------------------------------------------------------------- */
	phase = 33;
	iparm[7] = 2;         /* Max numbers of iterative refinement steps. */
						  /* Set right hand side to one. */
	for (i = 0; i < n; i++)
	{
		b[i] = 1;
	}
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
	if (error != 0)
	{
		printf("\nERROR during solution: %d", error);
		exit(3);
	}
	printf("\nSolve completed ... ");
	printf("\nThe solution of the system is: ");
	for (i = 0; i < n; i++)
	{
		printf("\n x [%d] = % f", i, x[i]);
	}
	printf("\n");
	/* -------------------------------------------------------------------- */
	/* .. Termination and release of memory. */
	/* -------------------------------------------------------------------- */
	phase = -1;           /* Release internal memory. */
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, &ddum, ia, ja, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);
}

void solveLDLT_MKLPardiso(const Eigen::SparseMatrix<double> &A, const Eigen::MatrixXd& B, Eigen::MatrixXd& Xf)
{
	/* Matrix data. */
	MKL_INT n = 8;
	MKL_INT ia[9] = { 1, 5, 8, 10, 12, 15, 17, 18, 19 };
	MKL_INT ja[18] =
	{ 1,    3,       6, 7,
		2, 3,    5,
		3,             8,
		4,       7,
		5, 6, 7,
		6,    8,
		7,
		8
	};
	double a[18] =
	{ 7.0,      1.0,           2.0, 7.0,
		-4.0, 8.0,      2.0,
		1.0,                     5.0,
		7.0,           9.0,
		5.0, 1.0, 5.0,
		-1.0,      5.0,
		11.0,
		5.0
	};
	MKL_INT mtype = -2;       /* Real symmetric matrix */
							  /* RHS and solution vectors. */
	double b[8], x[8];
	MKL_INT nrhs = 1;     /* Number of right hand sides. */
						  /* Internal solver memory pointer pt, */
						  /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
						  /* or void *pt[64] should be OK on both architectures */
	void *pt[64];
	/* Pardiso control parameters. */
	MKL_INT iparm[64];
	MKL_INT maxfct, mnum, phase, error, msglvl;
	/* Auxiliary variables. */
	MKL_INT i;
	double ddum;          /* Double dummy */
	MKL_INT idum;         /* Integer dummy. */
						  /* -------------------------------------------------------------------- */
						  /* .. Setup Pardiso control parameters. */
						  /* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++)
	{
		iparm[i] = 0;
	}
	iparm[0] = 1;         /* No solver default */
	iparm[1] = 2;         /* Fill-in reordering from METIS */
	iparm[3] = 0;         /* No iterative-direct algorithm */
	iparm[4] = 0;         /* No user fill-in reducing permutation */
	iparm[5] = 0;         /* Write solution into x */
	iparm[6] = 0;         /* Not in use */
	iparm[7] = 2;         /* Max numbers of iterative refinement steps */
	iparm[8] = 0;         /* Not in use */
	iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
	iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
	iparm[11] = 0;        /* Not in use */
	iparm[12] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
	iparm[13] = 0;        /* Output: Number of perturbed pivots */
	iparm[14] = 0;        /* Not in use */
	iparm[15] = 0;        /* Not in use */
	iparm[16] = 0;        /* Not in use */
	iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
	iparm[18] = -1;       /* Output: Mflops for LU factorization */
	iparm[19] = 0;        /* Output: Numbers of CG Iterations */
	maxfct = 1;           /* Maximum number of numerical factorizations. */
	mnum = 1;         /* Which factorization to use. */
	msglvl = 1;           /* Print statistical information in file */
	error = 0;            /* Initialize error flag */
						  /* -------------------------------------------------------------------- */
						  /* .. Initialize the internal solver memory pointer. This is only */
						  /* necessary for the FIRST call of the PARDISO solver. */
						  /* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++)
	{
		pt[i] = 0;
	}
	/* -------------------------------------------------------------------- */
	/* .. Reordering and Symbolic Factorization. This step also allocates */
	/* all memory that is necessary for the factorization. */
	/* -------------------------------------------------------------------- */
	phase = 11;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		printf("\nERROR during symbolic factorization: %d", error);
		exit(1);
	}
	printf("\nReordering completed ... ");
	printf("\nNumber of nonzeros in factors = %d", iparm[17]);
	printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
	/* -------------------------------------------------------------------- */
	/* .. Numerical factorization. */
	/* -------------------------------------------------------------------- */
	phase = 22;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		printf("\nERROR during numerical factorization: %d", error);
		exit(2);
	}
	printf("\nFactorization completed ... ");
	/* -------------------------------------------------------------------- */
	/* .. Back substitution and iterative refinement. */
	/* -------------------------------------------------------------------- */
	phase = 33;
	iparm[7] = 2;         /* Max numbers of iterative refinement steps. */
						  /* Set right hand side to one. */
	for (i = 0; i < n; i++)
	{
		b[i] = 1;
	}
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
	if (error != 0)
	{
		printf("\nERROR during solution: %d", error);
		exit(3);
	}
	printf("\nSolve completed ... ");
	printf("\nThe solution of the system is: ");
	for (i = 0; i < n; i++)
	{
		printf("\n x [%d] = % f", i, x[i]);
	}
	printf("\n");
	/* -------------------------------------------------------------------- */
	/* .. Termination and release of memory. */
	/* -------------------------------------------------------------------- */
	phase = -1;           /* Release internal memory. */
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, &ddum, ia, ja, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);
}

void solveLDLT_EigenPardiso(const Eigen::SparseMatrix<double> &A, const Eigen::MatrixXd& B, Eigen::MatrixXd& Xf)
{

}

/* ============================================================================
** ================================= CUDA =====================================
** ============================================================================*/

void printMatrix(int m, int n, const double*A, int lda, const char* name)
{
	for (int row = 0; row < m; row++) {
		for (int col = 0; col < n; col++) {
			double Areg = A[row + col*lda];
			printf("%s(%d,%d) = %f\n", name, row + 1, col + 1, Areg);
		}
	}
}

void testCUDA_LULinearSolver()
{
	// Defining necessary variables
	cusolverDnHandle_t cusolverH = NULL;
	cudaStream_t stream = NULL;

	cusolverStatus_t status = CUSOLVER_STATUS_SUCCESS;
	cudaError_t cudaStat1 = cudaSuccess;
	cudaError_t cudaStat2 = cudaSuccess;
	cudaError_t cudaStat3 = cudaSuccess;
	cudaError_t cudaStat4 = cudaSuccess;
	const int m = 3;
	const int lda = m;
	const int ldb = m;

	// Setting up the data
	double A[lda*m] = { 1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 10.0 };
	double B[m] = { 1.0, 2.0, 3.0 };
	double X[m]; /* X = A\B */
	double LU[lda*m]; /* L and U */
	int Ipiv[m];      /* host copy of pivoting sequence */
	int info = 0;     /* host copy of error info */

	double *d_A = NULL; /* device copy of A */
	double *d_B = NULL; /* device copy of B */
	int *d_Ipiv = NULL; /* pivoting sequence */
	int *d_info = NULL; /* error info */
	int  lwork = 0;     /* size of workspace */
	double *d_work = NULL; /* device workspace for getrf */

	const int pivot_on = 1;

	printf("example of getrf \n");

	if (pivot_on) {
		printf("pivot is on : compute P*A = L*U \n");
	}
	else {
		printf("pivot is off: compute A = L*U (not numerically stable)\n");
	}

	printf("A = (matlab base-1)\n");
	printMatrix(m, m, A, lda, "A");
	printf("=====\n");

	printf("B = (matlab base-1)\n");
	printMatrix(m, 1, B, ldb, "B");
	printf("=====\n");

	/* step 1: create cusolver handle, bind a stream */
	status = cusolverDnCreate(&cusolverH);
	assert(CUSOLVER_STATUS_SUCCESS == status);

	cudaStat1 = cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
	assert(cudaSuccess == cudaStat1);

	status = cusolverDnSetStream(cusolverH, stream);
	assert(CUSOLVER_STATUS_SUCCESS == status);

	/* step 2: copy A to device */
	cudaStat1 = cudaMalloc((void**)&d_A, sizeof(double) * lda * m);
	cudaStat2 = cudaMalloc((void**)&d_B, sizeof(double) * m);
	cudaStat2 = cudaMalloc((void**)&d_Ipiv, sizeof(int) * m);
	cudaStat4 = cudaMalloc((void**)&d_info, sizeof(int));
	assert(cudaSuccess == cudaStat1);
	assert(cudaSuccess == cudaStat2);
	assert(cudaSuccess == cudaStat3);
	assert(cudaSuccess == cudaStat4);

	cudaStat1 = cudaMemcpy(d_A, A, sizeof(double)*lda*m, cudaMemcpyHostToDevice);
	cudaStat2 = cudaMemcpy(d_B, B, sizeof(double)*m, cudaMemcpyHostToDevice);
	assert(cudaSuccess == cudaStat1);
	assert(cudaSuccess == cudaStat2);

	/* step 3: query working space of getrf */
	status = cusolverDnDgetrf_bufferSize(
		cusolverH,
		m,
		m,
		d_A,
		lda,
		&lwork);
	assert(CUSOLVER_STATUS_SUCCESS == status);

	cudaStat1 = cudaMalloc((void**)&d_work, sizeof(double)*lwork);
	assert(cudaSuccess == cudaStat1);

	/* step 4: LU factorization */
	if (pivot_on) {
		status = cusolverDnDgetrf(
			cusolverH,
			m,
			m,
			d_A,
			lda,
			d_work,
			d_Ipiv,
			d_info);
	}
	else {
		status = cusolverDnDgetrf(
			cusolverH,
			m,
			m,
			d_A,
			lda,
			d_work,
			NULL,
			d_info);
	}
	cudaStat1 = cudaDeviceSynchronize();
	assert(CUSOLVER_STATUS_SUCCESS == status);
	assert(cudaSuccess == cudaStat1);

	if (pivot_on) {
		cudaStat1 = cudaMemcpy(Ipiv, d_Ipiv, sizeof(int)*m, cudaMemcpyDeviceToHost);
	}
	cudaStat2 = cudaMemcpy(LU, d_A, sizeof(double)*lda*m, cudaMemcpyDeviceToHost);
	cudaStat3 = cudaMemcpy(&info, d_info, sizeof(int), cudaMemcpyDeviceToHost);
	assert(cudaSuccess == cudaStat1);
	assert(cudaSuccess == cudaStat2);
	assert(cudaSuccess == cudaStat3);

	if (0 > info) {
		printf("%d-th parameter is wrong \n", -info);
		exit(1);
	}
	if (pivot_on) {
		printf("pivoting sequence, matlab base-1\n");
		for (int j = 0; j < m; j++) {
			printf("Ipiv(%d) = %d\n", j + 1, Ipiv[j]);
		}
	}
	printf("L and U = (matlab base-1)\n");
	printMatrix(m, m, LU, lda, "LU");
	printf("=====\n");

	/*
	* step 5: solve A*X = B
	*       | 1 |       | -0.3333 |
	*   B = | 2 |,  X = |  0.6667 |
	*       | 3 |       |  0      |
	*
	*/
	if (pivot_on) {
		status = cusolverDnDgetrs(
			cusolverH,
			CUBLAS_OP_N,
			m,
			1, /* nrhs */
			d_A,
			lda,
			d_Ipiv,
			d_B,
			ldb,
			d_info);
	}
	else {
		status = cusolverDnDgetrs(
			cusolverH,
			CUBLAS_OP_N,
			m,
			1, /* nrhs */
			d_A,
			lda,
			NULL,
			d_B,
			ldb,
			d_info);
	}
	cudaStat1 = cudaDeviceSynchronize();
	assert(CUSOLVER_STATUS_SUCCESS == status);
	assert(cudaSuccess == cudaStat1);

	cudaStat1 = cudaMemcpy(X, d_B, sizeof(double)*m, cudaMemcpyDeviceToHost);
	assert(cudaSuccess == cudaStat1);

	printf("X = (matlab base-1)\n");
	printMatrix(m, 1, X, ldb, "X");
	printf("=====\n");

	/* free resources */
	if (d_A) cudaFree(d_A);
	if (d_B) cudaFree(d_B);
	if (d_Ipiv) cudaFree(d_Ipiv);
	if (d_info) cudaFree(d_info);
	if (d_work) cudaFree(d_work);

	if (cusolverH) cusolverDnDestroy(cusolverH);
	if (stream) cudaStreamDestroy(stream);

	cudaDeviceReset();
}

void testLDLTLinearSolver()
{
	// Defining necessary variables
	cusolverDnHandle_t cusolverH = NULL;
	cudaStream_t stream = NULL;

	cusolverStatus_t status = CUSOLVER_STATUS_SUCCESS;
	cudaError_t cudaStat1 = cudaSuccess;
	cudaError_t cudaStat2 = cudaSuccess;
	cudaError_t cudaStat3 = cudaSuccess;
	cudaError_t cudaStat4 = cudaSuccess;
	const int m = 3;
	const int lda = m;
	const int ldb = m;

	// Setting up the data
	double A[lda*m] = { 1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 10.0 };
	double B[m] = { 1.0, 2.0, 3.0 };
	double X[m]; /* X = A\B */
	double LDLT[lda*m]; /* LDLT */
	int Ipiv[m];      /* host copy of pivoting sequence */
	int info = 0;     /* host copy of error info */

	double *d_A = NULL; /* device copy of A */
	double *d_B = NULL; /* device copy of B */
	int *d_Ipiv = NULL; /* pivoting sequence */
	int *d_info = NULL; /* error info */
	int  lwork = 0;     /* size of workspace */
	double *d_work = NULL; /* device workspace for getrf */

	const int pivot_on = 1;

	printf("example of getrf \n");

	if (pivot_on) {
		printf("pivot is on : compute P*A = L*U \n");
	}
	else {
		printf("pivot is off: compute A = L*U (not numerically stable)\n");
	}

	printf("A = (matlab base-1)\n");
	//printMatrix(m, m, A, lda, "A");
	print_matrix("A", m, m, A, lda);
	printf("=====\n");

	printf("B = (matlab base-1)\n");
	//printMatrix(m, 1, B, ldb, "B");
	print_matrix("b", m, 1, B, ldb);
	printf("=====\n");

	/* step 1: create cusolver handle, bind a stream */
	status = cusolverDnCreate(&cusolverH);
	assert(CUSOLVER_STATUS_SUCCESS == status);

	cudaStat1 = cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
	assert(cudaSuccess == cudaStat1);

	status = cusolverDnSetStream(cusolverH, stream);
	assert(CUSOLVER_STATUS_SUCCESS == status);

	/* step 2: copy A to device */
	cudaStat1 = cudaMalloc((void**)&d_A, sizeof(double) * lda * m);
	cudaStat2 = cudaMalloc((void**)&d_B, sizeof(double) * m);
	cudaStat2 = cudaMalloc((void**)&d_Ipiv, sizeof(int) * m);
	cudaStat4 = cudaMalloc((void**)&d_info, sizeof(int));
	assert(cudaSuccess == cudaStat1);
	assert(cudaSuccess == cudaStat2);
	assert(cudaSuccess == cudaStat3);
	assert(cudaSuccess == cudaStat4);

	cudaStat1 = cudaMemcpy(d_A, A, sizeof(double)*lda*m, cudaMemcpyHostToDevice);
	cudaStat2 = cudaMemcpy(d_B, B, sizeof(double)*m, cudaMemcpyHostToDevice);
	assert(cudaSuccess == cudaStat1);
	assert(cudaSuccess == cudaStat2);

	/* step 3: query working space of getrf */	
	status = cusolverDnDsytrf_bufferSize(cusolverH,
			m,
			d_A,
			lda,
			&lwork);
	assert(CUSOLVER_STATUS_SUCCESS == status);

	cudaStat1 = cudaMalloc((void**)&d_work, sizeof(double)*lwork);
	assert(cudaSuccess == cudaStat1);

	/* Step 4: Factorization */
	status = cusolverDnDsytrf(cusolverH,
		CUBLAS_FILL_MODE_UPPER,
		m,
		d_A,
		lda,
		d_Ipiv,
		d_work,
		lwork,
		d_info);	
	cudaStat1 = cudaDeviceSynchronize();
	assert(CUSOLVER_STATUS_SUCCESS == status);
	assert(cudaSuccess == cudaStat1);

	
	cudaStat1 = cudaMemcpy(Ipiv, d_Ipiv, sizeof(int)*m, cudaMemcpyDeviceToHost);
	cudaStat2 = cudaMemcpy(LDLT, d_A, sizeof(double)*lda*m, cudaMemcpyDeviceToHost);
	cudaStat3 = cudaMemcpy(&info, d_info, sizeof(int), cudaMemcpyDeviceToHost);
	assert(cudaSuccess == cudaStat1);
	assert(cudaSuccess == cudaStat2);
	assert(cudaSuccess == cudaStat3);

	if (0 > info) {
		printf("%d-th parameter is wrong \n", -info);
		exit(1);
	}
	if (pivot_on) {
		printf("pivoting sequence, matlab base-1\n");
		for (int j = 0; j < m; j++) {
			printf("Ipiv(%d) = %d\n", j + 1, Ipiv[j]);
		}
	}
	printf("L and U = (matlab base-1)\n");
	//printMatrix(m, m, LU, lda, "LU");
	print_matrix("LDLT", m, m, LDLT, lda);
	printf("=====\n");

	/* Step 5: Solve Linear System */
	status = cusolverDnDgetrs(
		cusolverH,
		CUBLAS_OP_N,
		m,
		1, /* nrhs */
		d_A,
		lda,
		d_Ipiv,
		d_B,
		ldb,
		d_info);

	//status = cusolverDnDpotrsBatched(
	//		cusolverH,
	//		CUBLAS_FILL_MODE_UPPER,
	//		m,
	//		1,
	//		&d_A,
	//		lda,
	//		&d_B,
	//		ldb,
	//		d_info,
	//		m * (m+1) /2 );

	/*
	* step 5: solve A*X = B
	*       | 1 |       | -0.3333 |
	*   B = | 2 |,  X = |  0.6667 |
	*       | 3 |       |  0      |
	*
	*/
	
	cudaStat1 = cudaDeviceSynchronize();
	assert(CUSOLVER_STATUS_SUCCESS == status);
	assert(cudaSuccess == cudaStat1);

	cudaStat1 = cudaMemcpy(X, d_B, sizeof(double)*m, cudaMemcpyDeviceToHost);
	assert(cudaSuccess == cudaStat1);

	printf("X = (matlab base-1)\n");
	//printMatrix(m, 1, X, ldb, "X");
	print_matrix("X", m, 1, X, ldb);
	printf("=====\n");

	/* free resources */
	if (d_A) cudaFree(d_A);
	if (d_B) cudaFree(d_B);
	if (d_Ipiv) cudaFree(d_Ipiv);
	if (d_info) cudaFree(d_info);
	if (d_work) cudaFree(d_work);

	if (cusolverH) cusolverDnDestroy(cusolverH);
	if (stream) cudaStreamDestroy(stream);

	cudaDeviceReset();
}

void solveLDLTinCUDA(const Eigen::SparseMatrix<double> &MA, const Eigen::MatrixXd& MB, Eigen::MatrixXd& Xf)
{
	printf("This part is called \n");

	// Defining necessary variables
	cusolverDnHandle_t cusolverH = NULL;
	cudaStream_t stream = NULL;

	cusolverStatus_t status = CUSOLVER_STATUS_SUCCESS;
	cudaError_t cudaStat1 = cudaSuccess;
	cudaError_t cudaStat2 = cudaSuccess;
	cudaError_t cudaStat3 = cudaSuccess;
	cudaError_t cudaStat4 = cudaSuccess;
	const int m = MA.rows();
	const int lda = m;
	const int ldb = m;

	// Setting up the data
	double *A = (double*)std::malloc(m*lda * sizeof(double));
	Eigen::MatrixXd MA2(MA);
	A = MA2.data();
	double *B1 = (double*)std::malloc(m * sizeof(double));
	Eigen::MatrixXd BB(MB);
	B1 = BB.col(0).data();
	double *B2 = (double*)std::malloc(m * sizeof(double));
	B2 = BB.col(1).data();
	//double X[m]; /* X = A\B */
	double *X = (double*)std::malloc(m * sizeof(double));
	//double LDLT[lda*m]; /* LDLT */
	double *LDLT = (double*)std::malloc(m*lda * sizeof(double));
	//int Ipiv[m];      /* host copy of pivoting sequence */
	int *Ipiv = (int*)std::malloc(m * sizeof(int));
	int info = 0;     /* host copy of error info */

	double *d_A = NULL; /* device copy of A */
	double *d_B = NULL; /* device copy of B */
	int *d_Ipiv = NULL; /* pivoting sequence */
	int *d_info = NULL; /* error info */
	int  lwork = 0;     /* size of workspace */
	double *d_work = NULL; /* device workspace for getrf */

	const int pivot_on = 1;

	

	printf("A = (matlab base-1)\n");
	//printMatrix(m, m, A, lda, "A");
	//print_matrix("A", m, m, A, lda);
	printf("=====\n");

	printf("B = (matlab base-1)\n");
	//printMatrix(m, 1, B, ldb, "B");
	//print_matrix("b", m, 1, B1, ldb);
	printf("=====\n");

	/* step 1: create cusolver handle, bind a stream */
	status = cusolverDnCreate(&cusolverH);
	assert(CUSOLVER_STATUS_SUCCESS == status);

	cudaStat1 = cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
	assert(cudaSuccess == cudaStat1);

	status = cusolverDnSetStream(cusolverH, stream);
	assert(CUSOLVER_STATUS_SUCCESS == status);

	/* step 2: copy A to device */
	cout << "Copying to GPU " << std::endl; 
	cudaStat1 = cudaMalloc((void**)&d_A, sizeof(double) * lda * m);
	cudaStat2 = cudaMalloc((void**)&d_B, sizeof(double) * m);
	cudaStat2 = cudaMalloc((void**)&d_Ipiv, sizeof(int) * m);
	cudaStat4 = cudaMalloc((void**)&d_info, sizeof(int));
	assert(cudaSuccess == cudaStat1);
	assert(cudaSuccess == cudaStat2);
	assert(cudaSuccess == cudaStat3);
	assert(cudaSuccess == cudaStat4);

	cudaStat1 = cudaMemcpy(d_A, A, sizeof(double)*lda*m, cudaMemcpyHostToDevice);
	cudaStat2 = cudaMemcpy(d_B, B1, sizeof(double)*m, cudaMemcpyHostToDevice);
	assert(cudaSuccess == cudaStat1);
	assert(cudaSuccess == cudaStat2);

	cout << "Getting the buffer size." << std::endl;
	/* step 3: query working space of getrf */
	status = cusolverDnDsytrf_bufferSize(cusolverH,
		m,
		d_A,
		lda,
		&lwork);
	assert(CUSOLVER_STATUS_SUCCESS == status);

	cudaStat1 = cudaMalloc((void**)&d_work, sizeof(double)*lwork);
	assert(cudaSuccess == cudaStat1);


	
	/* Step 4: Factorization */
	cout << "Factorization" << std::endl;
	status = cusolverDnDsytrf(cusolverH,
		CUBLAS_FILL_MODE_UPPER,
		m,
		d_A,
		lda,
		d_Ipiv,
		d_work,
		lwork,
		d_info);
	cudaStat1 = cudaDeviceSynchronize();
	assert(CUSOLVER_STATUS_SUCCESS == status);
	assert(cudaSuccess == cudaStat1);

	cudaStat1 = cudaMemcpy(Ipiv, d_Ipiv, sizeof(int)*m, cudaMemcpyDeviceToHost);
	cudaStat2 = cudaMemcpy(LDLT, d_A, sizeof(double)*lda*m, cudaMemcpyDeviceToHost);
	cudaStat3 = cudaMemcpy(&info, d_info, sizeof(int), cudaMemcpyDeviceToHost);
	assert(cudaSuccess == cudaStat1);
	assert(cudaSuccess == cudaStat2);
	assert(cudaSuccess == cudaStat3);
		
	//printf("pivoting sequence, matlab base-1\n");
	//for (int j = 0; j < m; j++) {
	//	printf("Ipiv(%d) = %d\n", j + 1, Ipiv[j]);
	//}

	printf("LDLT = (matlab base-1)\n");
	//printMatrix(m, m, LU, lda, "LU");
	//print_matrix("LDLT", m, m, LDLT, lda);
	printf("=====\n");

	/* Step 5: Solve Linear System */
	cout << "Solving 1 system." << endl; 
	status = cusolverDnDgetrs(
		cusolverH,
		CUBLAS_OP_N,
		m,
		1, /* nrhs */
		d_A,
		lda,
		d_Ipiv,
		d_B,
		ldb,
		d_info);
	
	cudaStat1 = cudaDeviceSynchronize();
	assert(CUSOLVER_STATUS_SUCCESS == status);
	assert(cudaSuccess == cudaStat1);

	cudaStat1 = cudaMemcpy(X, d_B, sizeof(double)*m, cudaMemcpyDeviceToHost);
	assert(cudaSuccess == cudaStat1);
	
	printf("X = (matlab base-1)\n");
	//printMatrix(m, 1, X, ldb, "X");
	//print_matrix("X", m, 1, X, ldb);
	printf("=====\n");

	//Xf = Eigen::Map<Eigen::MatrixXd, 0, Eigen::OuterStride<>>(X, m, lda, Eigen::OuterStride<>(m));
	Xf.resize(MA.rows(), 2);
	Xf.col(0) = Eigen::Map<Eigen::VectorXd>(X, m);;

	// Solving second linear system
	cout << "Soving 2nd system" << endl; 
	cudaStat2 = cudaMemcpy(d_B, B2, sizeof(double)*m, cudaMemcpyHostToDevice);
	status = cusolverDnDgetrs(
		cusolverH,
		CUBLAS_OP_N,
		m,
		1, /* nrhs */
		d_A,
		lda,
		d_Ipiv,
		d_B,
		ldb,
		d_info);

	cudaStat1 = cudaDeviceSynchronize();
	assert(CUSOLVER_STATUS_SUCCESS == status);
	assert(cudaSuccess == cudaStat1);

	cudaStat1 = cudaMemcpy(X, d_B, sizeof(double)*m, cudaMemcpyDeviceToHost);
	assert(cudaSuccess == cudaStat1);

	printf("X = (matlab base-1)\n");
	//printMatrix(m, 1, X, ldb, "X");
	//print_matrix("X", m, 1, X, ldb);
	printf("=====\n");	
	cout << "Mpa to eigenformat" << endl; 
	Xf.col(1) = Eigen::Map<Eigen::VectorXd>(X, m);

	/* free resources */
	if (d_A) cudaFree(d_A);
	if (d_B) cudaFree(d_B);
	if (d_Ipiv) cudaFree(d_Ipiv);
	if (d_info) cudaFree(d_info);
	if (d_work) cudaFree(d_work);

	if (cusolverH) cusolverDnDestroy(cusolverH);
	if (stream) cudaStreamDestroy(stream);

	cudaDeviceReset();
}

void solveLUinCUDA(const Eigen::SparseMatrix<double> &AA, const Eigen::MatrixXd& BB, Eigen::MatrixXd& Xf)
{
	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					duration;
	t1 = chrono::high_resolution_clock::now();

	// Defining necessary variables
	cusolverDnHandle_t cusolverH = NULL;
	cudaStream_t stream = NULL;

	cusolverStatus_t status = CUSOLVER_STATUS_SUCCESS;
	cudaError_t cudaStat1 = cudaSuccess;
	cudaError_t cudaStat2 = cudaSuccess;
	cudaError_t cudaStat3 = cudaSuccess;
	cudaError_t cudaStat4 = cudaSuccess;
	const int m = AA.rows();
	const int lda = m;
	const int ldb = m;

	// Setting up the data
	//double A[lda*m] = { 1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 10.0 };
	double *A = (double*)std::malloc(m*lda * sizeof(double));
	Eigen::MatrixXd MA(AA);
	A = MA.data();
	//double B[m] = { 1.0, 2.0, 3.0 };
	double *B1 = (double*)std::malloc(m * sizeof(double));
	Eigen::MatrixXd BB1(BB);
	B1 = BB1.col(0).data();
	double *B2 = (double*)std::malloc(m * sizeof(double));
	B2 = BB1.col(1).data();
	//double X[m]; /* X = A\B */
	double *X = (double*)std::malloc(m * sizeof(double));		
	//double LU[lda*m]; /* L and U */
	double *LU = (double*)std::malloc(m*lda * sizeof(double));
	//int Ipiv[m];      /* host copy of pivoting sequence */
	int *Ipiv = (int*)std::malloc(m * sizeof(int));
	int info = 0;     /* host copy of error info */
	Xf.resize(m, 2);

	double *d_A = NULL; /* device copy of A */
	double *d_B = NULL; /* device copy of B */
	int *d_Ipiv = NULL; /* pivoting sequence */
	int *d_info = NULL; /* error info */
	int  lwork = 0;     /* size of workspace */
	double *d_work = NULL; /* device workspace for getrf */

	const int pivot_on = 1;

	printf("example of getrf \n");

	if (pivot_on) {
		printf("pivot is on : compute P*A = L*U \n");
	}
	else {
		printf("pivot is off: compute A = L*U (not numerically stable)\n");
	}

	printf("A = (matlab base-1)\n");
	//printMatrix(m, m, A, lda, "A");
	printf("=====\n");

	printf("B = (matlab base-1)\n");
	//printMatrix(m, 1, B, ldb, "B");
	printf("=====\n");

	/* step 1: create cusolver handle, bind a stream */
	status = cusolverDnCreate(&cusolverH);
	assert(CUSOLVER_STATUS_SUCCESS == status);

	cudaStat1 = cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
	assert(cudaSuccess == cudaStat1);

	status = cusolverDnSetStream(cusolverH, stream);
	assert(CUSOLVER_STATUS_SUCCESS == status);

	/* step 2: copy A to device */
	cudaStat1 = cudaMalloc((void**)&d_A, sizeof(double) * lda * m);
	cudaStat2 = cudaMalloc((void**)&d_B, sizeof(double) * m);
	cudaStat2 = cudaMalloc((void**)&d_Ipiv, sizeof(int) * m);
	cudaStat4 = cudaMalloc((void**)&d_info, sizeof(int));
	bool b1 = (cudaSuccess == cudaStat1);
	bool b2 = (cudaSuccess == cudaStat2);
	bool b3 = (cudaSuccess == cudaStat3);
	bool b4 = (cudaSuccess == cudaStat4);
	if (!b1 || !b2 || !b3 || !b4) {
		cout << "____Error! Cannot allocate memory at GPU" << endl; 
	}
	else {
		cout << "____Memory Allocation successful" << endl;
	}

	cudaStat1 = cudaMemcpy(d_A, A, sizeof(double)*lda*m, cudaMemcpyHostToDevice);
	cudaStat2 = cudaMemcpy(d_B, B1, sizeof(double)*m, cudaMemcpyHostToDevice);
	b1 = (cudaSuccess == cudaStat1);
	b2 = (cudaSuccess == cudaStat2);
	if (!b1 || !b2 ) {
		cout << "____Error! Cannot copy data to GPU" << endl;
	}
	else {
		cout << "____Data copying successful" << endl;
	}

	/* step 3: query working space of getrf */
	status = cusolverDnDgetrf_bufferSize(
		cusolverH,
		m,
		m,
		d_A,
		lda,
		&lwork);
	b1 = (CUSOLVER_STATUS_SUCCESS == status);
	if (!b1) {
		cout << "____Error! Buffering faileds" << endl;
	}
	else {
		cout << "____Buffering successful" << endl;
	}

	cudaStat1 = cudaMalloc((void**)&d_work, sizeof(double)*lwork);
	assert(cudaSuccess == cudaStat1);

	/* step 4: LU factorization */
	t1 = chrono::high_resolution_clock::now();
	cout << "> Factoring (LU in GPU)... \n";
		status = cusolverDnDgetrf(
			cusolverH,
			m,
			m,
			d_A,
			lda,
			d_work,
			d_Ipiv,
			d_info);
	
	cudaStat1 = cudaDeviceSynchronize();
	b1  = (CUSOLVER_STATUS_SUCCESS == status);
	b2  = (cudaSuccess == cudaStat1);
	if (!b1 || !b2 ) {
		cout << "____Error! LU Factorization fails" << endl;
	}
	else {
		cout << "____LU Factorization successful" << endl;
	}		
	t2 = chrono::high_resolution_clock::now();
	duration = t2 - t1;
	cout << "in " << duration.count() << " seconds" << endl;

	cudaStat1 = cudaMemcpy(Ipiv, d_Ipiv, sizeof(int)*m, cudaMemcpyDeviceToHost);

	cudaStat2 = cudaMemcpy(LU, d_A, sizeof(double)*lda*m, cudaMemcpyDeviceToHost);
	cudaStat3 = cudaMemcpy(&info, d_info, sizeof(int), cudaMemcpyDeviceToHost);
	b1 = (cudaSuccess == cudaStat1);
	b2 = (cudaSuccess == cudaStat2);
	b3 = (cudaSuccess == cudaStat3);
	if (!b1 || !b2 || !b3) {
		cout << "____Error! Copying factorization to GPU fails" << endl;
	}
	else {
		cout << "____Copying LU to GPU successful" << endl;
	}

	if (0 > info) {
		printf("%d-th parameter is wrong \n", -info);
		exit(1);
	}
	
		printf("pivoting sequence, matlab base-1\n");
		for (int j = 0; j < m; j++) {
			//printf("Ipiv(%d) = %d\n", j + 1, Ipiv[j]);
		}
	
	printf("L and U = (matlab base-1)\n");
	//printMatrix(m, m, LU, lda, "LU");
	printf("=====\n");
	
	
	status = cusolverDnDgetrs(
		cusolverH,
		CUBLAS_OP_N,
		m,
		1, /* nrhs */
		d_A,
		lda,
		d_Ipiv,
		d_B,
		ldb,
		d_info);
	
	cudaStat1 = cudaDeviceSynchronize();
	b1 = (CUSOLVER_STATUS_SUCCESS == status);
	b1 = (cudaSuccess == cudaStat1);
	if (!b1 || !b2) {
		cout << "____Error! Cannot Solve linear system" << endl;
	}
	else {
		cout << "____Solving linear system successful" << endl;
	}


	cudaStat1 = cudaMemcpy(X, d_B, sizeof(double)*m, cudaMemcpyDeviceToHost);
	b1 = (cudaSuccess == cudaStat1);
	if (!b1) {
		cout << "____Error! copy data back to CPU" << endl;
	}
	else {
		cout << "____Copying to CPU successful" << endl;
	}

	Xf.col(0) = Eigen::Map<Eigen::VectorXd>(X, m);


	// Second Basis
	cudaStat2 = cudaMemcpy(d_B, B1, sizeof(double)*m, cudaMemcpyHostToDevice);
	assert(cudaSuccess == cudaStat2);

	status = cusolverDnDgetrs(
		cusolverH,
		CUBLAS_OP_N,
		m,
		1, /* nrhs */
		d_A,
		lda,
		d_Ipiv,
		d_B,
		ldb,
		d_info);

	cudaStat1 = cudaDeviceSynchronize();
	assert(CUSOLVER_STATUS_SUCCESS == status);
	assert(cudaSuccess == cudaStat1);

	cudaStat1 = cudaMemcpy(X, d_B, sizeof(double)*m, cudaMemcpyDeviceToHost);
	assert(cudaSuccess == cudaStat1);

	Xf.col(1) = Eigen::Map<Eigen::VectorXd>(X, m);

	printf("X = (matlab base-1)\n");
	//printMatrix(m, 1, X, ldb, "X");
	printf("=====\n");

	/* free resources */
	if (d_A) cudaFree(d_A);
	if (d_B) cudaFree(d_B);
	if (d_Ipiv) cudaFree(d_Ipiv);
	if (d_info) cudaFree(d_info);
	if (d_work) cudaFree(d_work);

	if (cusolverH) cusolverDnDestroy(cusolverH);
	if (stream) cudaStreamDestroy(stream);

	cudaDeviceReset();
}

