#include "./include/global.hpp"
//#include<cula_sparse.h>
#include <cublas_v2.h>
#include <cusparse.h>



int main(int argc, char* argv[])
{
	const int bufsize = 512;
    	char buffer[bufsize];
	int m,n,S;
	double time_st,time_end,time_avg;
	//omp_set_num_threads(2);
//	printf("\n-----------------\nnumber of threads fired = %d\n-----------------\n",(int)omp_get_num_threads());
	if(argc!=2)
	{
		cout<<"Insufficient arguments"<<endl;
		return 1;
	}
	
	graph G;

	cerr<<"Start reading                    ";
//	time_st=dsecnd();
	G.create_graph(argv[1]);
//	time_end=dsecnd();
//	time_avg = (time_end-time_st);
//	cout<<"Success              "<<endl;
//	cerr<<"Reading time                     "<<time_avg<<endl;

	cerr<<"Constructing Matrices            ";
//	time_st=dsecnd();
	G.construct_MNA();
	G.construct_NA();
//	time_end=dsecnd();
//	time_avg = (time_end-time_st);
//	cerr<<"Done                 "<<time_avg<<endl;

//	G.construct_sparse_MNA();
	m=G.node_array.size()-1;
	n=G.voltage_edge_id.size();
	
	cout<<endl;
	cout<<"MATRIX STAT:"<<endl;
	cout<<"Nonzero elements:               "<<G.nonzero<<endl;
	cout<<"Number of Rows:                 "<<m+n<<endl;
	cout<<"Nonzero in G:			"<<G.Gnonzero<<endl;
	cout<<"Number of rows in G:		"<<m<<endl;
	cout<<"Nonzero in P: 			"<<G.Pnonzero<<endl;
	cout<<"Number of rows in P:		"<<m<<endl;


//	printf("\n Nonzero = %d", G.nonzero);
//	printf("\n Rows = %d", m+n);

	cout<<"MAT val:		       "<<endl;
	int i,j;

	G.Mat_val[0] += 100;
	G.Gmat[0] +=100;
/*
	for(i=0;i<G.Gnonzero;i++)
		cout<<" "<<G.Gmat[i];
	cout<<endl;
	for(i=0;i<G.Gnonzero;i++)
		cout<<" "<<G.Gcolumns[i];
	cout<<endl;	
	for(i=0;i<m+1;i++)
		cout<<" "<<G.GrowIndex[i];
	cout<<endl;
	for(i=0;i<m;i++)
		printf(" %.8f", G.b1[i]);
	cout<<endl;
	for(i=0;i<m;i++)
		printf(" %.8f", G.x1[i]);
	cout<<endl;	
	
	
*/	SuiteSparse_long *Gnz = (SuiteSparse_long*)calloc(m,sizeof(SuiteSparse_long));
	for(i=0;i<m;i++)
	{
	//	cout<<endl;
		SuiteSparse_long startindex=G.GrowIndex[i];
		SuiteSparse_long endindex=G.GrowIndex[i+1];
		Gnz[i] = endindex - startindex;

//		for(j=startindex;j<endindex;j++)
//			cout<<" "<<G.Gmat[j];
//		cout<<endl;
		
	}


/*	for(i=0;i<G.Pnonzero;i++)
		cout<<" "<<G.Pmat[i];
	cout<<endl;
	for(i=0;i<G.Pnonzero;i++)
		cout<<" "<<G.Pcolumns[i];
	cout<<endl;	
	for(i=0;i<m+1;i++)
		cout<<" "<<G.ProwIndex[i];
	cout<<endl;
/*	for(i=0;i<m;i++)
		printf(" %.8f", G.b1[i]);
	cout<<endl;
	for(i=0;i<m;i++)
		printf(" %.8f", G.x1[i]);
	cout<<endl;	
	
	
	for(i=0;i<m;i++)
	{
		cout<<endl;
		int startindex=G.ProwIndex[i];
		int endindex=G.ProwIndex[i+1];
		for(j=startindex;j<endindex;j++)
			cout<<" "<<G.Pmat[j];
		cout<<endl;
		
	}
	
/*	for(i=0;i<G.nonzero;i++)
		cout<<" "<<G.Mat_val[i];
	cout<<endl;
	for(i=0;i<G.nonzero;i++)
		cout<<" "<<G.columns[i];
	cout<<endl;	
	for(i=0;i<m+n+1;i++)
		cout<<" "<<G.rowIndex[i];
	cout<<endl;
	for(i=0;i<m+n;i++)
		printf(" %.8f", G.b[i]);
	cout<<endl;
	for(i=0;i<m+n;i++)
		printf(" %.8f", G.x[i]);
	cout<<endl;	
	
	
	for(i=0;i<m+n;i++)
	{
		cout<<endl;
		int startindex=G.rowIndex[i];
		int endindex=G.rowIndex[i+1];
		for(j=startindex;j<endindex;j++)
			cout<<" "<<G.Mat_val[j];
		cout<<endl;
		
	}
*/
/*	for (i=0;i<m+n+1;i++)
	{
		//cout<<endl;
		if(G.rowIndex[i]==G.rowIndex[i+1])
			break;
		
		for(j=G.rowIndex[i];j<G.rowIndex[i+1];j++)
		{
			if(G.Mat_val[j]>10)
				cout<<G.Mat_val[j]<<"\t";
		}
		//cout<<endl;
		/*for(j=G.rowIndex[i];j<G.rowIndex[i+1];j++)
		{
			cout<<G.columns[j]<<"\t";
		}
		//cout<<endl;
	}
	cout<<endl;
*/

//printing the matrix
	printf("\n Fine till here");
	printf("\n");
//	int* rowmIndex=(int*)calloc(m+1,sizeof(int));
	printf("\n Fine till here");
	printf("\n");
	//int rowmIndex[5]={1,2,3,4,5};
/*	for(i=0;i<m+1;i++)
	{
		rowmIndex[i]=G.rowIndex[i];
		printf(" %d", rowmIndex[i]);
	}
*/
	printf("\n Allocating GPU memory\n");
	cudaDeviceReset();
	size_t free, total;
	cudaMemGetInfo(&free, &total);
	printf("\n Free Mem = %lf MB, Total mem = %lf MB\n", (double)(free)/(1024*1024), (double)(total)/(1024*1024));


	double *dev_csrValA, *dev_b, *dev_x;
	int *dev_csrRowIdxA, *dev_csrColA;

	double *dev_GcsrVal, *dev_b1, *dev_x1;
	double *dev_PcsrVal, *dev_b2, *dev_x2;
	
	int *dev_GcsrRowIdx, *dev_PcsrRowIdx, *dev_GcsrCol, *dev_PcsrCol;
	

	
	cudaMalloc((void**)&dev_PcsrVal, G.Pnonzero*sizeof(double));
	cudaMalloc((void**)&dev_PcsrRowIdx, (m+1)*sizeof(int));
	cudaMalloc((void**)&dev_PcsrCol, G.Pnonzero*sizeof(int));


	cudaMalloc((void**)&dev_b1, (m)*sizeof(double));
	cudaMalloc((void**)&dev_b2, n*sizeof(double));
	cudaMalloc((void**)&dev_x1, m*sizeof(double));
	cudaMalloc((void**)&dev_x2, n*sizeof(double));

	cudaMemcpy(dev_b1, G.b1, (m)*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_x1, G.x1, (m)*sizeof(double), cudaMemcpyHostToDevice);

	cudaMemcpy(dev_PcsrVal, G.Pmat, G.Pnonzero*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_b2, G.b2, (n)*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_x2, G.x2, (n)*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_PcsrRowIdx, G.ProwIndex, (m+1)*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_PcsrCol, G.Pcolumns, (G.Pnonzero)*sizeof(int), cudaMemcpyHostToDevice);


	/* Matrix has been created and stored in CSR format.
	However, CHOLMOD requires CSC format. Since our matrix is symmetric positive definite, we can simply swap 
	csrColA with csrRowIdx and vice versa
	*/

	/* Starting the CHOLMOD routine now*/
	printf("\n Initiating CHOLMOD\n");
	cholmod_sparse *A, *P;
	cholmod_dense *x, *b, *r, *midvec;
	cholmod_factor *L;
	cholmod_common *Common, cm;
	Common = &cm;
	cholmod_l_start(Common);
//	&Common->useGPU=1;
	printf("\n m = %d, G.Gnonzero = %d\n", m, G.Gnonzero);	


	
	cholmod_sparse *C = cholmod_l_allocate_sparse((size_t)(m), (size_t)(m), (size_t)(G.Gnonzero), 1, 0, 1, 1, Common); 
//	P = cholmod_l_allocate_sparse((size_t)(m), (size_t)(n), (size_t)(G.Pnonzero), 1, 0, 0, 1, Common); 
//	printf("\n Allocated \n");

	C->itype = CHOLMOD_LONG;	
//	printf("\n Itype \n");
	C->p = &G.GrowIndex[0];
//	printf("\n Columns \n");
	C->nz = &Gnz[0];
//	printf("\n Rows \n");
	C->i = &G.Gcolumns[0];
	C->dtype = 0;
	C->x = &G.Gmat[0];
	 
/*	P->itype = CHOLMOD_LONG;
	P->p = &G.ProwIndex[0];
	P->nz = &Pnz[0];
	P->i = &G.Pcolumns[0];
	P->dtype = 0;
	P->x = &G.Pmat[0]; 
*/	 

	b = cholmod_l_allocate_dense((size_t)(m), 1, (size_t)(m), 1, Common);
	b->dtype=0;
	b->x = &G.b1[0];
	b->xtype = 1;

	printf("\n CHOLMOD manually set\n");
	cholmod_l_print_sparse(C, "A", Common);
	cholmod_l_print_dense(b, "b", Common);


	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	cudaEventRecord(start, 0);

	L = cholmod_l_analyze(C, Common);
	printf("\n Analysis: Flops: %g \t lnz: %g\n", Common->fl, Common->lnz);
	cholmod_l_factorize(C, L, Common);
	x = cholmod_l_solve(CHOLMOD_A, L, b, Common);
	
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	
	

	float elapsedTime;
	cudaEventElapsedTime(&elapsedTime, start, stop);
	printf("\n Time : %.6f secs :\n", elapsedTime);
	
	cholmod_l_print_dense(x, "X", Common);
	
	double *x1_mod = (double*)x->x;
	
	cudaMemcpy(dev_x1, x1_mod, m*sizeof(double), cudaMemcpyHostToDevice);
	
	cusparseStatus_t cuSparseStatus;
	cusparseHandle_t cuspHandle;
	cuSparseStatus = cusparseCreate(&cuspHandle);

	cusparseMatDescr_t descrP;
	cusparseCreateMatDescr(&descrP);

	cusparseSetMatType(descrP, CUSPARSE_MATRIX_TYPE_GENERAL);	
	cusparseSetMatIndexBase(descrP, CUSPARSE_INDEX_BASE_ZERO);
	
	
	double *dev_res1, *dev_simple;
	double *res1 = (double*)calloc(n,sizeof(double));
	cudaMalloc((void**)&dev_res1, n*sizeof(double));
	cudaMalloc((void**)&dev_simple, n*sizeof(double));
	
	const double alpha = 1.0, beta=0.0;
	//alpha = 1.0;
	//beta = 0.0;
	
	//solving P^T * G^-1 * b1 Result stored in dev_res1
	
	cuSparseStatus = cusparseDcsrmv(cuspHandle, CUSPARSE_OPERATION_TRANSPOSE, m, n, G.Pnonzero, &alpha, descrP, dev_PcsrVal, dev_PcsrRowIdx, dev_PcsrCol, dev_x1, &beta, dev_res1);
		
	if(cuSparseStatus == CUSPARSE_STATUS_SUCCESS)
	{
/*		cudaMemcpy(res1, dev_res1, n*sizeof(double), cudaMemcpyDeviceToHost);
		for(i=0;i<n;i++)
		{
			printf("\nres1[%d] = %.8f", i, res1[i]);
		}
		printf("\n P^T * G^-1 * b1 done! Vector stored in res1");
*/	}
	else
	{
		printf("\n P^T * G^-1 * b1 failed\n");
		exit(1);
	}
	
	const double alphaneg = -1.0;
		
		//Solving P^T * G^-1 * b1 - b2 ; Result stored in dev_res1
	
		
	cublasStatus_t cuBlasStatus;
	cublasHandle_t cubHandle;
	cuBlasStatus = cublasCreate(&cubHandle);
		
	cuBlasStatus = cublasDaxpy(cubHandle, n, &alphaneg, dev_b2, 1, dev_res1, 1);
	if(cuBlasStatus == CUBLAS_STATUS_SUCCESS)
	{
//		cudaMemcpy(res1, dev_res1, n*sizeof(double), cudaMemcpyDeviceToHost);
//		for(i=0;i<n;i++)
//		{
//			printf("\nres1[%d] = %.8f", i, res1[i]);
//		}
		printf("\n res1 = res1 - b2 done\n");
		 
	}
	else
	{
		printf("\n res1 = res1 - b2 failed\n");
	}
		
	
	
	
	///NOW COMPUTING G^-1 * P
	
	
	int k = 0;
	int breakloop=0;
	
	double **midMat = (double**)malloc(m*sizeof(double*));
	for(i=0;i<m;i++)
	{
		midMat[i] = (double*)calloc(n,sizeof(double));
	}
	
	cudaEventRecord(start, 0);
	
	for(i=0;i<n;i++)
	{
		breakloop = 0;
		double *vect = (double*)calloc(m,sizeof(double*));

		for(j=0;j<m;j++)
		{
			int startin = G.ProwIndex[j];
			int endin = G.ProwIndex[j+1];
			if(startin == endin)
				continue;
				
			k = startin;

			while(k<endin)
			{

				if(G.Pcolumns[k] == i)
				{	
					vect[j] = G.Pmat[k];
					breakloop=1;
					break;
				
				}
				k++;
			}
			if(breakloop == 1)
			{
				break;
			}
		}
		
		midvec = cholmod_l_allocate_dense((size_t)(m), 1, (size_t)(m), 1, Common);
		midvec->dtype=0;
		midvec->x=&vect[0];
		midvec->xtype = 1;
		
		cholmod_dense *res2;
		
		res2 = cholmod_l_solve(CHOLMOD_A, L, midvec, Common);

		
		double *re = (double*)res2->x;
		
//		printf("\n vector %d is:\n", i);
		int i1, j1, k1;
//		for(j1=0;j1<m;j1++)
//		{
//			midmat2flat[i+j1*n] = re[j1];
//			printf(" %lf", re[j1]);
//		}
//		printf("\n");
		for(i1=0;i1<m;i1++)
		{
			midMat[i1][i] = re[i1];
		}
		
		cholmod_l_free_dense(&midvec, Common);
			
	}
	
/*	printf("\n Midmat = \n");
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			printf(" %lf", midMat[i][j]);
		}
		printf("\n");
	} 
*/
	double *midMatflat = (double*)calloc((m*n),sizeof(double));
	double *dev_midMat;
	double *dev_solut;
	int counter = 0;
	for(i=0;i<n;i++)
	{
		for(j=0;j<m;j++)
		{
			midMatflat[counter] = midMat[j][i];
			counter++; 
		}
	}
	
	cudaMalloc((void**)&dev_midMat, m*n*sizeof(double));
	cudaMalloc((void**)&dev_solut, n*n*sizeof(double));
	
	cudaMemcpy(dev_midMat, midMatflat, m*n*sizeof(double), cudaMemcpyHostToDevice);
		
	//Solving P^T * midMat; Result stored in dev_solut
		
	cuSparseStatus = cusparseDcsrmm(cuspHandle, CUSPARSE_OPERATION_TRANSPOSE, m, n, n, G.Pnonzero, &alpha, descrP, dev_PcsrVal, dev_PcsrRowIdx, dev_PcsrCol, dev_midMat, m, &beta, dev_solut, n);
	
	if(cuSparseStatus == CUSPARSE_STATUS_SUCCESS)
	{
		printf("\n Solved P^T * G^-1 * P. Result stored in solut\n");
		
	}  
	else
	{
		printf("\n  Failed to Solve P^T * G^-1 * P \n");
		exit(1);
	}

/*	double *matGflat = (double*)calloc(n*n,sizeof(double));
	cudaMemcpy(matGflat, dev_solut, n*n*sizeof(double), cudaMemcpyDeviceToHost);
	counter = 0;
	printf("\nBefore LU starts\n");
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			printf(" %lf ", matGflat[counter]);
			counter++;
		}
		printf("\n");
	}
	printf("\n");
*/	
	cusolverStatus_t cuSolverStatus;
	
	
	cusolverDnHandle_t cudenHandle;
	cuSolverStatus = cusolverDnCreate(&cudenHandle);
	int Lwork = 0;
	cuSolverStatus = cusolverDnDgetrf_bufferSize(cudenHandle, n, n, dev_solut, n, &Lwork);
	
	if(cuSolverStatus == CUSOLVER_STATUS_SUCCESS)
	{
		printf("\n Buffer works\n Lwork = %d\n", Lwork);
	}
	else
	{
		exit(1);	
	}
	double *dev_Workspace;
	int *dev_Ipiv, *dev_Info;
	
	cudaMalloc((void**)&dev_Workspace, Lwork*sizeof(double));
	cudaMalloc((void**)&dev_Ipiv, n*sizeof(int));
	cudaMalloc((void**)&dev_Info, sizeof(int));
	
	//Calculating LU for dev_solut
//	double *nnmat = (double*)calloc(n*n,sizeof(double));
//	cudaMemcpy(nnmat, dev_solut, n*n*sizeof(double), cudaMemcpyDeviceToHost);
//	cuSolverStatus = cusolverDnDgetrfHost(cudenHandle, n, n, 

	cuSolverStatus = cusolverDnDgetrf(cudenHandle, n, n, dev_solut, n, dev_Workspace, dev_Ipiv, dev_Info);
	
	if(cuSolverStatus == CUSOLVER_STATUS_SUCCESS)
	{
		printf("\n solut has be defactorized into L and U. dev_Ipiv * solut = L * U\n");
	}
	else
	{
	
		printf("\n Unable to defactorize solut into LU\n");
		exit(1);	
	}
	
	//solving dev_solut * x = dev_res1. Result stored in dev_res1
	
	cuSolverStatus = cusolverDnDgetrs(cudenHandle, CUBLAS_OP_N, n, 1, dev_solut, n, dev_Ipiv, dev_res1, n, dev_Info);
	if(cuSolverStatus == CUSOLVER_STATUS_SUCCESS)
	{
		printf("\n Solution obtained for x2 \n");
	}
	else
	{
		printf("\n LU decomposition obtained by LU solver failed\n");
	}

/*	cudaMemcpy(G.x2, dev_res1, n*sizeof(double), cudaMemcpyDeviceToHost);
	printf("\n x2 = \n");
	for(i=0;i<n;i++)
	{
		printf("\n x2[%d] = %lf", i, G.x2[i]); 
	}
*/	
	
	double *dev_dummy;
	cudaMalloc((void**)&dev_dummy, m*sizeof(double));
	cudaMemset(dev_dummy, 0.0, m*sizeof(double));
	
	printf("\n Starting solving for x1 \n");
	//Solving for x1
	
		//Solving G^-1 * P * x2;  G^-1 * P is stored in midMat
	cuBlasStatus = cublasDgemv(cubHandle, CUBLAS_OP_N, m, n, &alpha, dev_midMat, m, dev_res1, 1, &beta, dev_dummy, 1);
	if(cuBlasStatus == CUBLAS_STATUS_SUCCESS)
	{
/*		double *toprint = (double*)calloc(m,sizeof(double));
		cudaMemcpy(toprint, dev_dummy, m*sizeof(double), cudaMemcpyDeviceToHost);
		printf("\n Intermediate vector :\n");
		for(i=0;i<m;i++)
		{
			printf("\ndummy[%d] = %lf", i, toprint[i]);
		}
*/		printf("\n midmat * x2 obtained. Stored in dummy\n");
	}
	else
	{
		printf("\n Failed to obtain midmat * x2\n");
	}
	
	cuBlasStatus = cublasDaxpy(cubHandle, m, &alphaneg, dev_dummy, 1, dev_x1, 1);
	if(cuBlasStatus == CUBLAS_STATUS_SUCCESS)
	{
		cudaMemcpy(G.x1, dev_x1, m*sizeof(double), cudaMemcpyDeviceToHost);
//		printf("\n x1 = \n");
//		for(i=0;i<m;i++)
//		{
//			printf("\n x1[%d] = %.15f", i, G.x1[i]);
//		}
		printf("\n x1 obtained");
	}
	else
	{
		printf("\n Failed to obtain x1");
	}

	printf("\n Solver finished its work\n");

/*		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);

		cudaEventElapsedTime(&elapsedTime, start, stop);
		printf("\n Time: %.6f msecs :\n", elapsedTime);
*/








	
	
	
	
	

	cholmod_l_finish(Common);
	return 0;

} 
