#include "../include/global.hpp"

#include <cuda_runtime.h>
#include <cusparse.h>
#include <cublas_v2.h>

#include "/usr/local/cuda-6.5/samples/common/inc/helper_functions.h"  // helper for shared functions common to CUDA Samples
#include "/usr/local/cuda-6.5/samples/common/inc/helper_cuda.h"       // helper function CUDA error checking and initialization


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
	printf("\n Nonzero = %ld", G.nonzero);
	printf("\n Rows = %d", m+n);

	cout<<"MAT val:		       "<<endl;
	int i,j;

//	G.Mat_val[0] +=100;
/*	
	for(i=0;i<G.nonzero;i++)
		cout<<" "<<G.Mat_val[i];
	cout<<endl;
	for(i=0;i<G.nonzero;i++)
		cout<<" "<<G.columns[i];
	cout<<endl;	
	for(i=0;i<m+n+1;i++)
		cout<<" "<<G.rowIndex[i];
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
	cerr<<"Solving Equations    "<<endl;


    double r1, b, alpha, alpham1, beta, r0, a, na;
    
    
    const double tol = 0.1;
    const int max_iter = 1000000;
    int *d_col, *d_row;
    double *d_val, *d_x, dot;
    double *d_r, *d_p, *d_Ax;
    int k;
    
    cublasHandle_t cublasHandle = 0;
    cublasStatus_t cublasStatus;
    cublasStatus = cublasCreate(&cublasHandle);

    checkCudaErrors(cublasStatus);

    /* Get handle to the CUSPARSE context */
    cusparseHandle_t cusparseHandle = 0;
    cusparseStatus_t cusparseStatus;
    cusparseStatus = cusparseCreate(&cusparseHandle);

    checkCudaErrors(cusparseStatus);

    cusparseMatDescr_t descr = 0;
    cusparseStatus = cusparseCreateMatDescr(&descr);

    checkCudaErrors(cusparseStatus);

    cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO);

    checkCudaErrors(cudaMalloc((void **)&d_col, G.nonzero*sizeof(int)));
    checkCudaErrors(cudaMalloc((void **)&d_row, (m+n+1)*sizeof(int)));
    checkCudaErrors(cudaMalloc((void **)&d_val, G.nonzero*sizeof(double)));
    checkCudaErrors(cudaMalloc((void **)&d_x, (m+n)*sizeof(double)));
    checkCudaErrors(cudaMalloc((void **)&d_r, (m+n)*sizeof(double)));
    checkCudaErrors(cudaMalloc((void **)&d_p, (m+n)*sizeof(double)));
    checkCudaErrors(cudaMalloc((void **)&d_Ax, (m+n)*sizeof(double)));

    cudaMemcpy(d_col, G.columns, G.nonzero*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_row, G.rowIndex, (m+n+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_val, G.Mat_val, G.nonzero*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_x, G.x, (m+n)*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_r, G.b, (m+n)*sizeof(double), cudaMemcpyHostToDevice);

    alpha = 1.0;
    alpham1 = -1.0;
    beta = 0.0;
    r0 = 0.;

    printf("\n Data transferred\n");

    	cudaEvent_t start,stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	
	cudaEventRecord(start, 0);

    cusparseDcsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE, (m+n), (m+n), G.nonzero, &alpha, descr, d_val, d_row, d_col, d_x, &beta, d_Ax);

    cublasDaxpy(cublasHandle, (m+n), &alpham1, d_Ax, 1, d_r, 1);
    cublasStatus = cublasDdot(cublasHandle, (m+n), d_r, 1, d_r, 1, &r1);

    k = 1;

    while (r1 > tol && k <= max_iter)
    {
        if (k > 1)
        {
            b = r1 / r0;
            cublasStatus = cublasDscal(cublasHandle, (m+n), &b, d_p, 1);
            cublasStatus = cublasDaxpy(cublasHandle, (m+n), &alpha, d_r, 1, d_p, 1);
        }
        else
        {
            cublasStatus = cublasDcopy(cublasHandle, (m+n), d_r, 1, d_p, 1);
        }

        cusparseDcsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, (m+n), (m+n), G.nonzero, &alpha, descr, d_val, d_row, d_col, d_p, &beta, d_Ax);
        cublasStatus = cublasDdot(cublasHandle, (m+n), d_p, 1, d_Ax, 1, &dot);
        a = r1 / dot;

        cublasStatus = cublasDaxpy(cublasHandle, (m+n), &a, d_p, 1, d_x, 1);
        na = -a;
        cublasStatus = cublasDaxpy(cublasHandle, (m+n), &na, d_Ax, 1, d_r, 1);

        r0 = r1;
        cublasStatus = cublasDdot(cublasHandle, (m+n), d_r, 1, d_r, 1, &r1);
//        cudaThreadSynchronize();
//        printf("iteration = %3d, residual = %e\n", k, sqrt(r1));
        k++;
    }
    
    	cudaEventRecord(stop, 0);
	
	cudaEventSynchronize(stop);
	float elapsedTime;
	cudaEventElapsedTime(&elapsedTime, start, stop);
	printf("Iterations = %3d\tTime : %.6f milli-seconds : \n", k, elapsedTime);

    cudaMemcpy(G.x, d_x, (m+n)*sizeof(double), cudaMemcpyDeviceToHost);

/*        printf("\n x = \n");
	for(i=0;i<(m+n);i++)
	{
		printf("\n x[%d] = %.8f", i, G.x[i]);
	}	
*/    float rsum, diff, err = 0.0;

    for (int i = 0; i < (m+n); i++)
    {
        rsum = 0.0;

        for (int j = G.rowIndex[i]; j < G.rowIndex[i+1]; j++)
        {
            rsum += G.Mat_val[j]*G.x[G.columns[j]];
        }

        diff = fabs(rsum - G.b[i]);

        if (diff > err)
        {
            err = diff;
        }
    }

    cusparseDestroy(cusparseHandle);
    cublasDestroy(cublasHandle);

/*    free(I);
    free(J);
    free(val);
    free(x);
    free(rhs);
*/    cudaFree(d_col);
    cudaFree(d_row);
    cudaFree(d_val);
    cudaFree(d_x);
    cudaFree(d_r);
    cudaFree(d_p);
    cudaFree(d_Ax);

    // cudaDeviceReset causes the driver to clean up all state. While
    // not mandatory in normal operation, it is good practice.  It is also
    // needed to ensure correct operation when the application is being
    // profiled. Calling cudaDeviceReset causes all profile data to be
    // flushed before the application exits
    cudaDeviceReset();
    
/*    printf("\n X is:\n");
    for(i=0;i<(m+n);i++)
    {
    	printf("\n X[%d] = %.4f", i, G.x[i]);
    }
*/    	
    printf("Test Summary:  Error amount = %f\n", err);
    exit((k <= max_iter) ? 0 : 1);




/*	
	culaSparseHandle handle;
    	if (culaSparseCreate(&handle) != culaSparseNoError)
    	{
        	// this should only fail under extreme conditions
        	std::cout << "fatal error: failed to create library handle!" << std::endl;
        	exit(EXIT_FAILURE);
    	}
	
	StatusChecker sc(handle);

	culaSparsePlan plan;
	sc = culaSparseCreatePlan(handle, &plan);
	
	culaSparseCsrOptions formatOpts;
	culaSparseCsrOptionsInit(handle, &formatOpts);
	//formatOpts.indexing=1;	
	sc = culaSparseSetDcsrData(handle, plan, &formatOpts, m+n, G.nonzero, &G.Mat_val[0], &G.rowIndex[0], &G.columns[0], &G.x[0], &G.b[0]);
	printf("\n Fine till here");
	printf("\n");	
	culaSparseConfig config;
    	sc = culaSparseConfigInit(handle, &config);
    	config.relativeTolerance = 1e-2;
    	config.maxIterations = 100000;
	config.divergenceTolerance = 20;
	culaSparseResult result;
    	
/*    	sc = culaSparseSetHostPlatform(handle, plan, 0);

    // set cg solver
    	sc = culaSparseSetCgSolver(handle, plan, 0);
	printf("\n Fine till here");
	printf("\n");    
    // perform solve (cg + no preconditioner on host)
	culaSparseResult result;
	sc = culaSparseExecutePlan(handle, plan, &config, &result);
	sc = culaSparseGetResultString(handle, &result, buffer, bufsize);
	std::cout << buffer << std::endl;
    	
    	if (culaSparsePreinitializeCuda(handle) == culaSparseNoError)
    {
        // change to cuda accelerated platform
        sc = culaSparseSetCudaPlatform(handle, plan, 0);

        // perform solve (cg + ilu0 on cuda)
        culaSparseGmresOptions solverOpts;
	culaSparseStatus status = culaSparseGmresOptionsInit(handle, &solverOpts);
	solverOpts.restart = 20;


	sc = culaSparseSetCgSolver(handle, plan, 0);

//        sc = culaSparseSetJacobiPreconditioner(handle, plan, 0);
//        sc = culaSparseSetIlu0Preconditioner(handle, plan, 0);
        sc = culaSparseExecutePlan(handle, plan, &config, &result);
        sc = culaSparseGetResultString(handle, &result, buffer, bufsize);
        std::cout << buffer << std::endl;

        // change preconditioner to fainv 
        // this avoids data transfer by using data cached by the plan
//        sc = culaSparseSetIlu0Preconditioner(handle, plan, 0);

        // perform solve (cg + fainv on cuda)
        // these timing results should indicate minimal overhead
/*        sc = culaSparseExecutePlan(handle, plan, &config, &result);
        sc = culaSparseGetResultString(handle, &result, buffer, bufsize);
        std::cout << buffer << std::endl;

        // change solver
        // this avoids preconditioner regeneration by using data cached by the plan
        sc = culaSparseSetBicgstabSolver(handle, plan, 0);

        // perform solve (bicgstab + fainv on cuda)
        // the timing results should indicate minimal overhead and preconditioner generation time
        sc = culaSparseExecutePlan(handle, plan, &config, &result);
        sc = culaSparseGetResultString(handle, &result, buffer, bufsize);
        std::cout << buffer << std::endl;
        
        sc = culaSparseSetGmresSolver(handle, plan, 0);

        // perform solve (bicgstab + fainv on cuda)
        // the timing results should indicate minimal overhead and preconditioner generation time
        sc = culaSparseExecutePlan(handle, plan, &config, &result);
        sc = culaSparseGetResultString(handle, &result, buffer, bufsize);
        std::cout << buffer << std::endl;
    }
    else
    {
        std::cout << "alert: no cuda capable gpu found" << std::endl;
    }

    // cleanup plan
    culaSparseDestroyPlan(plan);

    // cleanup handle
    culaSparseDestroy(handle);

    	FILE* myWriteFile;
  	myWriteFile=fopen("result.txt","w");
  	for (i = 0; i < n; i++)
    	{
		fprintf(myWriteFile,"%1f\n",G.x[i]);
    		//  printf ("\n x [%d] = % f", i, x[i]);
    	}
  	fprintf(myWriteFile,".end\n");
  	fclose(myWriteFile);
 
  	printf ("\n");
    	
    	
//	time_st=dsecnd();
//	solver(G.rowIndex,G.columns,G.Mat_val,G.b,G.x,m+n,G.nonzero);
//	time_end=dsecnd();
//	time_avg = (time_end-time_st);
//	printf("Successfully Solved in : %.6f secs\n",time_avg);
	cerr<<endl;

	cerr<<"Fillup Graph                    ";	
//	time_st=dsecnd();
	G.fillup_graph();
//	time_end=dsecnd();
//	time_avg = (time_end-time_st);
//	cerr<<"Done                 "<<time_avg<<endl;

	//G.output_graph_stdout();
	cerr<<"Matching KCL                    ";
//	time_st=dsecnd();
	G.check_kcl();
//	time_end=dsecnd();
//	time_avg = (time_end-time_st);
//	cerr<<"Done                 "<<time_avg<<endl;
	/*for (int i=0;i<m+n;i++)
	  {
	  cout<<"M"<<i<<endl;
	  for (int j=0;j<m+n;j++)
	  cout<<" "<<j<<"#"<<M[i][j]<<endl;
	  }*/
} 
