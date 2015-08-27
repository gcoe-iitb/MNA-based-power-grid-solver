#include "../include/global.hpp"

#include "klu.h"



int main(int argc, char* argv[])
{
	const int bufsize = 512;
    	char buffer[bufsize];
	int m,n,S;
	double time_st,time_end,time_avg;
	if(argc!=2)
	{
		cout<<"Insufficient arguments"<<endl;
		return 1;
	}
	
	graph G;

	cerr<<"Start reading                    ";

	G.create_graph(argv[1]);

	cerr<<"Constructing Matrices            ";
	G.construct_MNA();
	G.construct_NA();
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
	int i,j;

	
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	cerr<<"Solving Equations    "<<endl;
	klu_symbolic *Symbolic;
	klu_numeric *Numeric;
	klu_common Common;
	klu_defaults(&Common);
	cudaEventRecord(start, 0);

	Symbolic = klu_analyze((m+n), G.rowIndex, G.columns, &Common);
	Numeric = klu_factor(G.rowIndex, G.columns, G.Mat_val, Symbolic, &Common);
	klu_solve(Symbolic, Numeric, (m+n), 1, G.b, &Common);
	
	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);	
	float elapsedTime;
	cudaEventElapsedTime(&elapsedTime, start, stop);

	printf("\n Time: %.6f milliseconds \n", elapsedTime);

	klu_free_symbolic(&Symbolic, &Common);
	klu_free_numeric(&Numeric, &Common);
	
		
} 
