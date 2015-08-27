#include "./include/global.hpp"

void graph:: create_graph(char* filename)
{
	FILE * myReadFile;
	char type;
	int node1,node2,node3,i,k=0,j=0,l=0,max=-1,tot_edge=0,tem;
	double value;
	node node_temp;
	edge edge_temp;
	myReadFile = fopen (filename , "r");

	
	if(myReadFile!=NULL) 
	{
		while(1)
		{
			fscanf(myReadFile, "%c%*d %d %d %lf\n",&type,&node1,&node2,&value);
			//cout<<"type"<<type<<"node1"<<node1<<"node2"<<node2<<"value"<<value<<endl;
			if (type=='.')
				break;
		//	value=(condition)?1:0;
			node3=(node1>node2)?node1:node2;
		//	if (node2>node1)
		//		node3=node2;
			if (node3>max)
			{
				max=node3;
			}

			tot_edge++;
		}
	}
	fclose(myReadFile);
	for(k=0;k<=max;k++)
	{
		node_temp.node_id=k;
		node_array.push_back(node_temp);
	}

	myReadFile = fopen (filename, "r");
	if(myReadFile!=NULL)
	{
		while(1)
		{
			fscanf(myReadFile, "%c%*d %d %d %lf\n",&type,&node1,&node2,&value);
			if (type=='.')
				break;
		/*	node3=(node1>node2)?node1:node2;
			if(node3>max)
			{
				max=node3;
			}
			for(;k<=max;k++)
			{
				node_temp.node_id=k;
				node_array.push_back(node_temp);
			}*/
					

		/*	for (i=0;i<=max;i++)
				cout<<"##node"<<&node_array[i]<<"id"<<node_array[i].node_id<<endl;
			if (j>0)
	cout<<"#0edge0snode"<<(edge_array[0].s_node)->node_id<<endl;*/
			edge_temp.s_node=&(node_array[node1]);
			edge_temp.e_node=&(node_array[node2]);
		//	cout<<"snode"<<(edge_temp.s_node)->node_id<<"enode"<<(edge_temp.e_node)->node_id<<"snode"<<node_array[node1].node_id<<"enode"<<node_array[node2].node_id<<endl;
			edge_temp.type=type;
			edge_temp.value=value;
			//if(j>0)
//	cout<<"#1edge0snode"<<edge_array[0].s_node<<"id"<<(edge_array[0].s_node)->node_id<<endl;
			edge_array.push_back(edge_temp);
	//cout<<"#2edge0snode"<<edge_array[0].s_node<<"id"<<(edge_array[0].s_node)->node_id<<endl;
			if (type=='R')
			{
			/*	if(value<10^-6)
				{
					(node_array[node1].voltage_neighbor_id).push_back(l);
                               		(node_array[node2].voltage_neighbor_id).push_back(l);
                                	voltage_edge_id.push_back(j);
                                l++;
				}
				else*/ 
					(node_array[node1].neighbor_id).push_back(node2);
					(node_array[node1].incident_edge_id).push_back(j);	
				
				
					
					(node_array[node2].neighbor_id).push_back(node1);
					(node_array[node2].incident_edge_id).push_back(j);
				
				node_array[node1].conductance_sum=node_array[node1].conductance_sum+1.0/value;
				node_array[node2].conductance_sum=node_array[node2].conductance_sum+1.0/value;

//			cout<<"node#"<<node2<<" "<<node_array[node2].conductance_sum<<"value"<<value<<endl;
		
			}			
			else if(type=='V')
			{
				if((node_array[node1].voltage_neighbor_id).size() < 1)
				{
					(node_array[node1].voltage_neighbor_id).push_back(l);
					(node_array[node2].voltage_neighbor_id).push_back(l);
					voltage_edge_id.push_back(j);
					l++;
				}	
			}
			else
			{
				(node_array[node1].current_neighbor_id).push_back(j);
				(node_array[node2].current_neighbor_id).push_back(j);
				current_edge_id.push_back(j);
			}
			//cout<<"chk2 nbr"<<k;
			j++;
			tem=tot_edge/100;
			if(tem!=0&&j%tem==0)
				cerr<<" "<<j/tem<<"% ";

		}
	}

	fclose(myReadFile);
	cout<<endl;
	cout<<endl;
	cout<<"GRAPH STAT:"<<endl;
	cout<<"Total node:           "<<max<<endl;
	cout<<"Total edge:           "<<edge_array.size()<<endl;
	cout<<"Voltage source:       "<<l<<endl;
	cout<<"Current source:       "<<current_edge_id.size()<<endl;
	cout<<endl;
}

int graph:: construct_MNA()
{
	
	long m,n;	
	long i,j,k,ij,datatemp=0;
	nonzero=0;
	double value;
	m=node_array.size()-1;
	n=voltage_edge_id.size();
	long N=m;
	b=(double*)calloc(m+n,sizeof(double));
	x=(double*)calloc(m+n,sizeof(double));
	rowIndex=(int*)calloc(m+n+1,sizeof(int));
	
	printf("\n Initialization done");
	printf("\n");
	printf("\n");
	rowIndex[0] = 0;
/*	
	for (i=1;i<=m;i++)
	{
		if((node_array[i].neighbor_id).size()>4)
		{
			printf("\n Node %d has %d neighbors\n", i, (node_array[i].neighbor_id).size());
		}
	}
//	sleep(10);
*/	for (i=1;i<=m;i++)
	{
		rowIndex[i]=rowIndex[i-1]+(node_array[i].neighbor_id).size()+1+(node_array[i].voltage_neighbor_id).size();
	}
	nonzero=m*8;
//	rowIndex[m]-=1;
	for (i=m+1;i<=m+n;i++)
		rowIndex[i]=0;
	printf("\n nonzero = %ld", nonzero);	
//	printf("\n nonzero = %d", nonzero);
	printf("\n");	
	columns=(int*)calloc(nonzero,sizeof(int));
	Mat_val=(double*)calloc(nonzero,sizeof(double));
	
	printf("\n Entering 1st loop");
	printf("\n");	

//#pragma omp iparallel for private(i,j,ij)
	for (i = 1; i <= m; i++)
	{
		ij=rowIndex[i-1];
		Mat_val[ij]=node_array[i].conductance_sum;
		columns[ij]=i-1;
		ij++;
		for (j=0;j<(node_array[i].neighbor_id).size();j++)
		{
			Mat_val[ij]=-1/(edge_array[node_array[i].incident_edge_id[j]].value);
			columns[ij]=node_array[i].neighbor_id[j]-1;
			ij++;
		}
		for (j=0;j<(node_array[i].voltage_neighbor_id).size();j++)
		{
		//	printf("\n (node_array[%d].voltage_neighbor_id).size() = %d", i, (node_array[i].voltage_neighbor_id).size()); 
			int pos=node_array[i].voltage_neighbor_id[j];
			if((edge_array[voltage_edge_id[pos]].s_node)->node_id==i)
				Mat_val[ij]=1;
			else
				Mat_val[ij]=-1;
			columns[ij]=m+pos;
			rowIndex[columns[ij]+1]+=1;
			ij++;
		}
		for (j=0;j<(node_array[i].current_neighbor_id).size();j++)
		{
			if (edge_array[node_array[i].current_neighbor_id[j]].s_node->node_id==i)
				value=-edge_array[node_array[i].current_neighbor_id[j]].value;
			else
				value=edge_array[node_array[i].current_neighbor_id[j]].value;
			b[i-1]=b[i-1]+value;
		}

	}

	printf("\n Exiting first loop");
	printf("\n");
	
	for (i=m+1;i<=m+n;i++)
		rowIndex[i]+=rowIndex[i-1];
	
	nonzero=rowIndex[m+n];
	printf("\n Nonzero established = %ld\n", nonzero);
#pragma omp iparallel for private(i,j,ij)
	for (i = 1; i <= m; i++)
	{
		ij=rowIndex[i-1];
		for (j=0;j<(node_array[i].neighbor_id).size();j++)
		{
			ij++;
		}			
		for (j=0;j<(node_array[i].voltage_neighbor_id).size();j++)
		{
			int pos=node_array[i].voltage_neighbor_id[j];
			int temps=m+pos;
			columns[rowIndex[temps]]=i-1;
			Mat_val[rowIndex[temps]]=Mat_val[ij+1];
			ij++;
		}
	}

#pragma omp parallel for private(i)
	for(i=0;i<n;i++)
	{
		b[i+m]=edge_array[voltage_edge_id[i]].value;
//		x[i+m] = edge_array[current_edge_id[i]].value;
	}

#pragma omp iparallel for private(i,j)
	for(i=0;i<m;i++)
	{
		int startindex=rowIndex[i];
		int endindex=rowIndex[i+1]-1;
		int counter = 0;
		while(counter < endindex-startindex)
		{
			for(j=startindex;j<endindex;j++)
			{
				if(columns[j]>columns[j+1])
				{
					double tempdata=Mat_val[j];
					int tempcol=columns[j];
					columns[j]=columns[j+1];
					columns[j+1]=tempcol;
					Mat_val[j]=Mat_val[j+1];
					Mat_val[j+1]=tempdata;
				}
			}
			counter++;
		}
	}
//rowIndex[m]-=1;
//
	printf("\n Function done\n");
return(0);
}

int graph::construct_NA()
{
//	printf("\n NA function entered\n");
	int m,n;
	m = node_array.size()-1;
	n = voltage_edge_id.size();
	int i, j;
	b1 = (double*)calloc(m,sizeof(double));
	b2 = (double*)calloc(n,sizeof(double));
	GrowIndex = (SuiteSparse_long*)calloc((m+1),sizeof(SuiteSparse_long));
	ProwIndex = (int*)calloc((m+1),sizeof(int));
	Gnonzero = 0;
	Pnonzero = 0;

	GrowIndex[0] = 0;
	ProwIndex[0] = 0;
	for(i=0;i<m;i++)
	{
		int start = rowIndex[i];
		int end = rowIndex[i+1];	
		for(j=start;j<end;j++)
		{
			if(Mat_val[j]!=1)
			{
				Gnonzero++;
			}
			else
			{
				Pnonzero++;
			}
		}
	}
	//printf("\n First loop done\n");
	Gcolumns = (SuiteSparse_long*)calloc(Gnonzero,sizeof(SuiteSparse_long));
	Gmat = (double*)calloc(Gnonzero,sizeof(double));

	Pcolumns = (int*)calloc(Pnonzero, sizeof(int));
	Pmat = (double*)calloc(Pnonzero, sizeof(double));

	int Gcounter = 0;
	int Pcounter = 0;
	for(i=0;i<m;i++)
	{
		int Ginternalcounter = 0;
		int Pinternalcounter = 0;
		int start = rowIndex[i];
		int end = rowIndex[i+1];
	//	printf("\n i = %d\t start = %d\t end = %d\n", i,start,end);
		for(j=start;j<end;j++)
		{
			
			if(Mat_val[j] != 1)
			{
				Gcolumns[Gcounter] = columns[j];
				Gmat[Gcounter] = Mat_val[j];
				Ginternalcounter++;
				Gcounter++;
			}
			else
			{
				Pcolumns[Pcounter] = columns[j] - m;
				Pmat[Pcounter] = 1;
				Pinternalcounter++;
				Pcounter++;
			}
		}
		
		GrowIndex[i+1] = GrowIndex[i] + Ginternalcounter;
		ProwIndex[i+1] = ProwIndex[i] + Pinternalcounter;
	} 

	for(i=0;i<m;i++)
	{
		b1[i] = b[i];
	}
	for(i=0;i<n;i++)
	{	
		b2[i] = b[i+m];
	}
	
	x1 = (double*)calloc(m,sizeof(double));
	x2 = (double*)calloc(n,sizeof(double));
	
/*	
	cout<<endl;
	for(i=0;i<m;i++)
		printf(" %.8f", b1[i]);
	cout<<endl;
	for(i=0;i<n;i++)
		printf("%.8f", b2[i]);
	cout<<endl;
	for(i=0;i<m;i++)
		printf(" %.8f", x1[i]);
	cout<<endl;	
*/

}


int graph::construct_MNA3()
{	
	long m,n;	
	long i,j,k,ij,datatemp=0;
	nonzero=0;
	double value;
	m=node_array.size()-1;
	n=voltage_edge_id.size();
	long N=m;
	b=(double*)calloc(m,sizeof(double));
	x=(double*)calloc(m,sizeof(double));
	rowIndex=(int*)calloc(m+1,sizeof(int));
	
	printf("\n Initialization done");
	printf("\n");
	printf("\n");
	rowIndex[0] = 0;
/*	
	for (i=1;i<=m;i++)
	{
		if((node_array[i].neighbor_id).size()>4)
		{
			printf("\n Node %d has %d neighbors\n", i, (node_array[i].neighbor_id).size());
		}
	}
//	sleep(10);

*/	for (i=1;i<=m;i++)
	{
		rowIndex[i]=rowIndex[i-1]+(node_array[i].neighbor_id).size()+1;//+(node_array[i].voltage_neighbor_id).size();
	}
	nonzero=m*8;
//	rowIndex[m]-=1;
//	for (i=m+1;i<=m+n;i++)
//		rowIndex[i]=0;
	printf("\n nonzero = %ld", nonzero);	
//	printf("\n nonzero = %d", nonzero);
	printf("\n");	
	columns=(int*)calloc(nonzero,sizeof(int));
	Mat_val=(double*)calloc(nonzero,sizeof(double));
	
	printf("\n Entering 1st loop");
	printf("\n");	

//#pragma omp iparallel for private(i,j,ij)
	for (i = 1; i <= m; i++)
	{
		ij=rowIndex[i-1];
		Mat_val[ij]=node_array[i].conductance_sum;
		columns[ij]=i-1;
		ij++;
		for (j=0;j<(node_array[i].neighbor_id).size();j++)
		{
			Mat_val[ij]=-1/(edge_array[node_array[i].incident_edge_id[j]].value);
			columns[ij]=node_array[i].neighbor_id[j]-1;
			ij++;
		}
/*		for (j=0;j<(node_array[i].voltage_neighbor_id).size();j++)
		{
		//	printf("\n (node_array[%d].voltage_neighbor_id).size() = %d", i, (node_array[i].voltage_neighbor_id).size()); 
			int pos=node_array[i].voltage_neighbor_id[j];
			if((edge_array[voltage_edge_id[pos]].s_node)->node_id==i)
				Mat_val[ij]=1;
			else
				Mat_val[ij]=-1;
			columns[ij]=m+pos;
			rowIndex[columns[ij]+1]+=1;
			ij++;
		}
*/		for (j=0;j<(node_array[i].current_neighbor_id).size();j++)
		{
			if (edge_array[node_array[i].current_neighbor_id[j]].s_node->node_id==i)
				value=-edge_array[node_array[i].current_neighbor_id[j]].value;
			else
				value=edge_array[node_array[i].current_neighbor_id[j]].value;
			b[i-1]=b[i-1]+value;
		}

	}

	printf("\n Exiting first loop");
	printf("\n");
	
//	for (i=m+1;i<=m+n;i++)
//		rowIndex[i]+=rowIndex[i-1];
	
	nonzero=rowIndex[m];
	printf("\n Nonzero established = %ld\n", nonzero);
/*#pragma omp iparallel for private(i,j,ij)
	for (i = 1; i <= m; i++)
	{
		ij=rowIndex[i-1];
		for (j=0;j<(node_array[i].neighbor_id).size();j++)
		{
			ij++;
		}			
		for (j=0;j<(node_array[i].voltage_neighbor_id).size();j++)
		{
			int pos=node_array[i].voltage_neighbor_id[j];
			int temps=m+pos;
			columns[rowIndex[temps]]=i-1;
			Mat_val[rowIndex[temps]]=Mat_val[ij+1];
			ij++;
		}
	}

#pragma omp parallel for private(i)
	for(i=0;i<n;i++)
	{
		b[i+m]=edge_array[voltage_edge_id[i]].value;
//		x[i+m] = edge_array[current_edge_id[i]].value;
	}
*/
#pragma omp iparallel for private(i,j)
	for(i=0;i<m;i++)
	{
		int startindex=rowIndex[i];
		int endindex=rowIndex[i+1]-1;
		int counter = 0;
		while(counter < endindex-startindex)
		{
			for(j=startindex;j<endindex;j++)
			{
				if(columns[j]>columns[j+1])
				{
					double tempdata=Mat_val[j];
					int tempcol=columns[j];
					columns[j]=columns[j+1];
					columns[j+1]=tempcol;
					Mat_val[j]=Mat_val[j+1];
					Mat_val[j+1]=tempdata;
				}
			}
			counter++;
		}
	}
//rowIndex[m]-=1;
return(0);
}


int graph:: construct_MNA2()
{
	
	int m,n;	
	int i,j,k,ij,datatemp=0;
	nonzero=0;
	double value;
	m=node_array.size()-1;
	n=voltage_edge_id.size();
	int N=m+n;
	b=(double*)calloc(m+n,sizeof(double));
	x=(double*)calloc(m+n,sizeof(double));
	rowIndex=(int*)calloc(m+n+1,sizeof(int));
/*	for(i=0;i<m+n;i++)
		printf(" %.8f ", b[i]);
*/	printf("\n");
	//solValues=(_DOUBLE_PRECISION_t*)calloc(m+n,sizeof(_DOUBLE_PRECISION_t));
/*	for(i = 0; i < m+n; i++) 
	{    
		M[i] =(double*)calloc(m+n,sizeof(double));  	
	}	*/
	rowIndex[0] = 0;
	
	for (i=1;i<=m;i++)
	{
		rowIndex[i]=rowIndex[i-1]+(node_array[i].neighbor_id).size()+1+(node_array[i].voltage_neighbor_id).size();
	}
	
	for (i=m+1;i<=m+n;i++)
		rowIndex[i]=rowIndex[m];
	
	printf("\n");	
	nonzero=rowIndex[m]-1;
	columns=(int*)calloc(nonzero,sizeof(int));
	Mat_val=(double*)calloc(nonzero,sizeof(double));




#pragma omp iparallel for private(i,j,ij)
	for (i = 1; i <= m; i++)
	{
		ij=rowIndex[i-1];
//		printf("%d ",ij);
		Mat_val[ij]=node_array[i].conductance_sum;
		columns[ij]=i-1;
//		printf("%d ", columns[ij]);
		ij++;
//		printf("\n");

		
		for (j=0;j<(node_array[i].neighbor_id).size();j++)
		{
			//cout<<"j"<<j<<" "<<edge_array[node_array[i+1].incident_edge_id[j]].value<<endl;
		//	M=M+1/(edge_array[node_array[i].incident_edge_id[j]].value);

			Mat_val[ij]=-1/(edge_array[node_array[i].incident_edge_id[j]].value);
			columns[ij]=node_array[i].neighbor_id[j]-1;
//			printf("%d ",ij);
//			printf("%d ", columns[ij]);
			ij++;
			
//			printf("\n");
//			k=node_array[i+1].neighbor_id[j]-1;
/*			if(k>=0)
			{
				M[i][k]=M[i][k]-1/(edge_array[node_array[i+1].incident_edge_id[j]].value);
				datatemp++;
			}*/
		}
		//cout<<"Inner loop 1 done"<<endl;
		
		for (j=0;j<(node_array[i].voltage_neighbor_id).size();j++)
		{
			int pos=node_array[i].voltage_neighbor_id[j];
			if((edge_array[voltage_edge_id[pos]].s_node)->node_id==i)
				Mat_val[ij]=1;
			else
				Mat_val[ij]=-1;
			columns[ij]=m+pos;
			ij++;
		}
		for (j=0;j<(node_array[i].current_neighbor_id).size();j++)
		{
			if (edge_array[node_array[i].current_neighbor_id[j]].s_node->node_id==i)
				value=-edge_array[node_array[i].current_neighbor_id[j]].value;
			else
				value=edge_array[node_array[i].current_neighbor_id[j]].value;
			//int temp = i+1;
			b[i-1]=b[i-1]+value;
		}

	}
/*	for(i=0;i<m+n;i++)
		printf(" %.8f ", b[i]);
*/

#pragma omp parallel for private(i)
	for(i=0;i<n;i++)
		b[i+m]=edge_array[voltage_edge_id[i]].value;
			


return(0);
}







/*int graph :: construct_sparse_MNA()
{
	int i,j,m,n,row=1;
	m=node_array.size()-1;
	n=voltage_edge_id.size();
	for (i=0;i<m+n;i++)
	{
		rowIndex[i]=row;
		for (j=i;j<m+n;j++)
		{
			if(M[i][j]!=0)
			{
				Mat_val[row-1]=M[i][j];
				columns[row-1]=(j+1);
				row++;
			}
		}
	}
	rowIndex[m+n]=row;
	cout<<"nonzero"<<nonzero;
	for (i=0;i<nonzero;i++)
		cout<<"i"<<i<<"M"<<Mat_val[i]<<"col"<<columns[i]<<endl;
	cout<<"row"<<m+n<<"last_row_value"<<row<<"nonzero"<<nonzero<<endl;

	for (i=0;i<=m+n;i++)
		cout<<"i"<<i<<"row"<<rowIndex[i]<<endl;
	for (i=0;i<m+n;i++)
		cout<<"i"<<i<<"b"<<b[i]<<endl;
	ofstream arrayData("value.txt"); // File Creation
	ofstream rowindex("ia.txt");
	ofstream colindex("ja.txt");

	for(int k=0;k<nonzero;k++)
   	{
        	arrayData<<Mat_val[k]<<endl; //Outputs array to txtFile
		colindex<<columns[k]<<endl;
   	 }
	for(int k=0;k<m+n;k++)
	{
		rowindex<<rowIndex[k]<<endl;
	}
   
return(0);
}*/

int graph :: fillup_graph()
{
	int i;
	no_nodes=node_array.size();
	no_edges=edge_array.size();
#pragma omp parallel for private(i)
	for(i=1;i<no_nodes;i++)
	{
		node_array[i].node_potential=x[i-1];
	}
	node_array[0].node_potential=0;
#pragma omp parallel for private(i)
	for(i=0;i<no_edges;i++)
	{
		edge_array[i].branch_voltage=(edge_array[i].s_node)->node_potential-(edge_array[i].e_node)->node_potential;
		if(edge_array[i].type=='R')
			edge_array[i].branch_current=edge_array[i].branch_voltage/edge_array[i].value;
	}
#pragma omp parallel for private(i)
	for(i=0;i<voltage_edge_id.size();i++)
		edge_array[voltage_edge_id[i]].branch_current=x[i+no_nodes-1];
#pragma omp parallel for private(i)
	for(i=0;i<current_edge_id.size();i++)
		edge_array[current_edge_id[i]].branch_current=edge_array[current_edge_id[i]].value;
return(0);
}

void graph :: output_graph_stdout()
{
	int i,j;
	for(i=0;i<edge_array.size();i++)
	{
		cout<<""<<edge_array[i].type<<"     "<<edge_array[i].branch_voltage<<"     "<<edge_array[i].branch_current<<endl;
	}
/*	for(i=0;i<edge_array.size();i++)
	{
		cout<<"Edge#"<<i<<endl<<"snode "<<(edge_array[i].s_node)->node_id<<endl<<"enode "<<(edge_array[i].e_node)->node_id<<endl<<"type "<<edge_array[i].type<<endl<<"value "<<edge_array[i].value<<"voltage"<<edge_array[i].branch_voltage<<"current"<<edge_array[i].branch_current<<endl;
	}

	for(i=0;i<node_array.size();i++)
	{
		cout<<"NodeId"<<node_array[i].node_id<<"potential"<<node_array[i].node_potential<<endl;
	
	}*/
}

void graph :: check_kcl()
{
	int i,j,temp1,temp2,pos,temp;
	float *tot_current;
	tot_current=(float*)calloc(node_array.size(),sizeof(float));
#pragma omp parallel for private(i,j,temp1,temp2,pos) shared(tot_current)
	for(i=0;i<node_array.size();i++)
	{
		for(j=0;j<(node_array[i].incident_edge_id).size();j++)
		{

			temp1=(edge_array[node_array[i].incident_edge_id[j]].s_node)->node_id;
//			temp2=(edge_array[node_array[i].incident_edge_id[j]].e_node)->node_id;
			if(temp1==i)
				tot_current[i]=tot_current[i]+edge_array[node_array[i].incident_edge_id[j]].branch_current;
			else
				tot_current[i]=tot_current[i]-edge_array[node_array[i].incident_edge_id[j]].branch_current;
	//	cout<<" i "<<i<<" r "<<endl;
		}
		for(j=0;j<(node_array[i].voltage_neighbor_id).size();j++)
		{
			pos=node_array[i].voltage_neighbor_id[j];

			temp1=(edge_array[voltage_edge_id[pos]].s_node)->node_id;
	//		temp2=(edge_array[voltage_edge_id[pos]].e_node)->node_id;
			if(temp1==i)
				tot_current[i]=tot_current[i]+edge_array[voltage_edge_id[pos]].branch_current;
			else
				tot_current[i]=tot_current[i]-edge_array[voltage_edge_id[pos]].branch_current;
	//		cout<<" i "<<i<<" v "<<endl;

		}
		for(j=0;j<(node_array[i].current_neighbor_id).size();j++)
		{
			temp1=(edge_array[node_array[i].current_neighbor_id[j]].s_node)->node_id;
	//		temp2=(edge_array[node_array[i].current_neighbor_id[j]].e_node)->node_id;
			if(temp1==i)
				tot_current[i]=tot_current[i]+edge_array[node_array[i].current_neighbor_id[j]].branch_current;
			else
				tot_current[i]=tot_current[i]-edge_array[node_array[i].current_neighbor_id[j]].branch_current;
	//		cout<<" i "<<i<<" I "<<endl;

		}


	}
	for(i=0;i<node_array.size();i++)
	{
		for(j=0;j<(node_array[i].incident_edge_id).size();j++)
		{
			temp=(edge_array[node_array[i].incident_edge_id[j]].s_node)->node_id;
			if(temp==i)
			{
				temp=(edge_array[node_array[i].incident_edge_id[j]].e_node)->node_id;
				tot_current[temp]=tot_current[temp]-edge_array[node_array[i].incident_edge_id[j]].branch_current;
			}
			else
				tot_current[temp]=tot_current[temp]+edge_array[node_array[i].incident_edge_id[j]].branch_current;
		}
	}

	double max=0,min=0;

	for(i=1;i<node_array.size();i++)
	{
		//cout<<"total current"<<tot_current[i]<<endl;
		if(tot_current[i]>max)
			max=tot_current[i];
		if(tot_current[i]<min)
			min=tot_current[i];
	/*	if(tot_current[i]>.0001||tot_current[i]<-.0001)
		{
			cout<<"BAD CIRCUIT "<<tot_current[i]<<" FAULT IN NODE "<<i<<endl;
			break;
		}*/
	}
	cout<<"Error range  "<<min<<" < I < "<<max;
}
