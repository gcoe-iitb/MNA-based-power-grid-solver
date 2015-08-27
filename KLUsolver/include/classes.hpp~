class node
{
	public:
		int node_id,no_neighbor;
		vector<int> neighbor_id;
		vector<int> incident_edge_id;
		vector<int> voltage_neighbor_id;
		vector<int> current_neighbor_id;
		double node_potential,conductance_sum;
		node(): no_neighbor(0),node_potential(0),conductance_sum(0)
		{}
		~node(){}
};

class edge
{
	public:
		node* s_node;
		node* e_node;
		char type;
		//char type;
		double value;
		float branch_voltage;
		float branch_current;
	edge(): value(0),branch_voltage(0),branch_current(0)
	{}
	~edge(){}
};


class graph
{
	public:
		int no_nodes,no_edges;
		long nonzero;
//		double** M;
		//double* b;
		//double* x;
		//vector<float> Mat_values;
		//vector<int> Row_index,Col_index;
		int* rowIndex;
		int* columns;
		double* b;
		double* Mat_val;
		double* x;
		vector<node> node_array;
		vector<edge> edge_array;
		vector<int> voltage_edge_id;
		vector<int> current_edge_id;
		void create_graph(char*);
		void output_graph_stdout();
		int construct_MNA();
		int construct_MNA2();
//		int construct_sparse_MNA();
		int fillup_graph();
		void check_kcl();

	graph(): no_nodes(0),no_edges(0)
	{}
	~graph(){}
};
