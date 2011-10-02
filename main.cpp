#include "MyGraph.h"
#include <time.h>

int main()
{
	//graph types
	typedef adjacency_list<setS,vecS,undirectedS, MyNode, MyEdge> undirGraph;
	typedef adjacency_list<setS,vecS,directedS, MyNode, MyEdge> dirGraph;
	//typedef adjacency_list<vecS,vecS,undirectedS> U;

	//Time the Execution of the code
	clock_t start = clock();

	MyGraph<undirGraph> G_real("BrianInMap.txt",0);
	
	int n = num_vertices(G_real.g);
	int	e = num_edges(G_real.g);

	
	cout<<"Real Graph:"<<endl;
	cout<<num_vertices(G_real.g)<<endl;
	cout<<num_edges(G_real.g)<<endl;
	//cout<<G_real.computeClusteringCoefficient()<<endl<<endl;

	//G.printVE();

	//End doing stuff
	clock_t ends = clock();
	cout <<endl<< "Running Time : "<< (double) (ends - start) / CLOCKS_PER_SEC << endl;
	return 0;
}