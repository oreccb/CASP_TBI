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

	vector<int> R;
	for(int i=0;i<100;i++){
		R.push_back(99);
	}
	//MyGraph<undirGraph> G(100,1,.5); //need to look at this function again-not getting a correct CC
	MyGraph<undirGraph> G("BrianInMap.txt",0);

	cout<<num_vertices(G.g)<<endl;
	cout<<num_edges(G.g)<<endl;
	cout<<G.computeClusteringCoefficient()<<endl;

	//G.printVE();

	//End doing stuff
	clock_t ends = clock();
	cout <<endl<< "Running Time : "<< (double) (ends - start) / CLOCKS_PER_SEC << endl;
	return 0;
}