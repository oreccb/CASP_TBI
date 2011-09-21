#include "MyGraph.h"
#include <time.h>

int main()
{
	//graph types
	typedef adjacency_list<setS,vecS,undirectedS, MyNode, MyEdge> undirGraph;
	typedef adjacency_list<setS,vecS,directedS, MyNode, MyEdge> dirGraph;
	//typedef adjacency_list<vecS,vecS,undirectedS> U;

	//Time the Execution
	clock_t start = clock();

	vector<int> R;
	for(int i=0;i<100;i++){
		R.push_back(3);
	}
	MyGraph<undirGraph> G(R);

	cout<<num_vertices(G.g)<<endl;

	//End doing stuff
	clock_t ends = clock();
	cout <<endl<< "Running Time : "<< (double) (ends - start) / CLOCKS_PER_SEC << endl;
	return 0;
}