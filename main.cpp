#include "MyGraph.h"
#include <time.h>

int main()
{
	//graph types

	//use lists so we can remove edges
	typedef adjacency_list<listS,vecS,undirectedS, MyNode, MyEdge> undirGraph;

	typedef adjacency_list<setS,vecS,directedS, MyNode, MyEdge> dirGraph;
	//typedef adjacency_list<vecS,vecS,undirectedS> U;

	//Time the Execution of the code
	clock_t start = clock();

	MyGraph<undirGraph> G_real("BrianInMap.txt",0);
	
	
	cout<<"Real Graph:"<<endl;
	cout<<G_real.getNumVertices()<<endl;
	cout<<G_real.getNumEdges()<<endl<<endl;
	//cout<<G_real.computeClusteringCoefficient()<<endl<<endl;

	G_real.Infiltrate(1, .01, .4, .2, 1, 5, 0);
	//cout<<G_real.getNumVertices()<<endl;
	//cout<<G_real.getNumEdges()<<endl<<endl;

	//G.printVE();

	//End doing stuff
	clock_t ends = clock();
	cout <<endl<< "Running Time : "<< (double) (ends - start) / CLOCKS_PER_SEC << endl;
	return 0;
}