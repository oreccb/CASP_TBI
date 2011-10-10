#include "MyGraph.h"
#include <time.h>

int main()
{
	//graph types

	//use lists so we can remove edges
	typedef adjacency_list<setS,vecS,undirectedS, MyNode, MyEdge> undirGraph;

	typedef adjacency_list<setS,vecS,directedS, MyNode, MyEdge> dirGraph;
	//typedef adjacency_list<vecS,vecS,undirectedS> U;

	//Time the Execution of the code
	clock_t start = clock();

	MyGraph<undirGraph> G_real("BrianInMap.txt",0);
	
	
	cout<<"Real Graph:"<<endl;
	cout<<G_real.getNumVertices()<<endl;
	cout<<G_real.getNumEdges()<<endl<<endl;
	//cout<<G_real.computeClusteringCoefficient()<<endl<<endl;

	//////////////////////////////////////////////////TBI STUFF/////////////////////////
	//initialize seed
	srand( (int)time(NULL) );

	//need to add the SC node and set the SC_vertex value in class
	int SC_vert = G_real.TBI_add_vertex();
	G_real.setSC_vertex( SC_vert );
	cout<<"SC is vertex "<<G_real.getSC_vertex()<<endl;
	
	int total_itr = 200;
	for(int i=0; i<total_itr; i++)
	{
		G_real.Infiltrate(.1, //joion
						  .01, //leave
						  .4,  //pt
						  .2,  //po
						  1,   //alpha
						  5,   //budget
						  1);  //strategy
	}
	
	///////////////////////////////////////////////END TBI STUFF///////////////////////

	//cout<<"common neighbors: "<<G_real.numCommonNeighbors(144,145)<<endl;
	//cout<<G_real.getNumVertices()<<endl;
	//cout<<G_real.getNumEdges()<<endl<<endl;
	//G.printVE();

	//End doing stuff
	clock_t ends = clock();
	cout <<endl<< "Running Time : "<< (double) (ends - start) / CLOCKS_PER_SEC << endl;
	return 0;
}