#include "MyGraph.h"
#include <time.h>

int presentation2(double join,double leave, double pt, double po, int alpha, int budget, int strategy);

int main()
{
	//graph types

	//use lists so we can remove edges
	typedef adjacency_list<setS,vecS,undirectedS, MyNode, MyEdge> undirGraph;

	typedef adjacency_list<setS,vecS,directedS, MyNode, MyEdge> dirGraph;
	//typedef adjacency_list<vecS,vecS,undirectedS> U;

	//Time the Execution of the code
	clock_t start = clock();

	//initialize seed
	srand( (int)time(NULL) );

	//MyGraph<undirGraph> G_real("BrianInMap.txt",0);
	
	//MyGraph<undirGraph> G_real2(G_real);
	
	/*cout<<"Real Graph:"<<endl;
	cout<<G_real.getNumVertices()<<endl;
	cout<<G_real.getNumEdges()<<endl<<endl;*/
	//cout<<G_real.computeClusteringCoefficient()<<endl<<endl;

	//////////////////////////////////////////////////Presentation 2/////////////////////////
	ofstream out;
	out.open("thresholdtimeresults.csv");
	
	int timeforsuc;
	for(double pt = .05; pt<1.0; pt = pt + .05)
	{

		timeforsuc = presentation2(0,0, pt , .5, 1, 1, 1);
		cout<<pt<<","<<timeforsuc<<endl;
		out<<pt<<","<<timeforsuc<<endl;
	}

	///////////////////////////////////////////////END PRESENTATION 2///////////////////////


	
	//cout<<"common neighbors: "<<G_real.numCommonNeighbors(144,145)<<endl;
	//cout<<G_real.getNumVertices()<<endl;
	//cout<<G_real.getNumEdges()<<endl<<endl;
	//G.printVE();

	//End doing stuff
	clock_t ends = clock();
	cout <<endl<< "Running Time : "<< (double) (ends - start) / CLOCKS_PER_SEC << endl;
	return 0;
}