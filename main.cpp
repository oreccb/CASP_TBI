#include "MyGraph.h"
#include <time.h>

double presentation2(double join,double leave, double pt, double po, int alpha, int budget, int strategy);

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

	//MyGraph<undirGraph> G_real("Email-Enron.txt",0);
	
	//MyGraph<undirGraph> G_real("BrianInMap.txt",1);
	//MyGraph<undirGraph> G_real2(G_real);
	
	/*cout<<"Real Graph:"<<endl;
	cout<<G_real.getNumVertices()<<endl;
	cout<<G_real.getNumEdges()<<endl<<endl;
	cout<<G_real.computeClusteringCoefficient()<<endl<<endl;*/

	//////////////////////////////////////////////////Presentation 2/////////////////////////
	ofstream out;
	out.open("thresholdtimeresults_ENRONFIRSTTESTRANDOM.csv");
	
	int timeforsuc;
	double frac = 0;
	for(double pt = .05; pt<1.0; pt = pt + .05)
	{

		frac = presentation2(0.0,0.0, pt, .5, 2, 1, 1);
		cout<<pt<<","<<frac<<endl;
		out<<pt<<","<<frac<<endl;
	}

	//presentation2(0,0,.3,.5,2,1,1);

	///////////////////////////////////////////////END PRESENTATION 2///////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////
	////Find the average benchmark number of nodes that are infiltrated based on just 'ego'////////////
	//	//MyGraph<undirGraph> G_real("BrianInMap.txt",1);
	//	MyGraph<undirGraph> G_real("Email-Enron.txt",0);
	//	int SC_vert = G_real.TBI_add_vertex();
	//	G_real.setSC_vertex( SC_vert );

	//	ofstream outfile;
	//	outfile.open("benchmark.txt");
	//	int TrustValue = 0;
	//	int total_itr = 30000;
	//	
	//	cout<<"starting the simulation"<<endl;
	//	for(int i=0; i<total_itr; i++)
	//	{	
	//		cout<<"at iteration "<<i<<endl;
	//		TrustValue = G_real.Infiltrate_alt(0, //join
	//						  0, //leave
	//						  0,  //pt
	//						  .5,  //po
	//						  2,   //alpha
	//						  1,   //budget (number of nodes requested at each itr)
	//						  1);  //strategy

	//		cout<<"Trust Value at timestep "<<i<<": "<<TrustValue<<endl;
	//		outfile<<"Trust Value at timestep "<<i<<": "<<TrustValue<<endl;
	//		//results[i] = results[i] + TrustValue;
	//	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////

	
	//cout<<"common neighbors: "<<G_real.numCommonNeighbors(144,145)<<endl;
	//cout<<G_real.getNumVertices()<<endl;
	//cout<<G_real.getNumEdges()<<endl<<endl;
	//G.printVE();

	//End doing stuff
	clock_t ends = clock();
	cout <<endl<< "Running Time : "<< (double) (ends - start) / CLOCKS_PER_SEC << endl;
	return 0;
}