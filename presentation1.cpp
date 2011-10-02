#include "MyGraph.h"
#include <time.h>

int presentation1()
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

	
	//make R-type random graph based on deg distribution
	vector<int> R;
	R = G_real.degOfVertices();
	MyGraph<undirGraph> G_R(R);
	
	//G_real.computeDegreeDistribution();
	double p = .8;
	MyGraph<undirGraph> G_PA(n,e/n,p); //need to look at this function again-not getting a correct CC
	
	MyGraph<undirGraph> G_E(n,(double)((double)e/(((double)n*(double)(n-1))/2.0)) );
	//cout<<"E prob:"<<(double)((double)e/(((double)n*(double)(n-1))/2.0))<<endl;

	//cout<<(double)(2636.0/((147.0*146.0)/2.0))<<"blah"<<endl;
	cout<<"Real Graph:"<<endl;
	cout<<num_vertices(G_real.g)<<endl;
	cout<<num_edges(G_real.g)<<endl;
	cout<<G_real.computeClusteringCoefficient()<<endl<<endl;

	cout<<"R-type Random Graph:"<<endl;
	cout<<num_vertices(G_R.g)<<endl;
	cout<<num_edges(G_R.g)<<endl;
	cout<<G_R.computeClusteringCoefficient()<<endl<<endl;

	cout<<"E-type Random Graph:"<<endl;
	cout<<num_vertices(G_E.g)<<endl;
	cout<<num_edges(G_E.g)<<endl;
	cout<<G_E.computeClusteringCoefficient()<<endl<<endl;

	cout<<"PA-type Random Graph:"<<endl;
	cout<<num_vertices(G_PA.g)<<endl;
	cout<<num_edges(G_PA.g)<<endl;
	cout<<G_PA.computeClusteringCoefficient()<<endl<<endl;

	//G.printVE();

	//End doing stuff
	clock_t ends = clock();
	cout <<endl<< "Running Time : "<< (double) (ends - start) / CLOCKS_PER_SEC << endl;
	return 0;
}