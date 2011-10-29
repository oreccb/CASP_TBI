#include "MyGraph.h"
#include <time.h>


int presentation2(double join,double leave, double pt, double po, int alpha, int budget, int strategy);

double benchmark(double join,double leave, double pt, double po, int alpha, int budget, int strategy);


void testpq()
{
	vector< pair<int,double*> > thing;
	vector<double> graphP;
	graphP.push_back(.2);
	graphP.push_back(.34);
	graphP.push_back(.1);
	graphP.push_back(.7);
	graphP.push_back(.5);
	graphP.push_back(.05);
	graphP.push_back(.98);

	for(int i=0; i<7; i++)
	{
		double * ptr = &graphP[i];
		thing.push_back(pair<int,double*>(i,ptr) );
	}

	for(int i=0; i<7; i++)
	{
		
		cout<<thing[i].first<<","<<*thing[i].second<<endl;
	}
	cout<<endl<<endl;
	MyPriorityQueue PQ(thing);
	

	for(int i=0; i<7; i++)
	{
		
		cout<<PQ.getElement(i).first<<","<<*PQ.getElement(i).second<<endl;
	}
	cout<<endl<<endl;

	PQ.pop();

	for(int i=0; i<7; i++)
	{
		
		cout<<PQ.getElement(i).first<<","<<*PQ.getElement(i).second<<endl;
	}
	return;
}




int main()
{
	//graph types

	//use lists so we can remove edges
	typedef adjacency_list<setS,vecS,undirectedS, MyNode> undirGraph;

	typedef adjacency_list<setS,vecS,directedS, MyNode, MyEdge> dirGraph;
	//typedef adjacency_list<vecS,vecS,undirectedS> U;

	//Time the Execution of the code
	clock_t start = clock();

	//initialize seed
	srand( (int)time(NULL) );

	//MyGraph<undirGraph> G_real("Email-Enron.txt",0);
	
	MyGraph<undirGraph> G_real("BrianInMap.txt",1);
	//MyGraph<undirGraph> G_real2(G_real);
	
	/*cout<<"Real Graph:"<<endl;
	cout<<G_real.getNumVertices()<<endl;
	cout<<G_real.getNumEdges()<<endl<<endl;
	cout<<G_real.computeClusteringCoefficient()<<endl<<endl;*/

	//////////////////////////////////////////////////Presentation 2/////////////////////////
	//ofstream out;
	//out.open("thresholdtimeresults_ENRONFIRSTTESTGreedy.csv");
	//out.open("LinkedInTest_Etype_greedy.csv");

	/*int timeforsuc = 0;
	for(double pt = .05; pt<1.0; pt = pt + .05)
	{

		timeforsuc = presentation2(0.0,0.0, pt, .5, 2, 1, 2);
		cout<<pt<<","<<timeforsuc<<endl;
		out<<pt<<","<<timeforsuc<<endl;
	}*/

	//presentation2(0,0,.5,.5,1,1,2);

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

	//Real benchmark stuff
	//cout<<"frac due to ego: "<<benchmark(0.0,0.0,.0,.5,2,146,2)<<endl;


	//cout<<"common neighbors: "<<G_real.numCommonNeighbors(144,145)<<endl;
	//cout<<G_real.getNumVertices()<<endl;
	//cout<<G_real.getNumEdges()<<endl<<endl;
	//G.printVE();


	//testing priority queue
	//testpq();


	//End doing stuff
	clock_t ends = clock();
	cout <<endl<< "Running Time : "<< (double) (ends - start) / CLOCKS_PER_SEC << endl;
	return 0;
}