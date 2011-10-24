#include "MyGraph.h"
#include <time.h>

double presentation2(double join,double leave, double pt, double po, int alpha, int budget, int strategy)
{
	typedef adjacency_list<setS,vecS,undirectedS, MyNode, MyEdge> undirGraph;

	//will need to write out data to file
	ofstream outfile;
	outfile.open("results_ENRON.csv");


	//need to add the SC node and set the SC_vertex value in class
	/*int SC_vert = G_real.TBI_add_vertex();
	G_real.setSC_vertex( SC_vert );
	cout<<"SC is vertex "<<G_real.getSC_vertex()<<endl;*/
	
	double TrustValue = 0;
	int total_itr = 147;
	int num_sim = 10;

	vector<double> results(total_itr);
	//init results vector
	for(int i =0;i<total_itr;i++)
	{
		results[i] = 0.0;
	}

	clock_t start = clock();
	cout<<"starting to read in data"<<endl;

	
	//MyGraph<undirGraph> G_real("Email-Enron.txt",0);
	MyGraph<undirGraph> G_real("BrianInMap.txt",1);
	//MyGraph<undirGraph> G_real(146,(double)((double)2636.0/(((double)146.0*(double)145.0)/2.0)) );

	G_real.debugfile.open("debug.txt");


	clock_t ends = clock();
	cout <<endl<< "Running Time to read in data : "<< (double) (ends - start) / CLOCKS_PER_SEC << endl;

	int SC_vert = G_real.TBI_add_vertex();
	G_real.setSC_vertex( SC_vert );
	//cout<<"SC is vertex "<<G_real.getSC_vertex()<<endl;

	//do the simulation 100 times and average results
	for(int k=0; k<num_sim; k++)
	{
		//HACK - need a new graph at each itr OBVIOUSLY!
		//vector<int> R;
		//R = G_tempreal.degOfVertices();
		//MyGraph<undirGraph> G_real(R);
		
		
		//reinitialize graph
		G_real.reinit();

		//initialize the P attribute
		G_real.initP(pt,po,alpha);

		//init or reinitialize the priority queue that is used in the greedy strategy
		G_real.initq();

		//do 'total_itr' timesteps of simulation
		for(int i=0; i<total_itr; i++)
		{
			TrustValue = G_real.Infiltrate(join, //join
							  leave, //leave
							  pt,  //pt
							  po,  //po
							  alpha,   //alpha
							  budget,   //budget (number of nodes requested at each itr)
							  strategy);  //strategy

			cout<<"Simulation "<<k<<", Trust Value at timestep "<<i<<": "<<TrustValue<<endl;
			//outfile<<"Trust Value at timestep "<<i<<": "<<TrustValue<<endl;
			results[i] = results[i] + TrustValue;
		}

	}


	//take average of the simulations
	for(int i =0;i<total_itr;i++)
	{
		results[i] = results[i] / (double)num_sim;
		cout<<"Average Trust Value at timestep "<<i<<": "<<results[i]<<endl;
		outfile<<i<<","<<results[i]<<endl;
	}
	
	//int timetillsuccess = 200;
	//find the time to get 30% of network
	//for(int i =0;i<total_itr;i++)
	//{
	//	if(results[i] > .15)
	//	{
	//		cout<<"Time to get 10% of network: "<<i<<" "<<results[i]<<endl;
	//		timetillsuccess = i;
	//		break;
	//	}
	//	//outfile<<i<<","<<results[i]<<endl;
	//}

	//return timetillsuccess;
	return results[ total_itr-1 ];
	
}