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

	//MyGraph<undirGraph> G_real("BrianInMap.txt",0);
	
	//MyGraph<undirGraph> G_real2(G_real);
	
	/*cout<<"Real Graph:"<<endl;
	cout<<G_real.getNumVertices()<<endl;
	cout<<G_real.getNumEdges()<<endl<<endl;*/
	//cout<<G_real.computeClusteringCoefficient()<<endl<<endl;

	//////////////////////////////////////////////////Presentation 2/////////////////////////
	//will need to write out data to file
	ofstream outfile;
	outfile.open("results.csv");
	
	//initialize seed
	srand( (int)time(NULL) );

	//need to add the SC node and set the SC_vertex value in class
	/*int SC_vert = G_real.TBI_add_vertex();
	G_real.setSC_vertex( SC_vert );
	cout<<"SC is vertex "<<G_real.getSC_vertex()<<endl;*/
	
	double TrustValue = 0;
	int total_itr = 200;
	int num_sim = 500;

	vector<double> results(total_itr);
	//init results vector
	for(int i =0;i<total_itr;i++)
	{
		results[i] = 0.0;
	}

	//temporary
	MyGraph<undirGraph> G_tempreal("BrianInMap.txt",0);

	//do the simulation 100 times and average results
	for(int k=0; k<num_sim; k++)
	{
		//HACK - need a new graph at each itr OBVIOUSLY!
		vector<int> R;
		R = G_tempreal.degOfVertices();
		MyGraph<undirGraph> G_real(R);
		//MyGraph<undirGraph> G_real("BrianInMap.txt",0);
		//MyGraph<undirGraph> G_real(146,(double)((double)2636.0/(((double)146.0*(double)145.0)/2.0)) );
		int SC_vert = G_real.TBI_add_vertex();
		G_real.setSC_vertex( SC_vert );
		//cout<<"SC is vertex "<<G_real.getSC_vertex()<<endl;

		//do 'total_itr' timesteps of simulation
		for(int i=0; i<total_itr; i++)
		{
			TrustValue = G_real.Infiltrate(0, //join
							  0, //leave
							  .3,  //pt
							  .5,  //po
							  1,   //alpha
							  1,   //budget (number of nodes requested at each itr)
							  1);  //strategy

			//cout<<"Trust Value at timestep "<<i<<": "<<TrustValue<<endl;
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