#ifndef MYGRAPH_H
#define MYGRAPH_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <math.h>
#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/clustering_coefficient.hpp>
#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/exterior_property.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/random/linear_congruential.hpp>

//file with bundled proporties to use in graphs
#include "Bundled_Proporties.h"

using namespace std;
using namespace boost;

#define DEFAULT_GRAPH_SIZE 100

//Graph Class
template <class graphtype>
class MyGraph
{
private:
	int nodeLeft(double leave);						//removes vertices with probability 'leave'
	bool nodeJoin(double join);						//adds a node based on probability 'join', returns a bool if a vertex was added
	vector<int> genRequestList(int strategy, int budget);  //generate the list of nodes SC will request based on strat and budget

	graphtype g;									//graph
	int SC_vertex;									//index of stealth company node
public:
	
	MyGraph();										//create graph with DEFAULT_GRAPH_SIZE verticies and no edges
	MyGraph(int N, double p);						//create random erdos renyi graph
	MyGraph(vector<int> R);							//create graph with degree distribution R
	MyGraph(int N, int d, double p);				//create graph with mixed preferential attachment
	MyGraph(string datafile, int mode);				//create graph based on real datafile (mode tells what file type it is to read)

	float computeClusteringCoefficient();			//compute the clustering coefficient of the graph
	void computeDegreeDistribution();				//compute the degree distribution of the graph

	vector<int> degOfVertices();					//output vector<int> that is the degree of each vertex for each slot

	int getNumVertices();							//return the number of vertices of g
	int getNumEdges();								//return the number of edges of g
	void printVE();									//print vertices and edges

	int TBI_add_vertex();							//add a vertex to g and return the index of it
	void setSC_vertex( int v );						//set the SC vertex index
	int getSC_vertex();								//get the SC vertex index

	int Infiltrate(									//Run the trust based infiltration simulation
		double join,	//probability that an actor joins the network
		double leave,	//the probability that each node leaves the network during a timestep
		double pt,		//probability that actor will accept connection due to one shared connection with stealth company
		double po,		//base probability that an actor will accept connection due to ego
		double alpha,	//exponentiating factor in calculating ego probability
		int budget,	//Budget of stealth company (number of connection requests sent out at each timestep
		int stratagy //stratagy of stealth company is picking who they request
		);								
	
};


/////////////////////////////////CONSTRUCTORS///////////////////////////////////////////////

template <class graphtype>
MyGraph<graphtype>::MyGraph()
{
	g = graphtype(DEFAULT_GRAPH_SIZE);
	SC_vertex = -1;
}


//Make an erdos renyi graph with N nodes and edges with probability p
template <class graphtype>
MyGraph<graphtype>::MyGraph(int N, double p)
{
	typedef boost::erdos_renyi_iterator<boost::minstd_rand, graphtype> ERGen;
	boost::minstd_rand gen;
	g = graphtype(ERGen(gen, N, p), ERGen(), N);

	SC_vertex = -1;
}

//We are making a random graph with specified deg distribution R
//THIS MIGHT ONLY WORK FOR undirected graph right now, Fix if needed for directed graph
template <class graphtype>
MyGraph<graphtype>::MyGraph(vector<int> R)
{
	int i,j;
	graph_traits<graphtype>::vertex_iterator vi, vi_end;
    graph_traits<graphtype>::vertex_iterator vi2, vi2_end;
    graph_traits<graphtype>::edge_iterator ei, ei_end;

	typename graph_traits<graphtype>::edge_descriptor ed;
    bool inserted;

	//initialize seed
	srand( (int)time(NULL) );

	//make degree list(If R is [1 2 2 1] we create [0 1 1 2 2 3]
	vector< pair<int,int> > ordered_deg_list;
	for(i = 0; i<R.size(); i++)
	{
		for(j=0; j<R[i]; j++)
		{
			ordered_deg_list.push_back(pair<int,int>(i,0));
			//cout<< i <<" ";
		}
	}
	cout<<endl;
	
	//permute list
	random_shuffle( ordered_deg_list.begin(),ordered_deg_list.end() );
	//print out the new shuffled list
	/*for(i=0;i<ordered_deg_list.size();i++)
	{
		cout<<"("<<ordered_deg_list[i].first<<",";
		cout<<ordered_deg_list[i].second<<") ";
	}
	cout<<endl;*/

	g = graphtype(R.size());

	//sets the vertex id
	for (tie(vi, vi_end) = vertices(g), i=0 ; vi != vi_end; ++i, ++vi)
	{
		g[*vi].x = i;
	}
	
	//Add all the edges using the given algorithm described in the HW
	for (i=0; i<ordered_deg_list.size()-1; i++)
	{
		if(ordered_deg_list[i].second == 0)
		{
			for(j=i+1;j<ordered_deg_list.size(); j++)
			{
					if(ordered_deg_list[j].second != 1 && ordered_deg_list[i] != ordered_deg_list[j])
					{
					
						//add edge
						tie(ed, inserted) = add_edge(g[ordered_deg_list[i].first].x, g[ordered_deg_list[j].first].x, g);

						if(inserted){
							ordered_deg_list[j].second = 1;
							ordered_deg_list[i].second = 1;
							break;
						}
						else
						{
							continue;
						}
					}
			}
		}
    }
	

	//print out the verticies and edges
	/*for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
	{
		std::cout<<g[*vi].x <<" ";

	}
	std::cout<<std::endl;
	for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
        cout << "edge: " << source(*ei, g) << " to " << target(*ei, g) << endl;
    }
    cout<<endl;*/

	SC_vertex = -1;
	
}

//this constructor creates a graph with N nodes based on the mixed preferential attachment model for growth
//At each step, on node is added and the node makes d new edges
//each new edge will attach with probability p of being random and probability (1-p) of being preferentially attached
template <class graphtype>
MyGraph<graphtype>::MyGraph(int N, int d, double p)
{
	//working variables
	graph_traits<graphtype>::edge_descriptor ed;
	bool inserted;
	graph_traits<graphtype>::vertex_iterator vi2, vi2_end;
	graph_traits<graphtype>::edge_iterator ei, ei_end;
	graph_traits<graphtype>::vertex_descriptor vert1, vert2;

	double randnum = 0;
	int g_edges = 0;
	int g_verticies = 0;
	int ind_temp;

	vector<int> deg;

	//initialize seed
	srand( (int)time(NULL) );

	//need to initialize the graph with 2 nodes and an edge between them??
	vert1 = add_vertex(g);
	vert2 = add_vertex(g);
	add_edge(vert1,vert2,g);
	g_verticies = 2;
	g_edges = 1;
	
	//start tracking degrees of verticies
	deg.push_back(1);
	deg.push_back(1);

	//add a new vertex at every step in the growth process
	for(int i=2; i<N; i++)
	{
		vert1 = add_vertex(g);
		//loop through the d edges added
		for(int j=0; j<d; j++)
		{
			//decide whether to attach randomly or with pref attachment
			randnum = rand() / (RAND_MAX + 1);

			//attach randomly if randnum < p
			if(randnum < p)
			{
				//randomly find the vartex to attach an edge to
				ind_temp = rand() % g_verticies;
				tie(ed, inserted) = add_edge(i, ind_temp, g);
			}
			//attach with pref attachment if randnum >= p
			else
			{
				
				int r = rand() % g_edges;				//generate rand number between
				double s = 0;								//accumulation of degrees
				int k = 0;								//index degree vector
				while(r>s)
				{
					s = s + ((double)deg[k]/2.0);					//need to divide by two since we are doing an undirected graph
					k++;
				}
				deg[k]++;

			}
			g_edges++;
		}
		g_verticies++;
	}

	SC_vertex = -1;
}

//make graph based on real data from file
template <class graphtype>
MyGraph<graphtype>::MyGraph(string datafile, int mode)
{
	//working variables
	graph_traits<graphtype>::edge_descriptor ed;
	bool inserted;
	graph_traits<graphtype>::vertex_iterator vi2, vi2_end;
	graph_traits<graphtype>::edge_iterator ei, ei_end;
	graph_traits<graphtype>::vertex_descriptor vert1, vert2;

	int temp;
	string a_name;
	string line;
	vector< vector<int> > edgelist;
	vector<int> nextedges;

	//////////This mode is for the sample linkedIn data from Brian's InMap/////////////
	if(mode == 0)              
	{
		ifstream LI_inp;
		//We should be opening BrianInMap.txt which is the sample network data
		LI_inp.open(datafile.c_str());
		
		//loop through each line of file
		while ( getline(LI_inp, line) )
		{
			//make each line like an input stream so we can manipulate easier
			istringstream iss(line);
			
			//cout<<line<<endl; //dubugging

			//add the vertex and add attach the name to the vertex
			vert1 = add_vertex(g);
			iss >> temp;
			iss >> a_name;
			g[vert1].name = a_name;
			//cout<<g[vert1].name<<endl;  //debugging

			nextedges.clear();   //clear the temp vector!
			while (iss >> temp)
			{
				//make a vector of the vertex's edges
				nextedges.push_back(temp);
			}
			edgelist.push_back(nextedges);
			
		}
		
		//cout<<"size:"<<edgelist.size()<<endl;
		//put the edges in the graph from edgelist vector
		for(int i=0; i<edgelist.size(); i++)
		{
			for(int j=0; j<edgelist[i].size(); j++)
			{
				//add edge between node i and edgelist[i][j]
				tie(ed, inserted) = add_edge(i, edgelist[i][j], g);
			}

		}
	}


	else if(mode == 1)
	{
		//this will probably be for the real linkedIn data file if/when I get it
	}
	else
	{
		cout<<"You specified an unsupported mode....exiting"<<endl;
		exit(0);
	}

	SC_vertex = -1;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////END CONSTRUCTORS//////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

//This function computes the clustering coefficient of the graph g
template <class graphtype>
float MyGraph<graphtype>::computeClusteringCoefficient()
{
	typedef exterior_vertex_property<graphtype, float> ClusteringProperty;
	typedef ClusteringProperty::container_type ClusteringContainer;
	typedef ClusteringProperty::map_type ClusteringMap;

	ClusteringContainer coefs(num_vertices(g));
    ClusteringMap cm(coefs, g);
    float cc = all_clustering_coefficients(g, cm);
	//cout<<"clustering coefficient: "<<cc<<endl;

	return cc;
	
}

//output vector<int> that is the degree of each vertex for each slot
template <class graphtype>
vector<int> MyGraph<graphtype>::degOfVertices()
{
	vector<int> R;
	int deg;
	graph_traits<graphtype>::vertex_iterator u, u_end;
	//degree_size_type deg2;

	for (tie(u, u_end) = vertices(g); u != u_end; ++u)
	{
		deg = out_degree(*u,g);
		R.push_back(deg);
	}
	

	return R;
}


//This function computes the degree distribution of the graph g
//CAUTION:THIS MIGHT NOT BE WORKING RIGHT
template <class graphtype>
void MyGraph<graphtype>::computeDegreeDistribution()
{
	
	vector<int> R;
	int deg;
	int max_deg;
	graph_traits<graphtype>::vertex_iterator u, u_end;

	ofstream degOutfile;
	degOutfile.open("degOutfile.csv");

	for (tie(u, u_end) = vertices(g); u != u_end; ++u)
	{
		deg = out_degree(*u,g);
		R.push_back(deg);
	}
	
	sort( R.begin(), R.end() );
	//find max degree of graph
	max_deg = *max_element( R.begin(),R.end() );
	vector<int> deg_dist(max_deg);
	int temp;
	//loop through the different possible degrees and find the degree distribution 
	for(int i=0; i<max_deg; i++)
	{
		temp = 0;
		for(int j=max_deg; j>=0; j--)
		{
			if(R[j] >= i)
			{
				temp = temp + 1;
			}
		}
		deg_dist[i] = temp;

	}

	for(int i=0; i<R.size(); i++)
	{
		cout<<R[i]<<" ";
	}
	cout<<endl;
	cout<<max_deg<<endl;
	cout<<"start of deg dist: ";
	for(int i=0; i<deg_dist.size(); i++)
	{
		cout<<deg_dist[i]<<" ";
		degOutfile<<i<<","<<deg_dist[i]<<endl;
	}
	cout<<endl;

	return;
}

template <class graphtype>
int MyGraph<graphtype>::getNumVertices()
{
	return num_vertices(g);
}
template <class graphtype>
int MyGraph<graphtype>::getNumEdges()
{
	return num_edges(g);
}

//print out the verticies and edges
template <class graphtype>
void MyGraph<graphtype>::printVE()
{
	
	graph_traits<graphtype>::vertex_iterator vi, vi_end;
    graph_traits<graphtype>::edge_iterator ei, ei_end;

	for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
	{
		//std::cout<<g[*vi].x <<" ";
		cout<<*vi<<" ";
	}
	std::cout<<std::endl;
	for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
        cout << "edge: " << source(*ei, g) << " to " << target(*ei, g) << endl;
    }
    cout<<endl;

	return;
}

template <class graphtype>
int MyGraph<graphtype>::TBI_add_vertex() 
{
	return (int)add_vertex(g);
}

template <class graphtype>
void MyGraph<graphtype>::setSC_vertex( int v ) 
{	
	SC_vertex = v;
}
template <class graphtype>
int MyGraph<graphtype>::getSC_vertex() 
{	
	return SC_vertex;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////Trust Based Infiltration Simulation//////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <class graphtype>
int MyGraph<graphtype>::nodeLeft(double leave) 
{
	double randnum = 0;
	int num_removed = 0;
	graph_traits<graphtype>::vertex_iterator vi, vi_end;

	for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
	{
		randnum = (double)rand() / (double)(RAND_MAX + 1);
		if( randnum < leave && *vi != SC_vertex )  //check if the random number is less than the probability and if the vertex is not the SC
		{
			//cout<<leave<<"  "<<randnum<<endl;
			//cout<<"DEBUG: removing vertex "<<*vi<<endl;

			clear_vertex(*vi,g);
			remove_vertex(*vi,g);  //might be a bug here!  does removing mess with vertex iterator
			num_removed++;
		}
	}
	return num_removed;
}

template <class graphtype>
bool MyGraph<graphtype>::nodeJoin(double join)
{
	bool added = false;
	double randnum = 0;
	graph_traits<graphtype>::vertex_descriptor u;

	randnum = (double)rand() / (double)(RAND_MAX + 1);
	if( randnum < join )
	{
		u = add_vertex(g);
		added = true;
		cout<<"DEBUG: added vertex "<<u<<endl;
	}

	return added;
}

//return a vector of vertex indicies that the SC will request a connection based on the strategy
//Budget is an integer specifying how many requests the SC can send per timestep
//Strategy 1: random selection
//Strategy 2: preferential selection (send request the connections of people SC is currently connected to)
template <class graphtype>
vector<int> MyGraph<graphtype>::genRequestList(int strategy, int budget)
{
	vector<int> nodesRequested;
	int randvert = 0;
	
	switch(strategy)
	{
	case 1:
		{
			cout<<"DEBUG: Using Strategy 1 to pick requested nodes"<<endl;
			while(nodesRequested.size() < budget)
			{
				randvert = rand() % num_vertices(g);
				if(randvert != SC_vertex)
				{
					nodesRequested.push_back(randvert);
				}
			}
		}
		break;
	case 2:
		{
			cout<<"DEBUG: Using Strategy 2 to pick requested nodes"<<endl;

		}
		break;
	default:
		{
			cout<<"Invalid strategy, Exiting..."<<endl;
			exit(1);
		}
		break;
	}
	

	//DEBUG STUFF//////////////
	/*cout<<"DEBUG: Nodes to Request-- ";
	for(int i =0;i<nodesRequested.size();i++)
	{
		cout<<nodesRequested[i]<<" ";
	}
	cout<<endl;*/
	///////////////////////

	return nodesRequested;
}


template <class graphtype>
int MyGraph<graphtype>::Infiltrate(double join, double leave, double pt,double po, double alpha, int budget, int strategy)
{
	vector<int> nodesRequested;
	//initialize seed  -- do this in main
	//srand( (int)time(NULL) );

	//check each node if it leaves the network
	nodeLeft(leave);

	//check if a node joins the network
	nodeJoin(join);

	//find the nodes that the stealth company requested
	nodesRequested = genRequestList(strategy, budget);

	//find the nodes that accepted the connection requests and add them to the network

	//calculate the trust value and that is the return value


	return 0;

}









#endif
