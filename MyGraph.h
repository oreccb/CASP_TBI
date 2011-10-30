#ifndef MYGRAPH_H
#define MYGRAPH_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <math.h>
#include <algorithm>
#include <queue>
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/clustering_coefficient.hpp"
#include "boost/graph/undirected_graph.hpp"
#include "boost/graph/exterior_property.hpp"
#include "boost/graph/erdos_renyi_generator.hpp"
#include "boost/random/linear_congruential.hpp"

//#include "sortfile.h"
#include "MyPriorityQueue.h"

//file with bundled proporties to use in graphs
#include "Bundled_Proporties.h"

using namespace std;
using namespace boost;

#define DEFAULT_GRAPH_SIZE 100

//#define DEBUG

bool mysortt(pair<int,double*> pa, pair<int,double*> pb);

//Graph Class
template <class graphtype>
class MyGraph
{
private:
	int nodeLeft(double leave);						//removes vertices with probability 'leave'
	bool nodeJoin(double join);						//adds a node based on probability 'join', returns a bool if a vertex was added
	vector<int> genRequestList(double pt, double po, double alpha, int strategy, int budget);  //generate the list of nodes SC will request based on strat and budget
	int update(vector<int> nodesRequested, double pt, double po, double alpha);  //update network with connections that were accepted due to model
	double calcTrustVal();							//calculate the trust value
	double calcProbabilities(double pt, double po, double alpha, int nodeDeg, int numCN);   //calculates the probability that a certain node will accept the request based on the model

	graphtype g;									//graph
	int SC_vertex;									//stealth company vertex
	//vector<int> CNeigh;								//number of common neighbors of each vertex CAUTION: NOW IN bundled_proporties
	//priority_queue< pair<int,double>*, vector< pair<int,double>* >, CompareNode> q;  //this won't really work
	list< pair<int,double*> > q;

	
public:
	ofstream debugfile;

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

	int numCommonNeighbors(int SC, int v);			//Calc the number of common neighbors between 2 vertices of g
	int updateNumCommonNeighborsANDP(double pt,double po, double alpha,int u);  //update the CNeigh vector with the new common meighbor count for adjacent nodes to one that has been added

	void reinit();
	void initP(double pt, double po, double alpha);  //initialize theh P attribute for each node
	void updateP(double pt, double po, double alpha, int u);  //update the P attribute of a node
	void initq();    //initialize the priority queue that is used for the greedy strat(#2)

	double Infiltrate(									//Run the trust based infiltration simulation
		double join,	//probability that an actor joins the network
		double leave,	//the probability that each node leaves the network during a timestep
		double pt,		//probability that actor will accept connection due to one shared connection with stealth company
		double po,		//base probability that an actor will accept connection due to ego
		double alpha,	//exponentiating factor in calculating ego probability
		int budget,	//Budget of stealth company (number of connection requests sent out at each timestep
		int stratagy //stratagy of stealth company is picking who they request
		);	

	double Infiltrate_alt(double join, double leave, double pt,double po, double alpha, int budget, int strategy);
	int calcBenchVal();  //calc the benchmark value
	
};


/////////////////////////////////CONSTRUCTORS///////////////////////////////////////////////

template <class graphtype>
MyGraph<graphtype>::MyGraph()
{

	g = graphtype(DEFAULT_GRAPH_SIZE);
	SC_vertex = -1;
	
	typename graph_traits<graphtype>::vertex_iterator vi, vi_end;
	for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
	{
		g[*vi].covered = false;
		g[*vi].CNeigh = 0;
	}
}


//Make an erdos renyi graph with N nodes and edges with probability p
template <class graphtype>
MyGraph<graphtype>::MyGraph(int N, double p)
{
	typedef boost::erdos_renyi_iterator<boost::minstd_rand, graphtype> ERGen;
	boost::minstd_rand gen;
	g = graphtype(ERGen(gen, N, p), ERGen(), N);

	typename graph_traits<graphtype>::vertex_iterator vi, vi_end;
	for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
	{
		g[*vi].covered = false;
		g[*vi].CNeigh = 0;
	}

	
}

//We are making a random graph with specified deg distribution R
//THIS MIGHT ONLY WORK FOR undirected graph right now, Fix if needed for directed graph
template <class graphtype>
MyGraph<graphtype>::MyGraph(vector<int> R)
{
	unsigned int i;
	unsigned int j;
	typename graph_traits<graphtype>::vertex_iterator vi, vi_end;
   // typename graph_traits<graphtype>::vertex_iterator vi2, vi2_end;
    typename graph_traits<graphtype>::edge_iterator ei, ei_end;

	typename typename graph_traits<graphtype>::edge_descriptor ed;
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

	//set the 'covered' value of each vertex to false
	typename graph_traits<graphtype>::vertex_iterator vii, vii_end;
	for (tie(vii, vii_end) = vertices(g); vii != vii_end; ++vii)
	{
		g[*vii].covered = false;
		g[*vii].CNeigh = 0;
	}

	SC_vertex = -1;
	
}

//this constructor creates a graph with N nodes based on the mixed preferential attachment model for growth
//At each step, on node is added and the node makes d new edges
//each new edge will attach with probability p of being random and probability (1-p) of being preferentially attached
template <class graphtype>
MyGraph<graphtype>::MyGraph(int N, int d, double p)
{
	//working variables
	typename graph_traits<graphtype>::edge_descriptor ed;
	bool inserted;
	typename graph_traits<graphtype>::vertex_iterator vi2, vi2_end;
	typename graph_traits<graphtype>::edge_iterator ei, ei_end;
	typename graph_traits<graphtype>::vertex_descriptor vert1,vert2;

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
	for(unsigned int i=2; i<N; i++)
	{
		vert1 = add_vertex(g);
		//loop through the d edges added
		for(unsigned int j=0; j<d; j++)
		{
			//decide whether to attach randomly or with pref attachment
			randnum = (double)rand() / (double)(RAND_MAX + 1);

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

	//set the 'covered' value of each vertex to false
	typename graph_traits<graphtype>::vertex_iterator vi, vi_end;
	for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
	{
		g[*vi].covered = false;
		g[*vi].CNeigh = 0;
	}

	SC_vertex = -1;
}

//make graph based on real data from file
template <class graphtype>
MyGraph<graphtype>::MyGraph(string datafile, int mode)
{
	//working variables
	typename graph_traits<graphtype>::edge_descriptor ed;
	bool inserted;
	//typename graph_traits<graphtype>::vertex_iterator vi2, vi2_end;
	graph_traits<graphtype>::edge_iterator ei, ei_end;
	typename graph_traits<graphtype>::vertex_descriptor vert1, u, v;

	int a;
	int b;
	int temp;
	string a_name;
	string line;
	vector< vector<int> > edgelist;
	vector<int> nextedges;

	
	//This mode is for the Enron email data set
	//This data was obtained from http://snap.stanford.edu/data/email-Enron.html
	if(mode == 0)              
	{
		ifstream LI_inp;
		//We should be opening Email-Enron.txt which is the data of all the enron communication
		LI_inp.open(datafile.c_str());
		
		for(int i=0; i<36692; i++)
		{
			add_vertex(g);
			//cout<<"adding vertex "<<i<<endl;
		}


		//loop through each line of file, around 367662
		//while ( getline(LI_inp, line) )
		//{
		//	//make each line like an input stream so we can manipulate easier
		//	istringstream iss(line);
		//	
		//	//cout<<line<<endl; //dubugging
		//
		//	iss >> a;
		//	iss >> b;
		//
		//	//u = vertex(a,g);
		//	//v = vertex(b,g);
		//		 
		//	//add edge between nodes u and v
		//	tie(ed, inserted) = add_edge(a,b, g);
		//}

		int count = 0;
		while ( !LI_inp.eof() )
		{
			
			
			//cout<<"count: "<<++count<<endl; //dubugging
		
			LI_inp >> a;
			LI_inp >> b;
		
			//u = vertex(a,g);
			//v = vertex(b,g);
				
			//cout<<"DUBUGGING: adding edge :"<<a<<" "<<b<<endl;
			//add edge between nodes u and v
			tie(ed, inserted) = add_edge(a,b, g);
		}

	}
	///////////////////////////////////////////////////////////////////////////////////
	//////////This mode is for the sample linkedIn data from Brian's InMap/////////////
	else if(mode == 1)
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
				nextedges.push_back(temp - 1);
			}
			edgelist.push_back(nextedges);
			
		}
		
		//cout<<"size:"<<edgelist.size()<<endl;
		//put the edges in the graph from edgelist vector
		for(unsigned int i=0; i<edgelist.size(); i++)
		{
			for(unsigned int j=0; j<edgelist[i].size(); j++)
			{
				//add edge between node i and edgelist[i][j]
				tie(ed, inserted) = add_edge(i, edgelist[i][j], g);
			}

		}
	}
	else
	{
		cout<<"You specified an unsupported mode....exiting"<<endl;
		exit(0);
	}

	//set the 'covered' value of each vertex to false
	typename graph_traits<graphtype>::vertex_iterator vi, vi_end;
	for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
	{
		g[*vi].covered = false;
		g[*vi].CNeigh = 0;
	}

	SC_vertex = -1;
	
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////END CONSTRUCTORS//////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

template <class graphtype>
void MyGraph<graphtype>::reinit()
{
	typename graph_traits<graphtype>::vertex_descriptor SC;
	SC = vertex(SC_vertex,g);
	//clear edges of SC
	clear_vertex(SC,g);

	//make the covered and CNeigh attributes zero again
	typename graph_traits<graphtype>::vertex_iterator vi, vi_end;
	for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi) //loop through the adjacent vertices and update their CNeigh val
	{
		g[*vi].CNeigh = 0;
		g[*vi].covered = false;
	}
	

	return;
}



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
	typename graph_traits<graphtype>::vertex_iterator u, u_end;
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
	typename graph_traits<graphtype>::vertex_iterator u, u_end;

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
	for(unsigned int i=0; i<max_deg; i++)
	{
		temp = 0;
		for(unsigned int j=max_deg; j>=0; j--)
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
	
	typename graph_traits<graphtype>::vertex_iterator vi, vi_end;
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
	return add_vertex(g);
}

template <class graphtype>
void MyGraph<graphtype>::setSC_vertex( int v ) 
{	
	typename graph_traits<graphtype>::vertex_descriptor SC;
	SC_vertex = v;
	SC = vertex(SC_vertex, g);
	g[SC].name = "Stealth Company";
	g[SC].covered = true;
	
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
	typename graph_traits<graphtype>::vertex_iterator vi, vi_end;
	typename graph_traits<graphtype>::vertex_descriptor SC;

	SC = vertex(SC_vertex,g);

	for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
	{
		randnum = (double)rand() / (double)(RAND_MAX + 1);
		if( randnum < leave && *vi != SC )  //check if the random number is less than the probability and if the vertex is not the SC
		{
			//cout<<leave<<"  "<<randnum<<endl;
			//cout<<"DEBUG: removing vertex "<<*vi<<endl;

			clear_vertex(*vi,g);
			remove_vertex(*vi,g);  //might be a bug here!  does removing mess with vertex iterator
			//also need to update the SC_vertex!!!!!!!!
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
	typename graph_traits<graphtype>::vertex_descriptor u;

	//randnum = (double)rand() / (double)(RAND_MAX + 1);
	randnum = rand() % 100;
	cout<<"randnum: "<<randnum<<endl;
	if( randnum < (join*100) )  //use integer <join*100 instead 
	{
		u = add_vertex(g);
		g[u].name = "Random Guy";
		added = true;
		cout<<"DEBUG: added vertex "<<u<<endl;
	}

	return added;
}


template <class graphtype>
void MyGraph<graphtype>::initq()
{
	typename graph_traits<graphtype>::vertex_iterator vi, vi_end;
	int i =0;

	//for(int i=0; i< q.size();i++)
	//{
		q.clear();
	//}

	for (tie(vi, vi_end) = vertices(g),i=0; vi != vi_end; ++vi,++i) 
	{
		
		pair<int,double*> temp = pair<int,double*>(i,&g[*vi].P);
		q.push_back(temp);
	}
	q.sort(mysortt);
	
	return;
}

//return a vector of vertex indicies that the SC will request a connection based on the strategy
//Budget is an integer specifying how many requests the SC can send per timestep
//Strategy 1: random selection (use random number to find the requested nodes)
//Strategy 2: greedy selection (send request to the people who have the highest probability of accepting)
template <class graphtype>
vector<int> MyGraph<graphtype>::genRequestList(double pt, double po, double alpha, int budget,int strategy )
{
	vector<int> nodesRequested;
	int randvert = 0;
	typename graph_traits<graphtype>::vertex_descriptor randvertdescr;
	
	switch(strategy)
	{
	case 1:  //random strat
		{
			int max_itr = num_vertices(g) * 3;
			int count = 0;
			//cout<<"DEBUG: Using Strategy 1 to pick requested nodes"<<endl;
			while(nodesRequested.size() < (unsigned int)budget)
			{
				randvert = rand() % num_vertices(g);
				randvertdescr = vertex(randvert,g);
				if( randvert != SC_vertex  && g[randvertdescr].covered == false  ) 
				{
					nodesRequested.push_back(randvert);
					count = 0;  //reset the count of ineligible nodes
				}
				else
				{
					count++;	//increase the count of ineligible nodes
				}
				if(count > max_itr)  {
					//cout<<"exceeded max limit to find eligible node"<<endl;
					break;  //if we are looping for too long without finding an elegible node then break out of loop
				}
				else{
				}
			}
		}
		break;
	case 2:  //pick the nodes that have the highest trust value
		{
			//cout<<"DEBUG: Using Strategy 2 to pick requested nodes: greedy approach, pick nodes that have the highest chance to accept reuqest"<<endl;
		
			//this assumes that we have an up to date vector of probs for each node(stored in nodes itself
			//loop over the vector of probs and find the highest ones and request them
			
			int count = 0;
			//cout<<"DEBUG: Using Strategy 1 to pick requested nodes"<<endl;

			int i = 0;
			//vector< pair<int,double> > temp;
			
			q.sort(mysortt);
			
			pair<int,double*> highestPnode;
			
			while(nodesRequested.size() < (unsigned int)budget && !q.empty())
			{
				highestPnode = q.front();  //get the first element, should be the guy with highest prob
				q.pop_front();				//pop that guy from the list

				randvertdescr = vertex(highestPnode.first,g);
				if( highestPnode.first != SC_vertex  && g[randvertdescr].covered == false  ) 
				{
					nodesRequested.push_back(highestPnode.first);
					//cout<<"adding node "<<temp[j].first<<"to list cuz of prob"<<temp[j].second<<endl;
					
				}
				
				if(q.empty() == true)
				{
					break;
				}
			}

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

//find the number of common neighbors(adjacent verticies) between two nodes
//in real terms, the number of common 1st level connections
//this is decently slow
template <class graphtype>
int MyGraph<graphtype>::numCommonNeighbors(int u, int v)
{
	int count = 0;
	//typename graph_traits<graphtype>::vertex_descriptor U, V;
	//graph_traits<graphtype>::adjacency_iterator adjv, adjv_end;
	//graph_traits<graphtype>::adjacency_iterator adju, adju_end;

	//U = vertex(u,g);
	//V = vertex(v,g);

	//
	////cout<<"DEBUG: Common neighbors of "<<g[U].name<<" and "<<g[V].name<<endl;
	//for (tie(adju, adju_end) = adjacent_vertices(U,g); adju != adju_end; ++adju)
	//{
	//	for (tie(adjv, adjv_end) = adjacent_vertices(V,g); adjv != adjv_end; ++adjv) //for each vi, find the vertices adjacent to vi
	//	{
	//	
	//		if(*adju == *adjv)
	//		{
	//			//cout<<"DEBUG: "<<g[*adju].name<<" and  "<<g[*adjv].name<<endl;
	//			count++;
	//		}
	//	}
	//	

	//}
	cout<<"THIS IS AN OUTDATED FUNCTION"<<endl;
	exit(1);

	return count;
}



template <class graphtype>
int MyGraph<graphtype>::updateNumCommonNeighborsANDP(double pt, double po, double alpha, int u)
{
	typename graph_traits<graphtype>::vertex_descriptor U;
	typename graph_traits<graphtype>::adjacency_iterator adj, adj_end;

	U = vertex(u,g);

	for (tie(adj, adj_end) = adjacent_vertices(U,g); adj != adj_end; ++adj) //loop through the adjacent vertices and update their CNeigh val
	{
		//cout<<(int)*adj<<endl;
		//CNeigh[ (int)*adj ]++;  //old way
		//cout<<(int)*adj<<endl;
		g[*adj].CNeigh++;

		updateP(pt,po,alpha,*adj);

		
		
	}
	
	return 0;
}



template <class graphtype>
double MyGraph<graphtype>::calcProbabilities(double pt, double po, double alpha, int nodeDeg, int numCN)
{	
	double P = 0.0;
	double PE;
	double PT;
	int k = numCN;
	int d = nodeDeg;

	//calculate the prob that the node will accept due to trust
	PT = 1 - pow( (1 - pt), k );

	//calculate the prob that the node will accept due to ego
	PE = po / pow( (1+ d), alpha);

	//calculate the total probability that the node will accept the SC request
	P = 1 - ( (1 - PT)*(1 - PE) );
	
	//cout<<"DEBUG: Total probability is "<<P<<endl;

	return P;
}


//function that initializes everyone's P attribute to just their ego part since noone has linked with SC yet
template <class graphtype>
void MyGraph<graphtype>::initP(double pt, double po, double alpha)
{
	int k;
	int d;
	typename graph_traits<graphtype>::vertex_iterator vi, vi_end;
	//loop through vertices
	for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi) //loop through the adjacent vertices and update their CNeigh val
	{
		k = g[*vi].CNeigh;
		d = out_degree(*vi,g);
		g[*vi].P = calcProbabilities(pt,po,alpha,d,k);
	}
	

	return;
}

//updates the P of one node, called from updateNumCommonNeighborsANDP(..)
template <class graphtype>
void MyGraph<graphtype>::updateP(double pt, double po, double alpha, int u)
{
	int k,d;
	typename graph_traits<graphtype>::vertex_descriptor v;

	v = vertex(u,g);
	
	//cout<<"old P: "<<g[v].P<<endl;
	k = g[v].CNeigh;
	d = out_degree(v,g);
	g[v].P = calcProbabilities(pt,po,alpha,d,k);
	//cout<<"new P: "<<g[v].P<<endl<<endl;
	return;
}





//function to input a list of nodes to request and the probabilities associated with the model to find out which nodes accepted the connection
//request from the SC
//Arguments: 
//////nodesRequested -- vector of ints corresponding to node indecies to request
//////pt 
//////po
//////alpha
//Return: the number of nodes that accepted connection request
template <class graphtype>
int MyGraph<graphtype>::update(vector<int> nodesRequested, double pt, double po, double alpha)
{
	int total = 0; //total number of edges added (nodes infiltrated)
	//double PT;  //probability that the node will accept request due to trust
	//double PE;  //probability that the node will accept request due to ego
	double P;   //total probability that the node will accept the connection

	int nodeDeg;
	//int commonNeighbors;
	

	//for each requested vertices
	for(unsigned int i =0; i<nodesRequested.size(); i++)
	{
		//int numCN;
		double randnum;
		typename graph_traits<graphtype>::vertex_descriptor u, SC;

		//get the vertex desc of the requested node
		u = vertex(nodesRequested[i], g);
		//set the vertex's covered value to true since the node can't be requested again
		g[u].covered = true;
		
		SC = vertex(SC_vertex,g);

#ifdef DEBUG
		cout<<"Processing Request for node "<<nodesRequested[i]<<" who is "<<g[u].name<<endl;
		debugfile<<"Processing Request for node "<<nodesRequested[i]<<" who is "<<g[u].name<<endl;
#endif

		//find degree of node for ego
		//cout<<"out deg: "<<out_degree(nodesRequested[i],g)<<endl;
		//cout<<"in deg: "<<in_degree(nodesRequested[i],g)<<endl;
//		nodeDeg = out_degree( u ,g);
//#ifdef DEBUG
//		cout<<"Degree of node "<<nodesRequested[i]<<": "<<nodeDeg<<endl;
//#endif
//
//		//find neighbors between SC and node for trust
//		numCN = numCommonNeighbors(SC_vertex, nodesRequested[i]);  old way of finding neighbors
//		numCN = g[u].CNeigh;  //new way that is more efficient
//#ifdef DEBUG
//		cout<<"Number of neighbors in common with SC for node "<<nodesRequested[i]<<": "<<numCN<<endl;
//#endif
//
//		//calculate the PT and PE then P
//		P = calcProbabilities(pt, po, alpha, nodeDeg, numCN);
		P = g[u].P;
#ifdef DEBUG
		cout<<"Probability for node "<<nodesRequested[i]<<" to accept request is: "<<P<<endl;
		debugfile<<"Probability for node "<<nodesRequested[i]<<" to accept request is: "<<P<<endl;
#endif

		//randomly determine whether the request will be accepted
		randnum = (double)rand() / (double)(RAND_MAX + 1);
#ifdef DEBUG
		cout<<"DEBUG: randum num for adding edge: "<<randnum<<endl;
		debugfile<<"DEBUG: randum num for adding edge: "<<randnum<<endl;
#endif

		if(randnum < P)
		{
			//update the common neighbors vector before actually adding the edge due to the nature of the function 
			updateNumCommonNeighborsANDP(pt,po,alpha,nodesRequested[i]);

			//add function that updates the P attribute of nodes affected
			//updateP(pt, po, alpha, nodesRequested[i]);

			//if yes then add the edge and increment the total number of edges added
			add_edge(u , SC, g);
			total = total + 1;
			
#ifdef DEBUG
			cout<<"SC is now connected to "<<g[u].name<<" !!!"<<endl;
			debugfile<<"SC is now connected to "<<g[u].name<<" !!!"<<endl;
#endif
		}
		else
		{
			//if no then move on

		}
			
			
	}

	return total;
}


//calculate the trust value based on fraction of infected nodes (nodes connected to SC)
template <class graphtype>
double MyGraph<graphtype>::calcTrustVal()
{
	double TV = 0.0;
	int SC_deg = 0;
	int total_vertices = 0;

	typename graph_traits<graphtype>::vertex_descriptor SC;
	SC = vertex(SC_vertex,g);

	SC_deg = out_degree(SC,g);
	total_vertices = num_vertices(g);

	TV = (double)SC_deg / (double)total_vertices;
	
#ifdef DEBUG
	cout<<"TRUST VALUE (SC_deg/total vertices): "<<TV<<endl;
	cout<<"Trust value = "<<SC_deg<<"/"<<total_vertices<<endl;
#endif


	return TV;
}

//calculate just the # of connections that the SC node has
template <class graphtype>
int MyGraph<graphtype>::calcBenchVal()
{
	int TV = 0;
	int SC_deg = 0;

	typename graph_traits<graphtype>::vertex_descriptor SC;
	SC = vertex(SC_vertex,g);
	SC_deg = out_degree(SC,g);

	TV = SC_deg;
	

	return TV;
}


template <class graphtype>
double MyGraph<graphtype>::Infiltrate(double join, double leave, double pt,double po, double alpha, int budget, int strategy)
{
	vector<int> nodesRequested;
	double TrustValue;
	//initialize seed  -- do this in main
	//srand( (int)time(NULL) );

	//check each node if it leaves the network
	//nodeLeft(leave);  //don't call this yet since there are bugs

	//check if a node joins the network
	//nodeJoin(join);
	
	
	//find the nodes that the stealth company requested
	//clock_t start = clock();
	nodesRequested = genRequestList(pt,po,alpha,budget,strategy);
    //clock_t ends = clock();
	//cout <<endl<< "Running Time to genRequestList : "<< (double) (ends - start) / CLOCKS_PER_SEC << endl;
	//cout <<endl<< "Running Time to genRequestList : "<< (double) (ends - start)<< endl;
	
	//find the nodes that accepted the connection requests and add them to the network
	//start = clock();
	update(nodesRequested, pt, po, alpha);
	//ends = clock();
	//cout <<endl<< "Running Time to do the update : "<< (double) (ends - start) / CLOCKS_PER_SEC << endl;
	//cout <<endl<< "Running Time to do the update : "<< (double) (ends - start)<< endl;

	
	TrustValue = calcTrustVal();
	//calculate the trust value and that is the return value
	
#ifdef DEBUG
		cout<<endl<<endl;
#endif
	
	return TrustValue;

}

//do the infiltration and just return the # of connections of SC
template <class graphtype>
double MyGraph<graphtype>::Infiltrate_alt(double join, double leave, double pt,double po, double alpha, int budget, int strategy)
{
	vector<int> nodesRequested;
	int TrustValue;
	//initialize seed  -- do this in main
	//srand( (int)time(NULL) );

	//check each node if it leaves the network
	//nodeLeft(leave);  //don't call this yet since there are bugs

	//check if a node joins the network
	//nodeJoin(join);

	//find the nodes that the stealth company requested
	clock_t start = clock();
	nodesRequested = genRequestList(pt,po,alpha,budget,strategy);
    clock_t ends = clock();
	//cout <<endl<< "Running Time to genRequestList : "<< (double) (ends - start) / CLOCKS_PER_SEC << endl;
	cout <<endl<< "Running Time to genRequestList : "<< (double) (ends - start)<< endl;
	
	//find the nodes that accepted the connection requests and add them to the network
	start = clock();
	update(nodesRequested, pt, po, alpha);
	ends = clock();
	//cout <<endl<< "Running Time to genRequestList : "<< (double) (ends - start) / CLOCKS_PER_SEC << endl;
	cout <<endl<< "Running Time to genRequestList : "<< (double) (ends - start)<< endl;

	TrustValue = calcBenchVal();
	//calculate the trust value and that is the return value
	
#ifdef DEBUG
		cout<<endl<<endl;
#endif
	
	return TrustValue;

}



#endif
