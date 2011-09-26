#ifndef MYGRAPH_H
#define MYGRAPH_H

#include <iostream>
#include <vector>
#include <string>
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
	

public:
	graphtype g;

	MyGraph();										//create graph with DEFAULT_GRAPH_SIZE verticies and no edges
	MyGraph(int N, double p);						//create random erdos renyi graph
	MyGraph(vector<int> R);							//create graph with degree distribution R
	MyGraph(int N, int d, double p);				//create graph with mixed preferential attachment
	MyGraph(string datafile, int mode);				//create graph based on real datafile (mode tells what file type it is to read)

	float computeClusteringCoefficient();			//compute the clustering coefficient of the graph
	void computeDegreeDistribution();				//compute the degree distribution of the graph

	void printVE();									//print vertices and edges

	
};


/////////////////////////////////CONSTRUCTORS///////////////////////////////////////////////

template <class graphtype>
MyGraph<graphtype>::MyGraph()
{
	g = graphtype(DEFAULT_GRAPH_SIZE);
}


//Make an erdos renyi graph with N nodes and edges with probability p
template <class graphtype>
MyGraph<graphtype>::MyGraph(int N, double p)
{
	typedef boost::erdos_renyi_iterator<boost::minstd_rand, graphtype> ERGen;
	boost::minstd_rand gen;
	g = graphtype(ERGen(gen, N, p), ERGen(), N);
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
}


template <class graphtype>
MyGraph<graphtype>::MyGraph(string datafile, int mode)
{
	if(mode == 0)
	{
		ifstream LI_inp;
		//We should be opening BrianInMap.txt which is the sample network data
		LI_inp.open(datafile);

		LI


	}
	if(mode == 1)
	{
		//this will probably be for the real linkedIn data file if/when I get it
	}
	else
	{
		cout<<"You specified an unsupported mode....exiting"<<endl;
		exit();
	}
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

//This function computes the degree distribution of the graph g
template <class graphtype>
void MyGraph<graphtype>::computeDegreeDistribution()
{
	


	return;
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




#endif