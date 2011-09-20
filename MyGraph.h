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

using namespace std;
using namespace boost;

template <class graphtype>
class MyGraph
{
private:
	graphtype g;

public:
	MyGraph(int N, double p);						//create random erdos renyi graph
	MyGraph(vector<int> R);							//create graph with degree distribution R
	MyGraph(int N, int d, double p)					//create graph with mixed preferential attachment

	float computeClusteringCoefficient();			//compute the clustering coefficient of the graph
	void computeDegreeDistribution();				//compute the degree distribution of the graph


};


/////////////////////////////////CONSTRUCTORS///////////////////////////////////////////////

//Make an erdos renyi graph with N nodes and edges with probability p
template <class graphtype>
MyGraph<graphtype>::MyGraph(int N, double p)
{
	boost::minstd_rand gen;
	g = graphtype(ERGen(gen, N, p), ERGen(), N);
}

//We are making a random graph with specified deg distribution R
template <class graphtype>
MyGraph<graphtype>::MyGraph(vector<int> R)
{
	int i,j;
	graph_traits<MyGraph>::vertex_iterator vi, vi_end;
    graph_traits<MyGraph>::vertex_iterator vi2, vi2_end;
    graph_traits<MyGraph>::edge_iterator ei, ei_end;

	typename graph_traits<MyGraph>::edge_descriptor ed;
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
	for(i=0;i<ordered_deg_list.size();i++)
	{
		cout<<"("<<ordered_deg_list[i].first<<",";
		cout<<ordered_deg_list[i].second<<") ";
	}
	cout<<endl;

	//sets the vertex id
	//for (tie(vi, vi_end) = vertices(g), i=0 ; vi != vi_end; ++i, ++vi)
	//{
	//	g[*vi].x = i;
	//}

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
	for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
	{
		std::cout<<g[*vi].x <<" ";

	}
	std::cout<<std::endl;
	for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
        cout << "edge: " << source(*ei, g) << " to " << target(*ei, g) << endl;
    }
    cout<<endl;
	

	
	return(0);
}

//this constructor creates a graph with N nodes based on the mixed preferential attachment model for growth
//At each step, on node is added and the node makes d new edges
//each new edge will attach with probability p of being random and probability (1-p) of being preferentially attached
template <class graphtype>
MyGraph<graphtype>::MyGraph(int N, int d, double p)
{
	


}

/////////////////////////////////END CONSTRUCTORS//////////////////////////////////////////////////////////

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



#endif
