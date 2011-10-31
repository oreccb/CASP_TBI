#ifndef BUNDLED_PROPORTIES_H
#define BUNDLED_PROPORTIES_H
//This file contains the bundled proporties for the graphs

//Generic Graph Class///////////////////////////////////
struct MyNode {
	int x;
	std::string name;
	bool covered;
	int CNeigh;
	double P;
	int PQind;
};

struct MyEdge {
    int weight;
};
/////////////////////////////////////////////////////////

#endif
