#ifndef MYPRIORITYQUEUE_H
#define MYPRIORITYQUEUE_H


#include <iostream>
#include <vector>
//#include "MyGraph.h"
#include "Bundled_Proporties.h"

using namespace std;



//typedef adjacency_list<setS,vecS,undirectedS, MyNode> undirGraph;

class MyPriorityQueue
{
private:
	vector< struct MyNode * > pq;
	
	int lastind;  //index of the last real node in the queue (actual queue could be more cuz it needs to be balanced)


public:
	MyPriorityQueue();
	MyPriorityQueue( vector<struct MyNode*> nodeList );
	
	struct MyNode * pop();   //pops the first element off the top of the priority queue

	//void add(pair<int,double*> element);  //adds an element to the priority queue

	void heapify(int i);    //heapify the element at index i

	//void buildHeap(vector< pair<int,double*> >);

	void percolateUp(int index);    //does the percolate up function starting at the node that is at index, index
	//this will be used on each slot that is updated

	int relevantSize();  //# of nodes that matter
	int realSize();      //actual # of nodes in the pq vector obtained by pq.size() (could be a little higher than relevantSize()
	void setP(int index, double P); //sets the double value (second thing in pair) of the pair that the ptr is pointing at


	struct MyNode * getElement(int i);

	void pqClear();
	bool pqEmpty();

	~MyPriorityQueue();
};



#endif