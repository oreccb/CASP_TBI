#ifndef MYPRIORITYQUEUE_H
#define MYPRIORITYQUEUE_H

#include <iostream>
#include <vector>

using namespace std;

class MyPriorityQueue
{
private:
	vector< pair<int,double*> > pq;
	
	int lastind;  //index of the last real node in the queue (actual queue could be more cuz it needs to be balanced)


public:
	MyPriorityQueue();
	MyPriorityQueue( vector< pair<int,double*> > thing);
	
	pair<int,double*> pop();   //pops the first element off the top of the priority queue

	//void add(pair<int,double*> element);  //adds an element to the priority queue

	//void percDown(pair<int,double*> * eleptr);   //the element in pq pointed to by eleptr has changed to the pq needs to be updated by percolating down

	void heapify(int i);    //heapify the element at index i

	//void buildHeap(vector< pair<int,double*> >);


	int MyPriorityQueue::relevantSize();  //# of nodes that matter
	int MyPriorityQueue::realSize();      //actual # of nodes in the pq vector obtained by pq.size() (could be a little higher than relevantSize()

	pair<int,double*> getElement(int i);

};



#endif