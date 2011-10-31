#include "MyPriorityQueue.h"
#include <math.h>


#define DEFAULT_PQ_SIZE 8

double Log2( double n )  
{  
    return log( n ) / log( 2.0 );  
}


MyPriorityQueue::MyPriorityQueue()
{
	lastind = -1;
}

MyPriorityQueue::MyPriorityQueue( vector<struct MyNode*> nodeList )
{
	
	lastind = nodeList.size() - 1;

	for (int i =0; i<nodeList.size(); i++)
	{
		//add MyNode addresses that we got from graph to pq
		pq.push_back( nodeList[i] );
		nodeList[i]->PQind = i;
	}


	int balheapsize = pow(2, ceil(Log2( (double)(lastind+2) )) );
	pq.resize(balheapsize-1);
	int i = 0;
	for( i = (int)floor((double)pq.size()/2.0) - 1; i >= 0; i--)
	{
		
		heapify(i);
	}

}


void MyPriorityQueue::heapify(int i)
{
	
	int left = 2*i + 1;
	int right = 2*i + 2;
	int largest = i;
	if( left <= lastind && pq[left]->P > pq[i]->P )
	{
		largest = left;
	}
 
	if( right <= lastind && pq[right]->P > pq[largest]->P )
	{
		largest = right;
	}
 
	if( largest != i )
	{
		//swap pq[i] and pq[largest];
		struct MyNode * temp = pq[i];
		pq[i] = pq[largest];
		pq[i]->PQind = i;

		pq[largest] = temp;
		pq[largest]->PQind = largest;

	    heapify(largest);
	}
	
	return;
}


void MyPriorityQueue::percolateUp(int index)
{
	struct MyNode* node =  getElement(index);
	
	int parent = ceil( (double)((index - 1) / 2) );
	int left = 2*index + 1;
	int right = 2*index + 2;

	if( pq[index]->P > pq[parent]->P)
	{
		struct MyNode * temp = pq[index];
		pq[index] = pq[parent];
		pq[index]->PQind = index;
		pq[parent] = temp;
		pq[parent]->PQind = parent;

		percolateUp(parent);
	}
	
	return;
}


int MyPriorityQueue::relevantSize()
{
	return lastind+1;
}

int MyPriorityQueue::realSize()
{
	return pq.size();
}

struct MyNode* MyPriorityQueue::getElement(int i)
{
	return pq[i];
}

void MyPriorityQueue::setP(int index, double P)
{
	struct MyNode* node = getElement(index);
	node->P = P;
	percolateUp(index);
	return;
}

struct MyNode* MyPriorityQueue::pop()
{
	//get the value of top element
	struct MyNode* temp = pq[0];

	//re-heapify
	//cout<<"putting this on top"<<lastind<<endl;
	pq[0] = pq[lastind];
	pq[0]->PQind = 0;
	heapify(0);
	//setP(lastind,-1);
	lastind--;

	return temp;   //returning the root of the heap
}

void MyPriorityQueue::pqClear()
{
	pq.clear();
	return;
}

bool MyPriorityQueue::pqEmpty()
{
	if( pq.empty() )
	{
		return true;
	}
	else
	{
		return false;
	}

}

MyPriorityQueue::~MyPriorityQueue()
{
	pqClear();
	lastind = -1;
}