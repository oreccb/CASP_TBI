#include "MyPriorityQueue.h"
#include <math.h>

#define DEFAULT_PQ_SIZE 8

double Log2( double n )  
{  
    return log( n ) / log( 2.0 );  
}


MyPriorityQueue::MyPriorityQueue()
{
	pq.resize(DEFAULT_PQ_SIZE);
	lastind = 0;
}

MyPriorityQueue::MyPriorityQueue(vector< pair<int,double*> > thing)
{
	lastind = thing.size() - 1;

	pq = thing;
	int balheapsize = pow(2, ceil(Log2( (double)(thing.size()+1) )) );
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
	if( left <= lastind && *pq[left].second > *pq[i].second )
	{
		largest = left;
	}
 
	if( right <= lastind && *pq[right].second > *pq[largest].second )
	{
		largest = right;
	}
 
	if( largest != i )
	{
		//swap pq[i] and pq[largest];
		pair<int,double*> temp = pq[i];
		pq[i] = pq[largest];
		pq[largest] = temp;

	    heapify(largest);
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

pair<int,double*> MyPriorityQueue::getElement(int i)
{
	return pq[i];
}

pair<int, double*> MyPriorityQueue::pop()
{
	//get the value of top element
	pair<int, double*> temp = pq[0];

	//re-heapify
	pq[0] = pq[lastind];
	heapify(0);
	lastind--;

	return temp;   //returning the root of the heap
}
