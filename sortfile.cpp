#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

struct CompareNode                                                                                  
{
  bool operator()(pair<int,double*> pa, pair<int,double*> pb) const
  {
    return *(pa.second) > *(pb.second);
  }
};

bool mysortt(pair<int,double*> pa, pair<int,double*> pb)
{
	return *(pa.second) > *(pb.second);
}