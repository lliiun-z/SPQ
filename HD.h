#ifndef HD_H
#define HD_H

#include "Graph.h"

using namespace std;

//#define HD_DEBUG

// highway dimention;
class HD{
private:
	Graph& graph;
	
public:
	HD(Graph&);
	~HD(void);
	
	double run(VertexID, VertexID, vector<EdgeID>&, bool, vector<VertexID>&) const;
};


#endif
