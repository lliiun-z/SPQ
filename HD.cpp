#include "HD.h"

HD::HD(Graph& g):graph(g){}

HD::~HD(){}

double HD::run(VertexID source, VertexID target, vector<EdgeID>& path, bool unpack, vector<VertexID>& id_map) const{
	const LabelList& outlabel = graph.outLabelList[source];
	const LabelList& inlabel = graph.outLabelList[target];
	
	double dist = DBL_INFINITY;
	int meet_node = -1;
	int i = 0, j = 0;
	while( i < outlabel.size() && j < inlabel.size() ) {
		if( outlabel[i].id == inlabel[j].id ) {
			if ( outlabel[i].distance + inlabel[j].distance < dist ) {
				dist = outlabel[i].distance + inlabel[j].distance;
				meet_node = inlabel[j].id;
			}
			i++;
			j++;
		}else if (outlabel[i].id == target) {
			if ( outlabel[i].distance < dist ) {
				dist = outlabel[i].distance;
				meet_node = outlabel[i].id;
			}
			i++;
		}else if (inlabel[j].id == source) {
			if ( inlabel[j].distance < dist ) {
				dist = inlabel[j].distance;
				meet_node = inlabel[j].id;
			}
			j++;
		}else if ( outlabel[i].id < inlabel[j].id ) i++;
		else j++;
	}
	
	#ifdef HD_DEBUG
		cout << source << "->outlabellist: ";
		for ( int i = 0; i < outlabel.size(); i++ ) cout << "[" << outlabel[i].id << ", " << outlabel[i].distance << "] ";
		cout << endl;
		cout << target << "->inlabellist: ";	
		for ( int i = 0; i < outlabel.size(); i++ ) cout << "[" << inlabel[i].id << ", " << inlabel[i].distance << "] ";
		cout << endl;	
		cout << endl << "meet at node " << meet_node << endl;
	#endif
	
	// print path;
	int index_source = -1, index_target = -1;
	if (unpack) {
		for ( int i = 0; i < outlabel.size(); i++ ) {
			if ( outlabel[i].id == meet_node ) {
				index_source = i;
				break;
			}
		}
		
		for ( int i = 0; i < inlabel.size(); i++ ) {
			if ( inlabel[i].id == meet_node ) {
				index_target = i;
				break;
			}
		}
	}
	
	if (unpack && id_map.size() == 0) {
		cout << "Path from " << source << " to " << target << ": ";
		//if (meet_node != source) cout << source << "->";
		//cout << "index source = " << index_source << ", index_target = " << index_target << endl;
		//cout << graph.pathList[source].size() << ", " << graph.pathList[target].size() << endl;
		for (int begin = graph.pathIndex[source][index_source*2]; begin < graph.pathIndex[source][index_source*2+1]; begin++) {
			cout << "->" << graph.pathList[source][begin];
		}
		for (int begin = graph.pathIndex[target][index_target*2+1]-2; begin >= graph.pathIndex[target][index_target*2]; begin--) {
			cout << "->" << graph.pathList[target][begin] ;
		}
		cout << endl;
	}else if (unpack && id_map.size() != 0) {
		cout << "Path from " << id_map[source] << " to " << id_map[target] << ": ";
		for (int begin = graph.pathIndex[source][index_source*2]; begin < graph.pathIndex[source][index_source*2+1]; begin++) {
			cout << "->" << id_map[graph.pathList[source][begin]];
		}
		for (int begin = graph.pathIndex[target][index_target*2+1]-2; begin >= graph.pathIndex[target][index_target*2]; begin--) {
			cout << "->" << id_map[graph.pathList[target][begin]];
		}
		cout << endl;		
	}
	
	return dist;
}
