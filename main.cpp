#include "Contraction.h"
#include "Query.h"
#include "Dijkstra.h"
#include "PerformanceTimer.h"
#include "GraphUtil.h"
#include "HD.h"

#include <iostream>
#include <time.h>

using namespace std;

//#define TEST
//#define min(a, b) (a < b ? a : b)

void printUsage() {
	cerr << "Usage of CH:" << endl;
	cerr << "  coresearchd --arg1=val1 --arg2=val2" << endl << endl << endl;

	cerr << "Command line parameters:" << endl << endl;

	cerr << "--help" << endl;
	cerr << "--usage" << endl;
	cerr << "  prints this usage information and exits." << endl << endl;
	
	// cerr << "--debug" << endl;
	// cerr << "  Validate the results using BFS as benchmark." << endl << endl;
	
	// cerr << "--experiment" << endl;
	// cerr << "  Compare the performance of all methods." << endl << endl;

	cerr << "--query_num=1000" << endl;
	cerr << "  specifies the number of random shortest path queries to be 1000." << endl << endl;

	// cerr << "--resultfile=result.txt" << endl;
	// cerr << "  Specifies the file to output result." << endl << endl;
	
	cerr << "--filename=data.gra" << endl;
	cerr << "  use data.gra as input file." << endl << endl;
	
	cerr << "--unpack=true[default:false]" << endl;
	cerr << "  print out the path." << endl << endl;	
	
	cerr << "--experiment=1[default:0]" << endl;
	cerr << "  choose the experiment: 1: correctness; 2: performance; 0: both." << endl;
	cerr << "  this paramter is for debugging only." << endl << endl;		
	
	cerr << "--save_para=falsedefault:true]" << endl;
	cerr << "  true: use the parameters of last experiment if they are saved before." << endl; 
	cerr << "  false: generate new parameters for this experiment." << endl << endl;

	cerr << "--method=all[default:dijk]" << endl;
	cerr << "  dijk: only run the dijkstra method." << endl; 
	cerr << "  label: only run labeling method." << endl;
	cerr << "  all: run two methods." << endl;
	cerr << "  this parameter is for debugging only. " << endl << endl;
	
	cerr << "--idmap=yes[default:no]" << endl;
	cerr << "  yes: will read the .idmap file from the same directory of graph file." << endl; 
	cerr << "  no: otherwise." << endl << endl;
	
	cerr << "--path=yes[default:no]" << endl;
	cerr << "  yes: path information will be generated when produc labels." << endl; 
	cerr << "  no: otherwise." << endl << endl;	
}

void printParameters(const map<string, string>& cmdLineArgs) {
	map<string, string>::const_iterator iter;
	for (iter = cmdLineArgs.begin(); iter != cmdLineArgs.end(); iter++) {
		cout << iter->first << "=";
		if (iter->second=="")
			cout << "EMPTY";
		else 
			cout << iter->second;
		cout << " ";
	}
	cout << endl;
}

static void keepResult(const char* resultFileName, vector<string>& results) {
	ofstream out(resultFileName, ios_base::out|ios_base::app);
	for (int i = 0; i < results.size(); i++)
		out << results[i] << "\t";
	out << endl;
	out.close();
}

static void readFile(ifstream& in, vector<vector<double> >& data) {
	vector<vector<double> >(0).swap(data);
	string str;
	char* line;
	int count = 0;
	while (!in.eof()) {
		getline(in, str);
		if (str.size() == 0) break;
		line = new char[str.size()+1];
		strcpy(line, str.c_str());
		char* num = strtok(line, " \t#:");
		vector<double> tv;
		while(num != NULL) {
			tv.push_back(atof(num));
			num = strtok(NULL, " \t#:");
		}
		// vector<VertexID> tu;
		// for( int i = 0; i < tv.size(); i++ ) tu.push_back(tv[i]);
		data.push_back(tv);
		
		delete line;
	}
	in.close();
}

// check whether the queries have been generated or not;
// if yes, read the queries; otherwise, create new ones.
static void checkQuery(Graph& graph, string file_name, int num_query, vector<VertexID>& query, bool save_para){
	string query_file = file_name;
	size_t found = query_file.find_last_of('.');
	query_file.erase(query_file.begin()+found, query_file.end());
	
	char buff[20];
	sprintf(buff, "%d", num_query);
	query_file += "_";
	query_file += buff;
	query_file += ".query";
	
	ifstream infile(query_file.c_str());
	if(infile && save_para){
		vector<vector<double> > data;
		readFile(infile, data);
		for ( int i = 0; i < data.size(); i++ ) {
			for (int j = 0; j < min(2, data[i].size()); j++) {
				query.push_back(static_cast<VertexID>(data[i][j]));
			}	
		}
	}else{
		cout << "Generate new query file: " + query_file << ".";
		GraphUtil gu;
		ofstream out(query_file.c_str());
		gu.generateQuery(graph, graph.num_vertices(), num_query, out, query);
		cout << " Done!" << endl;
	}
}

static void readHierarchies(Graph& graph, string hierarchy_file){
	vector<ShortCuts> outshortcut, inshortcut;
	VertexList nodelist;
	int nodesize = graph.num_vertices();
	
	// outshortcutlist;
	vector<vector<double> > data;
	// cout << "file name:" << hierarchy_file << endl;
	size_t found = hierarchy_file.find_last_of(".");
	hierarchy_file.erase(hierarchy_file.begin()+found, hierarchy_file.end());
	hierarchy_file += ".hierarchy_out";
	ifstream in(hierarchy_file.c_str());
	//readFile(in, data);
	
	// for ( int i = 0; i < data.size(); i++ ) {
		// for ( int j = 0;  j < data[i].size(); j++ ) {
			// cout << data[i][j] << " ";
		// }
		// cout << endl;
	// }
	
	// exit(0);
	
	/*outshortcut.resize(nodesize);
	for ( int i = 0; i < data.size(); i++ ){
		cout << i << endl;
		int u = static_cast<int>(data[i][0]);
		int v = static_cast<int>(data[i][1]);
		if ( u > nodesize ) continue;
		// cout << i << " " << u << " " << v << endl;		
		ShortCutEdge tmp_sce;
		tmp_sce.target = v;
		tmp_sce.weight = data[i][2];
		for (int j = 3; j < data[i].size(); j++) tmp_sce.innerIDs.push_back(static_cast<int>(data[i][j]));

		outshortcut[u].push_back(tmp_sce);
	}*/
	
	// inshortcutlist;
	/*found = hierarchy_file.find_last_of("_");
	hierarchy_file.erase(hierarchy_file.begin()+found, hierarchy_file.end());
	hierarchy_file += "_in";
	in.open(hierarchy_file.c_str(), ifstream::in);
	readFile(in, data);*/
	
	
	// for ( int i = 0; i < data.size(); i++ ) {
		// for ( int j = 0;  j < data[i].size(); j++ ) {
			// cout << data[i][j] << " ";
		// }
		// cout << endl;
	// }	
	
	/*inshortcut.resize(nodesize);
	for ( int i = 0; i < data.size(); i++ ){
		int u = static_cast<int>(data[i][0]);
		int v = static_cast<int>(data[i][1]);
		if ( u > nodesize ) continue;
		//cout << i << " " << u << " " << v << endl;			
		ShortCutEdge tmp_sce;
		tmp_sce.target = v;
		tmp_sce.weight = data[i][2];
		for ( int j = 3; j < data[i].size(); j++ ) tmp_sce.innerIDs.push_back(static_cast<int>(data[i][j]));

		inshortcut[u].push_back(tmp_sce);
	}
	
	graph.insertShortcut(inshortcut, outshortcut);*/
	

	// nodelist;
	found = hierarchy_file.find_last_of("_");
	hierarchy_file.erase(hierarchy_file.begin()+found, hierarchy_file.end());
	hierarchy_file += "_node";
	cout << "file name:" << hierarchy_file << endl;	
	in.open(hierarchy_file.c_str(), ifstream::in);
	readFile(in, data);
	
	// for ( int i = 0; i < data.size(); i++ ) {
		// for ( int j = 0;  j < data[i].size(); j++ ) {
			// cout << data[i][j] << " ";
		// }
		// cout << endl;
	// }	
	
	nodelist.resize(nodesize);
	for ( int i = 0; i < data.size(); i++ ) {
		int u = static_cast<int>(data[i][0]);
		int v = static_cast<int>(data[i][1]);
		if (u > nodesize) continue;
		nodelist[u].rank = v;
		nodelist[u].id = u;
	}
	
	graph.insertNodeList(nodelist);
	// for ( int i = 0; i < graph.vertices().size(); i++) {
		// cout << i << " " << graph.vertices()[i].rank << endl;
	// }

	// read labels;
	// nodelist;
	found = hierarchy_file.find_last_of(".");
	hierarchy_file.erase(hierarchy_file.begin()+found, hierarchy_file.end());
	hierarchy_file += ".label_in";
	cout << "file name:" << hierarchy_file << endl;	
	in.open(hierarchy_file.c_str(), ifstream::in);
	readFile(in, data);	
	graph.insertInLabel(data);
	
	found = hierarchy_file.find_last_of(".");
	hierarchy_file.erase(hierarchy_file.begin()+found, hierarchy_file.end());
	hierarchy_file += ".label_out";
	cout << "file name:" << hierarchy_file << endl;	
	in.open(hierarchy_file.c_str(), ifstream::in);
	readFile(in, data);	
	graph.insertOutLabel(data);	
	
	// read path index and list;
	found = hierarchy_file.find_last_of(".");
	hierarchy_file.erase(hierarchy_file.begin()+found, hierarchy_file.end());
	hierarchy_file += ".path_index";
	cout << "file name:" << hierarchy_file << endl;	
	in.open(hierarchy_file.c_str(), ifstream::in);
	readFile(in, data);	
	graph.insertPathIndex(data);
	
	found = hierarchy_file.find_last_of(".");
	hierarchy_file.erase(hierarchy_file.begin()+found, hierarchy_file.end());
	hierarchy_file += ".path";
	cout << "file name:" << hierarchy_file << endl;	
	in.open(hierarchy_file.c_str(), ifstream::in);
	readFile(in, data);	
	graph.insertPathList(data);		
}

static void checkWritingCorrectness(Graph& graph, string hierarchy_file) {
	Graph tmp_g(graph.num_vertices());
	cout << "reading" << endl;
	readHierarchies(tmp_g, hierarchy_file);
	cout << "done" << endl;		
			
	// compare outshortcutlist;
	const ShortCuts& outshortcutOrig = graph.exportoutShortCutList();
	const ShortCuts& outshortcutDeri = tmp_g.exportoutShortCutList();
	if (outshortcutOrig.size() != outshortcutDeri.size()) { cout << "Size does not match!" << endl; exit(0); }
	for (int i = 0; i < outshortcutOrig.size(); i++) {
		// check id;
		if (outshortcutOrig[i].target != outshortcutDeri[i].target) {
			cout << "Node id does not match!" << endl;
			exit(0);
		}
		
		// check weight;
		if (fabs(outshortcutOrig[i].weight-outshortcutDeri[i].weight) > 1e-5) {
			cout << "Weight does not match!" << endl;
			exit(0);
		}
		
		// check inner ids;
		if (outshortcutOrig[i].innerIDs.size() != outshortcutDeri[i].innerIDs.size()) {
			cout << "InnerIDs.size does not match!" << endl;
			exit(0);
		}
		for (int j = 0; j < outshortcutOrig[i].innerIDs.size(); j++) {
			if (outshortcutOrig[i].innerIDs[j] != outshortcutDeri[i].innerIDs[j]){
				cout << "InnerIDs content does not match!" << endl;
				exit(0);
			}
		}
	}

	// compare shortcutindex;
	const vector<int>& outshortcutindexOrig = graph.exportoutShortCutIndex();
	const vector<int>& outshortcutindexDeri = tmp_g.exportoutShortCutIndex();
	if (outshortcutindexOrig.size() != outshortcutindexDeri.size()){
		cout << "index size does not match!" << endl;
		exit(0);
	}
	for (int i = 0; i < outshortcutindexOrig.size(); i++) {
		if (outshortcutindexOrig[i] != outshortcutindexDeri[i]) {
			cout << "index content does not match!" << endl;
			exit(0);
		}
	}
	
	
	// compare nodelist;
	int nodesize_o = graph.num_vertices();
	int nodesize_d = tmp_g.num_vertices();
	if (nodesize_o != nodesize_d) {
		cout << "node size does not match!" << endl;
		exit(0);
	}
	for (int i = 0; i < nodesize_o; i++) {
		if (graph[i].rank != tmp_g[i].rank){
			cout << "node rank does not match!" << endl;
			exit(0);
		}
		//cout << i << ": " << graph[i].id << "--" << tmp_g[i].id << endl;
		if (graph[i].id != tmp_g[i].id){
			cout << "node id does not match!" << endl;
			exit(0);
		}
	}		
}

// check whether hierachies have been built or not;
// if yes, read them in; otherwise, build them;
static void checkHierarchies(Graph& graph, string file_name, bool save_para, bool gp){
	string hierarchy_file = file_name;
	size_t found = hierarchy_file.find_last_of('.');
	hierarchy_file.erase(hierarchy_file.begin()+found, hierarchy_file.end());
	hierarchy_file += ".hierarchy_out";

	ifstream infile(hierarchy_file.c_str());
	if (infile && save_para){
		cout << "read ... ";
		readHierarchies(graph, hierarchy_file);
		cout << "done ... ";
	}else{
		cout << "Create new hierarchy file: " << hierarchy_file << ".";
		Contraction c(graph);
		cout << "Creating hierarchy ... ";
		c.computeShortcuts(gp);		
		
		//const vector<ShortCuts>& outshortcut = c.exportOutShortcut();
		//const vector<ShortCuts>& inshortcut  = c.exportInShortcut();
		const VertexList& nodelist = c.exportNodeList();
		
		// write them down;
		// cout << "file name:" << hierarchy_file << endl;
		ofstream out(hierarchy_file.c_str());		
		// for ( int i = 0; i < outshortcut.size(); i++ ) {
			// for ( int j = 0; j < outshortcut[i].size(); j++ ) {
				// out << i << "\t" << outshortcut[i][j].target << "\t" << outshortcut[i][j].weight << ": ";
				// for ( int k = 0; k < outshortcut[i][j].innerIDs.size(); k++ ) {
					// out << outshortcut[i][j].innerIDs[k] << " ";
				// }
				// out << endl;
			// }
		// }
 		out.close();
		
		// found = hierarchy_file.find_last_of("_");
		// hierarchy_file.erase(hierarchy_file.begin()+found, hierarchy_file.end());
		// hierarchy_file += "_in";
		out.open(hierarchy_file.c_str(), ofstream::out);		
		// for ( int i = 0; i < inshortcut.size(); i++ ) {
			// for ( int j = 0; j < inshortcut[i].size(); j++ ) {
				// out << i << "\t" << inshortcut[i][j].target << "\t" << outshortcut[i][j].weight << ": ";
				// for ( int k = 0; k < inshortcut[i][j].innerIDs.size(); k++ ) {
					// out << inshortcut[i][j].innerIDs[k] << " ";
				// }
				// out << endl;
			// }
		// }
		out.close();

		found = hierarchy_file.find_last_of("_");
		hierarchy_file.erase(hierarchy_file.begin()+found, hierarchy_file.end());
		hierarchy_file += "_node";
		// cout << "file name:" << hierarchy_file << endl;		
		out.open(hierarchy_file.c_str(), ofstream::out);	
		for ( int i = 0; i < nodelist.size(); i++ ) {
			out << i << ": " << nodelist[i].rank << endl;
		}
		out.close();
		// cout << "!!" << endl;
		
		found = hierarchy_file.find_last_of(".");
		hierarchy_file.erase(hierarchy_file.begin()+found, hierarchy_file.end());
		hierarchy_file += ".label_in";
		// cout << "file name:" << hierarchy_file << endl;		
		out.open(hierarchy_file.c_str(), ofstream::out);
		for ( int i = 0; i < nodelist.size(); i++ ) {
			LabelList& label = graph.exportInLabel(i);
			out << i << ": ";
			for ( int j = 0; j < label.size(); j++ ) {
				out << label[j].id << " " << label[j].distance << " ";
			}
			out << endl;
		}
		out.close();

		found = hierarchy_file.find_last_of(".");
		hierarchy_file.erase(hierarchy_file.begin()+found, hierarchy_file.end());
		hierarchy_file += ".label_out";
		// cout << "file name:" << hierarchy_file << endl;		
		out.open(hierarchy_file.c_str(), ofstream::out);
		for ( int i = 0; i < nodelist.size(); i++ ) {
			LabelList& label = graph.exportOutLabel(i);
			out << i << ": ";
			for ( int j = 0; j < label.size(); j++ ) {
				out << label[j].id << " " << label[j].distance << " ";
			}
			out << endl;
		}
		out.close();	

		found = hierarchy_file.find_last_of(".");
		hierarchy_file.erase(hierarchy_file.begin()+found, hierarchy_file.end());
		hierarchy_file += ".path_index";
		// cout << "file name:" << hierarchy_file << endl;		
		out.open(hierarchy_file.c_str(), ofstream::out);
		vector<vector<int> >& pathIndex = graph.exportPathIndex();
		for ( int i = 0; i < pathIndex.size(); i++ ) {
			out << i << ": ";
			for ( int j = 0; j < pathIndex[i].size(); j++ ) {
				out << pathIndex[i][j] << " ";
			}
			out << endl;
		}
		out.close();

		found = hierarchy_file.find_last_of(".");
		hierarchy_file.erase(hierarchy_file.begin()+found, hierarchy_file.end());
		hierarchy_file += ".path";
		// cout << "file name:" << hierarchy_file << endl;		
		out.open(hierarchy_file.c_str(), ofstream::out);
		vector<vector<VertexID> >& pathList = graph.exportPathList(); 
		for ( int i = 0; i < pathList.size(); i++ ) {
			out << i << ": ";
			for ( int j = 0; j < pathList[i].size(); j++ ) {
				out << pathList[i][j] << " ";
			}
			out << endl;
		}
		out.close();	

		// for testing
		// checkWritingCorrectness(graph, hierarchy_file);
	}
}

void readIDMap(string filename, vector<VertexID>& id_map){
	size_t found = filename.find_last_of(".");
	filename.erase(filename.begin()+found, filename.end());
	filename += ".idmap";
	
	ifstream in(filename.c_str());
	vector<vector<double> > data;
	readFile(in, data);
	
	vector<VertexID>(0).swap(id_map);
	for (int i = 0; i < data.size(); i++) {
		id_map.push_back(static_cast<VertexID>(data[i][1]));
	}
}

bool isNumeric(string input){
	if (input.size() == 0) return false;
	for (int i = 0; i < input.size(); i++) {
		if (input[i] < '0' || input[i] > '9') return false;
	}
	return true;
}

PerformanceTimer timer = PerformanceTimer::start();
double run_time = 0;

int main(int argc, char* argv[]) {
	//cout.setf(ios::fixed,ios::floatfield);
	cout.precision(4);

	// increase stack size;
	const rlim_t sSize = 1024 * 1024 * 512;
	struct rlimit rl;

	if (getrlimit(RLIMIT_STACK, &rl) == 0) {
		if (sSize > rl.rlim_cur){
			rl.rlim_cur = sSize;
			if (setrlimit(RLIMIT_STACK, &rl) != 0) {
				cerr << "Warning: could not increase stack size. Program might run into trouble..." << endl;
			}
		}
	}	
	srand( time(NULL) );

	// parse paramters;
	typedef map<string,string> mapType;
	map<string, string> cmdLineArgs;
	cmdLineArgs["filename"] = "";
	cmdLineArgs["query_num"] = "1000";
	cmdLineArgs["unpack"] = "true";
	cmdLineArgs["experiment"] = "0";
	cmdLineArgs["save_para"] = "false";
	cmdLineArgs["method"] = "dijk";
	cmdLineArgs["help"] = "unset";
	cmdLineArgs["usage"] = "unset";
	cmdLineArgs["debug"] = "unset";
	cmdLineArgs["idmap"] = "no";
	cmdLineArgs["path"] = "no";

	for (int i = 1; i < argc; i++) {
		string currArg = argv[i];
		size_t found = string::npos;
		
		for(mapType::const_iterator it = cmdLineArgs.begin(); it != cmdLineArgs.end(); ++it){
			string fullArg = "--" + it->first;
			found = currArg.find(fullArg);
			if (found != string::npos) {
				if (currArg.find(fullArg + "=") != string::npos) {
					cmdLineArgs[it->first] = currArg.substr(fullArg.size() + 1);
				} else {
					cmdLineArgs[it->first] = "";
				}
				break;
			}
		}

		if (found == string::npos) {
			cerr << "Invalid commandline argument: " << currArg << endl;
			printUsage();
			exit(1);
		}
	}

	if (cmdLineArgs["filename"] == "") {
		cerr << "No input file provided, please use command-line option --filename=data.gra" << endl;
		printUsage();
		exit(0);
	}
	if (cmdLineArgs["usage"] != "unset" || cmdLineArgs["help"] != "unset") {
		printUsage();
		exit(0);
	}	
	
	string filename = cmdLineArgs["filename"];
	int query_num = atoi(cmdLineArgs["query_num"].c_str());
	bool unpack = true;
	int experiment = atoi(cmdLineArgs["experiment"].c_str());
	bool save_para = false;
	if (cmdLineArgs["save_para"].compare("true") == 0) save_para = true;
	int method = 0;
	if (cmdLineArgs["method"].compare("label") == 0) method = 1;
	if (cmdLineArgs["method"].compare("all") == 0) method = 2;
	if (method == 0) save_para = true; // use dijkstra, we do not need proprocessing;
	bool gp = false;
	if (cmdLineArgs["path"].compare("yes") == 0) gp = true; // label includs path info; 
	if (!gp) unpack = false;
	
	ifstream infile(filename.c_str());
	if (!infile) {
		cout << "Error: Cannot open " << filename << endl;
		exit(0);
	}
	
	vector<VertexID> id_map;
	if (cmdLineArgs["idmap"].compare("yes") == 0) { // if has id map, we will map the id back to its original id;
		readIDMap(filename, id_map);
	}
	
	string query_file = filename;
	
	// main part;
	cout << "Building graph ... ";
	Graph g(infile);
	// int gsize = g.num_vertices();
	// g.printGraph();
	cout << "Done." << endl;
	
	// generate queries;
	cout << "Generating queries ... ";
	timer.reset();
	vector<VertexID> query;
	checkQuery(g, filename, query_num, query, save_para);	
	run_time = timer.reset();
	cout << "Done. It takes " << run_time << " ms." << endl;
	// print queries;
	// cout << "Print queries:" << endl;
	// for ( int i = 0; i < query.size(); i+=2 ) {
		// cout << query[i] << " " << query[i+1] << endl;
	// }
	// cout << endl;	
	
	// build contraction hierarchies;
	cout << "Generating hierarchies ... ";
	timer.reset();
	checkHierarchies(g, filename, save_para, gp);
	run_time = timer.reset();
	cout << "Done. It takes " << run_time << " ms." << endl;
	
	vector<EdgeID> path;	
	Query q(g); // contraction hierarchy method;
	Dijkstra dijk(g); // baseline method;
	HD h(g); // highway dimension method;

	// g.printRank();
	
	#ifdef TEST
	double d0 = 0, d1 = 0, d2 = 0;
	// Experiment 1: correctness check;
		if (experiment == 0 || experiment == 1) {
		cout << "Experiment 1: CORRECTNESS" << endl;
		int count = 0;	
		for ( int i = 0; i < query.size(); i=i+2 ) {
			d0 = h.run(query[i], query[i+1], path, unpack, id_map);
			//d1 = q.run(query[i], query[i+1], path);
			d2 = dijk.run(query[i], query[i+1], path);
			cout << d0 << ":" << d2;
			if ( fabs(d0-d2) > 1e-5 ) { count++; cout << ": wrong" << endl;}
			else cout << endl;
		}	
		cout << "query_num = " << query.size()/2 << ", wrong answer: " << count << endl;
	}
	// Experiment 2: performance;
	if (experiment == 0 || experiment == 2) {
		cout << "Experiment 2: PERFORMANCE" << endl;
		if (gp) {
			if (cmdLineArgs["unpack"].compare("true") == 0) unpack = true;
			else unpack=false;
		}
		// method: 1;
		if (method == 1 || method == 2) {
			cout << "Highway Dimension Method... ";
			timer.reset();	
			for ( int i = 0; i < query.size(); i=i+2 ) {
				d0 = h.run(query[i], query[i+1], path, unpack, id_map);
			}
			run_time = timer.reset();
			cout << "Done. It takes " << run_time << " ms." << endl;	
			// method: 2;
			/*cout << "Contraction Method... ";
			timer.reset();	
			for ( int i = 0; i < query.size(); i=i+2 ) {
				d1 = q.run(query[i], query[i+1], path);
			}
			run_time = timer.reset();
			cout << "Done. It takes " << run_time << " ms." << endl;*/
		}
		// method: 3;
		if (method == 0 || method == 2) {
			cout << "Dijkstra Method... ";
			timer.reset();
			for ( int i = 0; i < query.size(); i=i+2 ) {	
				d2 = dijk.run(query[i], query[i+1], path);
			}
			run_time = timer.reset();
			cout << "Done. It takes " << run_time << " ms." << endl;
		}
	}
	#endif

	while(true) {
		cout << "~~~~~~~~~~~~~~~~~~~NEW QUERY~~~~~~~~~~~~~~~~~~~" << endl;
		cout << "note: type \"stop\" to exit!" << endl;
		string input;
		
		int num_node = g.num_vertices();
		int source = 0, target = 0;
		cout << "Please input source node: ";
		while(true) {
			cin >> input;
			if (input.compare("stop") == 0) {
				cout << "program stops!" << endl;
				break;
			}
			if (!isNumeric(input)) {
				cerr << "input is invalid!" << endl;
				cerr << "please input again: ";
				continue;
			}			
			source = atoi(input.c_str());
			if (source < 0 || source >= num_node){
				cerr << "input is out of range!" << endl;
				cerr << "please input again: ";
			}else break;
		}
		
		if (input.compare("stop") == 0) return 0;		
		
		cout << "Please input target node: ";
		while(true) {
			cin >> input;
			if (input.compare("stop") == 0) {
				cout << "program stops!" << endl;
				break;
			}
			if (!isNumeric(input)) {
				cerr << "input is invalid!" << endl;
				cerr << "please input again: ";
				continue;
			}			
			target = atoi(input.c_str());
			if (target < 0 || target >= num_node){
				cerr << "input is out of range!" << endl;
				cerr << "please input again: ";
			}else break;
		}

		if (input.compare("stop") == 0) return 0;
		
		if (gp) {
			if (cmdLineArgs["unpack"].compare("true") == 0) unpack = true;
			else unpack = false;
		}else unpack = false;
		timer.reset();
		double distance = h.run(source, target, path, unpack, id_map);
		run_time = timer.reset();
		cout << "Query takes " << run_time << " ms." << endl;		
		
		cout << "Distance from " << source << " to " << target << ": " << distance << endl; 
	
	}
	return 0;
}

