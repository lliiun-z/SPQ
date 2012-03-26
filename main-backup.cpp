#include "GraphUtil.h"
#include "CoreSearch.h"

using namespace std;

//#define TEST

void printUsage() {
	cerr << "Usage of 6DegreeSearch:" << endl;
	cerr << "  coresearchd --arg1=val1 --arg2=val2" << endl << endl << endl;

	cerr << "Command line parameters:" << endl << endl;

	cerr << "--help" << endl;
	cerr << "--usage" << endl;
	cerr << "  Prints this usage information and exits." << endl << endl;
	
	cerr << "--debug" << endl;
	cerr << "  Validate the results using BFS as benchmark." << endl << endl;
	
	cerr << "--experiment" << endl;
	cerr << "  Compare the performance of all methods." << endl << endl;
	
	cerr << "--no-path" << endl;
	cerr << "  Compute shortest distance only." << endl << endl;
	
	cerr << "--write-corelabel=data_coresize.label" << endl;
	cerr << "  Output the core information and related labels to data_coresize.label." << endl << endl;

	cerr << "--query-num=1000" << endl;
	cerr << "  Specifies the number of random shortest path queries to be 1000." << endl << endl;

	cerr << "--coresize=1000" << endl;
	cerr << "  Specifies the number of core nodes to be 1000. (should be no more than 65535)" << endl << endl;

	cerr << "--radius=6" << endl;
	cerr << "  Specifies the maximal shortest path length to be 6." << endl << endl;

	cerr << "--method=coresearch" << endl;
	cerr << "  Specifies the algorithm to use. Valid implementations:" << endl;
	cerr << "  coresearch, bfs, bibfs (default=coresearch)." << endl << endl;

	cerr << "--load-corelabel=data_coresize.label" << endl;
	cerr << "  Specifies the intermediate file." << endl << endl;

	cerr << "--resultfile=result.txt" << endl;
	cerr << "  Specifies the file to output result." << endl << endl;
	
	cerr << "--filename=data.gra" << endl;
	cerr << "  Use data.gra as input file." << endl << endl;
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

double computeAvgError(vector<int>& distvec, vector<int>& truedist) {
	if (distvec.size()!=truedist.size())
		return MAX_VAL;
	double sum=0.0;
	int num=0, errnum=0;
	for (int i = 0; i < distvec.size(); i++) {
		if (distvec[i]>=truedist[i] && truedist[i]>0) {
			sum += ((distvec[i]-truedist[i])*1.0)/(truedist[i]*1.0);
			num++;
			if (distvec[i]>truedist[i]) {
				if (distvec[i]>12) distvec[i]=7;
				cout << distvec[i] << "->" << truedist[i] << endl;
				errnum++;
			}
		}
	}
	sum = (sum*1.0)/(num*1.0);
	cout << "Average relative error=" << sum << " errnum=" << errnum << " (" << num << ")" << endl;
	return sum;
}

#ifdef TEST
int main(int argc, char* argv[]) {
	if (argc == 1) {
		printUsage();
		return 1;
	}
	
	char* filename = argv[1];
	int coresize = atoi(argv[2]);
	ifstream infile(filename);
	if (!infile) {
		cout << "Error: Cannot open " << filename << endl;
		return -1;
	}
	Graph g(infile,false);
	int gsize = g.num_vertices();
	cout << "#vertex size:" << gsize << "\t#edges size:" << g.num_edges() << endl;
	
	g.printGraph();
	cout << "outlist" << endl;
	for (int i = 0; i < gsize; i++) {
		cout << i << ": ";
		forall_outneighbors(g,i,eit) {
			cout << *eit << " ";
		}
		cout << endl;
	}
	cout << "inlist" << endl;
	for (int i = 0; i < gsize; i++) {
		cout << i << ": ";
		forall_inneighbors(g,i,eit) {
			cout << *eit << " ";
		}
		cout << endl;
	}

//	bit_vector* corenodes = new bit_vector(gsize);
	vector<int> nodemap = vector<int>(gsize,0);
	for (int i = 0; i < nodemap.size(); i++)
		nodemap[i] = i;
	random_shuffle(nodemap.begin(),nodemap.end());
	cout << "nodemap: " << endl;
	for (int i = 0; i < nodemap.size(); i++)
		cout << i << ":" << nodemap[i] << " ";
	cout << endl;
//	g.reorderid(nodemap);
	g.printGraph();
	g.constructPortlist(3);
	for (int i = 0; i < gsize; i++) {
		cout << i << ": ";
		forall_outneighborport(g,i,eit,port_iter) {
			cout << *eit << "|" << *port_iter << " ";
		}
		cout << endl;
	}
//	if (true) return 0;
	
	
	
	CoreSearch cs(g,6,coresize);
	cs.createLabels();
	g.printGraph();
	vector<int> path;
	/*
	path = cs.shortestpath(2,0);
	cout << "path: ";
	for (int i = 0; i < path.size(); i++)
		cout << path[i] << " ";
	cout << endl;
	path = cs.shortestpath(9,1);
	cout << "path: ";
	for (int i = 0; i < path.size(); i++)
		cout << path[i] << " ";
	cout << endl;
	path = cs.shortestpath(9,2);
	cout << "path: ";
	for (int i = 0; i < path.size(); i++)
		cout << path[i] << " ";
	cout << endl;
	
	path = cs.shortestpath(8,9);
	cout << "path: ";
	for (int i = 0; i < path.size(); i++)
		cout << path[i] << " ";
	cout << endl;	
	path = cs.shortestpath(8,6);
	cout << "path: ";
	for (int i = 0; i < path.size(); i++)
		cout << path[i] << " ";
	cout << endl;
	path = cs.shortestpath(7,9);
	cout << "path: ";
	for (int i = 0; i < path.size(); i++)
		cout << path[i] << " ";
	cout << endl;	
	path = cs.shortestpath(3,4);
	cout << "path: ";
	for (int i = 0; i < path.size(); i++)
		cout << path[i] << " ";
	cout << endl;
	*/
//	if (true) return 0;
	
	
	int ref=0, searchspace=0, ss1, ss2, sdist, cdist; 
	long sum_ss1=0, sum_ss2=0, sum_ss=0;
	vector<vector<int> > distvec = vector<vector<int> >(2,vector<int>(gsize,0));
	vector<int> dists = vector<int>(gsize,0);
	vector<int> quevec = vector<int>(gsize,0);	
	/*
	for (int i = 0; i < gsize; i++) {
		for (int j = i+1; j < gsize; j++) { 
			cdist = cs.distance(i,j);
			sdist = GraphUtil::BiBFSDist(g,i,j,dists,quevec,ref,6,ss1);
			int tmp = GraphUtil::BFSDist(g,i,j,dists,quevec,ref,6,ss2);
			cout << i << "->" << j << " tmp=" << tmp << " ss2=" << ss2 << " sdist=" << sdist << " ss1=" << ss1 << " cdist=" << cdist<< endl;
			if (sdist!=tmp || sdist!=cdist) {
				cout << "Wrong: " << i << "->" << j << " sdist=" << sdist << " tmp=" << tmp << " cdist=" << cdist << endl;
				exit(0);
			}
			sum_ss1 += ss1;
			sum_ss2 += ss2;
			sum_ss  += cs.getSearchspace();
		}
	}
	cout << "sum_ss2BFS=" << sum_ss2 << " sum_ssCBiBFS=" << sum_ss << " sum_ss1BiBFS=" << sum_ss1 << endl;
	cout << "ss1=" << ss1 << " ss2=" << ss2 << " searchspace=" << cs.getSearchspace() << endl;
	if (true) return 0;
	*/
	
	
	vector<int> prev = vector<int>(gsize,-1);
	sum_ss1=0, sum_ss2=0, sum_ss=0;
	for (int i = 0; i < gsize; i++) {
		for (int j = 0; j < gsize; j++) {
			path.clear();
			path = cs.shortestpath(i,j);
			cdist = path.size()-1;
			cout << "Core path: ";
			for (int k = 0; k < path.size(); k++)
				cout << path[k] << " ";
			cout << endl;
			path.clear();
			path = GraphUtil::BiBFSPath(g,i,j,dists,quevec,prev,ref,6,ss1);
			cout << "path size=" << path.size() << endl;
			sdist = path.size()-1;
			cout << "BiBFS path: ";
			for (int k = 0; k < path.size(); k++)
				cout << path[k] << " ";
			cout << endl;
			path.clear();
			path = GraphUtil::BFSPath(g,i,j,dists,quevec,prev,ref,6,ss2);
			int tmp = path.size()-1;
			cout << "BFS path: ";
			for (int k = 0; k < path.size(); k++)
				cout << path[k] << " ";
			cout << endl;
			cout << i << "->" << j << " tmp=" << tmp << " ss2=" << ss2 << " sdist=" << sdist << " ss1=" << ss1 << endl;
			cout << "ss1=" << ss1 << " ss2=" << ss2 << " searchspace=" << cs.getSearchspace() << endl;
			if (sdist!=tmp || sdist!=cdist) {
				cout << "Wrong: " << i << "->" << j << " sdist=" << sdist << " tmp=" << tmp  << endl;
				exit(0);
			}
			sum_ss1 += ss1;
			sum_ss2 += ss2;
			sum_ss  += cs.getSearchspace();
		}
	}
	cout << "sum_ss2BFS=" << sum_ss2 << " sum_ssCBiBFS=" << sum_ss << " sum_ss1BiBFS=" << sum_ss1 << endl;
	
	double vm_usage, resident_set;
	Util::process_mem_usage(vm_usage, resident_set);
	if (true) return 0;
}
#else
int main(int argc, char* argv[]) {
	cout.setf(ios::fixed,ios::floatfield);
	cout.precision(4);

	// Increase stack size
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

	typedef map<string,string> mapType;
	map<string, string> cmdLineArgs;
	cmdLineArgs["filename"] = "";
	cmdLineArgs["help"] = "unset";
	cmdLineArgs["usage"] = "unset";
	cmdLineArgs["debug"] = "unset";
	cmdLineArgs["no-path"] = "unset";
	cmdLineArgs["no-bfs"] = "unset";
	cmdLineArgs["query-num"] = "10000";
	cmdLineArgs["query-distance"] = "-1";
	cmdLineArgs["nodeordertype"] = "1";
	cmdLineArgs["method"] = "coresearch";
	cmdLineArgs["coresize"] = "1000";
	cmdLineArgs["seednum"] = "1000";
	cmdLineArgs["radius"] = "6";
	cmdLineArgs["resultfile"] = "./result.txt";
	cmdLineArgs["load-corelabel"] = "";
	cmdLineArgs["write-corelabel"] = "";
	cmdLineArgs["experiment"] = "unset"; // compare all methods
	cmdLineArgs["degdecompexp"] = "unset";
	cmdLineArgs["estimate"] = "unset";
	cmdLineArgs["90effdia"] = "unset";
	cmdLineArgs["genqueryfile"] = "";
	cmdLineArgs["queryfile"] = "";

	for (int i = 1; i < argc; i++) {
		string currArg = argv[i];
		size_t found = string::npos;

		// Iterate over all possible keys in map
	    // Using a const_iterator since we are not going to change the values.
		for(mapType::const_iterator it = cmdLineArgs.begin(); it != cmdLineArgs.end(); ++it)
		{
			string fullArg = "--" + it->first;
			found = currArg.find(fullArg);
			if (found != string::npos) {
				if (currArg.find(fullArg + "=") != string::npos) {
					// record value
					cmdLineArgs[it->first] = currArg.substr(fullArg.size() + 1);
				} else {
					// no --argument=value, just --argument
					cmdLineArgs[it->first] = "";
				}
				break;
			}
		}

		if (found == string::npos) {
			// Commandline argument invalid
			cerr << "Invalid commandline argument: " << currArg << endl;
			printUsage();
			exit(1);
		}
	}

	if (cmdLineArgs["filename"] == "") {
		cerr << "No input file provided, please use command-line option --filename=data.gra" << endl;
		printUsage();
		exit (1);
	}

	if (cmdLineArgs["usage"] != "unset" || cmdLineArgs["help"] != "unset") {
		printUsage();
		exit(0);
	}

	 string method = cmdLineArgs["method"]; // options: coresearch, bfs, bibifs
	 int query_num = atoi(cmdLineArgs["query-num"].c_str());
	 int coresize = atoi(cmdLineArgs["coresize"].c_str());
	 int radius = atoi(cmdLineArgs["radius"].c_str());
	 int query_distance = atoi(cmdLineArgs["query-distance"].c_str());
	 int nodeordertype = atoi(cmdLineArgs["nodeordertype"].c_str());
	 int seednum = atoi(cmdLineArgs["seednum"].c_str());
	 string filename = cmdLineArgs["filename"];
	 string genqueryfile = cmdLineArgs["genqueryfile"];
	 bool debug = (cmdLineArgs["debug"] != "unset");
	 bool no_path = (cmdLineArgs["no-path"] != "unset");
	 bool no_bfs = (cmdLineArgs["no-bfs"] != "unset");
	 bool experiment = (cmdLineArgs["experiment"] != "unset");
	 bool degdecompexp = (cmdLineArgs["degdecompexp"] != "unset");
	 bool estimate = (cmdLineArgs["estimate"] != "unset");
	 bool effdia = (cmdLineArgs["90effdia"] != "unset");
	 string resultfile = cmdLineArgs["resultfile"];
	 string load_interfile = cmdLineArgs["load-corelabel"];
	 string write_interfile = cmdLineArgs["write-corelabel"];
	 string queryfile = cmdLineArgs["queryfile"];
	
	// output machine configuration
	vector<string> cpumeminfo;
	GraphUtil::readCpuMemInfo(cpumeminfo);
	for (int i = 0; i < cpumeminfo.size(); i++)
		cout << cpumeminfo[i] << endl;
	
	// write the results to resultfile
	vector<string> results;
	results.push_back(filename);
	if (experiment) results.push_back("exp");
	else results.push_back(method);
	if (no_path) results.push_back("DIST");
	else results.push_back("PATH");
	results.push_back(to_string(query_num));
	
	try {
		cout << "Directed Version" << endl;
		double run_time=0, label_time=0, query_time=0, peakmem=0, estquery_time=0;
		printParameters(cmdLineArgs);
		PerformanceTimer timer = PerformanceTimer::start();
		cout << "Reading graph file: " << filename << "..." << endl;
		ifstream infile(filename.c_str());
		if (!infile) {
			cout << "Error: Cannot open " << filename << endl;
			return -1;
		}	
		Graph g(infile,false);
		cout << "#V=" << g.num_vertices() << " #E=" << g.num_edges() << endl;
		int gsize=g.num_vertices();
		if (g.isDirectedGraph()) cout << "Graph is directed graph!" << endl;
		else cout << "Graph is undireceted graph!" << endl;
		
		if (genqueryfile!="") {
			if (query_distance==-1) {
				cout << "Generating random distance queries..." << endl;
				timer.reset();
				ofstream queryout(genqueryfile.c_str());
				GraphUtil::generateRandomQueries1(g,radius,query_num,queryout);
				queryout.close();
				run_time = timer.reset();
				cout << "Done. It took " << run_time << " (ms)" << endl;
			}
			else {
				cout << "Generating <=" << query_distance << " distance queries..." << endl;
				timer.reset();
				int qdist=max(query_distance,1);
				ofstream queryout(genqueryfile.c_str());
				GraphUtil::generateRandomBoundedQueries(g,qdist,query_num,queryout);
			//	GraphUtil::generateQueriesByDistanceAndOutput(g,qdist,query_num,queryout);
				queryout.close();
				run_time = timer.reset();
				cout << "Done. It took " << run_time << " (ms)" << endl;
			}
			return 0;
		}
		
		if (effdia) {
			cout << "Computing 90% effective diameter..." << endl;
			timer.reset();
			double diameter90 = GraphUtil::effective90diameter(g,seednum);
			run_time = timer.reset();
			cout << "Done. It took " << run_time << " (ms)" << endl;
			cout << "90% effective diameter=" << diameter90 << endl;
			return 0;
		}
		
		if (degdecompexp) {
			cout << "Performing graph decomposition after removing core nodes..." << endl;
			timer.reset();
			vector<int> components = GraphUtil::degreeDecomposition(g,coresize);
			int max_compsize=0, min_compsize=gsize;
			long sum=0, avg_compsize=0;
			sort(components.begin(),components.end());
			for (int i = 0; i < components.size(); i++) {
				if (max_compsize<components[i]) max_compsize=components[i];
				if (min_compsize>components[i]) min_compsize=components[i];
				sum += components[i];
				cout << i << "\t" << components[i] << endl;
			}
			avg_compsize = sum/components.size();
			cout << "Num_comps=" << components.size() << " max_compsize=" << max_compsize
					<< " min_compsize=" << min_compsize << " avg_compsize=" << avg_compsize << endl;
			run_time = timer.reset();
			cout << "Done. It took " << run_time << " (ms)" << endl;
			return 0;
		}
		
		int left=0, s, t;
		vector<int> src, trg, truedist, distvec;
		vector<int>::iterator sit, tit;
		cout << "Generating random queries (query_distance=" << query_distance << ")..." << endl;
		timer.reset();
		if (query_distance==-1) {
			if (queryfile=="") {
				srand48(time(NULL));
				while (left<query_num) {
					s = lrand48() % gsize;
					t = lrand48() % gsize;
					if (s==t) continue;
					src.push_back(s);
					trg.push_back(t);
					++left;
				}
			}
			else {
				ifstream queryin(queryfile.c_str());
				while (!queryin.eof()) {
					queryin >> s >> left >> t;
					src.push_back(s);
					trg.push_back(t);
					truedist.push_back(left);
				}
				queryin.close();
				query_num = src.size();
			}
		}
		else {
			int qdist=min(query_distance,radius);
			GraphUtil::generateQueriesByDistance(g,qdist,query_num,src,trg);
		}
		distvec = vector<int>(query_num,-1);
		run_time = timer.reset();
		cout << "Done. It took " << run_time << " (ms)" << endl;		

		int sdist=0, distance=100000, counter=0;
		vector<int> path;
		unsigned long ss_cs=0, ss_bibfs=0, ss_bfs=0, sum_ctimes=0, num_entries=0, num_corepathentries=0; // sum of search space
		if (experiment || method=="bibfs") {
			cout << "----------------------------------------BIBFS-----------------------------------------" << endl;
			label_time = 0;
			cout << "Initializing utility data structures..." << endl;
			timer.reset();
			vector<int> prev;
			if (!no_path) 
				prev = vector<int>(gsize,-1);
			vector<int> dist = vector<int>(gsize,0);
			vector<int> que = vector<int>(gsize,0);
			int ref=0, tmp_ss=0;
			run_time = timer.reset();
			cout << "Done. It took " << run_time << " (ms)" << endl;
			if (no_path) {
				cout << "Performing shortest distance queries..." << endl;
				timer.reset();
				for (sit = src.begin(), tit = trg.begin(); sit != src.end(); ++sit, ++tit) {
					sdist = GraphUtil::BiBFSDist(g,*sit,*tit,dist,que,ref,radius,tmp_ss);
					ss_bibfs += tmp_ss;
				}
				query_time = timer.reset();
			}
			else {
				cout << "Performing shortest path queries..." << endl;
				timer.reset();
				for (sit = src.begin(), tit = trg.begin(); sit != src.end(); ++sit, ++tit) {
					path = GraphUtil::BiBFSPath(g,*sit,*tit,dist,que,prev,ref,radius,tmp_ss);
					ss_bibfs += tmp_ss;
				}
				query_time = timer.reset();
			}
			cout << "BiBFSQuery time (query_num=" << query_num << "): " << query_time << "(ms)" << endl;
			cout << "Average query time=" << (query_time*1.0)/(query_num*1.0) << "(ms)" << endl;
			dist.clear();
			que.clear();
			prev.clear();
		}
		results.push_back(to_string(query_time));
		if (!no_bfs && (experiment||method=="bfs")) {
			cout << "-----------------------------------------BFS------------------------------------------" << endl;
			label_time = 0;
			cout << "Initializing utility data structures..." << endl;
			timer.reset();
			vector<int> prev;
			if (!no_path) 
				prev = vector<int>(gsize,-1);
			vector<int> dist = vector<int>(gsize,0);
			vector<int> que = vector<int>(gsize,0);
			int ref=0, tmp_ss=0;
			run_time = timer.reset();
			cout << "Done. It took " << run_time << " (ms)" << endl;
			if (no_path) {
				cout << "Performing shortest distance queries..." << endl;
				counter = 0;
				timer.reset();
				for (sit = src.begin(), tit = trg.begin(); sit != src.end() && counter<1000; ++sit, ++tit, ++counter) {
					sdist = GraphUtil::BFSDist(g,*sit,*tit,dist,que,ref,radius,tmp_ss);
					ss_bfs += tmp_ss;
				}
				query_time = timer.reset();
			}
			else {
				cout << "Performing shortest path queries..." << endl;
				timer.reset();
				/*
				for (sit = src.begin(), tit = trg.begin(); sit != src.end() && counter<5000; ++sit, ++tit, ++counter) {
					path = GraphUtil::BFSPath(g,*sit,*tit,dist,que,prev,ref,radius,tmp_ss);
					ss_bfs += tmp_ss;
				}
				*/
				int step=(int)(src.size()/2000);
				counter = 0;
				for (int i = 0; i < src.size(); i+=step, ++counter) {
					path = GraphUtil::BFSPath(g,src[i],trg[i],dist,que,prev,ref,radius,tmp_ss);
					ss_bfs += tmp_ss;
				}
				query_time = timer.reset();
			}
			cout << "BFSQuery time (query_num=" << counter << "): " << query_time << "(ms)" << endl;
			cout << "Average query time=" << (query_time*1.0)/(counter*1.0) << "(ms)" << endl;
			dist.clear();
			que.clear();
			prev.clear();	
		}
		results.push_back(to_string(query_time));
		// note that: make sure this is the last algorithm to be performed becuase queries will be mapped into new graph after node ordering
		if (experiment || method=="coresearch") {
			cout << "--------------------------------------CORESEARCH--------------------------------------" << endl;
			CoreSearch* cs = new CoreSearch(g,radius,coresize,nodeordertype,estimate);
			timer.reset();
			cout << "Performing core-based labeling..." << endl;
			cs->createLabels(src,trg,load_interfile);
			label_time = timer.reset();
			cout << "Done. Labeling took " << label_time << " (ms)" << endl; 
			if (write_interfile!="") {
				cout << "Writing core-label information to " << write_interfile << "..." << endl;
				timer.reset();
				ofstream outfile(write_interfile.c_str());
				cs->writeIntermediateResults(outfile);
				outfile.flush();
				outfile.close();
				run_time = timer.reset();
			}
			
			if (no_path) {
				cout << "Performing shortest distance queries..." << endl;
				if (debug) {
					cout << "Debug Mode" << endl;
					timer.reset();
					for (sit = src.begin(), tit = trg.begin(); sit != src.end(); ++sit, ++tit) {
						cs->test_distance(*sit,*tit);
						ss_cs += cs->getSearchspace();
						sum_ctimes += cs->getComparetimes();
					}
					query_time = timer.reset();
				}
				else {
					timer.reset();
					for (sit = src.begin(), tit = trg.begin(); sit != src.end(); ++sit, ++tit) {
						cs->distance(*sit,*tit);
						ss_cs += cs->getSearchspace();
						sum_ctimes += cs->getComparetimes();
					}
					query_time = timer.reset();
				}
			}
			else {
				cout << "Performing shortest path queries..." << endl;
				if (debug) {
					cout << "Debug Mode" << endl;
					timer.reset();
					for (sit = src.begin(), tit = trg.begin(); sit != src.end(); ++sit, ++tit) {
						cs->test_shortestpath(*sit,*tit);
						ss_cs += cs->getSearchspace();
						sum_ctimes += cs->getComparetimes();
					}
					query_time = timer.reset();
				}
				else {
					counter = 0;
					timer.reset();
					for (sit = src.begin(), tit = trg.begin(); sit != src.end(); ++sit, ++tit, ++counter) {
						cs->shortestpath(*sit,*tit,distance);
						ss_cs += cs->getSearchspace();
						sum_ctimes += cs->getComparetimes();
						distvec[counter] = distance;
					}
					query_time = timer.reset();
					
					
					// perform estimate version
					/*
					if (!estimate) {
						estimate = true;
						cout << "Performing estimated version of shortest path queries..." << endl;
						cs->setEstimate(true);
						counter = 0;
						timer.reset();
						for (sit = src.begin(), tit = trg.begin(); sit != src.end(); ++sit, ++tit, ++counter) {
							cs->shortestpath(*sit,*tit,distance);
						//	ss_cs += cs->getSearchspace();
						//	sum_ctimes += cs->getComparetimes();
							distvec[counter] = distance;
						}
						estquery_time = timer.reset();
					}
					*/
				}
			}
			num_corepathentries = cs->getCorepathsize();
			num_entries = cs->getNumRouteEntries();
			peakmem = cs->getPeakmemusage();
			delete cs;
			cout << "CSQuery time (query_num=" << query_num << "): " << query_time << "(ms) comparetimes=" 
					<< sum_ctimes << " num_RouteEntries=" << num_entries << " num_CorePathEntries=" << num_corepathentries 
					<< " peakmem=" << peakmem << "KB" << endl;
			cout << "Average query time=" << (query_time*1.0)/(query_num*1.0) << "(ms)" << endl;
			if (estimate) {
				double avgerror = computeAvgError(distvec,truedist);
				cout << "Est Average query time=" << (estquery_time*1.0)/(query_num*1.0) << "(ms)" << endl;
			}
		}
		results.push_back(to_string(label_time));
		results.push_back(to_string(query_time));
		cout << "searchspace_coresearch=" << ss_cs << " searchspace_bibfs=" << ss_bibfs << " searchspace_bfs=" << ss_bfs << endl;
		
		// append more output results
		results.push_back(to_string(ss_cs)); results.push_back(to_string(ss_bibfs)); results.push_back(to_string(ss_bfs));
		results.push_back(to_string(sum_ctimes));
		results.push_back(to_string(num_entries));
		results.push_back(to_string(num_corepathentries));
		results.push_back(to_string(peakmem)+"KB");
		for (int i = 0; i < cpumeminfo.size(); i++)
			results.push_back(cpumeminfo[i]);
		keepResult(resultfile.c_str(),results);
	} catch (string str) {
		cerr << "Exception: " << str << endl;
		cerr.flush();
	}
	
	return 0;
}
#endif
