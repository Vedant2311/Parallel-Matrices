#include "mpi.h"
#include <bits/stdc++.h>
using namespace std;

class MapReduce{

	public:
		vector<double> keyvalue;	
		MapReduce(int total_webpages);
		void map(vector<vector<int> > outgoing_links,vector<double> pageranks,int total_webpages);
		void reduce(int total_webpages);
	
};
