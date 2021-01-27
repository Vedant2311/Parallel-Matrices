/*  The compile and run commands are:

	Compile:- $ g++-7 -o mr-pr.o mr-pr.cpp
	Run:- $ ./mr-pr.o <input_file>.txt -o <output_file>.txt
		
	* Run command assumes the existence of a file <input_file>.txt in a directory named "test" in the cwd of the code
	* Run command also assumes the existence of a directory named "output" where it writes the outputs of the pagerank algorithm execution.  

*/

#include <bits/stdc++.h>
using namespace std;

// Functions for the proper usage of the code
const char *OUTPUT_ARG = "-o";

void usage() {
    cerr << "Use the below format " << endl
    	 << "pagerank <graph_file> -o <output_file>" << endl
         << " -o enable output " << endl;
}

// Global variables for the pagerank algorithm
vector<vector<double>> intermediate_page(100000);
vector<vector<int>> outgoing_links(100000);
vector<double> pageranks(100000);
vector<double> pageranks_calc(100000);

// Parameters of alpha and convergence rate
double alpha = 0.85;
double conv = 0.00001;
double dangling_pointer =0.0f;
int total_webpages=0;

// The mapper function. Uses the key and the outgoing links for the intermediate values
void emit(int key, double value){
	intermediate_page[key].push_back(value);
}

void mapper(int key, double page){
	int chilren = outgoing_links[key].size();
	if(chilren!=0){
		for(int i =0; i<chilren; i++){
			emit(outgoing_links[key][i],(double)(page/chilren));
		}
	}
	return;
}

// The reducer function. Gets the sum of all the divided probabilies coming from the parents of the nodes
void reducer(int key, vector<double> page){
	double pagerank_temp = 0.0f;
	for(int i =0; i<page.size(); i++){
		pagerank_temp += (double)page[i];
	}
	
	// Applying the pagerank update rule
	pageranks_calc[key] = pagerank_temp*alpha + (1-alpha)/total_webpages + alpha*dangling_pointer/total_webpages;
}

// Function to check the convergence
bool converging(vector<double> pageranks, vector<double> pageranks_calc, int total_webpages){
	bool isConverged= true;
	for (int i = 0; i < total_webpages; i++){
		if (abs(pageranks[i] - pageranks_calc[i])>conv){
			isConverged = false;
			break;
		}
	}
	return isConverged;
}

int main(int argc, char const *argv[])
{
	
	// File names
	string input_file;
	string output_file;
		
	// Reading the command line inputs		
	if (argc == 4){
		if (!strcmp(argv[2],OUTPUT_ARG)){
			input_file=argv[1];
			output_file=argv[3];
		}
		
		else{
			usage();
			exit(1);
		}
	}
	
	else{
		usage();
		exit(1);
	}
	
	int max_page_id = 0;
	int min_page_id = 0;	
	
	// Reading inputs from the specified file in the test folder of the standard pagerank repo
	ifstream fopen;
	fopen.open("test/" + input_file);
	while(!fopen.eof()){
		int a,b;
		fopen>>a>>b;
		outgoing_links[a].push_back(b);
		
		int max_temp = std::max(a,b);
		if (max_temp > max_page_id)
			max_page_id = max_temp;
		
		int min_temp = std::min(a,b);
		if (min_temp < min_page_id)
			min_page_id = min_temp;
	}
	
	// Total number of web-pages as: max_page_id - min_page_id + 1
	fopen.close();
	total_webpages = max_page_id - min_page_id + 1;

	cout << "File read successful" << endl;

	// Initializing the pageranks;
	double pgr = double(1.0f/total_webpages);
	for(int i=0; i<total_webpages; i++){
		pageranks[i] = pgr;
	}

	while(true){
		// Handling dangling pointers
		dangling_pointer=0.0;
		for(int i=0; i<total_webpages; i++){			
			if (outgoing_links[i].size()==0){
				dangling_pointer+=pageranks[i];
			}
		}		
		
		// Using the Map and Reduce
		for(int i=0; i<total_webpages; i++){
			mapper(i,pageranks[i]);
		}
		for(int i=0; i<total_webpages; i++){
			reducer(i,intermediate_page[i]);
			intermediate_page[i].clear();
		}	
		
		// Checking for convergence
		if(converging(pageranks, pageranks_calc, total_webpages)){
			break;
		}	
		else{
			for(int i=0; i<total_webpages; i++){
				pageranks[i] = pageranks_calc[i];
			}		
		}
	}

	// Writing to the file
	ofstream fout;
	fout.open("output/" + output_file);

	string str;
	double sumVal = 0.0;
	for(int i=0; i<total_webpages; i++){
		fout<<i<<" = "<<pageranks[i]<<endl;
		sumVal += pageranks[i];
	}	
	
	fout << "sum " << sumVal << endl;
	fout.close();

	return 0;
}
