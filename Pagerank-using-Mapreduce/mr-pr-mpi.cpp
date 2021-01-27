/*  The compile and run commands are:

	Compile:- $ mpic++ -o mr-pr-mpi.o mr-pr-mpi.cpp mapreduce.cpp
	Run:- $ mpirun -np <total_processors> ./mr-pr-mpi.o <input_file>.txt -o <output_file>.txt 
		
	* Run command assumes the existence of a file <input_file>.txt in a directory named "test" in the cwd of the code
	* Run command also assumes the existence of a directory named "output" where it writes the outputs of the pagerank algorithm execution.  

*/

#include <bits/stdc++.h>
#include "mapreduce.h"
#include "mpi.h"

using namespace std;

// Functions for the proper usage of the code
const char *OUTPUT_ARG = "-o";

void usage() {
    cerr << "Use the below format " << endl
    	 << "pagerank <graph_file> -o <output_file>" << endl
         << " -o enable output " << endl;
}

// Global variables for the pagerank algorithm
std::vector<vector<int> > outgoing_links(100000);
std::vector<double> pageranks(100000);
std::vector<double> pageranks_calc(100000);

// Parameters of alpha and convergence rate
double alpha = 0.85;
double conv = 0.00001;
double dangling_pointer =0.0f;
int total_webpages=0;


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


int main(int argc, char *argv[])
{

	MPI_Init(&argc, &argv);

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
	
	string temp = "test/" + input_file;
	fopen.open(temp.c_str());
	
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

	// Initializing the pageranks;
	double pgr = double(1.0f/total_webpages);
	for(int i=0; i<total_webpages; i++){
		pageranks[i] = pgr;
	}

	// Initializing the class 
	MapReduce *mr;

	while(true){
		// Handling dangling pointers
		dangling_pointer=0.0;
		for(int i=0; i<total_webpages; i++){			
			if (outgoing_links[i].size()==0){
				dangling_pointer+=pageranks[i];
			}
		}		

		MPI_Barrier(MPI_COMM_WORLD);

		mr = new MapReduce(total_webpages);
		
		mr->map(outgoing_links, pageranks, total_webpages);
		mr->reduce(total_webpages);

		// Getting the pageranks from the keyvalue list obtained from map-reduce
		for (int i = 0; i< total_webpages; i++){
			pageranks_calc[i] = (mr->keyvalue[i])*alpha + (1-alpha)/total_webpages + alpha*dangling_pointer/total_webpages;
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
	temp = "output/" + output_file;
	fout.open(temp.c_str());

	string str;
	double sumVal = 0.0;
	for(int i=0; i<total_webpages; i++){
		fout<<i<<" = "<<pageranks[i]<<endl;
		sumVal += pageranks[i];
	}	
	
	fout << "sum " << sumVal << endl;
	fout.close();
	
	MPI_Finalize();
	return 0;
}
