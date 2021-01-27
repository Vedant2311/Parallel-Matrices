/*  The compile and run commands are:

	Compile and Run
	mpic++ -c mr-pr-mpi-base.cpp -o mr-pr-mpi-base.o -I src
	mpic++ -g -O mr-pr-mpi-base.o src/libmrmpi_mpicc.a -o mr-pr-mpi-base
	mpirun -np 1 mr-pr-mpi-base bull.txt -o O_bull.txt
		
	* Run command assumes the existence of a file <input_file>.txt in a directory named "test" in the cwd of the code
	* Run command also assumes the existence of a directory named "output" where it writes the outputs of the pagerank algorithm execution.  

*/

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "sys/stat.h"
#include "mapreduce.h"
#include "keyvalue.h"
#include <bits/stdc++.h>

using namespace MAPREDUCE_NS;
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

// The function to get the Map value 
void mapper(int itask, KeyValue *kv, void *ptr){

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);

	for(int i =0; i<total_webpages/size; i++){
		int key = i + (rank * total_webpages/size);
		double page = pageranks[key];

		// Debug
	    double pgi = 0;
	    kv->add((char *) &key,sizeof(int),(char *) &pgi,sizeof(double));
		
		int children = outgoing_links[key].size();
		if (children!=0){
			for(int j=0; j<children; j++){
				int key_inner = outgoing_links[key][j];			
				double value = (double)(page/children);
				
				kv->add((char *) &key_inner,sizeof(int),(char *) &value,sizeof(double));		
			}
		}
	}

}

// Reduce function to get the sum of all the keys
void reduce(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr){
	char* values = multivalue;
	
	double pagerank_temp = 0.0f;
	for(int i =0; i<nvalues; i++){
		pagerank_temp += *(double*)values;	//Add the values
		values += valuebytes[i];	//TO the next value
	}
	
  	kv->add(key,keybytes,(char *) &pagerank_temp,sizeof(double));

}

// A special function to get the keys into the pages
void update(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr){

	int j = *(int*)(key);
	double pagerank_temp = *(double*)value;
	
	// Applying the pagerank update rule
	pageranks_calc[j] = pagerank_temp;

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


int main(int argc, char *argv[])
{

	MPI_Init(&argc, &argv);

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);

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

		mr = new MapReduce(MPI_COMM_WORLD);
		
		mr->map(size,mapper,NULL);

	    mr->gather(1);
	    mr->convert();
	    
	    mr->reduce(reduce,NULL);	    
	    mr->broadcast(0);
	    
	    MPI_Barrier(MPI_COMM_WORLD);
	    mr->map(mr,update,NULL);

		delete mr;
		
		// Getting the pageranks from the keyvalue list obtained from map-reduce
		for (int i = 0; i< total_webpages; i++){
			pageranks_calc[i] = (pageranks_calc[i])*alpha + (1-alpha)/total_webpages + alpha*dangling_pointer/total_webpages;
		}
		
		// Checking for convergence
		if(converging(pageranks, pageranks_calc, total_webpages)){
			break;
		}	
		else{
			for(int i=0; i<total_webpages; i++){
				pageranks[i] = pageranks_calc[i];
			}
			cout << pageranks[0] << endl;		
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(rank==0){
		// Writing to the file
		cout << "writing" << endl;
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
	}
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();
	return 0;
}
