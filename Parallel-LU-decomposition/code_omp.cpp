/*
				***OpenMP Implementation***
				
	Compile: 	g++ -std=c++11 -o code_omp -fopenmp code_omp.cpp
	Run:		./code_omp <Matrix Size> <No. of Threads>
*/
#include<iostream>
#include<vector>
#include<algorithm>
#include<string>
#include<cmath>
#include<ctime>
#include<chrono>
#include "omp.h"


using namespace std;

// Global Variables used for debugging and computing NORM
int THREADS =1;
bool DEBUG = false;
bool NORM = false;

// Generating the Random Matrix
void randomizeMatrix(vector<vector<double> > &a,int n){
//	a.resize(n);
	int tid = omp_get_thread_num();
	double x;
	struct drand48_data randBuffer;
	srand48_r(time(NULL)*tid, &randBuffer);
	# pragma omp parallel for num_threads(THREADS)
	for(int i=0;i<n;i++){
//		a[i].resize(n);
		for(int j=0;j<n;j++){
			drand48_r(&randBuffer, &x);
			a[i][j]=x;
		
		}
	}

}

//Basic Debugging functions
void printVector(vector<double> &v, int n){
	for(int i=0;i<n;i++)
		cout<<v[i]<<" ";
	cout<<endl;
}

void printMatrix(vector<vector<double> > &m, int n){
	for(int i=0;i<n;i++)
		printVector(m[i],n);
	cout<<endl;
}


//Main()
int main(int argc, char **argv){
	//input n:size of matrix, t: no. of threads

	//Program time: Start
	auto total_start = std::chrono::system_clock::now();
 
	//Argument parsing
	if(argc<3){
		
		cout<<"Missing arguments.\nPlease specify matrix size followed by number of threads"<<endl;
		return 0;
	}
	
	int n = stoi(argv[1]);	//size of matrix
	int t = stoi(argv[2]);	//Number of threads	
	THREADS =t;
	cout<<"Size: "<<n<<" Threads: "<<t<<endl;
    
	//Matrix initialization time:Start
	auto init_start = std::chrono::system_clock::now();
	vector<vector<double> > a(n,vector<double>(n,0));	
	vector<vector<double> > a_;

	randomizeMatrix(a,n);
	a_=a;
	
	//Matrix initialization time:Stop
	auto init_stop = std::chrono::system_clock::now();
	chrono::duration<double> init_time = init_stop-init_start;
	cout<<"Matrix initialization time: " << init_time.count() << "s\n";

	//To check initialization
	if(DEBUG) printMatrix(a,n);
	if(DEBUG) printMatrix(a_,n);

	//initialize pi as a vector of length n
	vector<double> pi(n);
	for(int i=0;i<n;i++) pi[i]=i;
	if(DEBUG) printVector(pi,n);
	
	//initialize u as an n x n matrix with 0s below the diagonal
	vector<vector<double> > u(n,vector<double>(n,0));	//all zeros

	//initialize l as an n x n matrix with 1s on the diagonal and 0s above the diagonal
	vector<vector<double> > l(n,vector<double>(n,0));	//all zeros	
	for (int i=0;i<n;i++)	l[i][i] =1;	//1s on the diagonal


/*----------------------------****Decomposition Begin****----------------------------*/

	cout<<"Decomposition starts"<<endl;
	//Decomposition time: Start
	auto dcmp_start = std::chrono::system_clock::now();
	
	for (int k =0;k<n;k++){

		//Finding row maximum
		double max = 0;
		int k_;

		# pragma omp parallel for num_threads(THREADS)
		for (int i=k;i<n;i++){
			if (max < abs(a[i][k])){
				# pragma omp critical
				max = abs(a[i][k]);
				k_ = i;
			}
		}
		
		if (max == 0){
			cout<<"error (singular matrix)"<<endl;
			return 0;
		}

		swap (pi[k], pi[k_]);	//swap values at k and k_
		swap (a[k], a[k_]);		//swap vectors at k and k_

		for(int i=0;i<k;i++)	//swap values at k and k_ for all i = 0 to k-1
			swap (l[k][i], l[k_][i]);

		//Computing L and U
		u[k][k] = a[k][k];
	        # pragma omp parallel for num_threads(THREADS)
		for(int i = k+1;i<n;i++){
			l[i][k] = a[i][k]/u[k][k];
			u[k][i] = a[k][i];
		}
		# pragma omp parallel for num_threads(THREADS) collapse(2)
		for(int i = k+1;i<n;i++)
			for(int j = k+1;j<n;j++)
				a[i][j] -= l[i][k]*u[k][j];
	}

	//Decomposition time: Stop
	auto dcmp_stop = std::chrono::system_clock::now();
	chrono::duration<double> dcmp_time = dcmp_stop-dcmp_start;
	cout<<"Decomposition time: " << dcmp_time.count() << "s\n";

/*----------------------------****Decomposition End****----------------------------*/

	if(!NORM) return 0;	//Exit if user doesn't want to print norm
	
	//Error Checking by computing L2,1 norm
	vector<vector<double> > res(n,vector<double>(n,0));	//all zeros
	
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++){
			res[i][j]=a_[pi[i]][j];
			for(int k=0;k<n;k++)
				res[i][j] -= l[i][k]*u[k][j];	
		}
	
	//Compute L2,1 norm of res
	double norm = 0;
	for(int j=0;j<n;j++){	//for each column
		double colsum =0;
		for(int i=0;i<n;i++)
			colsum+=res[i][j]*res[i][j];
		colsum = sqrt(colsum);
		norm+=colsum;
	}
	cout<<"Norm: "<<norm<<endl;
	
	//Program time: Stop
	auto total_stop = std::chrono::system_clock::now();
  	chrono::duration<double> total_time = total_stop-total_start;
	cout<<"Program time: " << total_time.count() << "s\n";
	return 0;
}
