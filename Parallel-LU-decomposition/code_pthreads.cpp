/*
				***Pthreads Implementation***
				
	Compile: 	g++ -std=c++11 -o code_pthreads code_pthreads.cpp -lpthread 
	Run:		./code_pthreads <Matrix Size> <No. of Threads>
*/
#include<iostream>
#include<vector>
#include<algorithm>
#include<string>
#include<cmath>
#include<ctime>
#include<chrono>
#include<pthread.h>

using namespace std;

// Global Variables used for the purposes of Debugging and/or Using Norm
int THREADS =1;
bool DEBUG = false;
bool NORM = false;

// Global Variables used for the pthread functions
int k_swap = 0;
int globalk_ = 0;
int n = 0;
int k_LU=0;

//initialize u as an n x n matrix with 0s below the diagonal
vector<vector<double> > u;

//initialize l as an n x n matrix with 1s on the diagonal and 0s above the diagonal
vector<vector<double> > l;	

//initialize the random input matrix as a global variable to be accessible by the pthread functions
vector<vector<double> > a;	
vector<vector<double> > a_;

// Helper Pthread functions
void* swap_thread(void* arg){

	int core = (long) arg;
	int min = (k_swap/THREADS)*(core);
	int max = (k_swap/THREADS)*(core+1);
	if (core==THREADS-1) max=k_swap;

	for(int i=min; i<max;i++){
		swap (l[k_swap][i], l[globalk_][i]);

	} 

}

void* LU_thread (void* arg){

	int core = (long) arg;
	int limit = n - k_LU - 1; 
	int min = (k_LU+1) + (limit/THREADS)*(core);
	int max = (k_LU+1) + (limit/THREADS)*(core+1);
	if (core==THREADS-1) max=n;

	for(int i=min; i<max;i++){

		l[i][k_LU] = a[i][k_LU]/u[k_LU][k_LU];
		u[k_LU][i] = a[k_LU][i];
	} 



}

void* final_thread (void* arg){

	int core = (long) arg;
	int limit = n - k_LU - 1; 
	int min = (k_LU+1) + (limit/THREADS)*(core);
	int max = (k_LU+1) + (limit/THREADS)*(core+1);
	if (core==THREADS-1) max=n;

	for(int i=min; i<max;i++){
		for(int j = k_LU+1;j<n;j++)
			a[i][j] -= l[i][k_LU]*u[k_LU][j];

	} 


}

// Random Number Generator
void randomizeMatrix(vector<vector<double> > &a,int n){
	a.resize(n);
	int tid = 1;
	double x;
	struct drand48_data randBuffer;
	srand48_r(time(NULL)*tid, &randBuffer);

	for(int i=0;i<n;i++){
		a[i].resize(n);
		for(int j=0;j<n;j++){
			drand48_r(&randBuffer, &x);
			a[i][j]=x;
		
		}
	}

}

// Basic Debugging functions
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


int main(int argc, char **argv){	
	
	//input n:size of matrix, t: no. of threads

	//Initialize clock
	auto total_start = std::chrono::system_clock::now();
 
	//Argument parsing
	if(argc<3){
		
		cout<<"Missing argumnets.\nPlease specify the matrix size followed by the number of threads"<<endl;
		return 0;
	}
	
	n = stoi(argv[1]);	//Size of matrix
	int t = stoi(argv[2]);	//Number of threads
	THREADS =t;
	cerr<<"Size: "<<n<<" Threads: "<<t<<endl;

	// Initializing L and U as n*n matrix with all zeros
	vector<vector<double> > u_temp(n,vector<double>(n,0));	
	u=u_temp;

	vector<vector<double> > l_temp(n,vector<double>(n,0));	
	l = l_temp;

	//Matrix initialization time:Start    
	auto rstart = std::chrono::system_clock::now();

	randomizeMatrix(a,n);
	a_=a;
	
	//Matrix initialization time:End
	auto rend = std::chrono::system_clock::now();
	chrono::duration<double> rtime = rend-rstart;
	cout<<"Initialization time: " << rtime.count() << "s\n";
	
	// Checking Initialization
	if(DEBUG) printMatrix(a,n);
	if(DEBUG) printMatrix(a_,n);

	//initialize pi as a vector of length n
	vector<double> pi(n);	
	for(int i=0;i<n;i++) pi[i]=i;
	if(DEBUG) printVector(pi,n);
	
	for (int i=0;i<n;i++)	l[i][i] =1;	//1s on the diagonal


/*----------------------------****Decomposition Begin****----------------------------*/

	cerr<<"Decomposition starts"<<endl;
	
	// Decomposition time: Start
	auto dcmp_start = std::chrono::system_clock::now();
	
	bool singular = false;
	for (int k =0;k<n;k++){

		// Finding row maximum
		double max = 0;
		int k_;

		for (int i=k;i<n;i++){
			if (max < abs(a[i][k])){
				max = abs(a[i][k]);
				k_ = i;
			}
		}
		
		if (max == 0){
			singular = true;
			cout<<"error (singular matrix)"<<endl;
			return 0;
		}

		// DEBUGGING: To find if the maximum value is obtained properly or not
		if(DEBUG) cerr<<"Max value found"<<endl;

		swap (pi[k], pi[k_]);	//swap values at k and k_
		swap (a[k], a[k_]);		//swap vectors at k and k_


		// Initializing the threads for the swapping function
		pthread_t threads_swap[THREADS];
		k_swap = k;
		globalk_ = k_;

                // Pthread optimization for swap (l[k][], l[k_][])	
		for (long i = 0; i < THREADS; i++) 
			pthread_create(&threads_swap[i], NULL, &swap_thread, (void*)i); 

		for (int i = 0; i < THREADS; i++) 
			pthread_join(threads_swap[i], NULL); 
		
		u[k][k] = a[k][k];

                // Pthread optimization for l[i][k] = a[i][k]/u[k][k] and u[k][i] = a[k][i]
		pthread_t threads_LU[THREADS];
		k_LU = k;

		for (long i = 0; i < THREADS; i++) 
			pthread_create(&threads_LU[i], NULL, &LU_thread, (void*)i); 

		for (int i = 0; i < THREADS; i++) 
			pthread_join(threads_LU[i], NULL); 

  
        	// Pthread optimization for the final computation of a[][]
		pthread_t threads_final[THREADS];

		for (long i = 0; i < THREADS; i++) 
			pthread_create(&threads_final[i], NULL, &final_thread, (void*)i); 

		for (int i = 0; i < THREADS; i++) 
			pthread_join(threads_final[i], NULL); 


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



