/*##TODO
				***Single Implementation***			
	Compile: 	mpicc serial.c -o code
	
	//V: Order of Inputs will be => A->m*n, B->n*p 
	Run (for general analysis):		mpirun -np <No of PROCESSORS> ./code <m> <n> <p>
	Run (for our requirement):		mpirun -np <No of PROCESSORS> ./code <N> <32> <N>
	
	ASSUMPTION HERE: The row size of A should be divisible by the number of processors
*/

#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>
#include<time.h>
#include"mpi.h"

// V: Added parts to debug the isEqual function since the values are doubles and not integers
#include<math.h>
double epsilon = 0.0001;

// Global Variables used for debugging and computing NORM
int THREADS =1;
bool DEBUG = false;
bool CHECK = true;
bool DEMO =false;

//Basic Debugging functions
void printVector(double *v, int n){
	for(int i=0;i<n;i++)
		printf("%f ",v[i]);
	printf("\n");
}

void printMatrix(double *C,int m, int n){	
	for(int i=0;i<m*n;i++){
		printf("%f ",C[i]);
		if((i+1)%n==0)
			printf("\n");
	}
	printf("--\n");
}


//Basic Output functions:for DEMO
void printVectorF(double *v, int n, char* file){
	FILE *f;
	f = fopen(file,"w");
	for(int i=0;i<n;i++)
		fprintf(f,"%f ",v[i]);
}

void printMatrixF(double *C,int m, int n,char* file){
	FILE *f;
	f = fopen(file,"w");
	for(int i=0;i<m*n;i++){
		fprintf(f,"%f ",C[i]);
		if((i+1)%n==0)
			fprintf(f,"\n");
	}
}


/*
// Generating the Random Vector

//To generate random numbers faster##
inline int fastrand() { 
  g_seed = (214013*g_seed+2531011); 
  return (g_seed>>16)&0x7FFF; 
} 
*/
void randomizeVector(double *A,int n){
	double a = 1000.0;	//RAND_MAX_FLOAT
    for (int i=0;i<n;i++)
        A[i]=((double)rand()/(double)(RAND_MAX)) * a;
}

// Return 0: A ≠ B, Return 1: A = B
int IsEqual(double *A, double *B, int m, int n){
	
	for(int i=0;i<m*n;i++)
		// V: Means that they are not equal
		if(fabs(A[i]-B[i]) > epsilon)
			return 0;
	return 1;
}

//Matrix Multiply: assuming that matrix elements are row-major stored in an array

void Multiply_serial(double *A, double *B, double *C, int m, int n, int p){
	for (int i = 0; i < m; i++){
		for (int j = 0; j < p; j++){
			C[i*p + j] = 0;
			for (int k = 0; k < n; k++)
				C[i*p + j] += A[i*n + k] * B[k*p + j];
		}
	}
}	

//Main()
int main(int argc, char **argv){

	// Initializing the MPI environment
	int rank, size;
	
	MPI_Init (&argc, &argv);      /* starts MPI */
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);        /* get current process id */
    MPI_Comm_size (MPI_COMM_WORLD, &size);        /* get number of processes */

	// Defining the offset and the number of rows
	int offset, rows;
	MPI_Status status;

	//input n:size of matrix, t: no. of threads
	srand(time(0));	//Seed rand
	
	//Argument parsing
	if(argc<4){
		
		printf("Missing arguments.\nPlease specify matrix sizes m n p(Output will be of m X p)\n");
		return 0;
	}
	
	int m = atoi(argv[1]);	//M
	int n = atoi(argv[2]);	//N
	int p = atoi(argv[3]);	//P
	
	//Just to avoid recompiling ##
	if(m+p+n>100)	DEBUG = false;
/*	
	string infile = "";	//For demo
	if(DEMO)
		infile = argv[4];	//Input matrix file
*/
    	
	//Declare Matrices: ##Use double* 2d different methods: https://www.geeksforgeeks.org/dynamically-allocate-2d-array-C/	
	double *A = (double *)malloc(m * n * sizeof(double)); 
	double *B = (double *)malloc(n * p * sizeof(double)); 
	double *C = (double *)malloc(m * p * sizeof(double)); 
	double *C_serial = (double *)malloc(m * p * sizeof(double)); 
		
	// Get the matrices for the different rows
	double *A_rows = (double *)malloc((m/size) * n * sizeof(double));
	double *C_rows = (double *)malloc((m/size) * p * sizeof(double));	

	if (rank==0)
		printf("Size A: %dx%d, Size B: %dx%d\n",m,n,n,p);

	//Program time: Start
	clock_t total_start = clock();

	//Generate Random Vectors
	randomizeVector(A,m*n);
	randomizeVector(B,n*p);

	//To check initialization
	if(DEBUG){printf("Matrix A: \n"); printMatrix(A,m,n);}
	if(DEBUG){printf("Matrix B: \n"); printMatrix(B,n,p);}

	double start;
	
    /*-----------------------------------Starting of Parallel Code--------------------------------------*/	
    if (rank==0){
		printf("\nParallel Multiplication starts\n");

		// Start time for the MPI code
	    start = MPI_Wtime();
	}
	
	// Scatter matrix A among all the processes
	MPI_Scatter(A, m*n/size, MPI_DOUBLE, A_rows, m*n/size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
	// Broadcast Second matrix to all the processes
	MPI_Bcast(B,n*p,MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	// Applying the Barrier
	MPI_Barrier(MPI_COMM_WORLD);
	
	// Doing the Matrix Multiplication Now
	for (int i=0; i<m/size; i++){
		for(int j=0; j<p; j++){
			C[i*p + j] = 0;
			for (int k=0; k<n; k++){
				C_rows[i*p+j] += A_rows[i*n + k]*B[k*p + j];
			}
		}
	}
	
	// Gathering the result back into the C matrix
	MPI_Gather(C_rows, m*p/size, MPI_DOUBLE, C, m*p/size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
	// Setting the barrier for all the processes to reach here
	MPI_Barrier	(MPI_COMM_WORLD);

	// Finalizing the MPI environment
	MPI_Finalize();
	
	if (rank==0){
		double finish = MPI_Wtime();
		printf("Parallel time: %fs\n", (finish - start));

		
		/*-----------------------------Starting of Serial Multiplication---------------------------------*/
		printf("\nSerial Multiplication starts\n");
		
		//Multiplication time: Start
		clock_t mul_start = clock();

		//multily A,B and get result in C
		Multiply_serial(A,B,C_serial,m,n,p);

		//Multiplication time: Stop
		clock_t mul_stop = clock();
		printf("Serial Multiplication time: %fs\n", (double)(mul_stop - mul_start)/CLOCKS_PER_SEC);

		
		/*-----------------------------Starting of comparision between the Serial and Parallel code----------------------------*/		
		printf("\nComparision starts\n");
		
		//Comparision time: Start
		clock_t comp_start = clock();

		//Error checking
		int eq =IsEqual(C,C_serial,m,p);

		// Printing the value accordingly 
		if(eq)	printf("Correct\n");
		else	printf("Error\n");

		//Comparision time: Stop
		clock_t comp_stop = clock();
		printf("Comparision time: %fs\n", (double)(comp_stop - comp_start)/CLOCKS_PER_SEC);
		
		//Program time: Stop
		clock_t total_stop = clock();
		printf("\nProgram time: %fs\n", ((double)(total_stop - total_start))/CLOCKS_PER_SEC); 
		

	}

	// Finalizing the MPI environment
//	MPI_Finalize();
	
	return 0;
}
