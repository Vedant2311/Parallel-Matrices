/*##TODO
				***Single Implementation***			
	Compile: 	mpicc serial.c -o code
	
	//V: Order of Inputs will be => A->m*n, B->n*p 
	Run (for general analysis):		mpirun -np <No of PROCESSORS> ./code <m> <n> <p>
	Run (for our requirement):		mpirun -np <No of PROCESSORS> ./code <N> <32> <N>
	
	
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

// Global Variable for the mtype 
int mtype;

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
	MPI_Request request;

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
	
	// Calculating the number of rows according to our assumption
	rows = m/(size-1);
	
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

	// Randomizing Matrix B
//	randomizeVector(B,n*p);
//	MPI_Bcast(B,n*p,MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
	// The Process 0: 
	if (rank==0){

		printf("Size A: %dx%d, Size B: %dx%d\n",m,n,n,p);

		//Program time: Start
		clock_t total_start = clock();
 
		//Generate Random Vectors
		randomizeVector(A,m*n);
		randomizeVector(B,n*p);

		//To check initialization
		if(DEBUG){printf("Matrix A: \n"); printMatrix(A,m,n);}
		if(DEBUG){printf("Matrix B: \n"); printMatrix(B,n,p);}

	
		/*---------------------------Starting of communications between the Process 0 and the other processes---------------------------*/
		printf("\nParallel Multiplication starts\n");

		// Start time for the MPI code
	    double start = MPI_Wtime();
		offset = 0;
		
		// Having a fixed mtype for the master to send the data to the other workers
		mtype = 1; // FROM_MASTER;
		
		// Sending the matrices A and B to all the workers
		for (int i=1; i<size; i++){
					
			// Sending the A and B matrices to the different processes
            MPI_Send(&A[offset*n + 0], rows*n, MPI_DOUBLE, i, mtype, MPI_COMM_WORLD);
            MPI_Send(B, n*p, MPI_DOUBLE, i, mtype, MPI_COMM_WORLD);
            
         	
         	// Adding the offset value for the next row. Assumed to be stored as row-major		
			offset = offset + rows;
		}
		
		// Having a fixed mtype for the master to receive data from the workers
		mtype = 2; //FROM_WORKER;
		
		for (int i=1; i<size; i++){
			
			offset = (i-1) * rows;			
            MPI_Recv(&C[offset*p + 0], rows*p, MPI_DOUBLE, i, mtype, MPI_COMM_WORLD,&status);			
			
		}
		
		
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
	
	
	// The other processes
	else{
	
		// Setting the mtype to receive
		mtype = 1; //FROM_MASTER
		offset = (rank-1)*rows;
		
		// Receiving the A and B matrix
		MPI_Irecv(A, rows*n, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD, &request);
//		MPI_Wait(&request,&status);
		MPI_Irecv(B, n*p, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD, &request);
		MPI_Wait(&request,&status);
		
		// Performing the Matrix Multiplication
		Multiply_serial(A,B,C,rows,n,p);	
				
		// Setting the mtype to send
		mtype = 2;  //FROM_WORKER
		
		// Sending the offset, number of rows, and the C matrix value
		MPI_Isend(C, rows*p, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD,&request); 
		MPI_Wait(&request,&status);
				
	}
	
	
	
	// Finalizing the MPI environment
	MPI_Finalize();
	return 0;
}
