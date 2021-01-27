/*##TODO
				***Single Implementation***			
	Compile: 	mpicc serial.c -o code
	Run:		mpirun -np <No of PROCESSORS> ./code <m> <n> <p>
*/

#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>
#include<time.h>

// Global Variables used for debugging and computing NORM
int THREADS =1;
bool DEBUG = true;
bool CHECK = true;
bool DEMO =false;

//Basic Debugging functions
void printVector(float *v, int n){
	for(int i=0;i<n;i++)
		printf("%f ",v[i]);
	printf("\n");
}

void printMatrix(float *C,int m, int n){	
	for(int i=0;i<m*n;i++){
		printf("%f ",C[i]);
		if((i+1)%n==0)
			printf("\n");
	}
	printf("--\n");
}


//Basic Output functions:for DEMO
void printVectorF(float *v, int n, char* file){
	FILE *f;
	f = fopen(file,"w");
	for(int i=0;i<n;i++)
		fprintf(f,"%f ",v[i]);
}

void printMatrixF(float *C,int m, int n,char* file){
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
void randomizeVector(float *A,int n){
	float a = 1000.0;	//RAND_MAX_FLOAT
    for (int i=0;i<n;i++)
        A[i]=((float)rand()/(float)(RAND_MAX)) * a;
}

// Return 0: A ≠ B, Return 1: A = B
int IsEqual(float *A, float *B, int m, int n){
	
	for(int i=0;i<m*n;i++)
		if(A[i]!=B[i])
			return 0;
	return 1;
}

//Matrix Multiply: assuming that matrix elements are row-major stored in an array

void Multiply_serial(float *A, float *B, float *C, int m, int n, int p){
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
	//input n:size of matrix, t: no. of threads
	
	srand(time(0));	//Seed rand
	
	//Program time: Start
	time_t total_start = time(NULL);
    printf("\nCurrent time is : %s\n",ctime(&total_start));
 
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
	printf("Size A: %dx%d, Size B: %dx%d\n",m,n,n,p);
    
    //Matrix initialization time:Start
	time_t init_start = time(NULL);
	
	//Declare Matrices: ##Use double* 2d different methods: https://www.geeksforgeeks.org/dynamically-allocate-2d-array-C/
    float *A = (float *)malloc(m * n * sizeof(float)); 
    float *B = (float *)malloc(n * p * sizeof(float)); 
    float *C = (float *)malloc(m * p * sizeof(float)); 
	float *C_serial = (float *)malloc(m * p * sizeof(float)); 
	//Generate Random Vectors
	randomizeVector(A,m*n);
	randomizeVector(B,n*p);
	
	//Matrix initialization time:Stop
	time_t init_stop = time(NULL);
	printf("Matrix initialization time: %ds\n", (init_stop - init_start));

	//To check initialization
	if(DEBUG){printf("Matrix A: \n"); printMatrix(A,m,n);}
	if(DEBUG){printf("Matrix B: \n"); printMatrix(B,n,p);}


/*----------------------------****Multiplication Begin****----------------------------*/

	printf("Multiplication starts\n");
	
	//Multiplication time: Start
	clock_t mul_start = clock();

	//multily A,B and get result in C
	Multiply_serial(A,B,C,m,n,p);

	//Multiplication time: Stop
	clock_t mul_stop = clock();
	printf("Serial Multiplication time: %fs\n", (double)(mul_stop - mul_start)/CLOCKS_PER_SEC);
	
	if(DEBUG){printf("Matrix C: \n"); printMatrix(C,m,p);}
/*----------------------------****Multiplication End****----------------------------*/
/*	if(DEMO){
		printMatrix(C,m,p,"./C"+itoa(m)+"x"+itoa(p)+".txt");
	}*/
/*----------------------------****C Written to files****----------------------------*/

	if(!CHECK) return 0;	//Exit if user doesn't want to check equality
	
	//Error Checking by comparing with serial implementation
	Multiply_serial(A,B,C_serial,m,n,p);
	int eq =IsEqual(C,C_serial,m,p);
	
	if(DEBUG){printf("Matrix C_serial: \n"); printMatrix(C_serial,m,p);}
	
	if(eq)	printf("Correct\n");
	else	printf("Error\n");
	
	//Program time: Stop
	time_t total_stop = time(NULL);
	printf("Program time: %ds\n", (total_stop - total_start));
	return 0;
}
