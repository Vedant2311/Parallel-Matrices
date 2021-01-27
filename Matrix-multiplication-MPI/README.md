# Matrix-multiplication-MPI

The standard matrix multiplication algorithm was parallelised here using different core-core communication strategies. These codes were progammed in C++ and MPI was used to implement the different communication strategies.

The function Multiply_serial() is used to perform multiplication without parallelism. Example code (assuming that matrix elements are row-major stored in an array):
```c
void Matrix_Multiply(float *A, float *B, float *C, int m, int n, int p){
	int i, j, k;
	for (i = 0; i < m; i++i){
		for (j = 0; j < p; j++){
			C[i*p + j] = 0;
			for (k = 0; k < n; k++)
				C[i*p + j] += A[i*n + k] * B[k*p + j];
		}
	}
}		
```

The function IsEqual() examines if two matrices are exactly the same. The function has the form:
```c
int IsEqual(float *A, float *B, int m, int n){
	for(int i=0;i<m*n;i++)
		if(A[i]!=B[i])
			return 0;
	return 1;
}
```

To implement a blocking P2P communication (MPI_Send/MPI_Recv) between processes, main() should at atleast perform the following steps. The similar details are followed for collective communication and non-blocking communication. 
```
Initialize and finalize MPI environment
Let Process #0 generate A (of size N*32) and B (of size 32*N) using random numbers in [0, 1].
Implement communications between Process 0 and other processes
Compute A*B using your parallel scheme (assuming the result is stored at C)
Let Process #0 compute A*B by Multiply_serial() function (assuming the result is stored at C_serial)
Let Process #0 check if C and C_serial are equal (using IsEqual() function) so as to verify your parallel algorithm
Let Process #0 get the running time of the two computations above
```

You can find the code usage and the different experimentations, observations, and explanations in the file **Report.pdf**
