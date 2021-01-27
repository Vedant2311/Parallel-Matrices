# Parallel-LU-decomposition

Two parallel implementations of LU decomposition that use Gaussian elimination to factor a dense N x N matrix into an upper-triangular one and a lower-triangular one are developed here. In matrix computations, pivoting involves finding the largest magnitude value in a row, column, or both and then interchanging rows and/or columns in the matrix for the next step in the algorithm. The purpose of pivoting is to reduce round-off error, which enhances numerical stability. Here row pivoting is used, a form of pivoting involves interchanging rows of a trailing submatrix based on the largest value in the current column. To perform LU decomposition with row pivoting, we will compute a permutation matrix P such that PA = LU. The permutation matrix keeps track of row exchanges performed.

Below is pseudocode for a sequential implementation of LU decomposition with row pivoting.
```
   inputs: a(n,n)
   outputs: π(n), l(n,n), and u(n,n)

   initialize π as a vector of length n
   initialize u as an n x n matrix with 0s below the diagonal
   initialize l as an n x n matrix with 1s on the diagonal and 0s above the diagonal
   for i = 1 to n
     π[i] = i
   for k = 1 to n
     max = 0
     for i = k to n
       if max < |a(i,k)|
         max = |a(i,k)|
         k' = i
     if max == 0
       error (singular matrix)
     swap π[k] and π[k']
     swap a(k,:) and a(k',:)
     swap l(k,1:k-1) and l(k',1:k-1)
     u(k,k) = a(k,k)
     for i = k+1 to n
       l(i,k) = a(i,k)/u(k,k)
       u(k,i) = a(k,i)
     for i = k+1 to n
       for j = k+1 to n
         a(i,j) = a(i,j) - l(i,k)*u(k,j)
          
   Here, the vector π is a compact representation of a permutation matrix p(n,n), 
   which is very sparse. For the ith row of p, π(i) stores the column index of
   the sole position that contains a 1.
```

One solution is developed using the Pthreads programming model and the other one using OpenMP. To check the answer, we also compute the sum of Euclidean norms of the columns of the residual matrix (this sum is known as the L2,1 norm) computed as PA-LU. This value should be very small.

You can find the key-setup, the program usage, and the different strategies incorporated in the **Report.pdf** document.

