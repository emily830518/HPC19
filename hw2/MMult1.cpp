// g++ -fopenmp -O3 -march=native MMult1.cpp && ./a.out

#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "utils.h"

#define BLOCK_SIZE 32
/* To index element (Row,Col) of a 2D array stored as 1D */
#define index(R, C, N)  ((C)*(N)) + (R)

using namespace std;
// Note: matrices are stored in column major order; i.e. the array elements in
// the (m x n) matrix C are stored in the sequence: {C_00, C_10, ..., C_m0,
// C_01, C_11, ..., C_m1, C_02, ..., C_0n, C_1n, ..., C_mn}
void MMult0(long m, long n, long k, double *a, double *b, double *c) {
  for (long j = 0; j < n; j++) {
    for (long p = 0; p < k; p++) {
      for (long i = 0; i < m; i++) {
        double A_ip = a[i+p*m];
        double B_pj = b[p+j*k];
        double C_ij = c[i+j*m];
        C_ij = C_ij + A_ip * B_pj;
        c[i+j*m] = C_ij;
      }
    }
  }
}

// void MMultOMP(long m, long n, long k, double *a, double *b, double *c) {
//   long i,p,j;
//   #pragma omp for collapse(3) private(i, p, j)
//   for (long i = 0; i < m; i++) {
//     for (long p = 0; p < k; p++) {
//       for (long j = 0; j < n; j++) {
//         double A_ip = a[i+p*m];
//         double B_pj = b[p+j*k];
//         double C_ij = c[i+j*m];
//         #pragma omp atomic
//         C_ij = C_ij + A_ip * B_pj;
//         c[i+j*m] = C_ij;
//       }
//     }
//   }
// }

// rearranged loop
void MMult1(long m, long n, long k, double *a, double *b, double *c) {
  // TODO: See instructions below
  for (long i = 0; i < m; i++) {
    for (long p = 0; p < k; p++) {
      for (long j = 0; j < n; j++) {
        double A_ip = a[i+p*m];
        double B_pj = b[p+j*k];
        double C_ij = c[i+j*m];
        C_ij = C_ij + A_ip * B_pj;
        c[i+j*m] = C_ij;
      }
    }
  }  
}

// read block(i,j) from Matrix M and store in Blk
void readintoBlock(int i, int j, long dim, double *M, double *Blk){
  for(int x=0;x<BLOCK_SIZE;x++){
    for(int y=0;y<BLOCK_SIZE;y++){
      // convert the block (x,y) to matrix (row,col);
      int row = i*BLOCK_SIZE+y;
      int col = j*BLOCK_SIZE+x;
      // convert matrxix (row, col) to actual array index and store into Blk
      Blk[index(y,x,BLOCK_SIZE)] = M[index(row,col,dim)];
    }
  }
}

// write Blk back to Matrix M
void writeBlockback(int i, int j, long dim, double *Blk, double *M){
  for(int x=0;x<BLOCK_SIZE;x++){
    for(int y=0;y<BLOCK_SIZE;y++){
      // convert the block (x,y) to matrix (row,col);
      int row = i*BLOCK_SIZE+y;
      int col = j*BLOCK_SIZE+x;
      // convert matrxix (row, col) to actual array index and write Blk back to M
      M[index(row,col,dim)] = Blk[index(y,x,BLOCK_SIZE)];
    }
  }  
}

// Matrix Multiplication Blocking version
void MMultBlocking(int numblks, long dim, double *a, double *b, double *c){
  // Initialize block A, B, C
    double* Ablck_ik = (double*) aligned_malloc(BLOCK_SIZE * BLOCK_SIZE * sizeof(double));
    double* Bblck_kj = (double*) aligned_malloc(BLOCK_SIZE * BLOCK_SIZE * sizeof(double));
    double* Cblck_ij = (double*) aligned_malloc(BLOCK_SIZE * BLOCK_SIZE * sizeof(double));
    // Initialize all to 0
    for(long i=0;i<BLOCK_SIZE*BLOCK_SIZE;i++){
      Ablck_ik[i]=0;
      Bblck_kj[i]=0;
      Cblck_ij[i]=0;
    }
    for(long i=0;i<numblks;i++){
      for(long j=0;j<numblks;j++){
        readintoBlock(i,j,dim,c,Cblck_ij); // read block C(i,j) into fast memory
        for(long k=0;k<numblks;k++){
          readintoBlock(i,k,dim,a,Ablck_ik); // read block A(i,k) into fast memory
          readintoBlock(k,j,dim,b,Bblck_kj); // read block B(k,j) into fast memory
          MMult0(BLOCK_SIZE,BLOCK_SIZE,BLOCK_SIZE,Ablck_ik,Bblck_kj,Cblck_ij); // do a matrix multiply on blocks
        }
        writeBlockback(i,j,dim,Cblck_ij,c); // write block C(i,j) back to slow memory
      }
    }
    aligned_free(Ablck_ik);
    aligned_free(Bblck_kj);
    aligned_free(Cblck_ij);  
}

// Matrix Multiplication OpenMP version
void MMultOpenMP(int numblks, long dim, double *a, double *b, double *c){
  // Initialize block A, B, C
  #pragma omp parallel
  {
    double* Ablck_ik = (double*) aligned_malloc(BLOCK_SIZE * BLOCK_SIZE * sizeof(double));
    double* Bblck_kj = (double*) aligned_malloc(BLOCK_SIZE * BLOCK_SIZE * sizeof(double));
    double* Cblck_ij = (double*) aligned_malloc(BLOCK_SIZE * BLOCK_SIZE * sizeof(double));
    // Initialize all to 0
    for(long i=0;i<BLOCK_SIZE*BLOCK_SIZE;i++){
      Ablck_ik[i]=0;
      Bblck_kj[i]=0;
      Cblck_ij[i]=0;
    }
    long i,j,k;
    #pragma omp for collapse(2) private(i,j,k)
    for(i=0;i<numblks;i++){
      for(j=0;j<numblks;j++){
        readintoBlock(i,j,dim,c,Cblck_ij); // read block C(i,j) into fast memory
        for(k=0;k<numblks;k++){
          readintoBlock(i,k,dim,a,Ablck_ik); // read block A(i,k) into fast memory
          readintoBlock(k,j,dim,b,Bblck_kj); // read block B(k,j) into fast memory
          MMult0(BLOCK_SIZE,BLOCK_SIZE,BLOCK_SIZE,Ablck_ik,Bblck_kj,Cblck_ij); // do a matrix multiply on blocks
        }
        writeBlockback(i,j,dim,Cblck_ij,c); // write block C(i,j) back to slow memory
      }
    }
    aligned_free(Ablck_ik);
    aligned_free(Bblck_kj);
    aligned_free(Cblck_ij);  
  }
}


int main(int argc, char** argv) {
  const long PFIRST = BLOCK_SIZE;
  const long PLAST = 2000;
  const long PINC = std::max(50/BLOCK_SIZE,1) * BLOCK_SIZE; // multiple of BLOCK_SIZE

  printf(" Dimension       Time    Gflop/s       GB/s        Error\n");
  for (long p = PFIRST; p < PLAST; p += PINC) {
    long m = p, n = p, k = p;
    // long NREPEATS = 1e9/(m*n*k)+1;
    long NREPEATS = 10;
    double* a = (double*) aligned_malloc(m * k * sizeof(double)); // m x k
    double* b = (double*) aligned_malloc(k * n * sizeof(double)); // k x n
    double* c = (double*) aligned_malloc(m * n * sizeof(double)); // m x n
    double* c_ref = (double*) aligned_malloc(m * n * sizeof(double)); // m x n

    // Initialize matrices
    for (long i = 0; i < m*k; i++) a[i] = drand48();
    for (long i = 0; i < k*n; i++) b[i] = drand48();
    for (long i = 0; i < m*n; i++) c_ref[i] = 0;
    for (long i = 0; i < m*n; i++) c[i] = 0;

    for (long rep = 0; rep < NREPEATS; rep++) { // Compute reference solution
      MMult0(m, n, k, a, b, c_ref);
    }
    int numblks = m/BLOCK_SIZE;
    Timer t;
    t.tic();
    for (long rep = 0; rep < NREPEATS; rep++) {
      // MMult0(m, n, k, a, b, c);
      // MMultBlocking(numblks, m, a, b, c);
      MMultOpenMP(numblks, m, a, b, c);
    }
    // for(int i=0;i<m*n;i++){
    //   if(c[i]!=c_ref[i])
    //     cout<<"NOT EQUAL"<<endl; break;
    // }
    // cout<<"EQUAL"<<endl;
    double time = t.toc();
    double flops = ((NREPEATS*2.0*m*n*k)/1e9)/time; // TODO: calculate from m, n, k, NREPEATS, time
    double bandwidth = ((NREPEATS*4.0*m*n*k*sizeof(double))/1e9)/time; // TODO: calculate from m, n, k, NREPEATS, time
    printf("%10d %10f %10f %10f", p, time, flops, bandwidth);

    double max_err = 0;
    for (long i = 0; i < m*n; i++) max_err = std::max(max_err, fabs(c[i] - c_ref[i]));
    printf(" %10e\n", max_err);


    aligned_free(a);
    aligned_free(b);
    aligned_free(c);
    aligned_free(c_ref);
  }

  return 0;
}

// * Using MMult0 as a reference, implement MMult1 and try to rearrange loops to
// maximize performance. Measure performance for different loop arrangements and
// try to reason why you get the best performance for a particular order?
//
//
// * You will notice that the performance degrades for larger matrix sizes that
// do not fit in the cache. To improve the performance for larger matrices,
// implement a one level blocking scheme by using BLOCK_SIZE macro as the block
// size. By partitioning big matrices into smaller blocks that fit in the cache
// and multiplying these blocks together at a time, we can reduce the number of
// accesses to main memory. This resolves the main memory bandwidth bottleneck
// for large matrices and improves performance.
//
// NOTE: You can assume that the matrix dimensions are multiples of BLOCK_SIZE.
//
//
// * Experiment with different values for BLOCK_SIZE (use multiples of 4) and
// measure performance.  What is the optimal value for BLOCK_SIZE?
//
//
// * Now parallelize your matrix-matrix multiplication code using OpenMP.
//
//
// * What percentage of the peak FLOP-rate do you achieve with your code?
//
//
// NOTE: Compile your code using the flag -march=native. This tells the compiler
// to generate the best output using the instruction set supported by your CPU
// architecture. Also, try using either of -O2 or -O3 optimization level flags.
// Be aware that -O2 can sometimes generate better output than using -O3 for
// programmer optimized code.
