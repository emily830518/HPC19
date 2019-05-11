/* This is OpenMP implementations of the Jacobi method in 2D to solve linear system Au = f*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"
#ifdef _OPENMP
	#include <omp.h>
#endif

/* To index element (Row,Col) of a 2D array stored as 1D and indexed by row */
#define index(R, C, N)  ((R)*(N)) + (C)

// get nearby value(left, right, up, down) of k-th u at point (i,j)
void getNearbyValue(int i,int j, int N, double& left, double& right, double& up, double& down, double* u){
    // first row
    if(i==0){
        if(j==0){
            left = 0;
            right = u[index(i,j+1,N)];
        }
        else if (j==N-1) {
            left = u[index(i,j-1,N)];
            right = 0;
        }
        else{
            left = u[index(i,j-1,N)];
            right = u[index(i,j+1,N)];
        }
        down = 0;
        up = u[index(i+1,j,N)];
    }
    // last row
    if (i==N-1){
        if(j==0){
            left = 0;
            right = u[index(i,j+1,N)];
        }
        else if (j==N-1) {
            left = u[index(i,j-1,N)];
            right = 0;
        }
        else{
            left = u[index(i,j-1,N)];
            right = u[index(i,j+1,N)];
        }
        up = 0;
        down = u[index(i-1,j,N)];        
    }
    // middle row
    else{
        if(j==0){
            left = 0;
            right = u[index(i,j+1,N)];
        }
        else if (j==N-1) {
            left = u[index(i,j-1,N)];
            right = 0;
        }
        else{
            left = u[index(i,j-1,N)];
            right = u[index(i,j+1,N)];
        }
        up = u[index(i+1,j,N)]; 
        down = u[index(i-1,j,N)];
    }
}

void jacobi2D(double* f, double* u, int N, int max_iteration, int threads){
    double h = 1.0/(N+1.0);
    double* u_next = (double*) malloc(N * N * sizeof(double)); // (k+1)-th u
    memcpy(u_next, u, N * N * sizeof(double));
    omp_set_num_threads(threads);
    for(int k=0;k<max_iteration;k++){
        // foreach point (i,j)
        #pragma omp parallel for collapse(2)
        for(int i=0;i<N;i++){
            for(int j=0;j<N;j++){
                double left,right,up,down;
                getNearbyValue(i,j,N,left,right,up,down,u);
                u_next[index(i,j,N)] = (pow(h,2)*f[index(i,j,N)]+up+down+left+right)/4.0;
            }
        }
        memcpy(u, u_next, N * N * sizeof(double));
    }
    free(u_next);
}

int main(int argc, char * argv[]){
    int maxite=5000;
    int N;
    int num_threads;
    if(argc!=3){
        fprintf(stderr, "usage: jacobi2D-omp <N> <num_threads>\n");
        exit(1);
    }
    N = atoi(argv[1]);
    num_threads = atoi(argv[2]);

    double* u = (double*) malloc(N * N * sizeof(double)); // solution
    double* f = (double*) malloc(N * N * sizeof(double)); // RHS

    for(int i=0;i<N*N;i++){
        u[i] = 0.0; // initial guess
        f[i] = 1.0; // right hand side equals 1
    }
    printf("Max thread number: %d\n", omp_get_max_threads()); // => 64
    Timer t;
    t.tic();
    jacobi2D(f,u,N,maxite,num_threads);
    double time = t.toc();
    printf("\n");
    printf("Total time taken: %10f seconds\n", time);
    
    // for(int i=0;i<N*N;i++){
    //     printf("%lf\n", u[i]);
    // }

    free(u);
    free(f);
    return 0;
}