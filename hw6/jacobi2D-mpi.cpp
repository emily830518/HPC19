/* This is MPI implementations of the Jacobi method in 2D to solve linear system Au = f*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

/* To index element (Row,Col) of a 2D array stored as 1D and indexed by row */
#define index(R, C, N)  ((R)*(N)) + (C)

// get nearby value(left, right, up, down) of k-th u at point (i,j)
void getNearbyValue(int i,int j, int N, double& left, double& right, double& up, double& down, double* u){
    left = u[index(i,j-1,N)];
    right = u[index(i,j+1,N)];
    up = u[index(i+1,j,N)];
    down = u[index(i-1,j,N)];
}

void top(int i, int j, double* lu, double* buffer, int Nl, int rank, int perrow){
    if(i>0){
        int send = rank;
        int receive = index(i-1,j,perrow);
        for(int k = 0; k<Nl+2; k++)
            buffer[k] = lu[Nl+2+k];
        MPI_Send(buffer, Nl+2, MPI_DOUBLE, receive, 0, MPI_COMM_WORLD);
    }
    if(i<perrow-1){
        int receive = rank;
        int send = index(i+1,j,perrow);  
        MPI_Recv(buffer,Nl+2,MPI_DOUBLE,send,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        for(int k = 0; k<Nl+2; k++)
            lu[(Nl+1)*(Nl+2)+k] = buffer[k];
    }
}

void bottom(int i, int j, double* lu, double* buffer, int Nl, int rank, int perrow){
    if(i<perrow-1){
        int send = rank;
        int receive = index(i+1,j,perrow);
        for(int k=0;k<Nl+2;k++)
            buffer[k] = lu[Nl*(Nl+2)+k];
        MPI_Send(buffer,Nl+2,MPI_DOUBLE,receive,0,MPI_COMM_WORLD);
    }
    if(i>0){
        int receive = rank;
        int send = index(i-1,j,perrow);
        MPI_Recv(buffer,Nl+2,MPI_DOUBLE,send,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        for(int k=0;k<Nl+2;k++)
            lu[k] = buffer[k];
    }
}

void left(int i, int j, double* lu, double* buffer, int Nl, int rank, int perrow){
    if(j<perrow-1){
        int send = rank;
        int receive = index(i,j+1,perrow);
        for(int k=0;k<Nl+2;k++)
            buffer[k] = lu[(k+1)*(Nl+2)-2];
        MPI_Send(buffer,Nl+2,MPI_DOUBLE,receive,0,MPI_COMM_WORLD);
    }
    if(j>0){
        int receive = rank;
        int send = index(i,j-1,perrow);
        MPI_Recv(buffer,Nl+2,MPI_DOUBLE,send,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        for(int k=0;k<Nl+2;k++)
            lu[k*(Nl+2)] = buffer[k];
    }
}

void right(int i, int j, double* lu, double* buffer, int Nl, int rank, int perrow){
    if(j>0){
        int send = rank;
        int receive = index(i,j-1,perrow);
        for(int k=0;k<Nl+2;k++)
            buffer[k] = lu[k*(Nl+2)+1];
        MPI_Send(buffer,Nl+2,MPI_DOUBLE,receive,0,MPI_COMM_WORLD);
    }
    if(j<perrow-1){
        int receive = rank;
        int send = index(i,j+1,perrow);
        MPI_Recv(buffer,Nl+2,MPI_DOUBLE,send,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        for(int k=0;k<Nl+2;k++)
            lu[(k+1)*(Nl+2)-1] = buffer[k];
    }
}

void communicateGhost(double* lu, double* buffer, int Nl, int rank, int perrow){
    int pi, pj;
    pi = rank/perrow;
    pj = rank%perrow;
    top(pi,pj,lu,buffer,Nl,rank,perrow);
    bottom(pi,pj,lu,buffer,Nl,rank,perrow);
    left(pi,pj,lu,buffer,Nl,rank,perrow);
    right(pi,pj,lu,buffer,Nl,rank,perrow);
}

void jacobi2D(double* lf, double* lu, int Nl, int N, int max_iteration, int rank, int pnum){
    double h = 1.0/(N+1.0);
    double* lu_next = (double*) malloc((Nl+2) * (Nl+2) * sizeof(double)); // (k+1)-th u
    memcpy(lu_next, lu, (Nl+2) * (Nl+2) * sizeof(double));
    double* buffer = (double*) malloc((Nl+2)*sizeof(double));
    int perrow = sqrt(pnum+0.001);
    for(int k=0;k<max_iteration;k++){
        // foreach point (i,j)
        for(int i=0;i<(Nl+2);i++){
            for(int j=0;j<(Nl+2);j++){
                if (i>0 && i<Nl+1 && j>0 && j<Nl+1){
                    double left,right,up,down;
                    getNearbyValue(i,j,Nl+2,left,right,up,down,lu);
                    lu_next[index(i,j,Nl+2)] = (pow(h,2)*lf[index(i,j,Nl+2)]+up+down+left+right)/4.0;
                }
            }
        }
        communicateGhost(lu_next,buffer,Nl,rank,perrow);
        memcpy(lu, lu_next, (Nl+2) * (Nl+2) * sizeof(double));
    }
    free(lu_next);
    free(buffer);
}

int main(int argc, char * argv[]){
    int maxite = 500;
    int N;
    int rank, pnum;
    if(argc!=2){
        fprintf(stderr, "usage: jacobi2D-mpi <N>\n");
        exit(1);
    }
    N = atoi(argv[1]);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &pnum);

    int j = (int)(log(pnum)/log(4));
    if(pow(4,j)!=pnum){
        fprintf(stderr,"process number is not divisible by 4^j!");
        MPI_Abort(MPI_COMM_WORLD, 0);
    }
    int Nl = N/pow(2,j);
    if(pow(2,j)*Nl != N){
        fprintf(stderr,"N is not divisible by 2^j!");
        MPI_Abort(MPI_COMM_WORLD, 0);        
    }
    if (rank == 0) printf("Nl is %d\n", Nl);

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    printf("Rank %d/%d running on %s.\n", rank, pnum, processor_name);

    double* lu = (double*) malloc((Nl+2) * (Nl+2) * sizeof(double)); // solution
    double* lf = (double*) malloc((Nl+2) * (Nl+2) * sizeof(double)); // RHS

    for(int i=0;i<(Nl+2)*(Nl+2);i++){
        lu[i] = 0.0; // initial guess
        lf[i] = 1.0; // right hand side equals 1
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double tt = MPI_Wtime();

    jacobi2D(lf,lu,Nl,N,maxite,rank,pnum);

    MPI_Barrier(MPI_COMM_WORLD);
    double elapsedtime = MPI_Wtime()-tt;
    if(rank==0) printf("Total time taken: %10f seconds\n", elapsedtime);
    
    free(lu);
    free(lf);
    MPI_Finalize();
    return 0;
}