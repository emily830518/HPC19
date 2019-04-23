#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define N 100000

int main(int argc, char** argv) {
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);
  // Find out rank, size
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  int token=0; //initialize token as 0 at first
  // Receive from the lower process and send to the higher process. Take care
  // of the special case when you are the first process to prevent deadlock.
  MPI_Barrier(MPI_COMM_WORLD);
  double tt = MPI_Wtime();
  for(int n=0;n<N;n++){
    if (world_rank != 0) {
      MPI_Recv(&token, 1, MPI_INT, world_rank - 1, 0, MPI_COMM_WORLD,
              MPI_STATUS_IGNORE);
      if(n==N-1){
        printf("Process %d received token %d from process %d\n", world_rank, token,
            world_rank - 1);
      }
    } 
    token = token + world_rank; // add rank to token before sending
    MPI_Send(&token, 1, MPI_INT, (world_rank + 1) % world_size, 0,
            MPI_COMM_WORLD);
    // Now process 0 can receive from the last process. This makes sure that at
    // least one MPI_Send is initialized before all MPI_Recvs (again, to prevent
    // deadlock)
    if (world_rank == 0) {
      MPI_Recv(&token, 1, MPI_INT, world_size - 1, 0, MPI_COMM_WORLD,
              MPI_STATUS_IGNORE);
      if(n==N-1){
        printf("Process %d received token %d from process %d\n", world_rank, token,
            world_size - 1);
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  tt = MPI_Wtime() - tt;
  if (!world_rank) printf("int_ring latency: %e ms\n", tt/N * 1000);
  MPI_Finalize();
}