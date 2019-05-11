// Parallel sample sort
#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include <stdlib.h>
#include <algorithm>

int main( int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  int rank, p;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  // Number of random numbers per processor (this should be increased
  // for actual tests or could be passed in through the command line
  int N = 10000;

  int* vec = (int*)malloc(N*sizeof(int));
  // seed random number generator differently on every core
  srand((unsigned int) (rank + 393919));

  // fill vector with random integers
  for (int i = 0; i < N; ++i) {
    vec[i] = rand();
  }
  printf("rank: %d, first entry: %d\n", rank, vec[0]);

  // timing
  MPI_Barrier(MPI_COMM_WORLD);
  double tt = MPI_Wtime();

  // sort locally
  std::sort(vec, vec+N);

  // sample p-1 entries from vector as the local splitters, i.e.,
  // every N/P-th entry of the sorted vector
  int *lsplitters, *gsplitters, *sdispls, *rdispls, *sendcnts, *recvcnts, *recvarr;

  lsplitters = (int*)malloc((p-1)*sizeof(int));
  sdispls = (int*)malloc(p*sizeof(int));
  rdispls = (int*)malloc(p*sizeof(int));
  sendcnts = (int*)malloc(p*sizeof(int));
  recvcnts = (int*)malloc(p*sizeof(int));

  for(int i=0;i<p-1;i++){
    lsplitters[i] = vec[(i+1)*N/p];
  }

  // every process communicates the selected entries to the root
  // process; use for instance an MPI_Gather
  if(rank==0){
    gsplitters = (int*)malloc(p*(p-1)*sizeof(int));
  }
  MPI_Gather(lsplitters, p-1, MPI_INT, gsplitters, p-1, MPI_INT, 0, MPI_COMM_WORLD);

  // root process does a sort and picks (p-1) splitters (from the
  // p(p-1) received elements)
  if(rank==0){
    std::sort(gsplitters, gsplitters+p*(p-1));
    for(int i=0; i<p-1; i++)
      lsplitters[i] = gsplitters[(i+1)*(p-1)];
  }

  // root process broadcasts splitters to all other processes
  MPI_Bcast(lsplitters, p-1, MPI_INT, 0, MPI_COMM_WORLD);

  // every process uses the obtained splitters to decide which
  // integers need to be sent to which other process (local bins).
  // Note that the vector is already locally sorted and so are the
  // splitters; therefore, we can use std::lower_bound function to
  // determine the bins efficiently.
  sdispls[0] = 0;
  for(int i=1; i<p; i++){
    sdispls[i] = std::lower_bound(vec,vec+N,lsplitters[i-1])-vec;
    sendcnts[i-1] = sdispls[i]-sdispls[i-1];
  }
  sendcnts[p-1] = N-sdispls[p-1];

  // Hint: the MPI_Alltoallv exchange in the next step requires
  // send-counts and send-displacements to each process. Determining the
  // bins for an already sorted array just means to determine these
  // counts and displacements. For a splitter s[i], the corresponding
  // send-displacement for the message to process (i+1) is then given by,
  // sdispls[i+1] = std::lower_bound(vec, vec+N, s[i]) - vec;

  // send and receive: first use an MPI_Alltoall to share with every
  // process how many integers it should expect, and then use
  // MPI_Alltoallv to exchange the data
  MPI_Alltoall(sendcnts,1,MPI_INT,recvcnts,1,MPI_INT,MPI_COMM_WORLD);
  int totalcnts = 0;
  for(int i=0; i<p; i++){
    totalcnts += recvcnts[i];
  }
  totalcnts+=recvcnts[p-1];

  recvarr = (int*)malloc(totalcnts*sizeof(int));

  int rS = 0;
  for(int i=0;i<p;i++){
    rdispls[i] = rS;
    rS += recvcnts[i];
  }

  MPI_Alltoallv(vec, sendcnts, sdispls, MPI_INT, recvarr, recvcnts, rdispls, MPI_INT, MPI_COMM_WORLD);

  // do a local sort of the received data
  std::sort(recvarr, recvarr+totalcnts);
  
  // timing
  MPI_Barrier(MPI_COMM_WORLD);
  double elapsedtime = MPI_Wtime() - tt;
  if (rank == 0) {
    printf("Time: %f seconds.\n", elapsedtime);
  }

  // every process writes its result to a file
  { // Write output to a file
    FILE* fd = NULL;
    char filename[256];
    snprintf(filename, 256, "output%02d.txt", rank);
    fd = fopen(filename,"w+");

    if(NULL == fd) {
      printf("Error opening file \n");
      return 1;
    }

    for(int n = 0; n < totalcnts; ++n)
      fprintf(fd, "%d\n", recvarr[n]);

    fclose(fd);
  }

  // free all memory allocation
  if(rank==0) free(gsplitters);
  free(lsplitters);
  free(sdispls);
  free(rdispls);
  free(sendcnts);
  free(recvcnts);
  free(recvarr);
  free(vec);
  MPI_Finalize();
  return 0;
}
