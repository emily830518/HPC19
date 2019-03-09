/******************************************************************************
* FILE: omp_bug2.c
* DESCRIPTION:
*   Another OpenMP program with a bug. 
* AUTHOR: Blaise Barney 
* LAST REVISED: 04/06/05 
******************************************************************************/
/******************************************************************************
* Comment:
*   There are several bugs in this program.
*   (1) tid and i declaration should be moved from top into the section of #pragma omp parallel.
*   (2) modifying 'float total' to 'long int total' to solve the float point calculation problem.
*   (3) There is a race condition, because two threads try to modify the memory address of total variable at the same time.
*       To solve this, I have to add #pragma omp atomic to make it to be an atomic operation.
******************************************************************************/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[]) 
{
int nthreads; // removing tid and i decalaration 
long int total=0; // modifying from float to long int.

/*** Spawn parallel region ***/
#pragma omp parallel 
  {
  /* Obtain thread number */
  int tid = omp_get_thread_num(); // declare tid here so that tid will become private in each thread.
  /* Only master thread does this */
  if (tid == 0) {
    nthreads = omp_get_num_threads();
    printf("Number of threads = %d\n", nthreads);
    }
  printf("Thread %d is starting...\n",tid);

  #pragma omp barrier

  /* do some work */
  #pragma omp for schedule(dynamic,10)
  for (int i=0; i<1000000; i++){ // declare i here so that i will become private in each thread
    #pragma omp atomic // add this line to solve race condition, make 'total = total + i*1.0;' an atomic operation
     total = total + i*1.0;
  }

  printf ("Thread %d is done! Total= %ld\n",tid,total); // modifying from %e to %ld to match long int

  } /*** End of parallel region ***/
}