/******************************************************************************
* FILE: omp_bug5.c
* DESCRIPTION:
*   Using SECTIONS, two threads initialize their own array and then add
*   it to the other's array, however a deadlock occurs.
* AUTHOR: Blaise Barney  01/29/04
* LAST REVISED: 04/06/05
******************************************************************************/
/******************************************************************************
* Comment:
* The deadlock happens because in section one, it tries to lock a first and then lock b.
* In section two, it tries to lock b first adn then lock a. 
* This is the factor to cause deadlock, because when the thread(running section 1) locks a and the other thread(running section 2) locks b.
* Then, the thread(running section 1) wants to lock b, but b is locked by the other thread(running on section 2), so the thread waiting for b to unlock.
* Similarly, the thread(running section 2) wants to lock a, but a is locked by the other thread(runnign on section 1), so the thread waiting for a to unlock.
* This is a deadlock. To solve this, we just need to adjust the order of lock in section 2.
******************************************************************************/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#define N 1000000
#define PI 3.1415926535
#define DELTA .01415926535

int main (int argc, char *argv[]) 
{
int nthreads, tid, i;
float a[N], b[N];
omp_lock_t locka, lockb;

/* Initialize the locks */
omp_init_lock(&locka);
omp_init_lock(&lockb);

/* Fork a team of threads giving them their own copies of variables */
#pragma omp parallel shared(a, b, nthreads, locka, lockb) private(tid)
  {

  /* Obtain thread number and number of threads */
  tid = omp_get_thread_num();
  #pragma omp master
    {
    nthreads = omp_get_num_threads();
    printf("Number of threads = %d\n", nthreads);
    }
  printf("Thread %d starting...\n", tid);
  #pragma omp barrier

  #pragma omp sections nowait
    {
    #pragma omp section
      {
      printf("Thread %d initializing a[]\n",tid);
      omp_set_lock(&locka);
      for (i=0; i<N; i++)
        a[i] = i * DELTA;
      omp_set_lock(&lockb);
      printf("Thread %d adding a[] to b[]\n",tid);
      for (i=0; i<N; i++)
        b[i] += a[i];
      omp_unset_lock(&lockb);
      omp_unset_lock(&locka);
      }

    #pragma omp section
      {
      printf("Thread %d initializing b[]\n",tid);
      // omp_set_lock(&lockb);
      omp_set_lock(&locka); // lock a first as section 1
      for (i=0; i<N; i++)
        b[i] = i * PI;
      // omp_set_lock(&locka);
      omp_set_lock(&lockb); // then, lock b
      printf("Thread %d adding b[] to a[]\n",tid);
      for (i=0; i<N; i++)
        a[i] += b[i];
      // omp_unset_lock(&locka);
      // omp_unset_lock(&lockb);
      omp_unset_lock(&lockb); // to maintain proper lock ordering, unlocking in reverse order. So release b first
      omp_unset_lock(&locka); // then, release a
      }
    }  /* end of sections */
  }  /* end of parallel region */

}

