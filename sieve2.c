/*
 *   Sieve of Eratosthenes
 *
 *   Programmed by Michael J. Quinn
 *
 *   Last modification: 7 September 2001
 */

#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MIN(a, b)  ((a)<(b)?(a):(b))

#define BLOCK_LOW(id, p, n)  ((id)*(n)/(p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW((id)+1, p, n) -1)
#define BLOCK_SIZE(id, p, n) (BLOCK_LOW((id)+1, p, n) - BLOCK_LOW((id), p, n))

int main(int argc, char *argv[]) {
    unsigned long int count;        /* Local prime count */
    double elapsed_time; /* Parallel execution time */
    unsigned long int first;        /* Index of first multiple */
    unsigned long int global_count = 0; /* Global prime count */
    unsigned long long int high_value;   /* Highest value on this proc */
    unsigned long int i;
    int id;           /* Process ID number */
    unsigned long int index;        /* Index of current prime */
    unsigned long long int low_value;    /* Lowest value on this proc */
    char *marked;       /* Portion of 2,...,'n' */
    char *marked_sqrt;  /* Portion of 2,...,'sqrt(n)' */
    unsigned long long int n;            /* Sieving from 2, ..., 'n' */
    int p;            /* Number of processes */
    unsigned long int proc0_size;   /* Size of proc 0's subarray */
    unsigned long int prime;        /* Current prime */
    unsigned long int size;         /* Elements in 'marked' */

    MPI_Init(&argc, &argv);

    MPI_Barrier(MPI_COMM_WORLD);
    /* Start the timer */

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    elapsed_time = -MPI_Wtime();

    if (argc != 2) {
        if (!id) printf("Command line: %s <m>\n", argv[0]);
        MPI_Finalize();
        exit(1);
    }

    n = atoll(argv[1]);

    /* Figure out this process's share of the array, as
       well as the integers represented by the first and
       last array elements */

    low_value = 2 + id * (n - 1) / p;
    high_value = 1 + (id + 1) * (n - 1) / p;
//    size = high_value - low_value + 1;
    size = BLOCK_SIZE(id, p, n-1); // local size per process
    /* Bail out if all the primes used for sieving are
       not all held by process 0 */

    proc0_size = (n - 1) / p;

    if ((2 + proc0_size) < (int) sqrt((double) n)) {
        if (!id) printf("Too many processes\n");
        MPI_Finalize();
        exit(1);
    }

    /* Allocate this process's share of the array. */

    marked      = (char *) malloc(size/2);
    // each process need to do this but avoid broadcast
    marked_sqrt = (char *) malloc(sqrt(n));

    if (marked == NULL || marked_sqrt == NULL) {
        printf("Cannot allocate enough memory\n");
        MPI_Finalize();
        exit(1);
    }

    for (i = 0; i < size/2;  i++)        marked[i]      = 0;
    for (i = 0; i < (int)(sqrt(n)); i++) marked_sqrt[i] = 0;
    
    //if (!id) index = 0;
    index = 0;
    prime = 3;
    do {
        if (prime * prime > low_value)
            first = prime * prime - low_value;
        else {
            if (!(low_value % prime)) first = 0;
            else first = prime - (low_value % prime);
        }
        for (i = first; i < size; i += prime){
            if((i+low_value)%2 == 1)  marked[i/2] = 1;
        }
        // each process do this part to avoid broadcast
        for(i = (prime*prime-2) ; i < (int)(sqrt(n)) ; i+= prime){
            if((prime*prime)%2 == 1) marked_sqrt[i/2] = 1;
        }
        //if (!id) {
        while (marked_sqrt[++index]);
        prime = (index * 2) + 3;
        //}
        //MPI_Bcast(&prime, 1, MPI_INT, 0, MPI_COMM_WORLD);
    } while (prime * prime <= n);
    count = 0;
    for (i = 0; i < size; i++)
        if ( ((i+low_value)%2 == 1) && !marked[i/2]) count++;
    //if (p > 1)
        MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM,
                   0, MPI_COMM_WORLD);

    /* Stop the timer */

    elapsed_time += MPI_Wtime();


    /* Print the results */
    for(i = 0 ; i < size ; i++){
      if(!marked[i]) printf("prime: %d\n", 2*i+3);
    }
    if (!id) {
        printf("The total number of prime: %ld, total time: %10.6f, total node %d\n", global_count+1, elapsed_time, p);
    }
    MPI_Finalize();
    return 0;

}


