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
    unsigned long int i, m;
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

    unsigned long int cache_size    = 16*1024; // unit: Byte
    unsigned long int line_size     = 64;      // unit: Byte
    unsigned long int num_of_lines  = (unsigned long int)(cache_size/ line_size); //  256
    unsigned long int item_per_line = (unsigned long int)(line_size / sizeof(int));  //16 int
    unsigned long int item_per_cap  = (unsigned long int)(cache_size / sizeof(int));  //16 int
    unsigned long int iter = 0;
    unsigned long int cache_low;
    unsigned long int cache_high;
    int top_prime;
/* Tardis cache info
 *  L1 D$ size:      16KB -> (4K int)
 *  L1 D$ assoc:     4 
 *  L1 D$ LINE size: 64B -> (16 int) 
 *  # of L1 D$ lines: 256 lines
 * */
/*
 * line 0:  3,   5,   7,   9,  11,  13,  15,  17,  19,  21,  23,  25,  27,  29,  31,  33 (3, 5)
 * line 1: 35,  37,  39,  41,  43,  45,  47,  49,  51,  53,  55,  57,  59,  61,  63,  65 (3, 5, 7)
 * line 2: 67,  69,  71,  73,  75,  77,  79,  81,  83,  85,  87,  89,  91,  93,  95,  97 (3, 5, 7)
 * line 3: 99, 101, 103, 105, 107, 109, 111, 113, 115, 117, 119, 121, 123, 125, 127, 129 (3, 5, 7, 11)  
 *
 * line  : 2*0 +3,..., 2*(n-1)+3 
 *         2*n+3,..., 2*(2*n-1)+3
 *         2*(2*n)+3,..., 2*(3*n-1)+3
 *
 * line i: 2*(n*i)+3,...,2*(n*(i+1)-1)+3
 *      =  2*i*n+3,...,2*(i*n+n-1)+3 = 2*i*n+2*n+1, i: id of line within cap, n: id of item within cache line
 *      =  32*i+3,...,32*i+33, for 0 <= i < 256
 * line 255: 32*255+3,..., 32*255+33 (3, 5, 7, 11, 13, 17, 19,...,89)
 *        =  8163,...,8193  (sqrt(8193) = 90.5...)
 *
 * for kth cap: 
 *
 * */

    MPI_Init(&argc, &argv);

    MPI_Barrier(MPI_COMM_WORLD);
    /* Start the timer */

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    elapsed_time = -MPI_Wtime();

    if (argc != 3) {
        if (!id) printf("Command line: %s <m> <item_per_cap>\n", argv[0]);
        MPI_Finalize();
        exit(1);
    }

    n = atoll(argv[1]);
    item_per_cap = atoll(argv[2]);
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
// find mark_sqrt first
    prime = 3;
    do {
        for(i = ((prime*prime)-2) ; i < (int)(sqrt(n)) ; i+= prime){
            if((i+2)%2 == 1) marked_sqrt[i/2] = 1;
        }
        while (marked_sqrt[++index]);
        prime = (index * 2) + 3;
    } while (prime * prime <= n);
//   =  2*i*n+3,...,2*(i*n+n-1)+3 = 2*i*n+2*n+1, i: id of line within cap, n: id of item within cache line
//   k*(M)+2*i*n+3, ..., k*(M)+2*(i+1)*n+1 
    
    do{
        // for each cache cap section, know the highest prime we need to check    
        iter+=1;   
        cache_low  = low_value + (iter-1)*item_per_cap;//+2*(num_of_lines)*(item_per_line-1)+1;  
        cache_high = low_value + iter    *item_per_cap;//+2*(num_of_lines)*(item_per_line-1)+1;  
        top_prime  = (int)sqrt(cache_high); // the top prime of this cache cap at current iter
//printf("iter: %d, top_prime ceil: %d, item_per_cap: %d\n", iter, top_prime, item_per_cap);
        prime = 3;
        index = 0;
        do{        
            if (prime * prime > cache_low)
                first = prime * prime - cache_low;
            else {
                if (!(cache_low % prime)) first = 0;
                else first = prime - (cache_low % prime);
            }
            for (i = first+(iter-1)*item_per_cap; i < MIN(first+iter*item_per_cap, size); i += prime){
                if((i+cache_low)%2 == 1){
                    marked[i/2] = 1;
                }
            }
            while (marked_sqrt[++index]);
            prime = (index * 2) + 3;
        }while(prime <= top_prime);
    }while(iter*item_per_cap <= size);
//    do {
//        if (prime * prime > low_value)
//            first = prime * prime - low_value;
//        else {
//            if (!(low_value % prime)) first = 0;
//            else first = prime - (low_value % prime);
//        }
//        for (i = first; i < size; i += prime){
//            if((i+low_value)%2 == 1)  marked[i/2] = 1;
//        }
//        // each process do this part to avoid broadcast
//        for(i = ((prime*prime)-2) ; i < (int)(sqrt(n)) ; i+= prime){
//            if((i+2)%2 == 1) marked_sqrt[i/2] = 1;
//        }
//        //if (!id) {
//        
//        while (marked_sqrt[++index]);
//        prime = (index * 2) + 3;
//        //}
//        //MPI_Bcast(&prime, 1, MPI_INT, 0, MPI_COMM_WORLD);
//    } while (prime * prime <= n);
    count = 0;
    for (i = 0; i < size; i++)
        if ( ((i+low_value)%2 == 1) && !marked[i/2]) count++;
    //if (p > 1)
        MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM,
                   0, MPI_COMM_WORLD);

    /* Stop the timer */

    elapsed_time += MPI_Wtime();


    /* Print the results */
    if (!id) {
        printf("The total number of prime: %ld, total time: %10.6f, total node %d\n", global_count+1, elapsed_time, p);
    }
    MPI_Finalize();
    return 0;

}


