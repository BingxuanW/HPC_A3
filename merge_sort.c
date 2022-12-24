#include <stdio.h>      // Printf
#include <time.h>       // Timer
#include <math.h>       // Logarithm
#include <stdlib.h>     // Malloc
#include "mpi.h"        // MPI Library
#include <string.h>
#include <random>

void merge(int *, int *, int, int, int);
void mergeSort(int *, int *, int, int);

int main(int argc, char** argv[]) {
    /********** Initialize MPI **********/
    int world_rank;
    int world_size;

    MPI_INIT(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    /* ------------------ Generate random my_array[] ------------------*/
    /* N, base_seed will be passed as command line arguments */
    /********** Divide the array in equal-sized chunks **********/
    int N = (int)atoi(argv[1])/world_size;              // 1st argument after ./sort.cc, size of my_array[]
    uint32_t base_seed = (uint32_t)atoi(argv[2]);      // 2nd argument after ./sort.cc
    /********** Send each subarray to each process **********/
    uint32_t* my_array = (uint32_t *) malloc(sizeof(uint32_t) * N);
    /* Initialize the random number generator for the given base_seed
    * plus an offset for the MPI rank of the node, such that on every
    * node different numbers are generated.
    */
    std::mt19937 generator(base_seed + mpi_rank);
    /* Generate N pseudo-random uint32_t’s */
    std::uniform_int_distribution<uint32_t> distribution(0, std::numeric_limits<uint32_t>::max());
    for (int idx = 0; idx < N; idx++)
    {
        my_array[idx] = distribution(generator);
    }

    /********** Divide the array in equal-sized chunks **********/
    // int size = n/world_size;
    /********** Send each subarray to each process **********/
    // int *sub_array = malloc(size * sizeof(int));
    // MPI_Scatter(original_array, size, MPI_INT, sub_array, size, MPI_INT, 0, MPI_COMM_WORLD);
    
    /********** Perform the mergesort on each process **********/
    int *tmp_array = (uint32_t *) malloc(sizeof(uint32_t) * N);
    mergeSort(my_array, tmp_array, 0, (N - 1));
    
    /********** Gather the sorted subarrays into one **********/
    int *sorted = NULL;
    if(world_rank == 0) {
        
        sorted = (uint32_t *) malloc(sizeof(uint32_t) * N * world_size);

        }
    
    MPI_Gather(MY_array, N, MPI_INT, sorted, N, MPI_INT, 0, MPI_COMM_WORLD);
    
    /********** Make the final mergeSort call **********/
    if(world_rank == 0) {
        
        int *other_array = (uint32_t *) malloc(sizeof(uint32_t) * N * world_size);
        mergeSort(sorted, other_array, 0, ( N * world_size - 1));
        
        /********** Display the sorted array **********/
        printf("This is the sorted array: ");
        for(int i = 0; i < N * world_size; i++) {
            
            printf("%d ", sorted[i]);
            
            }
            
        printf("\n");
        printf("\n");
            
        /********** Clean up root **********/
        free(sorted);
        free(other_array);
            
        }
    
    /********** Clean up rest **********/
    free(original_array);
    free(sub_array);
    free(tmp_array);
    
    /********** Finalize MPI **********/
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    
    }

/********** Merge Function **********/
void merge(int *a, int *b, int l, int m, int r) {
    
    int h, i, j, k;
    h = l;
    i = l;
    j = m + 1;
    
    while((h <= m) && (j <= r)) {
        
        if(a[h] <= a[j]) {
            
            b[i] = a[h];
            h++;
            
            }
            
        else {
            
            b[i] = a[j];
            j++;
            
            }
            
        i++;
        
        }
        
    if(m < h) {
        
        for(k = j; k <= r; k++) {
            
            b[i] = a[k];
            i++;
            
            }
            
        }
        
    else {
        
        for(k = h; k <= m; k++) {
            
            b[i] = a[k];
            i++;
            
            }
            
        }
        
    for(k = l; k <= r; k++) {
        
        a[k] = b[k];
        
        }
        
    }

/********** Recursive Merge Function **********/
void mergeSort(int *a, int *b, int l, int r) {
    
    int m;
    
    if(l < r) {
        
        m = (l + r)/2;
        
        mergeSort(a, b, l, m);
        mergeSort(a, b, (m + 1), r);
        merge(a, b, l, m, r);
        
        }
        
    }
