#include <stdio.h>      // Printf
#include <time.h>       // Timer
#include <math.h>       // Logarithm
#include <stdlib.h>     // Malloc
#include "mpi.h"        // MPI Library
#include <string.h>
#include <random>

#define MASTER 0        // Who should do the final processing?
#define OUTPUT_NUM 10   // Number of elements to display in output

/* Variables */
double timer_start;
double timer_end;
// uint32_t process_rank;
// uint32_t num_processes;
// uint32_t *array;
// uint32_t array_size;

/* Functions */ 
void CompareLow(int bit, int N, int mpi_rank, uint32_t my_array[]);
void CompareHigh(int bit, int N, int mpi_rank, uint32_t my_array[]);
int ComparisonFunc(const void * a, const void * b);

unsigned int Log2n(unsigned int n)
{
	return (n>1) ? 1+Log2n(n/2) : 0;
}

int main(int argc, char * argv[]) {
    /* ------------------ For MPI ------------------*/
    int mpi_rank;
    int P;  // number of processors

    // MPI_Init(&argc, &argv);
    // MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
    // MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    
    /* ------------------ Generate random my_array[] ------------------*/
    /* N, base_seed will be passed as command line arguments */
    int N = (int)atoi(argv[1]);              // 1st argument after ./sort.cc, size of my_array[]
    uint32_t base_seed = (uint32_t)atoi(argv[2]);      // 2nd argument after ./sort.cc
    // uint32_t base_seed = strtoul(atoi(argv[2], NULL, 10); // 该函数待研究
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

    // array_size = atoi(argv[1]) / num_processes;
    // array = (uint32_t *) malloc(array_size * sizeof(uint32_t));

    // srand(time(NULL)+process_rank*num_processes);  
    // for (i = 0; i < array_size; i++) {
    //     array[i] = rand() % (atoi(argv[1]));
    // }

    MPI_Barrier(MPI_COMM_WORLD);                // wait for everything to be ready

    int dimensions = Log2n(P);

    if (mpi_rank == MASTER)
    {
        printf("Number of Processes spawned: %d\n", P);
        timer_start = MPI_Wtime();
    }

    qsort(my_array, N, sizeof(uint32_t), ComparisonFunc);

   // printf("My rank is %d \t", process_rank);
   // for(i=0;i<array_size;i++)
	//{
	//	printf("%d\t", array[i]);
		//uint32_t s=0;	
		//for(s=0;s<=process_rank;s++)
		//	printf(".");
//	}
	
   // MPI_Barrier(MPI_COMM_WORLD);
    int i, j;
    for (i = 0; i <dimensions; i++)
    {
        for (j = i; j>=0; j--) {
		//MPI_Barrier(MPI_COMM_WORLD);
		
            if (( (mpi_rank >> (i + 1)) % 2 == 0 && 
                  (mpi_rank >> j) % 2 == 0) || ((mpi_rank >> (i + 1)) % 2 != 0 && 
                  (mpi_rank >> j) % 2 != 0 ))
            {
                CompareLow(j, N, mpi_rank, my_array);
            }
            else
            {
                CompareHigh(j, N, mpi_rank, my_array);
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (mpi_rank == MASTER)
    {
        timer_end = MPI_Wtime();
        printf("Time Elapsed (Sec): %f\n", timer_end - timer_start);
    }

      //printf("Displaying sorted array (only 10 elements for quick verification), My rank is %d\n", process_rank);

    for (i = 0; i < N; i++) {
        printf("My rank is %d, no is %d\n",process_rank, array[i]);
    	for(j = 0; j <= process_rank; j++)
    		printf(".");
           }
     }
    printf("\n\n");

    free(my_array);

    MPI_Finalize();
    return 0;
}


int ComparisonFunc(const void * a, const void * b) {
    return ( * (uint32_t *)a - * (uint32_t *)b );
}


void CompareLow(int j, int N, int mpi_rank, uint32_t my_array[]) {
    int i;
    uint32_t min;

   // printf("My rank is %d Pairing with %d in CL\n", process_rank, process_rank^(1<<j));
    
    int send_counter = 0;
    uint32_t *buffer_send = (uint32_t *) malloc( (N + 1) * sizeof(uint32_t) );
   // printf("Trying to send local max in CL:%d\n", array[array_size-1]);
    MPI_Send(
        &my_array[N - 1],     
        1,                          
        MPI_UINT32_T,                    
        mpi_rank ^ (1 << j),    
        0,
        MPI_COMM_WORLD              
    );

    int recv_counter;
    uint32_t *buffer_recieve = (uint32_t *) malloc( (N + 1) * sizeof(uint32_t) );
    MPI_Recv(
        &min,                       
        1,                          
        MPI_UINT32_T,                    
        mpi_rank ^ (1 << j),    
        0,                          
        MPI_COMM_WORLD,             
        MPI_STATUS_IGNORE           
    );
    
    for (i = N-1; i >= 0; i--)
    {
        if ( my_array[i] > min)
        {
	        send_counter++;
            buffer_send[send_counter] = my_array[i];        
        }
        else
        {
            break;      
        }
    }

    buffer_send[0] = (uint32_t)send_counter;

    MPI_Send(
        buffer_send,                
        send_counter+1,               
        MPI_UINT32_T,                    
        mpi_rank ^ (1 << j),    
        0,                          
        MPI_COMM_WORLD              
    );

    MPI_Recv(
        buffer_recieve,             
        N+1,                 
        MPI_UINT32_T,                    
        mpi_rank ^ (1 << j),    
        0,                          
        MPI_COMM_WORLD,             
        MPI_STATUS_IGNORE           
    );

    uint32_t *temp_array = (uint32_t *)malloc( N * sizeof(uint32_t) );
    for( i=0 ; i < N ; i++ )
    {
	    temp_array[i]= my_array[i];
    }
   
    int buffer_size = (int)buffer_recieve[0] ;
    int k = 1 ;
    int m = 0 ;   

 
    k = 1 ;
    for ( i = 0; i < N ; i++ )
    {
        if( temp_array[m] <= buffer_recieve[k] )
        {
            my_array[i] = temp_array[m] ;
            m++ ;
        }
        else if( k <= buffer_size )
        {
            my_array[i] = buffer_recieve[k] ;
            k++ ;
        }
    }

    qsort( my_array, N, sizeof(uint32_t), ComparisonFunc );

    free(buffer_send);
    free(buffer_recieve); 
    return;
}


void CompareHigh(int j, int N, int mpi_rank, uint32_t my_array[]) {
    int i;
    uint32_t max;
    int recv_counter;
    uint32_t *buffer_recieve = (uint32_t *) malloc( (N + 1) * sizeof(uint32_t) );
    MPI_Recv(
        &max,                       
        1,                          
        MPI_UINT32_T,                    
        mpi_rank ^ (1 << j),    
        0,                          
        MPI_COMM_WORLD,             
        MPI_STATUS_IGNORE           
    );

   // printf("Received max from pair in CH:%d\n",max);
    int send_counter = 0;
    uint32_t *buffer_send = (uint32_t *) malloc(( N + 1) * sizeof(uint32_t));
    MPI_Send(
        &my_array[0],                  
        1,                          
        MPI_UINT32_T,                    
        mpi_rank ^ (1 << j),    
        0,                          
        MPI_COMM_WORLD              
    );

   // printf("Sending min to my pair from CH:%d\n", array[0]);
    for (i = 0; i < N ; i++)
    {
        if ( my_array[i] < max)
        {
	   // printf("Buffer sending in CH: %d\n", array[i]);
            send_counter++;
		    buffer_send[send_counter] = my_array[i];
        }
        else
        {
            break;      
        }
    }
    
    MPI_Recv(
        buffer_recieve,             
        N+1,                 
        MPI_UINT32_T,                    
        mpi_rank ^ (1 << j),    
        0,                          
        MPI_COMM_WORLD,             
        MPI_STATUS_IGNORE           
    );
    recv_counter = (int)buffer_recieve[0];
    buffer_send[0] = (uint32_t)send_counter;
    //printf("Send counter in CH: %d\n", send_counter);

  //  for(i=0;i<=send_counter;i++)
//	printf(" %d>> ", buffer_send[i]);
    MPI_Send(
        buffer_send,               
        send_counter+1,               
        MPI_UINT32_T,                    
        mpi_rank ^ (1 << j),    
        0,                          
        MPI_COMM_WORLD              
    );
    uint32_t *temp_array = (uint32_t *) malloc( N * sizeof(uint32_t));
    //memcpy(temp_array, array, array_size * sizeof(uint32_t));
    for( i=0; i<N ; i++ )
    {
	    temp_array[i] = my_array[i];
    }   

    int k = 1 ;
    int m = N-1 ;
    int buffer_size = (int)buffer_recieve[0];
		
    //for(i=0;i<=buffer_size;i++)
//	printf(" %d> ", buffer_recieve[i]);
    for (i = N-1 ; i >= 0 ; i--)
    {
	//printf("Buffer receive ele in CH: %d\n", buffer_recieve[i]);
        //if (buffer_recieve[i] > array[0]) {
          //  array[0] = buffer_recieve[i];
        //} else {
          //  break;      
        //}
      //  printf("buffer_rec[k] is %d, temp_array[m] is %d\n",buffer_recieve[k], temp_array[m]);
//	printf("M is %d k is %d i is %d\n",m, k, i);
        if( temp_array[m] >= buffer_recieve[k] )
        {
            my_array[i]=temp_array[m];
            m--;
        }
        else if( k <= buffer_size )
        {
            my_array[i] = buffer_recieve[k] ;
            k++;
        }
    }

    qsort( my_array, N, sizeof(uint32_t), ComparisonFunc );
	
    free( buffer_send );
    free( buffer_recieve );
    return;
}
