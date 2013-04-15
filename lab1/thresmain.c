#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "ppmio.h"
#include "thresfilter.h"
#include "mpi.h"

int main (int argc, char ** argv) {
    int xsize, ysize, colmax;
    pixel src[MAX_PIXELS];
    struct timespec stime, etime;
    
    
    // Composite type
    buf_t item;
    // Element of the type
    MPI_Datatype buf_t_mpi;
    // MPI type to commit
    int
    block_lengths [] = { 1, 10 };
    // Lengths of type elements
    MPI_Datatype block_types [] = { MPI_INT, MPI_DOUBLE }
    ;
    //Set types
    MPI_Aint start, displ[2];
    MPI_Address( &item, &start );
    MPI_Address( &item.id, &displ[0] );
    MPI_Address( &item.data[0], &displ[1] );
    displ[0] -= start;
    // Displacement relative to address of start
    displ[1] -= start;
    // Displacement relative to address of start
    MPI_Type_struct( 2, block_lengths, displ, block_typ
                    es, &buf_t_mpi );
    MPI_Type_commit( &buf_t_mpi );
    

    /* Take care of the arguments */

    if (argc != 3) {
        fprintf(stderr, "Usage: %s infile outfile\n", argv[0]);
        exit(1);
    }

    /* read file */
    if(read_ppm (argv[1], &xsize, &ysize, &colmax, (char *) src) != 0)
        exit(1);

    if (colmax > 255) {
        fprintf(stderr, "Too large maximum color-component value\n");
        exit(1);
    }
    
    int rank, numprocs,left,right;
    const int Message_size = 128;
    double Message[Message_size];
    double sum;
    
    //-- Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    

    if (rank == 0) {
        MPI_Scatterv(void *sendbuf, int *sendcnts, int *displs,
                     MPI_Datatype sendtype, void *recvbuf, int recvcnt,
                     MPI_Datatype recvtype,
                     int root, MPI_Comm comm)
    }
    
    
    printf("Has read the image, calling filter\n");

    clock_gettime(CLOCK_REALTIME, &stime);

    thresfilter(xsize, ysize, src);

    clock_gettime(CLOCK_REALTIME, &etime);
    
    
    //-- Terminate MPI.
    MPI_Finalize();

    printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
	   1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;

    /* write result */
    printf("Writing output file\n");
    
    if(write_ppm (argv[2], xsize, ysize, (char *)src) != 0)
      exit(1);


    return(0);
}
