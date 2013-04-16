#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "ppmio.h"
#include "thresfilter.h"
#include "mpi.h"

#define ROOT 0

int main (int argc, char ** argv) {
    int rank, numprocs;

    //-- Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int xsize, ysize, colmax;
    struct timespec stime, etime;

    pixel pix;
    MPI_Datatype pix_mpi;
    // MPI type to commit
    int block_lengths [] = { 1, 1, 1 };
    // Lengths of type elements
    MPI_Datatype block_types [] = { MPI_UNSIGNED_CHAR, MPI_UNSIGNED_CHAR, MPI_UNSIGNED_CHAR };

    //Set types
    MPI_Aint start, displ[3];

    MPI_Address( &pix, &start );
    MPI_Address( &pix.r, &displ[0] );
    MPI_Address( &pix.g, &displ[1] );
    MPI_Address( &pix.b, &displ[2] );

    displ[0] -= start;
    displ[1] -= start;
    displ[2] -= start;

    MPI_Type_struct( 3, block_lengths, displ, block_types, &pix_mpi );
    MPI_Type_commit( &pix_mpi );

    pixel src[MAX_PIXELS];

    /* Take care of the arguments */

    if (argc != 3) {
        fprintf(stderr, "Usage: %s infile outfile\n", argv[0]);
        exit(1);
    }

    int num[numprocs];
    int distances[numprocs];

    if (rank == ROOT) {
        /* read file */
        if(read_ppm (argv[1], &xsize, &ysize, &colmax, (char *) src) != 0)
            exit(1);

        if (colmax > 255) {
            fprintf(stderr, "Too large maximum color-component value\n");
            exit(1);
        }

        printf("File initialized\n");

        int i;
        for (i = 0; i < numprocs; i++) {
            num[i] = 1000;
            distances[i] = i * numprocs;
        }

        printf("Generated problem sizes\n");
    }

    pixel *recvbuffer = (pixel*) malloc(1000 * sizeof(pixel));

    if(recvbuffer == 0) {
       printf("Allocation failed");
    } else {
       printf("Mallocated %i bytes\n", (int) (1000 * sizeof(pixel)));
    }

    printf("Scattering");

    MPI_Scatter( src, 1000, pix_mpi, rbuf, 1000, pix_mpi, ROOT, MPI_COMM_WORLD);
    //MPI_Scatterv(src, num, distances, pix_mpi, recvbuffer, numprocs,
    //             pix_mpi, ROOT, MPI_COMM_WORLD);


    printf("Has read the image, calling filter\n");

    clock_gettime(CLOCK_REALTIME, &stime);

    thresfilter(100, 10, recvbuffer);

    clock_gettime(CLOCK_REALTIME, &etime);


    printf("Gathering");


    pixel* rbuf = (int *)malloc(numprocs * 1000 * sizeof(pix_mpi));
    MPI_Gather( src, 1000, pix_mpi, rbuf, 1000, pix_mpi, ROOT, MPI_COMM_WORLD);

    //MPI_Gatherv(recvbuffer, 1000, pix_mpi, src, num, distances,
    //            pix_mpi, ROOT, MPI_COMM_WORLD);

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
