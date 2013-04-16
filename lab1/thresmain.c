#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <time.h>
#include "ppmio.h"
#include "thresfilter.h"
#include "mpi.h"

#define ROOT 0

int main (int argc, char ** argv) {
    int xsize, ysize, colmax, taskid, ntasks;
    pixel src[MAX_PIXELS];
    struct timespec stime, etime;

    MPI_Datatype pixel_mpi;
    MPI_Datatype type[3] = { MPI_UNSIGNED_CHAR, MPI_UNSIGNED_CHAR, MPI_UNSIGNED_CHAR };
    int blocklen[3] = { 1, 1, 1 };
    MPI_Aint start, disp[3];
 
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    MPI_Address( &src[0], &start );
    MPI_Address( &src[0].r, &disp[0] );
    MPI_Address( &src[0].g, &disp[1] );
    MPI_Address( &src[0].b, &disp[2] );
 
    disp[0] -= start;
    disp[1] -= start;
    disp[2] -= start;

    MPI_Type_create_struct(3, blocklen, disp, type, &pixel_mpi);
    MPI_Type_commit(&pixel_mpi);

    int buffsize = 10;
    
    /* Take care of the arguments */
    if (taskid == ROOT) {
        printf("I am root =)\n");
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

        buffsize = xsize * ysize / ntasks + 1;
    }
    printf("buffsize: %i, xsize: %i, ysize: %i, ntasks: %i\n", buffsize, xsize, ysize, ntasks);

    pixel recvbuff[MAX_PIXELS];


    MPI_Scatter(src, buffsize, pixel_mpi, 
                recvbuff, buffsize, pixel_mpi,
                0, MPI_COMM_WORLD);

//MPI_Scatter(sendbuff[0],buffsize,MPI_DOUBLE,
//                   recvbuff,buffsize,MPI_DOUBLE,
//                   0,MPI_COMM_WORLD);

//MPI_Gather(&buffsum,1,MPI_DOUBLE,
//                   buffsums,1, MPI_DOUBLE,
//                   0,MPI_COMM_WORLD);



    printf("Has read the image, calling filter\n");

    clock_gettime(CLOCK_REALTIME, &stime);

    thresfilter(buffsize, src);

    clock_gettime(CLOCK_REALTIME, &etime);

    MPI_Gather(src, buffsize, pixel_mpi,
               recvbuff, buffsize, pixel_mpi,
               0, MPI_COMM_WORLD);

    if (taskid == ROOT) {
        printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
    	   1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;

        /* write result */
        printf("Writing output file\n");
        
        if(write_ppm (argv[2], xsize, ysize, (char *)recvbuff) != 0)
          exit(1);
    } 


    MPI_Finalize();

    return(0);
}
