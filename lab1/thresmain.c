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
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    pixel src[MAX_PIXELS];
    struct timespec stime, etime;
    struct timespec tstime, tetime;

    // Create a custom MPI datatype for pixel
    pixel item;
    MPI_Datatype pixel_mpi;
    MPI_Datatype type[3] = { MPI_UNSIGNED_CHAR, MPI_UNSIGNED_CHAR, MPI_UNSIGNED_CHAR };
    int blocklen[] = { 1, 1, 1 };
    MPI_Aint start, disp[3];

    MPI_Address( &item, &start );
    MPI_Address( &item.r, &disp[0] );
    MPI_Address( &item.g, &disp[1] );
    MPI_Address( &item.b, &disp[2] );

    disp[0] -= start;
    disp[1] -= start;
    disp[2] -= start;

    MPI_Type_struct(3, blocklen, disp, type, &pixel_mpi);
    MPI_Type_commit(&pixel_mpi);

    unsigned int buffsize, threshold_level;

    /* Take care of the arguments */
    if (taskid == ROOT) {
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

        /* Calculate chunk size for nodes. Every node gets a little bit "too much"
           data if the image size is not dividable by the number of processes */
        buffsize = ceil((float)ysize / (float)ntasks) * xsize;

        clock_gettime(CLOCK_REALTIME, &stime);
        // Calculate medium value of all pixels in the image
        int i;
        for(i = 0, threshold_level = 0; i < buffsize; i++) {
            threshold_level += (uint)src[i].r + (uint)src[i].g + (uint)src[i].b;
        }
        threshold_level /= (xsize * ysize);

        clock_gettime(CLOCK_REALTIME, &etime);
        printf("Threshold level calculation took: %g secs\n", taskid, (etime.tv_sec  - stime.tv_sec) +
                1e-9*(etime.tv_nsec  - stime.tv_nsec));

        clock_gettime(CLOCK_REALTIME, &tstime);
        clock_gettime(CLOCK_REALTIME, &stime);
    }

    // Broadcast chunk size to alla nodes
    MPI_Bcast(&buffsize, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    pixel recvbuff[MAX_PIXELS];

    // Send the image in chunks to all nodes
    MPI_Scatter(src, buffsize, pixel_mpi,
                recvbuff, buffsize, pixel_mpi,
                ROOT, MPI_COMM_WORLD);

    // Broadcast threshold level
    MPI_Bcast(&threshold_level, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    if (taskid == ROOT) {
        clock_gettime(CLOCK_REALTIME, &etime);
        printf("Data shuffeling took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
           1e-9*(etime.tv_nsec  - stime.tv_nsec));
    }

    clock_gettime(CLOCK_REALTIME, &stime);

    // Run the filter on the recieved chunk
    thresfilter(buffsize, recvbuff, threshold_level);

    clock_gettime(CLOCK_REALTIME, &etime);
    printf("Filtering at %i took: %g secs\n", taskid, (etime.tv_sec  - stime.tv_sec) +
        1e-9*(etime.tv_nsec  - stime.tv_nsec));

    if (taskid == ROOT) {
        clock_gettime(CLOCK_REALTIME, &stime);
    }

    // Collect the data from the nodes
    MPI_Gather(recvbuff, buffsize, pixel_mpi,
               recvbuff, buffsize, pixel_mpi,
               ROOT, MPI_COMM_WORLD);

    if (taskid == ROOT) {
        clock_gettime(CLOCK_REALTIME, &etime);
        printf("Gather took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
           1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;

        clock_gettime(CLOCK_REALTIME, &tetime);
        printf("Everything took: %g secs\n", (tetime.tv_sec  - tstime.tv_sec) +
           1e-9*(tetime.tv_nsec  - tstime.tv_nsec));

        /* write result */
        printf("Writing output file\n");

        if(write_ppm (argv[2], xsize, ysize, (char *)recvbuff) != 0)
          exit(1);
    }


    MPI_Finalize();

    return(0);
}
