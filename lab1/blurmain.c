#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "ppmio.h"
#include "blurfilter.h"
#include "gaussw.h"
#include "mpi.h"

#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )

#define ROOT 0
#define MAX_RAD 1000

int main (int argc, char ** argv) {
   int taskid, ntasks;

    int xsize, ysize, colmax;
    pixel src[MAX_PIXELS];
    double w[MAX_RAD];

    struct timespec stime, etime;
    struct timespec tstime, tetime;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

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

    int buffsize, radius, startY, endY;

    /* Take care of the arguments */
    if (argc != 4) {
        fprintf(stderr, "Usage: %s radius infile outfile\n", argv[0]);
        exit(1);
    }
    radius = atoi(argv[1]);
    if((radius > MAX_RAD) || (radius < 1)) {
        fprintf(stderr, "Radius (%d) must be greater than zero and less then %d\n", radius, MAX_RAD);
        exit(1);
    }

    if (taskid == ROOT) {
        /* read file */
        if(read_ppm (argv[2], &xsize, &ysize, &colmax, (char *) src) != 0)
            exit(1);

        if (colmax > 255) {
            fprintf(stderr, "Too large maximum color-component value\n");
            exit(1);
        }

        /* filter */
        printf("Has read the image, generating coefficients\n");
        get_gauss_weights(radius, w);
    }

    // Broadcast the gaussian weight vector
    MPI_Bcast(w, MAX_RAD, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    // Broadcast image dimensions
    MPI_Bcast(&xsize, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&ysize, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    // Calculate chunk size
    buffsize = ceil((float)ysize / (float)ntasks) * xsize;
    pixel recvbuff[MAX_PIXELS];

    int sendcnts[ntasks], displs[ntasks], result_write_starts[ntasks], recievecounts[ntasks];
    int i;
    // Generate sendcount and displacement vectors for Scatterv
    for (i = 0; i < ntasks; i++) {
        // Send enought neighbors to make it possible to also calculate
        // blur in the edges of the chunk
        sendcnts[i] = buffsize + 2 * radius * xsize;
        displs[i] = max(0, i * buffsize);
    }

    clock_gettime(CLOCK_REALTIME, &tstime);

    // Send the image in chunks to all nodes
    MPI_Scatterv(src, sendcnts, displs,
                 pixel_mpi, recvbuff, buffsize + 2 * radius * xsize,
                 pixel_mpi, ROOT, MPI_COMM_WORLD);

    clock_gettime(CLOCK_REALTIME, &stime);

    // Run the filter on the recieved chunk
    blurfilter(xsize, (ysize / ntasks) + 2 * radius, recvbuff, radius, w, taskid);

    clock_gettime(CLOCK_REALTIME, &etime);
    printf("Filtering at %i took: %g secs\n", taskid, (etime.tv_sec  - stime.tv_sec) +
        1e-9*(etime.tv_nsec  - stime.tv_nsec));

    // Generate sendcount and displacement vectors for Scatterv
    for (i = 0; i < ntasks; i++) {
        result_write_starts[i] = i * buffsize + xsize * radius;
        // Only send as much of the chunk that is really useful data
        recievecounts[i] = buffsize;
    }

    // Start writing from the beginning of the buffer if root
    result_write_starts[0] = 0;

    // Since the root node has no overlap in the beginning, we need to
    // send a little bit more from that node than from the rest.
    recievecounts[0] = buffsize + xsize * radius;

    pixel* result_read_start;
    if(taskid==ROOT) {
        // Root-node has no duplicated data in the beginning
        result_read_start = recvbuff;
    } else {
        // Jump over the duplicated data in the beginning of each chunk
        result_read_start = recvbuff + xsize * radius;
    }

    MPI_Gatherv(result_read_start, recievecounts[taskid], pixel_mpi,
                src, recievecounts, result_write_starts,
                pixel_mpi, ROOT, MPI_COMM_WORLD);

    clock_gettime(CLOCK_REALTIME, &tetime);

    MPI_Finalize();


    /* write result */
    if (taskid == ROOT) {
        printf("Everything took: %g secs\n", (tetime.tv_sec  - tstime.tv_sec) +
           1e-9*(tetime.tv_nsec  - tstime.tv_nsec));


        printf("Writing output file\n");

        if(write_ppm (argv[3], xsize, ysize, (char *)src) != 0)
          exit(1);
    }

    return(0);
}
