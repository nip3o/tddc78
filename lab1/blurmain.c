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


    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

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

    printf("Calling filter\n");

    MPI_Bcast(w, MAX_RAD, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&xsize, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&ysize, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    buffsize = xsize * ceil((float)ysize / (float)ntasks);
    pixel recvbuff[MAX_PIXELS];

    int sendcnts[ntasks], displs[ntasks], result_write_starts[ntasks], recievecounts[ntasks];
    int i;
    for (i = 0; i < ntasks; i++) {
        //printf("disp: %i, sendcnts: %i, buffsize: %i, taskid: %i\n", i*buffsize-radius*xsize, buffsize + radius*xsize, buffsize, taskid);
        sendcnts[i] = buffsize + 2 * radius * xsize;
        displs[i] = max(0, i * buffsize); // + 0 * radius * xsize);
    }

    MPI_Scatterv(src, sendcnts, displs,
                 pixel_mpi, recvbuff, buffsize + 2 * radius * xsize,
                 pixel_mpi, ROOT, MPI_COMM_WORLD);

    clock_gettime(CLOCK_REALTIME, &stime);

    blurfilter(xsize, (ysize / ntasks) + 2 * radius, recvbuff, radius, w, taskid);

    clock_gettime(CLOCK_REALTIME, &etime);

    for (i = 0; i < ntasks; i++) {
        recievecounts[i] = buffsize;
        result_write_starts[i] = i * buffsize + xsize * radius;
    }

    recievecounts[0] = buffsize + xsize * radius;
    recievecounts[ntasks-1] = buffsize + xsize * radius;

    result_write_starts[0] = 0;

    char fname[15];
    sprintf(fname, "%d.ppm", taskid);
    write_ppm (fname, xsize, ysize, (char *)recvbuff);

    pixel* result_read_start;
    if(taskid==ROOT) {
        result_read_start = recvbuff;
    } else {
        result_read_start = recvbuff + xsize * radius;
    }

    int sendcount = recievecounts[taskid];
    MPI_Gatherv(result_read_start, sendcount, pixel_mpi,
                src, recievecounts, result_write_starts,
                pixel_mpi, ROOT, MPI_COMM_WORLD);

    MPI_Finalize();

    printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
	   1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;

    /* write result */
    if (taskid == ROOT) {
        printf("Writing output file\n");

        if(write_ppm (argv[3], xsize, ysize, (char *)src) != 0)
          exit(1);
    }

    return(0);
}
