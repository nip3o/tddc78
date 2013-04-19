#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "ppmio.h"
#include "blurfilter.h"
#include "gaussw.h"
#include "mpi.h"

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

    buffsize = xsize * ysize / ntasks + 1;

    pixel recvbuff[MAX_PIXELS];

    MPI_Scatter(src, buffsize, pixel_mpi, 
                recvbuff, buffsize, pixel_mpi,
                ROOT, MPI_COMM_WORLD);

    clock_gettime(CLOCK_REALTIME, &stime);

    startY = (taskid) * (ysize / ntasks);
    endY = (taskid + 1) * (ysize / ntasks);

    printf("x = %i, ysize: %i, buffsize: %i, startY: %i, endY: %i, from: %i, radius: %i\n", xsize, ysize, buffsize, startY, endY, taskid, radius);

    blurfilter(xsize, 0, (ysize / ntasks)+1, recvbuff, radius, w, taskid);

    clock_gettime(CLOCK_REALTIME, &etime);

    MPI_Gather(recvbuff, buffsize, pixel_mpi,
               recvbuff, buffsize, pixel_mpi,
               ROOT, MPI_COMM_WORLD);

    printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
	   1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;

    /* write result */
    if (taskid == ROOT) {
        printf("Writing output file\n");
    
        if(write_ppm (argv[3], xsize, ysize, (char *)recvbuff) != 0)
          exit(1);
    }

    MPI_Finalize();


    return(0);
}
