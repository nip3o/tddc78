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

    int buffsize;
    
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

        buffsize = xsize * ysize / ntasks + 1;
    }
    
    MPI_Bcast(&buffsize, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    pixel recvbuff[MAX_PIXELS];

    MPI_Scatter(src, buffsize, pixel_mpi, 
                recvbuff, buffsize, pixel_mpi,
                ROOT, MPI_COMM_WORLD);

    printf("Has read the image, calling filter\n");

    clock_gettime(CLOCK_REALTIME, &stime);

    int i;
    unsigned int threshold_level;
    if (taskid == ROOT) {

      for(i = 0, threshold_level = 0; i < buffsize; i++) {
        threshold_level += (uint)src[i].r + (uint)src[i].g + (uint)src[i].b;
      }

      threshold_level /= (xsize * ysize);
    }

    MPI_Bcast(&threshold_level, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    thresfilter(buffsize, recvbuff, threshold_level);

    clock_gettime(CLOCK_REALTIME, &etime);

    printf("buffsize: %i, xsize: %i, ysize: %i, ntasks: %i\n", buffsize, xsize, ysize, ntasks);

    MPI_Gather(recvbuff, buffsize, pixel_mpi,
               recvbuff, buffsize, pixel_mpi,
               ROOT, MPI_COMM_WORLD);

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
