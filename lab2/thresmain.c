#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "ppmio.h"
#include "thresfilter.h"

#define NUM_THREADS 2


struct thread_data {
   int thread_id;
   int start;
   int end;
   int threshold_level;
   pixel *chunk;
};

void *haaard_workwork(void *targs)
{
    printf("Haaard working\n");
    struct thread_data *data = (struct thread_data *) targs;

 //   pixel *src = malloc(sizeof(pixel) * MAX_PIXELS);
 //   memcpy( &data->chunk, &src, sizeof(pixel) * MAX_PIXELS);


    printf("Starting filtering on %ld with start %i, end %i, threshold_level: %i...\n", data->thread_id, data->start, data->end, data->threshold_level);

    thresfilter(data->start, data->end, data->chunk, data->threshold_level);

    char s[25];
    sprintf(s, "images/%d.ppm", data->thread_id);
    if(write_ppm (s, 676, 763, (char *)data->chunk) != 0)
      exit(1);

    pthread_exit(NULL);
}

int main (int argc, char ** argv) {
    pthread_t threads[NUM_THREADS];

    int xsize, ysize, colmax;
    pixel src[MAX_PIXELS];
    struct timespec stime, etime;
    unsigned int threshold_level;

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

    int i;
    for(i = 0, threshold_level = 0; i < (xsize * ysize); i++) {
        threshold_level += (int)src[i].r + (int)src[i].g + (int)src[i].b;
    }
    threshold_level /= (xsize * ysize);

//    clock_gettime(CLOCK_REALTIME, &stime);

    struct thread_data data[NUM_THREADS];

    printf("Setting up thread data...\n");

//    int chunksize = ceil(MAX_PIXELS / NUM_THREADS);

    int t, rc;
    for(t = 0; t < NUM_THREADS; t++) {
        data[t].thread_id = t;

        data[t].start = t * ceil(ysize / NUM_THREADS) * xsize;
        data[t].end = (t + 1) * ceil(ysize / NUM_THREADS) * xsize;

        data[t].threshold_level = threshold_level;
        data[t].chunk = src;

        printf("Creating thread...%i\n", t);

        rc = pthread_create(&threads[t], NULL, haaard_workwork, (void *) &data[t]);

        if (rc){
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }

//    clock_gettime(CLOCK_REALTIME, &etime);

//    printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
//       1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;

    /* write result */
    printf("Writing output file\n");

    if(write_ppm (argv[2], xsize, ysize, (char *)src) != 0)
      exit(1);


    return(0);
}
