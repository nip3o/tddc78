#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "ppmio.h"
#include "thresfilter.h"

#define NUM_THREADS 4


struct thread_data {
   int thread_id;
   int start;
   int end;
   int threshold_level;
   pixel *image;
};

void *haaard_workwork(void *targs)
{
    struct thread_data *data = (struct thread_data *) targs;

    printf("Starting filtering on %d with start %i, end %i, threshold_level: %i...\n", data->thread_id, data->start, data->end, data->threshold_level);

    thresfilter(data->start, data->end, data->image, data->threshold_level);

    pthread_exit(NULL);
}

int main (int argc, char ** argv) {
    pthread_t threads[NUM_THREADS];

    int xsize, ysize, colmax;
    struct timespec stime, etime;
    unsigned int threshold_level;

    pixel* src = (pixel*)malloc(MAX_PIXELS * sizeof(pixel));
    if(!src) {
        printf("Could not mallocate memory\n");
        exit(1);
    }

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
    // 5 * P FLOPS
    for(i = 0, threshold_level = 0; i < (xsize * ysize); i++) {
        threshold_level += (int)src[i].r + (int)src[i].g + (int)src[i].b;
    }
    threshold_level /= (xsize * ysize);

    struct thread_data data[NUM_THREADS];
    printf("Setting up thread data...\n");

    int t, rc;
    for(t = 0; t < NUM_THREADS; t++) {
        data[t].thread_id = t;

        data[t].start = t * ceil(ysize / NUM_THREADS) * xsize;
        if (t + 1 == NUM_THREADS)
            data[t].end = ysize * xsize;
        else
            data[t].end = (t + 1) * ceil(ysize / NUM_THREADS) * xsize;

        data[t].threshold_level = threshold_level;
        data[t].image = src;

        printf("Creating thread...%i\n", t);

       // clock_gettime(CLOCK_REALTIME, &stime);

        rc = pthread_create(&threads[t], NULL, haaard_workwork, (void *) &data[t]);

        if (rc){
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }

    for(t = 0; t < NUM_THREADS; t++) {
        pthread_join(threads[t], NULL);
    }

   // clock_gettime(CLOCK_REALTIME, &etime);
    printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
	       1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;

    /* write result */
    printf("Writing output file\n");

    if(write_ppm (argv[2], xsize, ysize, (char *)src) != 0)
      exit(1);

    free(src);
    free(out);

    return(0);
}
