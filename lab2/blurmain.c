#include <pthread.h>
#include <semaphore.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <errno.h>

#include "ppmio.h"
#include "blurfilter.h"
#include "gaussw.h"


#define NUM_THREADS 2


#define ROOT 0
#define MAX_RAD 1000


struct thread_data {
    int thread_id;

    int xsize;
    int startY;
    int endY;
    int radius;

    sem_t* semaphore;
    pthread_mutex_t* zone_lock;

    double *w;
    pixel *image;
};

void *unsynced_blur(void *targs)
{
    struct thread_data *data = (struct thread_data *) targs;

    printf("Starting filtering on %d with startY %i, endY %i\n", data->thread_id, data->startY, data->endY);
    blurfilter(data->xsize, data->startY, data->endY, data->image, data->radius,
               data->w, data->thread_id,
               data->semaphore, data->zone_lock);

    pthread_exit(NULL);
}


int main (int argc, char ** argv) {
    int radius;
    int xsize, ysize, colmax;
    struct timespec stime, etime;

    double w[MAX_RAD];

    pixel* src = (pixel*)malloc(MAX_PIXELS * sizeof(pixel));
    if(!src) {
        printf("Could not mallocate memory\n");
        exit(1);
    }

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

    /* read file */
    if(read_ppm (argv[2], &xsize, &ysize, &colmax, (char *) src) != 0)
        exit(1);

    if (colmax > 255) {
        fprintf(stderr, "Too large maximum color-component value\n");
        exit(1);
    }

    printf("Has read the image, generating coefficients\n");


    /* filter */
    get_gauss_weights(radius, w);

    printf("Calling filter\n");

//    clock_gettime(CLOCK_REALTIME, &stime);

    struct thread_data data[NUM_THREADS];
    printf("Setting up thread data...\n");

    pthread_t threads[NUM_THREADS];
    int t, rc;

    pthread_mutex_t zone_lock;
    pthread_mutex_init( &zone_lock, NULL);

    sem_t semaphore;
    sem_init(&semaphore, 0, 0);

    for(t = 0; t < NUM_THREADS; t++) {
        data[t].thread_id = t;

        data[t].xsize = xsize;

        if(t == 0) {
            data[t].startY = 0;
        } else {
            data[t].startY = t * ceil(ysize / NUM_THREADS);
        }
        if(t == NUM_THREADS - 1) {
            data[t].endY = ysize;
        } else {
            data[t].endY = (t + 1) * ceil(ysize / NUM_THREADS);
        }

        data[t].radius = radius;
        data[t].w = w;
        data[t].image = src;

        data[t].zone_lock = &zone_lock;
        data[t].semaphore = &semaphore;

        printf("Creating thread...%i\n", t);
//        clock_gettime(CLOCK_REALTIME, &stime);

        pthread_create(&threads[t], NULL, unsynced_blur, (void *) &data[t]);
    }


    for(t = 0; t < NUM_THREADS; t++) {
        pthread_join(threads[t], NULL);
    }

//    blurfilter(xsize, 226, 286, src, radius, w);


//    clock_gettime(CLOCK_REALTIME, &etime);

//    printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
//       1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;

    /* write result */
    printf("Writing output file\n");

    if(write_ppm (argv[3], xsize, ysize, (char *)src) != 0)
      exit(1);


    return(0);
}
