/*
  File: blurfilter.c

  Implementation of blurfilter function.

 */
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>
#include <assert.h>
#include <stdbool.h>

#include "blurfilter.h"
#include "ppmio.h"

#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )


pixel* pix(pixel* image, const int xx, const int yy, const int xsize)
{
  register int off = xsize*yy + xx;

#ifdef DBG
  if(off >= MAX_PIXELS) {
    fprintf(stderr, "\n Terribly wrong: %d %d %d\n",xx,yy,xsize);
  }
#endif
  return (image + off);
}

void blurfilter(const int xsize, const int startY, const int endY,
                pixel* src, pixel* out,
                const int radius, const double *w, const int thread_id, const int ysize) {
    int x,y,x2,y2, wi;
    double r,g,b,n, wc;
    bool unlocked = false;

    pixel* dst = (pixel*)malloc(MAX_PIXELS * sizeof(pixel));
    if(!src) {
        printf("Could not mallocate dst memory\n");
        exit(1);
    }


    // 2 * P * (22 * R + 25) FLOPS

    // ysize times
    for (y = max(0, startY - radius); y < min(ysize, endY + radius); y++) {
        // xsize times
        for (x = 0; x < xsize; x++) {
            // 9 FLOP
            r = w[0] * pix(src, x, y, xsize)->r;
            g = w[0] * pix(src, x, y, xsize)->g;
            b = w[0] * pix(src, x, y, xsize)->b;
            n = w[0];
            for ( wi = 1; wi <= radius; wi++) {
                // 22 FLOP * radius
                wc = w[wi];
                x2 = x - wi;
                if(x2 >= 0) {
                    r += wc * pix(src, x2, y, xsize)->r;
                    g += wc * pix(src, x2, y, xsize)->g;
                    b += wc * pix(src, x2, y, xsize)->b;
                    n += wc;
                }
                x2 = x + wi;
                if(x2 < xsize) {
                    r += wc * pix(src, x2, y, xsize)->r;
                    g += wc * pix(src, x2, y, xsize)->g;
                    b += wc * pix(src, x2, y, xsize)->b;
                    n += wc;
                }
            }
            // 16 FLOP
            pix(dst,x,y, xsize)->r = r/n;
            pix(dst,x,y, xsize)->g = g/n;
            pix(dst,x,y, xsize)->b = b/n;
        }
    }

    // ysize times
    for (y = startY; y < endY; y++) {
        // xsize times
        for (x = 0; x < xsize; x++) {
            r = w[0] * pix(dst, x, y, xsize)->r;
            g = w[0] * pix(dst, x, y, xsize)->g;
            b = w[0] * pix(dst, x, y, xsize)->b;
            n = w[0];
            for ( wi = 1; wi <= radius; wi++) {
                wc = w[wi];
                y2 = y - wi;
                if(y2 >= 0) {
                    r += wc * pix(dst, x, y2, xsize)->r;
                    g += wc * pix(dst, x, y2, xsize)->g;
                    b += wc * pix(dst, x, y2, xsize)->b;
                    n += wc;
                }
                y2 = y + wi;
                if(y2 <= endY + radius) {
                    r += wc * pix(dst, x, y2, xsize)->r;
                    g += wc * pix(dst, x, y2, xsize)->g;
                    b += wc * pix(dst, x, y2, xsize)->b;
                    n += wc;
                }
            }

            pix(out,x,y, xsize)->r = r/n;
            pix(out,x,y, xsize)->g = g/n;
            pix(out,x,y, xsize)->b = b/n;
        }
    }
}
