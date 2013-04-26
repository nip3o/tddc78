#include "thresfilter.h"

void thresfilter(const int nump, pixel* src, const unsigned int threshold_level){
#define uint unsigned int

  uint i, psum;

  for(i = 0; i < nump; i++) {
    psum = (uint)src[i].r + (uint)src[i].g + (uint)src[i].b;

    if(threshold_level > psum) {
      src[i].r = src[i].g = src[i].b = 0;
    }
    else {
      src[i].r = src[i].g = src[i].b = 255;
    }
  }
}
