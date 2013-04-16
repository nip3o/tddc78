#!/bin/sh
rm *.out
rm images/thresh.ppm
icc -Nmpi -o thresc thresmain.c ppmio.c thresfilter.c
if [ $? = 0 ] ; then
   sbatch script.sh
   squeue -u $USER
fi

