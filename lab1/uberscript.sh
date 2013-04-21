#!/bin/sh
rm *.out

if [ $1 = 0 ] ; then
  icc -Nmpi -o thresc thresmain.c ppmio.c thresfilter.c
  if [ $? = 0 ] ; then
    sbatch script.sh
    squeue -u $USER
  fi
elif [ $1 = 1 ] ; then
  icc -Nmpi -o thresc thresmain.c ppmio.c thresfilter.c
  if [ $? = 0 ] ; then
    mpirun -np 4 ./thresc images/im1.ppm images/thresh.ppm
  fi 
elif [ $1 = 2 ] ; then
  icc -Nmpi -o blurc blurmain.c ppmio.c blurfilter.c gaussw.c
  if [ $? = 0 ] ; then
    mpirun -np 4 ./blurc 15 images/im1.ppm images/blur.ppm
  fi
elif [ $1 = 3 ] ; then
  icc -Nmpi -o blurc blurmain.c ppmio.c blurfilter.c gaussw.c
  if [ $? = 0 ] ; then
    sbatch scriptBlur.sh
    squeue -u $USER
  fi
fi

