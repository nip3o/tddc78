#!/bin/bash
#SBATCH -N 4
#SBATCH -t 00:01:00
mpprun ./thresc images/im2.ppm images/thresh.ppm
