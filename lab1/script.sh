#!/bin/bash
#SBATCH -N 2
#SBATCH -t 00:00:10
mpprun ./thresc images/im4.ppm images/thresh.ppm
