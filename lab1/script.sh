#!/bin/bash
#SBATCH -N 1
#SBATCH -t 00:00:30
mpprun ./thresc images/im4.ppm images/thresh.ppm
