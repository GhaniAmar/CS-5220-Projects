#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export OMP_NUM_THREADS=4
./sph.x -s 3e-2
exit 0