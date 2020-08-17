#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N colaptes_show_coords_gallus_reordered
#$ -q omni
#$ -pe sm 12
#$ -P quanah


/home/jmanthey/mummer/bin/show-coords -l -q -T colaptes_reordered_filtered.delta > colaptes_reordered_filtered.txt


