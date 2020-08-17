#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N colaptes_filter_gallus_reordered
#$ -q omni
#$ -pe sm 12
#$ -P quanah


/home/jmanthey/mummer/bin/delta-filter -m -i 70 colaptes_reordered.delta > colaptes_reordered_filtered.delta
