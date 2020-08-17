#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N colaptes_hi_c_juicer
#$ -q omni
#$ -pe sm 8
#$ -P quanah

/lustre/work/johruska/juicer/scripts/juicer.sh -z /lustre/work/johruska/juicer/references/c05_colaptes_arcs3.fasta \
-p /lustre/work/johruska/juicer/restriction_sites/c05_colaptes_arcs3.chrom.sizes \
-y /lustre/work/johruska/juicer/restriction_sites/c05_colaptes_arcs3_Sau3AI.txt \
-d /lustre/scratch/johruska/colaptes/colaptes_reference_genome/Colaptes_HI-C_2 \
-q omni -l omni -s Sau3AI -a 'Colaptes genome assembly Hi-C' \
-D /lustre/work/johruska/juicer -Q 48:00:00 -L 48:00:00 -t 8
