#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N colaptes_repeatmodeler_build_database
#$ -q omni
#$ -pe sm 8
#$ -P quanah

/home/jmanthey/RepeatModeler-open-1.0.11/BuildDatabase -engine ncbi -name colaptes_hic c05_colaptes_arcs3.FINAL.reordered_2.fasta  
