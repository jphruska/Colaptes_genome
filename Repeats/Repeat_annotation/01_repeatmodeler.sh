#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N colaptes_repeatmodeler_hic 
#$ -q omni
#$ -pe sm 8
#$ -P quanah

/home/jmanthey/RepeatModeler-open-1.0.11/RepeatModeler -database colaptes_hic -engine ncbi -pa 11 

