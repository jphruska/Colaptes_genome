#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N colaptes_repeatmodeler 
#$ -q omni
#$ -pe sm 8
#$ -P quanah

/home/jmanthey/RepeatModeler-open-1.0.11/RepeatModeler -database colaptes -engine ncbi -pa 11 

