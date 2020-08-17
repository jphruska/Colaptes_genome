#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N colaptes_repeatmasker_hic_2
#$ -q omni
#$ -pe sm 24
#$ -P quanah


/home/jmanthey/RepeatMasker/RepeatMasker -pa 24 -s -lib /home/jmanthey/RepeatMasker/Libraries/custom_library_certhia_colaptes.fa c05_colaptes_arcs3.FINAL.reordered_2.fasta -a -dir ./repeatmasker_2 
