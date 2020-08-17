#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N colaptes_nucmer_gallus
#$ -q omni
#$ -pe sm 12
#$ -P quanah

/home/jmanthey/mummer/bin/nucmer -p colaptes_reordered -t 12 GCF_000002315.6_GRCg6a_genomic.fna c05_colaptes_arcs3.FINAL.reordered_2.fasta   
