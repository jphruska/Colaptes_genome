#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N colaptes_pilon_1
#$ -q Yoda 
#$ -P communitycluster  
#$ -pe sm 20

java -jar -Xmx500G /lustre/work/johruska/pilon-1.22.jar --genome /lustre/scratch/johruska/colaptes/canu_assembly/colaptes.contigs.fasta --frags /lustre/scratch/johruska/colaptes/colaptes_scripts/aln_illumina_pacbio_colaptes.bam --fix all --output colaptes_pilon_1 --outdir ../colaptes_pilon_1 --changes --tracks --threads 12  
