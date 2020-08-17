#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N colaptes_genome_index
#$ -q omni
#$ -pe sm 2
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=24G

module load intel bwa samtools java 

cd /lustre/scratch/johruska/colaptes/canu_assembly/
bwa index colaptes.contigs.fasta
samtools faidx colaptes.contigs.fasta
/lustre/work/johruska/gatk-4.0.8.1/gatk --java-options "-Xmx10g" CreateSequenceDictionary -R colaptes.contigs.fasta
