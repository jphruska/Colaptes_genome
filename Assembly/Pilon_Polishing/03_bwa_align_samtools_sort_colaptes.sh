#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N colaptes_bwa_align_samtools_sort
#$ -q omni
#$ -pe sm 8
#$ -P xlquanah 
#$ -l h_rt=120:00:00
#$ -l h_vmem=24G

module load intel samtools bwa 

bwa mem -t 8 /lustre/scratch/johruska/colaptes/canu_assembly/colaptes.contigs.fasta /lustre/scratch/johruska/colaptes/chromium_10x/cleaned/5469-JDM-0003_S1_L001_R1_001.fastq /lustre/scratch/johruska/colaptes/chromium_10x/cleaned/5469-JDM-0003_S1_L001_R2_002.fastq | samtools sort > aln_illumina_pacbio_colaptes.bam

