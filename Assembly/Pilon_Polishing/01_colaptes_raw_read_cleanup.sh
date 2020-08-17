#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N colaptes_raw_read_cleanup
#$ -q omni
#$ -pe sm 2
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=24G

module load java 

/lustre/work/johruska/bbmap/bbduk.sh in1=/lustre/scratch/johruska/colaptes/chromium_10x/raw/5469-JDM-0003_S1_L001_R1_001.fastq.gz \
in2=/lustre/scratch/johruska/colaptes/chromium_10x/raw/5469-JDM-0003_S1_L001_R2_001.fastq.gz out1=/lustre/scratch/johruska/colaptes/chromium_10x/cleaned/5469-JDM-0003_S1_L001_R1_001.fastq.gz \
out2=/lustre/scratch/johruska/colaptes/chromium_10x/cleaned/5469-JDM-0003_S1_L001_R2_002.fastq.gz minlen=50 ftl=10 qtrim=rl trimq=10 ktrim=r k=25 mink=7 ref=/lustre/work/johruska/bbmap/resources/adapters.fa \ 
hdist=1 tbo tpe
