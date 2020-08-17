#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N c_codeml
#$ -q omni
#$ -pe sm 1
#$ -P quanah
#$ -l h_rt=12:00:00
#$ -l h_vmem=12G
#$ -t 1-9813

# create a new directory for this alignment and move to it
mkdir ${SGE_TASK_ID}_work
cd ${SGE_TASK_ID}_work/

# copy a control file to this directory
cp ../../codeml.ctl .

# copy the alignment to this directory
cp ../${SGE_TASK_ID}_aligned_trimmed.fasta .

# change the control file input and output names 
sed -i "s/blank.fasta/${SGE_TASK_ID}_aligned_trimmed.fasta/g" codeml.ctl
sed -i "s/blank_output.txt/${SGE_TASK_ID}_output.txt/g" codeml.ctl

# run 1st iteration of codeml
codeml codeml.ctl


