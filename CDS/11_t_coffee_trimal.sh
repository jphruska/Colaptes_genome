#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N c_t-coffe_align
#$ -q omni
#$ -pe sm 1
#$ -P quanah
#$ -l h_rt=12:00:00
#$ -l h_vmem=12G
#$ -t 1-9813

# translate the sequences in the fasta file
t_coffee -other_pg seq_reformat -in ${SGE_TASK_ID}.fasta -action +translate -output fasta_seq > ${SGE_TASK_ID}_protein.fasta

# align the sequences
t_coffee ${SGE_TASK_ID}_protein.fasta -mode mcoffee

# back translate to nucleotides
t_coffee -other_pg seq_reformat -in ${SGE_TASK_ID}.fasta -in2 ${SGE_TASK_ID}_protein.aln -action +thread_dna_on_prot_aln \
-output fasta_aln > ${SGE_TASK_ID}_aligned.fasta

# remove gaps from alignment
trimal -in ${SGE_TASK_ID}_aligned.fasta -out ${SGE_TASK_ID}_aligned_trimmed.fasta -nogaps -fasta

