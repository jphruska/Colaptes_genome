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


# move to appropriate directory
cd ${SGE_TASK_ID}_work/

# copy certhia control file to this directory
cp ../../codeml_2.ctl .

# change the control file input and output names for colaptes
sed -i "s/blank.fasta/${SGE_TASK_ID}_aligned_trimmed.fasta/g" codeml_2.ctl
sed -i "s/blank_output.txt/${SGE_TASK_ID}_output_colaptes.txt/g" codeml_2.ctl
sed -i 's/unrooted_tree.tre/unrooted_tree_colaptes.tre/g' codeml_2.ctl

# run colaptes iteration of codeml
codeml codeml_2.ctl

# change the control file input and output names for apaloderma
sed -i "s/${SGE_TASK_ID}_output_colaptes.txt/${SGE_TASK_ID}_output_apaloderma.txt/g" codeml_2.ctl
sed -i 's/unrooted_tree_colaptes.tre/unrooted_tree_apaloderma.tre/g' codeml_2.ctl

# run apaloderma iteration of codeml
codeml codeml_2.ctl

# change the control file input and output names for buceros
sed -i "s/${SGE_TASK_ID}_output_apaloderma.txt/${SGE_TASK_ID}_output_buceros.txt/g" codeml_2.ctl
sed -i 's/unrooted_tree_apaloderma.tre/unrooted_tree_buceros.tre/g' codeml_2.ctl

# run buceros iteration of codeml
codeml codeml_2.ctl

# change the control file input and output names for merops
sed -i "s/${SGE_TASK_ID}_output_buceros.txt/${SGE_TASK_ID}_output_merops.txt/g" codeml_2.ctl
sed -i 's/unrooted_tree_buceros.tre/unrooted_tree_merops.tre/g' codeml_2.ctl

# run merops iteration of codeml
codeml codeml_2.ctl

# add all likelihoods and omegas to a single file
grep '^lnL' ${SGE_TASK_ID}_output.txt >> /lustre/scratch/johruska/colaptes/colaptes_reference_genome/colaptes_cds/${SGE_TASK_ID}_total_output.txt
grep '^omega' ${SGE_TASK_ID}_output.txt >> /lustre/scratch/johruska/colaptes/colaptes_reference_genome/colaptes_cds/${SGE_TASK_ID}_total_output.txt
grep '^lnL' ${SGE_TASK_ID}_output_colaptes.txt >> /lustre/scratch/johruska/colaptes/colaptes_reference_genome/colaptes_cds/${SGE_TASK_ID}_total_output.txt
grep '^w (dN/dS)' ${SGE_TASK_ID}_output_colaptes.txt >> /lustre/scratch/johruska/colaptes/colaptes_reference_genome/colaptes_cds/${SGE_TASK_ID}_total_output.txt
grep '^lnL' ${SGE_TASK_ID}_output_apaloderma.txt >> /lustre/scratch/johruska/colaptes/colaptes_reference_genome/colaptes_cds/${SGE_TASK_ID}_total_output.txt
grep '^w (dN/dS)' ${SGE_TASK_ID}_output_apaloderma.txt >> /lustre/scratch/johruska/colaptes/colaptes_reference_genome/colaptes_cds/${SGE_TASK_ID}_total_output.txt
grep '^lnL' ${SGE_TASK_ID}_output_buceros.txt >> /lustre/scratch/johruska/colaptes/colaptes_reference_genome/colaptes_cds/${SGE_TASK_ID}_total_output.txt
grep '^w (dN/dS)' ${SGE_TASK_ID}_output_buceros.txt >> /lustre/scratch/johruska/colaptes/colaptes_reference_genome/colaptes_cds/${SGE_TASK_ID}_total_output.txt
grep '^lnL' ${SGE_TASK_ID}_output_merops.txt >> /lustre/scratch/johruska/colaptes/colaptes_reference_genome/colaptes_cds/${SGE_TASK_ID}_total_output.txt
grep '^w (dN/dS)' ${SGE_TASK_ID}_output_merops.txt >> /lustre/scratch/johruska/colaptes/colaptes_reference_genome/colaptes_cds/${SGE_TASK_ID}_total_output.txt

# move the alignment fasta to the output directory 
mv ${SGE_TASK_ID}_aligned_trimmed.fasta ../../output
