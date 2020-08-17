#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N colap_scaf1
#$ -q omni
#$ -pe sm 24
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G

module load gnu

cd working1

/lustre/work/johruska/3d-dna/run-asm-pipeline.sh -r 2 -i 25000 --sort-output \
--editor-coarse-resolution 100000 --editor-coarse-region 300000 --editor-fine-resolution 25000 \
--editor-repeat-coverage 3 \
--polisher-coarse-resolution 100000 --polisher-coarse-region 300000 --polisher-fine-resolution 25000 \
--splitter-coarse-resolution 100000 --splitter-coarse-region 300000 --splitter-fine-resolution 25000 \
/lustre/work/johruska/juicer/references/c05_colaptes_arcs3.fasta /lustre/scratch/johruska/colaptes/colaptes_reference_genome/Colaptes_HI-C_2/aligned/merged_nodups.txt
