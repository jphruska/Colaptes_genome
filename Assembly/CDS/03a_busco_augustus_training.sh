#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N busco_augustus_training_take_3_3.0.2
#$ -q omni
#$ -pe sm 36
#$ -P xlquanah
#$ -l h_rt=120:00:00

source activate busco.3.0.2

export BUSCO_CONFIG_FILE="/lustre/scratch/johruska/colaptes/colaptes_reference_genome/busco/busco_3.0.2/config.ini"

run_busco -i Colaptes_rnd1.all.maker.transcripts1000.fasta -c 36 -o busco_output_round_1_3.0.2 -l /home/jmanthey/busco/tetrapoda_odb9/ --mode genome -sp human --long --augustus_parameters='--progress=true'

