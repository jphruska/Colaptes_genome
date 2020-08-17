#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N colaptes_round_1_maker_retry_2
#$ -q omni
#$ -pe mpi 108
#$ -P quanah
#$ -l h_rt=48:00:00

mpiexec -n 108 /home/jmanthey/maker/bin/maker -base Caur_rnd1.2 maker_round1_opts.ctl maker_bopts.ctl maker_exe.ctl -tries 5
