#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N colaptes_round_2_maker
#$ -q omni
#$ -pe mpi 108
#$ -P quanah
#$ -l h_rt=48:00:00

mpiexec -n 108 /home/jmanthey/maker/bin/maker -base Caur_rnd2 maker_round2_opts.ctl maker_bopts.ctl maker_exe_round2.ctl
