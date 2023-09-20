#!/bin/bash

#SBATCH --partition=full

#SBATCH --job-name=T2_IMT2112
#SBATCH --output=log.out

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

mpic++ main.cpp
time mpirun a.out
