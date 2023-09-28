#!/bin/bash

#SBATCH --partition=full

#SBATCH --job-name=debug
#SBATCH --output=debug.out

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

mpic++ main.cpp
time mpirun a.out
