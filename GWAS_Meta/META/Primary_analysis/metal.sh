#!/bin/bash

#SBATCH --job-name=meta
#SBATCH --error=%x-%j.error
#SBATCH --out=%x-%j.out
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-0:0:0


/home/lyul1/generic-metal/metal      metal_cd.txt
/home/lyul1/generic-metal/metal      metal_uc.txt
/home/lyul1/generic-metal/metal      metal_ibd.txt
