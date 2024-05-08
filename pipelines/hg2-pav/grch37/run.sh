#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=3G
#SBATCH --partition=batch
#SBATCH --time=36:0:0

singularity run --bind "$(pwd):/mnt" --pwd /mnt library://becklab/pav/pav:2.3.4 -c 32

