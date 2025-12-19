#!/bin/bash
#SBATCH --job-name="My first job"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=10:00
#SBATCH --output=firstjob.out

cd /home/ulg/info0939/qbinstok/HPC-EM-wave-propagation
make run
