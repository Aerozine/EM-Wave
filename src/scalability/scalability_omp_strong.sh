#!/bin/bash
#SBATCH --job-name="OMP strong scalability"
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=04:00:00
#SBATCH --output="OMP_strong_scaling_slurm.log"

module load Info0939Tools
numCpuCores=$SLURM_CPUS_ON_NODE
startDir=$(pwd)
executable="$HOME/EM-Wave/hpc_project"
args="3 20000"

for binding in close spread; do
    export OMP_PROC_BIND=${binding}
    
    numThreads=1
    while [[ $numThreads -le $numCpuCores ]]; do
        export OMP_NUM_THREADS=${numThreads}

        workDir="omp_scalability_strong"
        output="${binding}_${numThreads}_threads.log"

        mkdir -p ${workDir} && cd ${workDir}
        ${executable} ${args} > ${output}
        cd ${startDir}

        (( numThreads *= 2 ))
    done
done