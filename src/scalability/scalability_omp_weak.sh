#!/bin/bash
#SBATCH --job-name="OMP weak scalability"
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --output="OMP_weak_scaling_slurm.log"

module load Info0939Tools
numCpuCores=$SLURM_CPUS_ON_NODE
startDir=$(pwd)
executable="$HOME/EM-Wave/hpc_project"
arg="3"

for binding in close spread; do
    export OMP_PROC_BIND=${binding}
    
    numThreads=1
    grid_size=2000
    while [[ $numThreads -le $numCpuCores ]]; do
        export OMP_NUM_THREADS=${numThreads}

        workDir="omp_scalability_weak"
        output="${binding}_${numThreads}_threads_${grid_size}_size.log"

        mkdir -p ${workDir} && cd ${workDir}
        ${executable} ${arg} ${grid_size} > ${output}
        cd ${startDir}

        (( numThreads *= 4 ))
        (( grid_size *= 2 ))
    done
done