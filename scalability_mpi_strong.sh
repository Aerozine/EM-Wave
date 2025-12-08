#!/bin/bash
#SBATCH --job-name="MPI strong scalability"
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=02:00:00

module load Info0939Tools
numCpuCores=$SLURM_CPUS_ON_NODE
startDir=$(pwd)
executable="/home/ulg/info0939/qbinstok/EM-Wave/hpc_project"
args="3 5000"

for binding in close spread; do
    export OMP_PROC_BIND=${binding}
    
    numThreads=1
    while [[ $numThreads -le $numCpuCores ]]; do
        export OMP_NUM_THREADS=${numThreads}

        workDir="omp_${binding}_${numThreads}"
        output="${numThreads}_threads.out"

        mkdir -p ${workDir} && cd ${workDir}
        ${executable} ${args} > ${output}
        cd ${startDir}

        (( numThreads *= 2 ))
    done
done