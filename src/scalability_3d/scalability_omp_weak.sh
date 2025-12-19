#!/bin/bash
#SBATCH --job-name="OMP weak scalability"
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=08:00:00
#SBATCH --output="OMP_weak_scaling_3d_slurm.log"

module load Info0939Tools
numCpuCores=$SLURM_CPUS_ON_NODE
startDir=$(pwd)
executable="$HOME/EM-Wave/hpc_project_omp_3d"
arg="3"

for binding in close spread; do
    export OMP_PROC_BIND=${binding}
    
    numThreads=1
    grid_size=160
    while [[ $numThreads -le $numCpuCores ]]; do
        export OMP_NUM_THREADS=${numThreads}

        workDir="omp_scalability_weak_3d"
        output="${binding}_${numThreads}_threads_${grid_size}_size.log"

        mkdir -p ${workDir} && cd ${workDir}
        ${executable} ${arg} ${grid_size} > ${output}
        cd ${startDir}

        (( numThreads *= 2 ))

        grid_size=$(awk -v g="$grid_size" 'BEGIN { printf "%.0f\n", g * (2)^(1/3) }')
    done
done