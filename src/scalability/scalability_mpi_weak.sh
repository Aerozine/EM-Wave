#!/bin/bash -l
#SBATCH --job-name="MPI weak scaling"
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=04:00:00
#SBATCH --output="mpi_weak_scaling_slurm.log"

module purge
module load releases/2021b
module load Info0939Tools

numCpuCores=$SLURM_CPUS_ON_NODE
startDir=$(pwd)

executable="$HOME/EM-Wave/hpc_project_mpi"
args="3"

echo "Job info"
echo "--------"
echo
echo "    Job ID:" $SLURM_JOB_ID
echo " Node list:" $SLURM_JOB_NODELIST
echo "Start time:" $(date +"%d-%m-%Y %H:%M:%S")
echo

echo "Scaling experiment"
echo "MPI - weak scale"
echo "------------------"
echo

export SLURM_DISTRIBUTION=block:block

for binding in close spread;
do
  export OMP_PROC_BIND=${binding}

  numRanks=$SLURM_JOB_NUM_NODES
  grid_size=2000

  while [[ $numRanks -le $((numCpuCores * SLURM_JOB_NUM_NODES)) ]];
  do
    export OMP_NUM_THREADS=1

    if [[ ${binding} == "spread" ]]; then
      srunBind=$(seq -s ',' 0 $((numCpuCores / (numRanks / SLURM_JOB_NUM_NODES))) $((numCpuCores - 1)))
    else
      srunBind=$(seq -s ',' 0 $((numRanks / SLURM_JOB_NUM_NODES - 1)))
    fi

    workDir="mpi_scalability_weak"
    output="${binding}_${SLURM_JOB_NUM_NODES}nodes_${numRanks}ranks_${grid_size}.log"
    mkdir -p ${workDir}

    echo "[$(date +"%d-%m-%Y %H:%M:%S")] Running with binding=${binding}, numranks=${numRanks}, workdir=${workDir}"

    start_time=$(date +%s.%N)
    cd ${workDir}
    srun --nodes=$SLURM_JOB_NUM_NODES     \
         --ntasks=${numRanks}             \
         --cpu-bind="map_cpu:${srunBind}" \
         ${executable} ${args} > ${output}
    cd ${startDir}
    end_time=$(date +%s.%N)
    elapsed=$(printf "%.6f" "$(echo "$end_time - $start_time" | bc)")

    echo "[$(date +"%d-%m-%Y %H:%M:%S")] Done in ${elapsed} secs"
    echo

    (( numRanks *= 4 ))
    (( grid_size *= 2 ))
  done
done

echo "End time:" $(date +"%d-%m-%Y %H:%M:%S")
