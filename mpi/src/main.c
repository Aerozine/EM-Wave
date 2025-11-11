#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "params.h"
#include "solver.h"

#if defined(_OPENMP)
#include <omp.h>
#define GET_TIME() (omp_get_wtime()) // wall time
#else
#define GET_TIME() ((double)clock() / CLOCKS_PER_SEC) // cpu time
#endif

#define GET(data, i, j) ((data)->values[(data)->nx * (j) + (i)])
#define SET(data, i, j, val) ((data)->values[(data)->nx * (j) + (i)] = (val))

void print_perf_data(struct PerformanceData *perf_data,
                     struct MpiParams *mpi_params) {
  fprintf(stdout, "Performance : \n");
  fprintf(stdout, "\t- Time : %g\n", perf_data->time);
  fprintf(stdout, "\t- MUpdates/s : %g\n",
          perf_data->MUps_per_sec * mpi_params->num_ranks);
}

void print_sim_params(struct SimulationParams *sim_params) {
  printf("Solving problem %d:\n", sim_params->problem_id);
  printf(" - space %gm x %gm (dx=%g, dy=%g; "
         "nx=%d, ny=%d)\n",
         sim_params->steps[0] * sim_params->size_of_space[0],
         sim_params->steps[1] * sim_params->size_of_space[1],
         sim_params->steps[0], sim_params->steps[1],
         sim_params->size_of_space[0], sim_params->size_of_space[1]);
  printf(" - time %gs (dt=%g, nt=%d)\n",
         sim_params->steps[sim_params->ndim] *
             sim_params->size_of_space[sim_params->ndim],
         sim_params->steps[sim_params->ndim],
         sim_params->size_of_space[sim_params->ndim]);
}

int main(int argc, char **argv) {
  if (argc != 2 && argc != 3) {
    printf("Usage: %s <problem_id> [--use-mpi]\n", argv[0]);
    return EXIT_FAILURE;
  }

  struct SimulationParams sim_params;
  struct PerformanceData perf_data;
  struct PhysicalParams phys_params;
  struct MpiParams mpi_params;

  init_params(&sim_params, &perf_data, &phys_params, &mpi_params);
  sim_params.problem_id = atoi(argv[1]);

  if (argc == 3 && !strcmp(argv[2], "--use-mpi"))
    mpi_params.use_mpi = true;

  if (mpi_params.use_mpi) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &(mpi_params.rank));
    MPI_Comm_size(MPI_COMM_WORLD, &(mpi_params.num_ranks));
  }

  if (set_params(&sim_params)) {
    if (mpi_params.use_mpi)
      MPI_Finalize();

    free_params(&sim_params, &mpi_params);
    return EXIT_FAILURE;
  }

  if (mpi_params.use_mpi && mpi_params.rank == 0)
    print_sim_params(&sim_params);

  if (solve(&sim_params, &phys_params, &perf_data, &mpi_params)) {
    if (mpi_params.use_mpi)
      MPI_Finalize();

    free_params(&sim_params, &mpi_params);
    return EXIT_FAILURE;
  }

  if (mpi_params.use_mpi && mpi_params.rank == 0)
    print_perf_data(&perf_data, &mpi_params);

  if (mpi_params.use_mpi)
    MPI_Finalize();

  free_params(&sim_params, &mpi_params);

  return EXIT_SUCCESS;
}
