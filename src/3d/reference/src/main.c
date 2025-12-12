#include "main.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
// #include <string.h>
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

#define GET(data, i, j, k)                                                     \
  ((data)->values[(data->nx * data->ny) * k + (data)->nx * (j) + (i)])
#define SET(data, i, j, k, val)                                                \
  ((data)->values[((data)->nx * (data)->ny) * k + (data)->nx * (j) + (i)] =    \
       (val))

void print_perf_data(struct PerformanceData *perf_data) {
  fprintf(stdout, "Performance : \n");
  fprintf(stdout, "\t- Time : %g\n", perf_data->time);
  fprintf(stdout, "\t- MUpdates/s : %g\n", perf_data->MUps_per_sec);
}

void print_sim_params(struct SimulationParams *sim_params) {
  printf("Solving problem %d:\n", sim_params->problem_id);
  printf(" - space %gm x %gm x %gm (dx=%g, dy=%g, dz=%g; "
         "nx=%d, ny=%d)\n",
         sim_params->steps[0] * sim_params->size_of_space[0],
         sim_params->steps[1] * sim_params->size_of_space[1],
         sim_params->steps[2] * sim_params->size_of_space[2],
         sim_params->steps[0], sim_params->steps[1], sim_params->steps[2],
         sim_params->size_of_space[0], sim_params->size_of_space[1]);
  printf(" - time %gs (dt=%g, nt=%d)\n",
         sim_params->steps[sim_params->ndim] *
             sim_params->size_of_space[sim_params->ndim],
         sim_params->steps[sim_params->ndim],
         sim_params->size_of_space[sim_params->ndim]);
}

int main(int argc, char **argv) {
  if ((argc != 2 && argc != 3) ||
      (argc == 2 && (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help") ||
                     !strcmp(argv[1], "help")))) {
    printf("Usage: %s <problem_id> [size of cube grid, eg 500] \n", argv[0]);
    return EXIT_FAILURE;
  }

  struct SimulationParams sim_params;
  struct PerformanceData perf_data;
  struct PhysicalParams phys_params;

  init_params(&sim_params, &perf_data, &phys_params);
  sim_params.problem_id = atoi(argv[1]);

  if (set_params(&sim_params)) {
    free_params(&sim_params);
    return EXIT_FAILURE;
  }

  if (argc == 3) {
    for (int i = 0; i < sim_params.ndim; i++)
      sim_params.size_of_space[i] = atoi(argv[2]);
  }

  print_sim_params(&sim_params);

  if (solve(&sim_params, &phys_params, &perf_data)) {
    free_params(&sim_params);
    return EXIT_FAILURE;
  }

  print_perf_data(&perf_data);

  free_params(&sim_params);

  DEBUG_PRINT("Finished the computation.\n");

  return EXIT_SUCCESS;
}
