#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "params.h"

int init_params(struct SimulationParams *sim_params,
                struct PerformanceData *perf_data,
                struct PhysicalParams *phys_params) {
  sim_params->sampling_rate = 1;
  sim_params->problem_id = 1;
  sim_params->ndim = 3;
  sim_params->size_of_space = NULL;
  sim_params->steps = NULL;
  sim_params->size_of_space = malloc(sizeof(int) * (sim_params->ndim + 1));
  sim_params->steps = malloc(sizeof(double) * (sim_params->ndim + 1));
  if (sim_params->size_of_space == NULL || sim_params->steps == NULL) {
    free(sim_params->size_of_space);
    free(sim_params->steps);
    return EXIT_FAILURE;
  }

  for (int i = 0; i < sim_params->ndim + 1; i++) {
    sim_params->size_of_space[i] = 1;
    sim_params->steps[i] = 1.;
  }

  perf_data->time = -1;
  perf_data->MUps_per_sec = -1;

  phys_params->eps = 8.854187817e-12;
  phys_params->mu = 1.2566370614359173e-06;

  return EXIT_SUCCESS;
}

int set_params(struct SimulationParams *sim_params) {
  switch (sim_params->problem_id) {
  case 1: // small size problem
    for (int i = 0; i < sim_params->ndim; i++) {
      sim_params->size_of_space[i] = 500;
      sim_params->steps[i] = (3.e8 / 2.4e9) / 20.;
    }
    // To redo according to rediscretized scheme
    // TODO
    sim_params->steps[sim_params->ndim] =
        0.5 / (3.e8 * sqrt(1. / (sim_params->steps[0] * sim_params->steps[0]) +
                           1. / (sim_params->steps[1] * sim_params->steps[1]) +
                           1. / (sim_params->steps[2] *
                                 sim_params->steps[2]))); // cfl / 2
    sim_params->size_of_space[sim_params->ndim] = 500;
    sim_params->sampling_rate = 5; // save 1 step out of 5
    break;
  case 2: // larger size problem, usable for initial scaling tests
    for (int i = 0; i < sim_params->ndim; i++) {
      sim_params->size_of_space[i] = 16000;
      sim_params->steps[i] = (3.e8 / 2.4e9) / 20.;
    }
    sim_params->steps[sim_params->ndim] =
        0.5 / (3.e8 * sqrt(1. / (sim_params->steps[0] * sim_params->steps[0]) +
                           1. / (sim_params->steps[1] * sim_params->steps[1]) +
                           1. / (sim_params->steps[2] *
                                 sim_params->steps[2]))); // cfl / 2
    sim_params->size_of_space[sim_params->ndim] = 500;
    sim_params->sampling_rate = 0; // don't save results
    break;
  case 3: // larger size problem, usable for initial scaling tests
    for (int i = 0; i < sim_params->ndim; i++) {
      sim_params->size_of_space[i] = 100;
      sim_params->steps[i] = (3.e8 / 2.4e9) / 40.;
    }
    sim_params->steps[sim_params->ndim] =
        0.5 / (3.e8 * sqrt(1. / (sim_params->steps[0] * sim_params->steps[0]) +
                           1. / (sim_params->steps[1] * sim_params->steps[1]) +
                           1. / (sim_params->steps[2] *
                                 sim_params->steps[2]))); // cfl / 2
    sim_params->size_of_space[sim_params->ndim] = 500;
    sim_params->sampling_rate = 0;
    break;
  default:
    printf("Error: unknown problem id %d\n", sim_params->problem_id);
    return EXIT_FAILURE;
    break;
  }

  return EXIT_SUCCESS;
}

void free_params(struct SimulationParams *sim_params) {
  free(sim_params->size_of_space);
  free(sim_params->steps);
}