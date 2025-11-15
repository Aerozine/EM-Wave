#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "params.h"

int init_params(struct SimulationParams *sim_params,
                struct PerformanceData *perf_data,
                struct PhysicalParams *phys_params,
                struct MpiParams *mpi_params) {
  sim_params->sampling_rate = 1;
  sim_params->problem_id = 1;
  sim_params->ndim = 2;
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

  mpi_params->use_mpi = false;
  mpi_params->rank = 0;
  mpi_params->num_ranks = 1;
  mpi_params->procs_per_dim = NULL;
  mpi_params->periods = NULL;
  mpi_params->coords = NULL;
  mpi_params->neighbours = NULL;
  mpi_params->sizes = NULL;
  mpi_params->send_sizes = NULL;
  mpi_params->sizes = malloc(sizeof(int) * sim_params->ndim);
  mpi_params->send_sizes = malloc(sizeof(int) * sim_params->ndim);
  if (!mpi_params->sizes || !mpi_params->send_sizes) {
    free(sim_params->size_of_space);
    free(sim_params->steps);
    free(mpi_params->sizes);
    free(mpi_params->send_sizes);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

int set_params(struct SimulationParams *sim_params) {
  switch (sim_params->problem_id) {
  case 1: // small size problem
    sim_params->steps[0] = sim_params->steps[1] =
        (3.e8 / 2.4e9) / 20.; // wavelength / 20
    sim_params->size_of_space[0] = sim_params->size_of_space[1] = 500;
    sim_params->steps[sim_params->ndim] =
        0.5 / (3.e8 * sqrt(1. / (sim_params->steps[0] * sim_params->steps[0]) +
                           1. / (sim_params->steps[1] *
                                 sim_params->steps[1]))); // cfl / 2
    sim_params->size_of_space[sim_params->ndim] = 500;
    sim_params->sampling_rate = 5; // save 1 step out of 5
    break;
  case 2: // larger size problem, usable for initial scaling tests
    sim_params->steps[0] = sim_params->steps[1] =
        (3.e8 / 2.4e9) / 40.; // wavelength / 40
    sim_params->size_of_space[0] = sim_params->size_of_space[1] = 16000;
    sim_params->steps[sim_params->ndim] =
        0.5 / (3.e8 * sqrt(1. / (sim_params->steps[0] * sim_params->steps[0]) +
                           1. / (sim_params->steps[1] *
                                 sim_params->steps[1]))); // cfl / 2
    sim_params->size_of_space[sim_params->ndim] = 500;
    sim_params->sampling_rate = 0; // don't save results
    break;
  case 3: // larger size problem, usable for initial scaling tests
    sim_params->steps[0] = sim_params->steps[1] =
        (3.e8 / 2.4e9) / 40.; // wavelength / 40
    sim_params->size_of_space[0] = sim_params->size_of_space[1] = 16000;
    sim_params->steps[sim_params->ndim] =
        0.1 / (3.e8 * sqrt(1. / (sim_params->steps[0] * sim_params->steps[0]) +
                           1. / (sim_params->steps[1] *
                                 sim_params->steps[1]))); // cfl / 2
    sim_params->size_of_space[sim_params->ndim] = 50;
    sim_params->sampling_rate = 2;
    break;
  default:
    printf("Error: unknown problem id %d\n", sim_params->problem_id);
    return EXIT_FAILURE;
    break;
  }

  return EXIT_SUCCESS;
}

void free_params(struct SimulationParams *sim_params,
                 struct MpiParams *mpi_params) {
  free(sim_params->size_of_space);
  free(sim_params->steps);

  free(mpi_params->procs_per_dim);
  free(mpi_params->periods);
  free(mpi_params->coords);
  free(mpi_params->neighbours);
  free(mpi_params->sizes);
  free(mpi_params->send_sizes);
}