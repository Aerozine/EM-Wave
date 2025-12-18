#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "main.h"
#include "solver.h"

#if defined(_OPENMP)
#include <omp.h>
#define GET_TIME() (omp_get_wtime()) // wall time
#else
#define GET_TIME() ((double)clock() / CLOCKS_PER_SEC) // cpu time
#endif

#define GET(data, i, j) ((data)->values[(data)->nx * (j) + (i)])
#define SET(data, i, j, val) ((data)->values[(data)->nx * (j) + (i)] = (val))

void init_params(struct SimulationParams *params) {
  params->dx = 1.;
  params->dy = 1.;
  params->dt = 1.;
  params->nx = 1;
  params->ny = 1;
  params->nt = 1;
  params->sampling_rate = 1;
}

void set_params(struct SimulationParams *params, int problem_id) {
  switch (problem_id) {
  case 1:                                           // small size problem
    params->dx = params->dy = (3.e8 / 2.4e9) / 20.; // wavelength / 20
    params->nx = params->ny = 500;
    params->dt = 0.5 / (3.e8 * sqrt(1. / (params->dx * params->dx) +
                                    1. / (params->dy * params->dy))); // cfl / 2
    params->nt = 500;
    params->sampling_rate = 5; // save 1 step out of 5
    break;
  case 2: // larger size problem, usable for initial scaling tests
    params->dx = params->dy = (3.e8 / 2.4e9) / 40.; // wavelength / 40
    params->nx = params->ny = 16000;
    params->dt = 0.5 / (3.e8 * sqrt(1. / (params->dx * params->dx) +
                                    1. / (params->dy * params->dy))); // cfl / 2
    params->nt = 500;
    params->sampling_rate = 0; // don't save results
    break;
  default:
    printf("Error: unknown problem id %d\n", problem_id);
    exit(EXIT_FAILURE);
  }
}

int main(int argc, char **argv) {
#ifndef STABILITY_STUDY

  if (argc != 2) {
    printf("Usage: %s problem_id\n", argv[0]);
    return 1;
  }
#else

  if (argc != 3) {
    printf("Usage: %s problem_id dt \n", argv[0]);
    return 1;
  }
#endif

  struct SimulationParams sim_params;
  init_params(&sim_params);

  struct PhysicalParams phys_params;
  phys_params.eps = 8.854187817e-12;
  phys_params.mu = 1.2566370614359173e-06;

  int problem_id = atoi(argv[1]);
  set_params(&sim_params, problem_id);

#ifdef STABILITY_STUDY

  sim_params.dt =
      ((float)atof(argv[2])) /
      (3.e8 * sqrt(1. / (sim_params.dx * sim_params.dx) +
                   1. / (sim_params.dy * sim_params.dy))); // cfl / 2
#endif /* ifndef STABILITY_STUDY */
#ifndef STABILITY_STUDY
  printf("Solving problem %d:\n", problem_id);
  printf(" - space %gm x %gm (sim_params.dx=%g, sim_params.dy=%g; "
         "sim_params.nx=%d, sim_params.ny=%d)\n",
         sim_params.dx * sim_params.nx, sim_params.dy * sim_params.ny,
         sim_params.dx, sim_params.dy, sim_params.nx, sim_params.ny);
  printf(" - time %gs (sim_params.dt=%g, sim_params.nt=%d)\n",
         sim_params.dt * sim_params.nt, sim_params.dt, sim_params.nt);
#endif
  if (solve(&sim_params, &phys_params, problem_id))
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
