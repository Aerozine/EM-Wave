#include "solver.h"
#include "data.h"
#include "params.h"

#include <math.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

int solve(struct SimulationParams *sim_params,
          struct PhysicalParams *phys_params,
          struct PerformanceData *perf_data) {

  struct data ex, ey, ez, hx, hy, hz;
  float val = 0.;
  // nx, ny
  if (init_data(&ez, "ez", sim_params->size_of_space[0],
                sim_params->size_of_space[1], sim_params->size_of_space[2],
                sim_params->steps[0], sim_params->steps[1],
                sim_params->steps[2], val) ||
      init_data(&hx, "hx", sim_params->size_of_space[0],
                sim_params->size_of_space[1], sim_params->size_of_space[2],
                sim_params->steps[0], sim_params->steps[1],
                sim_params->steps[2], val) ||
      init_data(&hy, "hy", sim_params->size_of_space[0],
                sim_params->size_of_space[1], sim_params->size_of_space[2],
                sim_params->steps[0], sim_params->steps[1],
                sim_params->steps[2], val) ||
      init_data(&ex, "ex", sim_params->size_of_space[0],
                sim_params->size_of_space[1], sim_params->size_of_space[2],
                sim_params->steps[0], sim_params->steps[1],
                sim_params->steps[2], val) ||
      init_data(&ey, "ey", sim_params->size_of_space[0],
                sim_params->size_of_space[1], sim_params->size_of_space[2],
                sim_params->steps[0], sim_params->steps[1],
                sim_params->steps[2], val) ||
      init_data(&hz, "hz", sim_params->size_of_space[0],
                sim_params->size_of_space[1], sim_params->size_of_space[2],
                sim_params->steps[0], sim_params->steps[1],
                sim_params->steps[2], val)) {
    printf("Error: could not allocate data\n");
    return EXIT_FAILURE;
  }

  double start = GET_TIME();

  for (int n = 0; n < sim_params->size_of_space[sim_params->ndim]; n++) {
    if (n && (n % (sim_params->size_of_space[sim_params->ndim] / 10)) == 0) {
      double time_sofar = GET_TIME() - start;
      double eta =
          (sim_params->size_of_space[sim_params->ndim] - n) * time_sofar / n;
#ifndef STABILITY_STUDY
      printf("Computing time step %d/%d (ETA: %g seconds)     \r", n,
             sim_params->size_of_space[sim_params->ndim], eta);
      fflush(stdout);
#endif
    }

    // update hx and hy
    double chx1 = sim_params->steps[sim_params->ndim] /
                  (sim_params->steps[2] * phys_params->mu);
    double chx2 = sim_params->steps[sim_params->ndim] /
                  (sim_params->steps[1] * phys_params->mu);
    for (int k = 0; k < sim_params->size_of_space[2] - 1; k++) {
      for (int j = 0; j < sim_params->size_of_space[1] - 1; j++) {
        for (int i = 0; i < sim_params->size_of_space[0]; i++) {
          double hx_ij = GET(&hx, i, j, k) +
                         chx1 * (GET(&ey, i, j, k + 1) - GET(&ey, i, j, k)) -
                         chx2 * (GET(&ez, i, j + 1, k) - GET(&ez, i, j, k));
          SET(&hx, i, j, k, hx_ij);
        }
      }
    }

    double chy1 = sim_params->steps[sim_params->ndim] /
                  (sim_params->steps[0] * phys_params->mu);
    double chy2 = sim_params->steps[sim_params->ndim] /
                  (phys_params->mu * sim_params->steps[2]);
    for (int k = 0; k < sim_params->size_of_space[2] - 1; k++) {
      for (int j = 0; j < sim_params->size_of_space[1]; j++) {
        for (int i = 0; i < sim_params->size_of_space[0] - 1; i++) {
          double hy_ij = GET(&hy, i, j, k) +
                         chy1 * (GET(&ez, i + 1, j, k) - GET(&ez, i, j, k)) -
                         chy2 * (GET(&ex, i, j, k + 1) - GET(&ex, i, j, k));
          SET(&hy, i, j, k, hy_ij);
        }
      }
    }

    double chz1 = sim_params->steps[sim_params->ndim] /
                  (phys_params->mu * sim_params->steps[1]);
    double chz2 = sim_params->steps[sim_params->ndim] /
                  (phys_params->mu * sim_params->steps[0]);

    for (int k = 0; k < sim_params->size_of_space[2]; k++) {
      for (int j = 0; j < sim_params->size_of_space[1] - 1; j++) {
        for (int i = 0; i < sim_params->size_of_space[0] - 1; i++) {
          double hz_ij = GET(&hz, i, j, k) +
                         chz1 * (GET(&ex, i, j + 1, k) - GET(&ex, i, j, k)) -
                         chz2 * (GET(&ey, i + 1, j, k) - GET(&ey, i, j, k));
          SET(&hz, i, j, k, hz_ij);
        }
      }
    }

    double cex1 = sim_params->steps[sim_params->ndim] /
                  (phys_params->eps * sim_params->steps[1]);
    double cex2 = sim_params->steps[sim_params->ndim] /
                  (phys_params->eps * sim_params->steps[2]);
    for (int k = 1; k < sim_params->size_of_space[2]; k++) {
      for (int j = 1; j < sim_params->size_of_space[1]; j++) {
        for (int i = 0; i < sim_params->size_of_space[0]; i++) {
          double ex_ij = GET(&ex, i, j, k) +
                         cex1 * (GET(&hz, i, j, k) - GET(&hz, i, j - 1, k)) -
                         cex2 * (GET(&hy, i, j, k) - GET(&hy, i, j, k - 1));
          SET(&ex, i, j, k, ex_ij);
        }
      }
    }

    double cey1 = sim_params->steps[sim_params->ndim] /
                  (phys_params->eps * sim_params->steps[2]);
    double cey2 = sim_params->steps[sim_params->ndim] /
                  (phys_params->eps * sim_params->steps[0]);
    for (int k = 1; k < sim_params->size_of_space[2]; k++) {
      for (int j = 0; j < sim_params->size_of_space[1]; j++) {
        for (int i = 1; i < sim_params->size_of_space[0]; i++) {
          double ey_ij = GET(&ey, i, j, k) +
                         cey1 * (GET(&hx, i, j, k) - GET(&hx, i, j, k - 1)) -
                         cey2 * (GET(&hz, i, j, k) - GET(&hz, i - 1, j, k));
          SET(&ey, i, j, k, ey_ij);
        }
      }
    }

    // update ez
    double cez1 = sim_params->steps[sim_params->ndim] /
                  (sim_params->steps[0] * phys_params->eps),
           cez2 = sim_params->steps[sim_params->ndim] /
                  (sim_params->steps[1] * phys_params->eps);
    for (int k = 0; k < sim_params->size_of_space[2]; k++) {
      for (int j = 1; j < sim_params->size_of_space[1]; j++) {
        for (int i = 1; i < sim_params->size_of_space[0]; i++) {
          double ez_ij = GET(&ez, i, j, k) +
                         cez1 * (GET(&hy, i, j, k) - GET(&hy, i - 1, j, k)) -
                         cez2 * (GET(&hx, i, j, k) - GET(&hx, i, j - 1, k));
          SET(&ez, i, j, k, ez_ij);
        }
      }
    }

#ifdef STABILITY_STUDY
    // +42 is to decentrate from the source
    // maybe using FLT_DIG ?
    printf("%.12f \n", GET(&ez, (sim_params->size_of_space[0] >> 1) + 42,
                           (sim_params->size_of_space[1] >> 1) + 42));
#endif
    // impose source
    double t = n * sim_params->steps[sim_params->ndim];
    switch (sim_params->problem_id) {
    case 1:
    case 2:
      // sinusoidal excitation at 2.4 GHz in the middle of the domain
      SET(&ez, sim_params->size_of_space[0] / 2,
          sim_params->size_of_space[1] / 2, sim_params->size_of_space[2] / 2,
          sin(2. * M_PI * 2.4e9 * t));
      break;
    default:
      printf("Error: unknown source\n");
      break;
    }

    // output step data in VTK format
    if (sim_params->sampling_rate && !(n % sim_params->sampling_rate)) {
      write_data_vtk(&ex, n, 0);
      write_data_vtk(&ey, n, 0);
      write_data_vtk(&ez, n, 0);
      // write_data_vtk(&hx, n, 0);
      // write_data_vtk(&hy, n, 0);
    }
  }

  // write VTK manifest, linking to individual step data files
  write_manifest_vtk("ex", sim_params->steps[sim_params->ndim],
                     sim_params->size_of_space[sim_params->ndim],
                     sim_params->sampling_rate, 1);
  write_manifest_vtk("ey", sim_params->steps[sim_params->ndim],
                     sim_params->size_of_space[sim_params->ndim],
                     sim_params->sampling_rate, 1);
  write_manifest_vtk("ez", sim_params->steps[sim_params->ndim],
                     sim_params->size_of_space[sim_params->ndim],
                     sim_params->sampling_rate, 1);
  // write_manifest_vtk("hx", sim_params->steps[sim_params->ndim],
  // sim_params->size_of_space[sim_params->ndim], sampling_rate, 1);
  // write_manifest_vtk("hy", sim_params->steps[sim_params->ndim],
  // sim_params->size_of_space[sim_params->ndim], sampling_rate, 1);

  double time = GET_TIME() - start;
#ifndef STABILITY_STUDY

  printf("\nDone: %g seconds (%g MUpdates/s)\n", time,
         1.e-6 * (double)sim_params->size_of_space[0] *
             (double)sim_params->size_of_space[1] *
             (double)sim_params->size_of_space[2] *
             (double)sim_params->size_of_space[sim_params->ndim] / time);

#endif
  free_data(&ex);
  free_data(&ey);
  free_data(&ez);
  free_data(&hx);
  free_data(&hy);
  free_data(&hz);

  return EXIT_SUCCESS;
}