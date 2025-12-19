#include "solver.h"
#include "data.h"
#include "main.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>

int solve(struct SimulationParams *sim_params,
          struct PhysicalParams *phys_params, int problem_id) {

  DEBUG_PRINT("Starting computation on process %d, with :\n\tnx from %d to "
              "%d\n\tny from %d to %d\n",
              0, 0, sim_params->nx - 1, 0, sim_params->ny - 1);

  struct data ez, hx, hy;
  if (init_data(&ez, "ez", sim_params->nx, sim_params->ny, sim_params->dx,
                sim_params->dy, 0.) ||
      init_data(&hx, "hx", sim_params->nx, sim_params->ny - 1, sim_params->dx,
                sim_params->dy, 0.) ||
      init_data(&hy, "hy", sim_params->nx - 1, sim_params->ny, sim_params->dx,
                sim_params->dy, 0.)) {
    printf("Error: could not allocate data\n");
    return EXIT_FAILURE;
  }

  double start = GET_TIME();

  for (int n = 0; n < sim_params->nt; n++) {
    if (n && (n % (sim_params->nt / 10)) == 0) {
      double time_sofar = GET_TIME() - start;
      double eta = (sim_params->nt - n) * time_sofar / n;
#ifndef STABILITY_STUDY
      printf("Computing time step %d/%d (ETA: %g seconds)     \r", n,
             sim_params->nt, eta);
      fflush(stdout);
#endif
    }

    // update hx and hy
    double chy = sim_params->dt / (sim_params->dy * phys_params->mu);
    for (int j = 0; j < sim_params->ny - 1; j++) {
      for (int i = 0; i < sim_params->nx; i++) {
        double hx_ij =
            GET(&hx, i, j) - chy * (GET(&ez, i, j + 1) - GET(&ez, i, j));
        SET(&hx, i, j, hx_ij);
      }
    }
    double chx = sim_params->dt / (sim_params->dx * phys_params->mu);
    for (int j = 0; j < sim_params->ny; j++) {
      for (int i = 0; i < sim_params->nx - 1; i++) {
        double hy_ij =
            GET(&hy, i, j) + chx * (GET(&ez, i + 1, j) - GET(&ez, i, j));
        SET(&hy, i, j, hy_ij);
      }
    }

    // update ez
    double cex = sim_params->dt / (sim_params->dx * phys_params->eps),
           cey = sim_params->dt / (sim_params->dy * phys_params->eps);
    for (int j = 1; j < sim_params->ny - 1; j++) {
      for (int i = 1; i < sim_params->nx - 1; i++) {
        double ez_ij = GET(&ez, i, j) +
                       cex * (GET(&hy, i, j) - GET(&hy, i - 1, j)) -
                       cey * (GET(&hx, i, j) - GET(&hx, i, j - 1));
        SET(&ez, i, j, ez_ij);
      }
    }

#ifdef STABILITY_STUDY
    // +42 is to decentrate from the source
    // maybe using FLT_DIG ?
    printf("%.12f \n",
           GET(&ez, (sim_params->nx >> 1) + 42, (sim_params->ny >> 1) + 42));
#endif
    // impose source
    double t = n * sim_params->dt;
    switch (problem_id) {
    case 1:
    case 2:
      // sinusoidal excitation at 2.4 GHz in the middle of the domain
      SET(&ez, sim_params->nx / 2, sim_params->ny / 2,
          sin(2. * M_PI * 2.4e9 * t));
      break;
    default:
      printf("Error: unknown source\n");
      break;
    }

    // output step data in VTK format
    if (sim_params->sampling_rate && !(n % sim_params->sampling_rate)) {
      write_data_vtk(&ez, n, 0);
      // write_data_vtk(&hx, n, 0);
      // write_data_vtk(&hy, n, 0);
    }
  }

  // write VTK manifest, linking to individual step data files
  write_manifest_vtk("ez", sim_params->dt, sim_params->nt,
                     sim_params->sampling_rate, 1);
  // write_manifest_vtk("hx", sim_params->dt, sim_params->nt, sampling_rate, 1);
  // write_manifest_vtk("hy", sim_params->dt, sim_params->nt, sampling_rate, 1);

  double time = GET_TIME() - start;
#ifndef STABILITY_STUDY

  printf("\nDone: %g seconds (%g MUpdates/s)\n", time,
         1.e-6 * (double)sim_params->nx * (double)sim_params->ny *
             (double)sim_params->nt / time);

#endif
  free_data(&ez);
  free_data(&hx);
  free_data(&hy);

  return EXIT_SUCCESS;
}
