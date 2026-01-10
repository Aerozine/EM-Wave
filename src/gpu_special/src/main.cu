#include "data.h"
#include "main.h"
#include "solver.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#if defined(_OPENMP)
#include <omp.h>
#define GET_TIME() (omp_get_wtime()) // wall time
#else
#define GET_TIME() ((real)clock() / CLOCKS_PER_SEC) // cpu time
#endif

#define GET(data, i, j) ((data)->values[(data)->nx * (j) + (i)])
#define SET(data, i, j, val) ((data)->values[(data)->nx * (j) + (i)] = (val))

// Uniform permittivity (vacuum/air)
real eps_uniform(real x, real y, real eps0) { return eps0; }

// Vertical barrier in the middle of the domain
real eps_vertical_barrier(real x, real y, real eps0) {
  // Domain size comes from actual simulation parameters
  // For problem 3: dx=0.00625, nx=500 → domain = 3.125m
  real domain_size = 3.125; // meters
  real domain_center_x = domain_size / 2.0;
  real barrier_width = domain_size * 0.04; // 4% of domain width

  // Higher permittivity in the barrier region (offset to the right of source)
  real barrier_x = domain_center_x + domain_size * 0.2;
  if (fabs(x - barrier_x) < barrier_width / 2.0) {
    return eps0 * 25.0; // Relative permittivity of 25 (very high contrast!)
  }
  return eps0;
}

// TWO parallel barriers - creates a waveguide/cavity effect
real eps_two_barriers(real x, real y, real eps0) {
  real domain_size = 3.125;
  real domain_center_x = domain_size / 2.0;
  real barrier_width = domain_size * 0.04;

  // First barrier to the left of source
  real barrier1_x = domain_center_x - domain_size * 0.15;
  // Second barrier to the right of source
  real barrier2_x = domain_center_x + domain_size * 0.25;

  if (fabs(x - barrier1_x) < barrier_width / 2.0 ||
      fabs(x - barrier2_x) < barrier_width / 2.0) {
    return eps0 * 25.0; // Very high permittivity
  }
  return eps0;
}

// Cross-shaped obstacle
real eps_cross_obstacle(real x, real y, real eps0) {
  real domain_size = 3.125;
  real domain_center_x = domain_size / 2.0;
  real domain_center_y = domain_size / 2.0;

  real obstacle_x = domain_center_x + domain_size * 0.2;
  real arm_width = domain_size * 0.08;
  real arm_length = domain_size * 0.3;

  // Vertical arm
  bool in_vertical = (fabs(x - obstacle_x) < arm_width / 2.0) &&
                     (fabs(y - domain_center_y) < arm_length / 2.0);

  // Horizontal arm
  bool in_horizontal = (fabs(y - domain_center_y) < arm_width / 2.0) &&
                       (fabs(x - obstacle_x) < arm_length / 2.0);

  if (in_vertical || in_horizontal) {
    return eps0 * 20.0;
  }
  return eps0;
}

// Circular dielectric obstacle
real eps_circular_obstacle(real x, real y, real eps0) {
  real domain_size = 3.125;
  real domain_center_x = domain_size / 2.0;
  real domain_center_y = domain_size / 2.0;
  real obstacle_radius = domain_size * 0.12;

  // Offset to the right of source
  real obstacle_x = domain_center_x + domain_size * 0.25;

  real dx = x - obstacle_x;
  real dy = y - domain_center_y;
  real r = sqrt(dx * dx + dy * dy);

  if (r < obstacle_radius) {
    return eps0 * 20.0; // Very high permittivity - strong reflection!
  }
  return eps0;
}

// Multiple circles - scattering array
real eps_circle_array(real x, real y, real eps0) {
  real domain_size = 3.125;
  real domain_center_x = domain_size / 2.0;
  real domain_center_y = domain_size / 2.0;
  real circle_radius = domain_size * 0.08;

  // Array of 4 circles
  real positions[4][2] = {{domain_center_x + domain_size * 0.15,
                           domain_center_y + domain_size * 0.15},
                          {domain_center_x + domain_size * 0.15,
                           domain_center_y - domain_size * 0.15},
                          {domain_center_x + domain_size * 0.30,
                           domain_center_y + domain_size * 0.15},
                          {domain_center_x + domain_size * 0.30,
                           domain_center_y - domain_size * 0.15}};

  for (int i = 0; i < 4; i++) {
    real dx = x - positions[i][0];
    real dy = y - positions[i][1];
    real r = sqrt(dx * dx + dy * dy);

    if (r < circle_radius) {
      return eps0 * 15.0;
    }
  }
  return eps0;
}

// Multiple barriers (double slit)
real eps_double_slit(real x, real y, real eps0) {
  real domain_size = 3.125;
  real domain_center_x = domain_size / 2.0;
  real domain_center_y = domain_size / 2.0;

  real barrier_x =
      domain_center_x + domain_size * 0.25; // To the right of source
  real barrier_width = domain_size * 0.04;
  real slit_width = domain_size * 0.10;
  real slit_separation = domain_size * 0.35;

  // In the barrier region
  if (fabs(x - barrier_x) < barrier_width / 2.0) {
    // Check if we're in a slit
    real y_from_center = y - domain_center_y;
    real slit_center_offset = slit_separation / 2.0;

    bool in_slit1 = fabs(y_from_center - slit_center_offset) < slit_width / 2.0;
    bool in_slit2 = fabs(y_from_center + slit_center_offset) < slit_width / 2.0;

    if (in_slit1 || in_slit2) {
      return eps0; // Air in slits
    }
    return eps0 * 30.0; // Very high permittivity barrier - strong contrast!
  }
  return eps0;
}

// Checkerboard pattern - creates interesting diffraction
real eps_checkerboard(real x, real y, real eps0) {
  real domain_size = 3.125;
  real domain_center_x = domain_size / 2.0;
  real domain_center_y = domain_size / 2.0;

  // Start checkerboard to the right of source
  real start_x = domain_center_x + domain_size * 0.15;
  real end_x = domain_center_x + domain_size * 0.40;

  if (x >= start_x && x <= end_x) {
    real cell_size = domain_size * 0.05;
    int cell_i = (int)((x - start_x) / cell_size);
    int cell_j = (int)(y / cell_size);

    // Checkerboard pattern
    if ((cell_i + cell_j) % 2 == 0) {
      return eps0 * 16.0;
    }
  }
  return eps0;
}

// Gradient permittivity (smoothly varying)
real eps_gradient(real x, real y, real eps0) {
  real domain_size = 3.125;
  // Strong exponential gradient from left to right
  real rel_eps = 1.0 + 9.0 * (x / domain_size);
  return eps0 * rel_eps;
}

// Lens shape - focusing effect!
real eps_lens(real x, real y, real eps0) {
  real domain_size = 3.125;
  real domain_center_x = domain_size / 2.0;
  real domain_center_y = domain_size / 2.0;

  // Biconvex lens positioned to the right of source
  real lens_center_x = domain_center_x + domain_size * 0.25;
  real lens_center_y = domain_center_y;

  // Lens parameters
  real lens_max_radius = domain_size * 0.25; // Maximum radius at center
  real lens_thickness = domain_size * 0.12;  // Thickness at center

  // Distance from lens center
  real dx = x - lens_center_x;
  real dy = y - lens_center_y;
  real r = sqrt(dy * dy); // Radial distance from lens axis

  // Create biconvex lens shape
  // Lens is thicker in the center, thinner at edges
  if (r < lens_max_radius) {
    // Radius of curvature for the lens surfaces
    real R = lens_max_radius * 1.2; // Radius of curvature

    // Calculate the lens thickness at this radial position
    // Using circular arc formula: thickness varies as sqrt(R^2 - r^2)
    real edge_offset = sqrt(R * R - lens_max_radius * lens_max_radius);
    real local_offset = sqrt(R * R - r * r);
    real half_thickness = (local_offset - edge_offset) * 0.5;

    // Check if we're inside the lens volume
    if (fabs(dx) < half_thickness) {
      return eps0 * 16.0; // High permittivity for focusing
    }
  }

  return eps0;
}

// Alternative: Graded-index (GRIN) lens - smoother focusing
real eps_grin_lens(real x, real y, real eps0) {
  real domain_size = 3.125;
  real domain_center_x = domain_size / 2.0;
  real domain_center_y = domain_size / 2.0;

  real lens_center_x = domain_center_x + domain_size * 0.25;
  real lens_center_y = domain_center_y;
  real lens_radius = domain_size * 0.2;

  real dx = x - lens_center_x;
  real dy = y - lens_center_y;
  real r = sqrt(dx * dx + dy * dy);

  // Graded index: higher permittivity in center, lower at edges
  // This creates a parabolic index profile
  if (r < lens_radius) {
    real normalized_r = r / lens_radius;
    // Parabolic profile: n^2 = n0^2 * (1 - (r/a)^2)
    real rel_eps = 16.0 * (1.0 - 0.6 * normalized_r * normalized_r);
    if (rel_eps < 1.0)
      rel_eps = 1.0;
    return eps0 * rel_eps;
  }

  return eps0;
}

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
  case 3: // Problem with vertical barrier
    params->dx = params->dy = (3.e8 / 2.4e9) / 20.; // wavelength / 20
    params->nx = params->ny = 500;
    params->dt = 0.5 / (3.e8 * sqrt(1. / (params->dx * params->dx) +
                                    1. / (params->dy * params->dy))); // cfl / 2
    params->nt = 1500;
    params->sampling_rate = 5; // save 1 step out of 5
    break;
  case 4: // Circular obstacle
    params->dx = params->dy = (3.e8 / 2.4e9) / 20.;
    params->nx = params->ny = 500;
    params->dt = 0.5 / (3.e8 * sqrt(1. / (params->dx * params->dx) +
                                    1. / (params->dy * params->dy)));
    params->nt = 1500;
    params->sampling_rate = 5;
    break;
  case 5: // Double slit
    params->dx = params->dy = (3.e8 / 2.4e9) / 20.;
    params->nx = params->ny = 500;
    params->dt = 0.5 / (3.e8 * sqrt(1. / (params->dx * params->dx) +
                                    1. / (params->dy * params->dy)));
    params->nt = 1500;
    params->sampling_rate = 5;
    break;
  case 6: // Gradient permittivity
    params->dx = params->dy = (3.e8 / 2.4e9) / 20.;
    params->nx = params->ny = 500;
    params->dt = 0.5 / (3.e8 * sqrt(1. / (params->dx * params->dx) +
                                    1. / (params->dy * params->dy)));
    params->nt = 1500;
    params->sampling_rate = 5;
    break;
  case 7: // Two barriers (waveguide/cavity)
    params->dx = params->dy = (3.e8 / 2.4e9) / 20.;
    params->nx = params->ny = 500;
    params->dt = 0.5 / (3.e8 * sqrt(1. / (params->dx * params->dx) +
                                    1. / (params->dy * params->dy)));
    params->nt = 2000;
    params->sampling_rate = 5;
    break;
  case 8: // Circle array
    params->dx = params->dy = (3.e8 / 2.4e9) / 20.;
    params->nx = params->ny = 500;
    params->dt = 0.5 / (3.e8 * sqrt(1. / (params->dx * params->dx) +
                                    1. / (params->dy * params->dy)));
    params->nt = 1500;
    params->sampling_rate = 5;
    break;
  case 9: // Checkerboard
    params->dx = params->dy = (3.e8 / 2.4e9) / 20.;
    params->nx = params->ny = 500;
    params->dt = 0.5 / (3.e8 * sqrt(1. / (params->dx * params->dx) +
                                    1. / (params->dy * params->dy)));
    params->nt = 1500;
    params->sampling_rate = 5;
    break;
  case 10: // Cross obstacle
    params->dx = params->dy = (3.e8 / 2.4e9) / 20.;
    params->nx = params->ny = 500;
    params->dt = 0.5 / (3.e8 * sqrt(1. / (params->dx * params->dx) +
                                    1. / (params->dy * params->dy)));
    params->nt = 1500;
    params->sampling_rate = 5;
    break;
  case 11: // Lens
    params->dx = params->dy = (3.e8 / 2.4e9) / 20.;
    params->nx = params->ny = 500;
    params->dt = 0.5 / (3.e8 * sqrt(1. / (params->dx * params->dx) +
                                    1. / (params->dy * params->dy)));
    params->nt = 2000; // Longer to see focusing effect
    params->sampling_rate = 5;
    break;
  case 12: // GRIN lens (graded index)
    params->dx = params->dy = (3.e8 / 2.4e9) / 20.;
    params->nx = params->ny = 500;
    params->dt = 0.5 / (3.e8 * sqrt(1. / (params->dx * params->dx) +
                                    1. / (params->dy * params->dy)));
    params->nt = 6000;
    params->sampling_rate = 5;
    break;
  default:
    printf("Error: unknown problem id %d\n", problem_id);
    exit(EXIT_FAILURE);
  }
}

PermittivityFunc get_permittivity_func(int problem_id) {
  switch (problem_id) {
  case 1:
  case 2:
    return eps_uniform;
  case 3:
    return eps_vertical_barrier;
  case 4:
    return eps_circular_obstacle;
  case 5:
    return eps_double_slit;
  case 6:
    return eps_gradient;
  case 7:
    return eps_two_barriers;
  case 8:
    return eps_circle_array;
  case 9:
    return eps_checkerboard;
  case 10:
    return eps_cross_obstacle;
  case 11:
    return eps_lens;
  case 12:
    return eps_grin_lens;
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

  // Get the appropriate permittivity function
  PermittivityFunc eps_func = get_permittivity_func(problem_id);

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
  if (solve(&sim_params, &phys_params, problem_id, eps_func))
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}