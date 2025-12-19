#include "data.h"
#include "main.h"
#include "solver.h"

#include <cuda_runtime.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>

// Block dimensions
#define BLOCK_X 32
#define BLOCK_Y 16
__constant__ real c_chx, c_chy, c_dt_over_dx, c_dt_over_dy, c_dt_over_mu;

__global__ void apply_source_kern(real *__restrict__ ez, size_t pitch,
                                  int src_i, int src_j, real value) {
  if (threadIdx.x != 0 || blockIdx.x != 0)
    return;
  ez[PIDX(src_i, src_j, pitch)] = value;
}

__global__ void upd_h_kern(real *__restrict__ hx, size_t hx_pitch,
                           real *__restrict__ hy, size_t hy_pitch,
                           const real *__restrict__ ez, size_t ez_pitch, int nx,
                           int ny) {

  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  __shared__ real s_ez[BLOCK_Y + 1][BLOCK_X + 1];

  // Load i,j tile
  if (i < nx && j < ny) {
    s_ez[threadIdx.y][threadIdx.x] = ez[PIDX(i, j, ez_pitch)];
  }

  // Load j+1 for hx
  if (threadIdx.y == blockDim.y - 1 && i < nx && j + 1 < ny) {
    s_ez[threadIdx.y + 1][threadIdx.x] = ez[PIDX(i, j + 1, ez_pitch)];
  }

  // Load i+1 for hy
  if (threadIdx.x == blockDim.x - 1 && i + 1 < nx && j < ny) {
    s_ez[threadIdx.y][threadIdx.x + 1] = ez[PIDX(i + 1, j, ez_pitch)];
  }
  __syncthreads();

  real ez_curr = s_ez[threadIdx.y][threadIdx.x];

  // update hx: hx[i,j] -= chy * (ez[i,j+1] - ez[i,j])
  if (i < nx && j < ny - 1) {
    real ez_next_j = s_ez[threadIdx.y + 1][threadIdx.x];
    hx[PIDX(i, j, hx_pitch)] -= c_chy * (ez_next_j - ez_curr);
  }

  // update hy: hy[i,j] += chx * (ez[i+1,j] - ez[i,j])
  if (i < nx - 1 && j < ny) {
    real ez_next_i = s_ez[threadIdx.y][threadIdx.x + 1];
    hy[PIDX(i, j, hy_pitch)] += c_chx * (ez_next_i - ez_curr);
  }
}

__global__ void upd_ez_kern(real *__restrict__ ez, size_t ez_pitch,
                            const real *__restrict__ hx, size_t hx_pitch,
                            const real *__restrict__ hy, size_t hy_pitch,
                            const real *__restrict__ eps_inv, size_t eps_pitch,
                            int nx, int ny) {

  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  __shared__ real s_hx[BLOCK_Y + 1][BLOCK_X + 1];
  __shared__ real s_hy[BLOCK_Y + 1][BLOCK_X + 1];

  // Load main Hx tile
  if (i < nx && j > 0 && j < ny) {
    s_hx[threadIdx.y][threadIdx.x] = hx[PIDX(i, j - 1, hx_pitch)];
  } else {
    s_hx[threadIdx.y][threadIdx.x] = 0.0;
  }

  // Load halo: hx[i, j]
  if (threadIdx.y == blockDim.y - 1 && i < nx && j < ny - 1) {
    s_hx[threadIdx.y + 1][threadIdx.x] = hx[PIDX(i, j, hx_pitch)];
  } else if (threadIdx.y == blockDim.y - 1) {
    s_hx[threadIdx.y + 1][threadIdx.x] = 0.0;
  }

  // Load Hy into shared memory
  if (i > 0 && i < nx && j < ny) {
    s_hy[threadIdx.y][threadIdx.x] = hy[PIDX(i - 1, j, hy_pitch)];
  } else {
    s_hy[threadIdx.y][threadIdx.x] = 0.0;
  }

  // Load halo: hy
  if (threadIdx.x == blockDim.x - 1 && i < nx - 1 && j < ny) {
    s_hy[threadIdx.y][threadIdx.x + 1] = hy[PIDX(i, j, hy_pitch)];
  } else if (threadIdx.x == blockDim.x - 1) {
    s_hy[threadIdx.y][threadIdx.x + 1] = 0.0;
  }

  __syncthreads();

  // early exit
  if (i <= 0 || i >= nx - 1 || j <= 0 || j >= ny - 1)
    return;

  real hx_curr = s_hx[threadIdx.y + 1][threadIdx.x];
  real hx_prev = s_hx[threadIdx.y][threadIdx.x];
  real hy_curr = s_hy[threadIdx.y][threadIdx.x + 1];
  real hy_prev = s_hy[threadIdx.y][threadIdx.x];

  real dhx = hx_curr - hx_prev;
  real dhy = hy_curr - hy_prev;

  // Load local permittivity (inverse for efficiency)
  int idx_ez = PIDX(i, j, ez_pitch);
  int idx_eps = PIDX(i, j, eps_pitch);
  real local_eps_inv = eps_inv[idx_eps];

  // Update Ez with spatially varying permittivity
  // Ez += dt/eps * (d Hy/dx - d Hx/dy)
  ez[idx_ez] += local_eps_inv * (c_dt_over_dx * dhy - c_dt_over_dy * dhx);
}

int solve(struct SimulationParams *sim_params,
          struct PhysicalParams *phys_params, int problem_id,
          PermittivityFunc eps_func) {

  DEBUG_PRINT("Starting computation on process %d, with :\n\tnx from %d to "
              "%d\n\tny from %d to %d\n",
              0, 0, sim_params->nx - 1, 0, sim_params->ny - 1);

  struct data ez, hx, hy, eps_inv;
  if (init_data(&ez, "ez", sim_params->nx, sim_params->ny, sim_params->dx,
                sim_params->dy, 0.) ||
      init_data(&hx, "hx", sim_params->nx, sim_params->ny - 1, sim_params->dx,
                sim_params->dy, 0.) ||
      init_data(&hy, "hy", sim_params->nx - 1, sim_params->ny, sim_params->dx,
                sim_params->dy, 0.) ||
      init_data(&eps_inv, "eps_inv", sim_params->nx, sim_params->ny,
                sim_params->dx, sim_params->dy, 0.)) {
    printf("Error: could not allocate data\n");
    return EXIT_FAILURE;
  }

  // Initialize spatially varying permittivity
  real *eps_inv_host =
      (real *)malloc(sim_params->nx * sim_params->ny * sizeof(real));
  if (!eps_inv_host) {
    printf("Error: could not allocate host memory for permittivity\n");
    return EXIT_FAILURE;
  }

  printf("Initializing permittivity field...\n");
  // Compute permittivity at each grid point
  int count_modified = 0;
  real min_eps = 1e30, max_eps = -1e30;
  for (int j = 0; j < sim_params->ny; j++) {
    for (int i = 0; i < sim_params->nx; i++) {
      real x = i * sim_params->dx;
      real y = j * sim_params->dy;
      real local_eps = eps_func(x, y, phys_params->eps);
      // Store inverse for computational efficiency
      eps_inv_host[j * sim_params->nx + i] = 1.0 / local_eps;

      // Track statistics
      real rel_eps = local_eps / phys_params->eps;
      if (rel_eps < min_eps)
        min_eps = rel_eps;
      if (rel_eps > max_eps)
        max_eps = rel_eps;

      // Count non-vacuum points
      if (fabs(local_eps - phys_params->eps) > 1e-15) {
        count_modified++;
      }
    }
  }
  printf("Permittivity field initialized:\n");
  printf("  - %d/%d points have modified permittivity\n", count_modified,
         sim_params->nx * sim_params->ny);
  printf("  - Relative permittivity range: %.2f to %.2f\n", min_eps, max_eps);
  printf("  - This should show strong wave interaction if > 1.0!\n");

  // Copy to device
  CUDA_CHECK(cudaMemcpy2D(eps_inv.values, eps_inv.pitch, eps_inv_host,
                          sim_params->nx * sizeof(real),
                          sim_params->nx * sizeof(real), sim_params->ny,
                          cudaMemcpyHostToDevice));

  // Verify copy worked by reading back a few values
  real test_val;
  CUDA_CHECK(cudaMemcpy2D(&test_val, sizeof(real), eps_inv.values,
                          eps_inv.pitch, sizeof(real), 1,
                          cudaMemcpyDeviceToHost));
  printf("  - Device verification: eps_inv[0,0] = %.6e (should be %.6e)\n",
         test_val, eps_inv_host[0]);

  // Test a point in the middle where barrier should be
  int test_i = sim_params->nx / 2 + sim_params->nx / 5; // Right of center
  int test_j = sim_params->ny / 2;
  CUDA_CHECK(cudaMemcpy2D(
      &test_val, sizeof(real),
      (char *)eps_inv.values + test_j * eps_inv.pitch + test_i * sizeof(real),
      eps_inv.pitch, sizeof(real), 1, cudaMemcpyDeviceToHost));
  printf("  - Device verification: eps_inv[%d,%d] = %.6e (should be %.6e)\n",
         test_i, test_j, test_val,
         eps_inv_host[test_j * sim_params->nx + test_i]);

  // Optionally save permittivity field for visualization
  if (sim_params->sampling_rate) {
    struct data eps_viz;
    eps_viz.name = "eps_rel";
    eps_viz.nx = sim_params->nx;
    eps_viz.ny = sim_params->ny;
    eps_viz.dx = sim_params->dx;
    eps_viz.dy = sim_params->dy;
    // Convert back to relative permittivity for visualization
    real *eps_rel_host =
        (real *)malloc(sim_params->nx * sim_params->ny * sizeof(real));
    for (int i = 0; i < sim_params->nx * sim_params->ny; i++) {
      eps_rel_host[i] = (1.0 / eps_inv_host[i]) / phys_params->eps;
    }
    eps_viz.values = eps_rel_host;
    write_data_vtk(&eps_viz, 0, 0);
    free(eps_rel_host);
  }

  free(eps_inv_host);

  // Precomputed constants
  real chy = sim_params->dt / (sim_params->dy * phys_params->mu);
  real chx = sim_params->dt / (sim_params->dx * phys_params->mu);
  real dt_over_dx = sim_params->dt / sim_params->dx;
  real dt_over_dy = sim_params->dt / sim_params->dy;

  int src_i = sim_params->nx / 2;
  int src_j = sim_params->ny / 2;

  // Copy constants to constant memory
  cudaMemcpyToSymbol(c_chx, &chx, sizeof(real));
  cudaMemcpyToSymbol(c_chy, &chy, sizeof(real));
  cudaMemcpyToSymbol(c_dt_over_dx, &dt_over_dx, sizeof(real));
  cudaMemcpyToSymbol(c_dt_over_dy, &dt_over_dy, sizeof(real));

  dim3 blockDim(BLOCK_X, BLOCK_Y);

  dim3 gridDim_h((sim_params->nx + blockDim.x - 1) / blockDim.x,
                 (sim_params->ny + blockDim.y - 1) / blockDim.y);

  dim3 gridDim_ez((sim_params->nx + blockDim.x - 1) / blockDim.x,
                  (sim_params->ny + blockDim.y - 1) / blockDim.y);

  // Pre-allocate host memory for sampling
  real *ez_host_values = NULL;
  if (sim_params->sampling_rate) {
    ez_host_values = (real *)malloc(ez.nx * ez.ny * sizeof(real));
    if (!ez_host_values) {
      printf("Error: could not allocate host memory for sampling\n");
      cudaFree(ez.values);
      cudaFree(hx.values);
      cudaFree(hy.values);
      cudaFree(eps_inv.values);
      return EXIT_FAILURE;
    }
  }

  real start = GET_TIME();

  for (int n = 0; n < sim_params->nt; n++) {
    if (n && (n % (sim_params->nt / 10)) == 0) {
      real time_sofar = GET_TIME() - start;
      real eta = (sim_params->nt - n) * time_sofar / n;
      printf("Computing time step %d/%d (ETA: %g seconds)     \r", n,
             sim_params->nt, eta);
      fflush(stdout);
    }

    upd_h_kern<<<gridDim_h, blockDim>>>(hx.values, hx.pitch, hy.values,
                                        hy.pitch, ez.values, ez.pitch,
                                        sim_params->nx, sim_params->ny);

    upd_ez_kern<<<gridDim_ez, blockDim>>>(
        ez.values, ez.pitch, hx.values, hx.pitch, hy.values, hy.pitch,
        eps_inv.values, eps_inv.pitch, sim_params->nx, sim_params->ny);

    // Impose source
    switch (problem_id) {
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
    case 8:
    case 9:
    case 10:
    case 11:
    case 12: {
      real t = n * sim_params->dt;
      real source_val = sin(2.0 * M_PI * 2.4e9 * t);
      apply_source_kern<<<1, 1>>>(ez.values, ez.pitch, src_i, src_j,
                                  source_val);
      break;
    }
    default:
      printf("Error: unknown source\n");
      break;
    }

    // Output VTK data
    if (sim_params->sampling_rate && !(n % sim_params->sampling_rate)) {
      cudaMemcpy2D(ez_host_values, ez.nx * sizeof(real), ez.values, ez.pitch,
                   ez.nx * sizeof(real), ez.ny, cudaMemcpyDeviceToHost);

      struct data ez_host;
      ez_host.name = ez.name;
      ez_host.nx = ez.nx;
      ez_host.ny = ez.ny;
      ez_host.dx = ez.dx;
      ez_host.dy = ez.dy;
      ez_host.values = ez_host_values;

      write_data_vtk(&ez_host, n, 0);
    }
  }

  cudaDeviceSynchronize();

  write_manifest_vtk("ez", sim_params->dt, sim_params->nt,
                     sim_params->sampling_rate, 1);

  real time = GET_TIME() - start;
  printf("\nDone: %g seconds (%g MUpdates/s)\n", time,
         1.e-6 * (real)sim_params->nx * (real)sim_params->ny *
             (real)sim_params->nt / time);

  cudaFree(ez.values);
  cudaFree(hx.values);
  cudaFree(hy.values);
  cudaFree(eps_inv.values);
  if (ez_host_values)
    free(ez_host_values);

  return EXIT_SUCCESS;
}