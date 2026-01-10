#include "data.h"
#include "main.h"
#include "solver.h"

#include <cuda_runtime.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
// to be refactored
// Block dimensions
#define BLOCK_X 32
#define BLOCK_Y 16
__constant__ real c_chx, c_chy, c_cex, c_cey;

__global__ void  __launch_bounds__(BLOCK_X * BLOCK_Y) apply_source_kern(real *__restrict__ ez, size_t pitch,
                                  int src_i, int src_j, real value) {
  if (threadIdx.x != 0 || blockIdx.x != 0)
    return;
  ez[PIDX(src_i, src_j, pitch)] = value;
}
__global__ void  __launch_bounds__(BLOCK_X * BLOCK_Y) upd_h_kern(real *__restrict__ hx, size_t hx_pitch,
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
  // reduce shared memory reads
  // increase reg per threads
  //real ez_curr = s_ez[threadIdx.y][threadIdx.x];

  // update hx
  // hx[i,j] -= chy * (ez[i,j+1] - ez[i,j])
  if (i < nx && j < ny - 1) {
    //real ez_next_j = s_ez[threadIdx.y + 1][threadIdx.x];
    hx[PIDX(i, j, hx_pitch)] -= c_chy * (s_ez[threadIdx.y + 1][threadIdx.x] - s_ez[threadIdx.y][threadIdx.x]);
  }

  // update hy
  // hy[i,j] += chx * (ez[i+1,j] - ez[i,j])
  if (i < nx - 1 && j < ny) {
    //real ez_next_i = s_ez[threadIdx.y][threadIdx.x + 1];
    hy[PIDX(i, j, hy_pitch)] += c_chx * (s_ez[threadIdx.y][threadIdx.x + 1] - s_ez[threadIdx.y][threadIdx.x]);
  }
}

__global__ void  __launch_bounds__(BLOCK_X * BLOCK_Y) upd_ez_kern(real *__restrict__ ez, size_t ez_pitch,
                            const real *__restrict__ hx, size_t hx_pitch,
                            const real *__restrict__ hy, size_t hy_pitch,
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
  // Load main Hy tile (offset by -1 in i to get i-1 values)
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

  // Update Ez with coalesced write
  int idx_ez = PIDX(i, j, ez_pitch);
  ez[idx_ez] += c_cex * dhy - c_cey * dhx;
}

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

  // Precomputed constants
  real chy = sim_params->dt / (sim_params->dy * phys_params->mu);
  real chx = sim_params->dt / (sim_params->dx * phys_params->mu);
  real cex = sim_params->dt / (sim_params->dx * phys_params->eps);
  real cey = sim_params->dt / (sim_params->dy * phys_params->eps);
  int src_i = sim_params->nx / 2;
  int src_j = sim_params->ny / 2;
  // Copy constants to constant memory (64KB constant cache, very fast)
  cudaMemcpyToSymbol(c_chx, &chx, sizeof(real));
  cudaMemcpyToSymbol(c_chy, &chy, sizeof(real));
  cudaMemcpyToSymbol(c_cex, &cex, sizeof(real));
  cudaMemcpyToSymbol(c_cey, &cey, sizeof(real));

  dim3 blockDim(BLOCK_X, BLOCK_Y);

  // Grid for full domain h
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

    upd_ez_kern<<<gridDim_ez, blockDim>>>(ez.values, ez.pitch, hx.values,
                                          hx.pitch, hy.values, hy.pitch,
                                          sim_params->nx, sim_params->ny);

    // Impose source
    switch (problem_id) {
    case 1:
    case 2: {
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
  if (ez_host_values)
    free(ez_host_values);

  return EXIT_SUCCESS;
}