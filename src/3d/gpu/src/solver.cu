#include "data.h"
#include "main.h"
#include "solver.h"

#include <cuda_runtime.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>

#define CUDA_CHECK(call)                                                       \
  do {                                                                         \
    cudaError_t err = call;                                                    \
    if (err != cudaSuccess) {                                                  \
      fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(err),     \
              __FILE__, __LINE__);                                             \
      exit(EXIT_FAILURE);                                                      \
    }                                                                          \
  } while (0)

#define GETK(data, nx, ny, i, j, k)                                            \
  ((data)[(nx) * (ny) * (k) + (nx) * (j) + (i)])
#define SETK(data, nx, ny, i, j, k, val)                                       \
  ((data)[(nx) * (ny) * (k) + (nx) * (j) + (i)] = (val))

// Update Ex: dEx/dt = (1/eps) * (dHz/dy - dHy/dz)
__global__ void upd_ex_kern(double *ex, const double *hy, const double *hz,
                            int nx, int ny, int nz, double cey, double cez) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  if (!(i < nx - 1 && j > 0 && j < ny - 1 && k > 0 && k < nz - 1))
    return;

  double ex_ijk =
      GETK(ex, nx, ny, i, j, k) +
      cey * (GETK(hz, nx - 1, ny - 1, i, j, k) -
             GETK(hz, nx - 1, ny - 1, i, j - 1, k)) -
      cez * (GETK(hy, nx - 1, ny, i, j, k) - GETK(hy, nx - 1, ny, i, j, k - 1));
  SETK(ex, nx, ny, i, j, k, ex_ijk);
}

// Update Ey: dEy/dt = (1/eps) * (dHx/dz - dHz/dx)
__global__ void upd_ey_kern(double *ey, const double *hx, const double *hz,
                            int nx, int ny, int nz, double cex, double cez) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;

  if (!(i > 0 && i < nx - 1 && j < ny - 1 && k > 0 && k < nz - 1))
    return;
  double ey_ijk = GETK(ey, nx, ny, i, j, k) +
                  cez * (GETK(hx, nx, ny - 1, i, j, k) -
                         GETK(hx, nx, ny - 1, i, j, k - 1)) -
                  cex * (GETK(hz, nx - 1, ny - 1, i, j, k) -
                         GETK(hz, nx - 1, ny - 1, i - 1, j, k));
  SETK(ey, nx, ny, i, j, k, ey_ijk);
}

// Update Ez: dEz/dt = (1/eps) * (dHy/dx - dHx/dy)
__global__ void upd_ez_kern(double *ez, const double *hx, const double *hy,
                            int nx, int ny, int nz, double cex, double cey) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;

  // Ez is at (i, j, k+0.5) and has dimensions [nx][ny][nz]
  if (!(i > 0 && i < nx - 1 && j > 0 && j < ny - 1 && k < nz - 1))
    return;
  double ez_ijk =
      GETK(ez, nx, ny, i, j, k) +
      cex *
          (GETK(hy, nx - 1, ny, i, j, k) - GETK(hy, nx - 1, ny, i - 1, j, k)) -
      cey * (GETK(hx, nx, ny - 1, i, j, k) - GETK(hx, nx, ny - 1, i, j - 1, k));
  SETK(ez, nx, ny, i, j, k, ez_ijk);
}

// Update Hx: dHx/dt = (1/mu) * (dEy/dz - dEz/dy)
__global__ void upd_hx_kern(double *hx, const double *ey, const double *ez,
                            int nx, int ny, int nz, double chy, double chz) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  if (!(i < nx && j < ny - 1 && k < nz - 1))
    return;
  double hx_ijk =
      GETK(hx, nx, ny - 1, i, j, k) +
      chz * (GETK(ey, nx, ny, i, j, k + 1) - GETK(ey, nx, ny, i, j, k)) -
      chy * (GETK(ez, nx, ny, i, j + 1, k) - GETK(ez, nx, ny, i, j, k));
  SETK(hx, nx, ny - 1, i, j, k, hx_ijk);
}

// Update Hy: dHy/dt = (1/mu) * (dEz/dx - dEx/dz)
__global__ void upd_hy_kern(double *hy, const double *ex, const double *ez,
                            int nx, int ny, int nz, double chx, double chz) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  if (!(i < nx - 1 && j < ny && k < nz - 1))
    return;
  double hy_ijk =
      GETK(hy, nx - 1, ny, i, j, k) +
      chx * (GETK(ez, nx, ny, i + 1, j, k) - GETK(ez, nx, ny, i, j, k)) -
      chz * (GETK(ex, nx, ny, i, j, k + 1) - GETK(ex, nx, ny, i, j, k));
  SETK(hy, nx - 1, ny, i, j, k, hy_ijk);
}

// Update Hz: dHz/dt = (1/mu) * (dEx/dy - dEy/dx)
__global__ void upd_hz_kern(double *hz, const double *ex, const double *ey,
                            int nx, int ny, int nz, double chx, double chy) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;

  if (!(i < nx - 1 && j < ny - 1 && k < nz))
    return;
  double hz_ijk =
      GETK(hz, nx - 1, ny - 1, i, j, k) +
      chy * (GETK(ex, nx, ny, i, j + 1, k) - GETK(ex, nx, ny, i, j, k)) -
      chx * (GETK(ey, nx, ny, i + 1, j, k) - GETK(ey, nx, ny, i, j, k));
  SETK(hz, nx - 1, ny - 1, i, j, k, hz_ijk);
}

int solve(struct SimulationParams *sim_params,
          struct PhysicalParams *phys_params, int problem_id) {

  DEBUG_PRINT("Starting computation on process %d, with:\n\tnx from %d to "
              "%d\n\tny from %d to %d\n\tnz from %d to %d\n",
              0, 0, sim_params->nx - 1, 0, sim_params->ny - 1, 0,
              sim_params->nz - 1);

  // Allocate field components
  struct data ex, ey, ez, hx, hy, hz;

  if (init_data(&ex, "ex", sim_params->nx, sim_params->ny, sim_params->nz,
                sim_params->dx, sim_params->dy, sim_params->dz, 0.) ||
      init_data(&ey, "ey", sim_params->nx, sim_params->ny, sim_params->nz,
                sim_params->dx, sim_params->dy, sim_params->dz, 0.) ||
      init_data(&ez, "ez", sim_params->nx, sim_params->ny, sim_params->nz,
                sim_params->dx, sim_params->dy, sim_params->dz, 0.) ||
      init_data(&hx, "hx", sim_params->nx, sim_params->ny - 1,
                sim_params->nz - 1, sim_params->dx, sim_params->dy,
                sim_params->dz, 0.) ||
      init_data(&hy, "hy", sim_params->nx - 1, sim_params->ny,
                sim_params->nz - 1, sim_params->dx, sim_params->dy,
                sim_params->dz, 0.) ||
      init_data(&hz, "hz", sim_params->nx - 1, sim_params->ny - 1,
                sim_params->nz, sim_params->dx, sim_params->dy, sim_params->dz,
                0.)) {
    printf("Error: could not allocate data\n");
    return EXIT_FAILURE;
  }

  // Precompute constants
  double chx = sim_params->dt / (sim_params->dx * phys_params->mu);
  double chy = sim_params->dt / (sim_params->dy * phys_params->mu);
  double chz = sim_params->dt / (sim_params->dz * phys_params->mu);
  double cex = sim_params->dt / (sim_params->dx * phys_params->eps);
  double cey = sim_params->dt / (sim_params->dy * phys_params->eps);
  double cez = sim_params->dt / (sim_params->dz * phys_params->eps);

  // Setup 3D grid dimensions for each field component
  dim3 blockDim(8, 8, 8);
  // Grid for E fields (full nx, ny, nz)
  dim3 gridDim_e((sim_params->nx + blockDim.x - 1) / blockDim.x,
                 (sim_params->ny + blockDim.y - 1) / blockDim.y,
                 (sim_params->nz + blockDim.z - 1) / blockDim.z);

  // Grid for Hx (nx, ny-1, nz-1)
  dim3 gridDim_hx((sim_params->nx + blockDim.x - 1) / blockDim.x,
                  (sim_params->ny - 1 + blockDim.y - 1) / blockDim.y,
                  (sim_params->nz - 1 + blockDim.z - 1) / blockDim.z);

  // Grid for Hy (nx-1, ny, nz-1)
  dim3 gridDim_hy((sim_params->nx - 1 + blockDim.x - 1) / blockDim.x,
                  (sim_params->ny + blockDim.y - 1) / blockDim.y,
                  (sim_params->nz - 1 + blockDim.z - 1) / blockDim.z);

  // Grid for Hz (nx-1, ny-1, nz)
  dim3 gridDim_hz((sim_params->nx - 1 + blockDim.x - 1) / blockDim.x,
                  (sim_params->ny - 1 + blockDim.y - 1) / blockDim.y,
                  (sim_params->nz + blockDim.z - 1) / blockDim.z);

  double start = GET_TIME();

  for (int n = 0; n < sim_params->nt; n++) {
    if (n && (n % (sim_params->nt / 10)) == 0) {
      double time_sofar = GET_TIME() - start;
      double eta = (sim_params->nt - n) * time_sofar / n;
      printf("Computing time step %d/%d (ETA: %g seconds)     \r", n,
             sim_params->nt, eta);
      fflush(stdout);
    }

    // Update H fields
    upd_hx_kern<<<gridDim_hx, blockDim>>>(hx.values, ey.values, ez.values,
                                          sim_params->nx, sim_params->ny,
                                          sim_params->nz, chy, chz);

    upd_hy_kern<<<gridDim_hy, blockDim>>>(hy.values, ex.values, ez.values,
                                          sim_params->nx, sim_params->ny,
                                          sim_params->nz, chx, chz);

    upd_hz_kern<<<gridDim_hz, blockDim>>>(hz.values, ex.values, ey.values,
                                          sim_params->nx, sim_params->ny,
                                          sim_params->nz, chx, chy);

    // Update E fields
    upd_ex_kern<<<gridDim_e, blockDim>>>(ex.values, hy.values, hz.values,
                                         sim_params->nx, sim_params->ny,
                                         sim_params->nz, cey, cez);
    CUDA_CHECK(cudaGetLastError());

    upd_ey_kern<<<gridDim_e, blockDim>>>(ey.values, hx.values, hz.values,
                                         sim_params->nx, sim_params->ny,
                                         sim_params->nz, cex, cez);
    CUDA_CHECK(cudaGetLastError());

    upd_ez_kern<<<gridDim_e, blockDim>>>(ez.values, hx.values, hy.values,
                                         sim_params->nx, sim_params->ny,
                                         sim_params->nz, cex, cey);

    // Impose source
    double t = n * sim_params->dt;
    switch (problem_id) {
    case 1:
    case 2: {
      // sinusoidal excitation at 2.4 GHz in the middle of the domain
      double source_val = sin(2. * M_PI * 2.4e9 * t);
      int src_idx = (sim_params->nz / 2) * sim_params->nx * sim_params->ny +
                    (sim_params->ny / 2) * sim_params->nx +
                    (sim_params->nx / 2);
      cudaMemcpy(&ez.values[src_idx], &source_val, sizeof(double),
                 cudaMemcpyHostToDevice);
      break;
    }
    default:
      printf("Error: unknown source\n");
      break;
    }

    // Output step data in VTK format
    if (sim_params->sampling_rate && !(n % sim_params->sampling_rate)) {
      // Write Ex
      struct data ex_host;
      ex_host.name = ex.name;
      ex_host.nx = ex.nx;
      ex_host.ny = ex.ny;
      ex_host.nz = ex.nz;
      ex_host.dx = ex.dx;
      ex_host.dy = ex.dy;
      ex_host.dz = ex.dz;
      ex_host.values = (double *)malloc(ex.nx * ex.ny * ex.nz * sizeof(double));
      cudaMemcpy(ex_host.values, ex.values,
                 ex.nx * ex.ny * ex.nz * sizeof(double),
                 cudaMemcpyDeviceToHost);
      write_data_vtk(&ex_host, n, 0);
      free(ex_host.values);

      // Write Ey
      struct data ey_host;
      ey_host.name = ey.name;
      ey_host.nx = ey.nx;
      ey_host.ny = ey.ny;
      ey_host.nz = ey.nz;
      ey_host.dx = ey.dx;
      ey_host.dy = ey.dy;
      ey_host.dz = ey.dz;
      ey_host.values = (double *)malloc(ey.nx * ey.ny * ey.nz * sizeof(double));
      cudaMemcpy(ey_host.values, ey.values,
                 ey.nx * ey.ny * ey.nz * sizeof(double),
                 cudaMemcpyDeviceToHost);
      write_data_vtk(&ey_host, n, 0);
      free(ey_host.values);

      // Write Ez
      struct data ez_host;
      ez_host.name = ez.name;
      ez_host.nx = ez.nx;
      ez_host.ny = ez.ny;
      ez_host.nz = ez.nz;
      ez_host.dx = ez.dx;
      ez_host.dy = ez.dy;
      ez_host.dz = ez.dz;
      ez_host.values = (double *)malloc(ez.nx * ez.ny * ez.nz * sizeof(double));
      cudaMemcpy(ez_host.values, ez.values,
                 ez.nx * ez.ny * ez.nz * sizeof(double),
                 cudaMemcpyDeviceToHost);
      write_data_vtk(&ez_host, n, 0);
      free(ez_host.values);
    }
  }

  // Write VTK manifests
  write_manifest_vtk("ex", sim_params->dt, sim_params->nt,
                     sim_params->sampling_rate, 1);
  write_manifest_vtk("ey", sim_params->dt, sim_params->nt,
                     sim_params->sampling_rate, 1);
  write_manifest_vtk("ez", sim_params->dt, sim_params->nt,
                     sim_params->sampling_rate, 1);

  double time = GET_TIME() - start;
#ifndef STABILITY_STUDY
  printf("\nDone: %g seconds (%g MUpdates/s)\n", time,
         1.e-6 * (double)sim_params->nx * (double)sim_params->ny *
             (double)sim_params->nz * (double)sim_params->nt / time);
#endif

  free_data(&ex);
  free_data(&ey);
  free_data(&ez);
  free_data(&hx);
  free_data(&hy);
  free_data(&hz);

  return EXIT_SUCCESS;
}