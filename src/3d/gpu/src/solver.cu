#include "data.h"
#include "main.h"
#include "solver.h"

#include <cuda_runtime.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>

__global__ void apply_source_kern(double *ez, int src_idx, double value) {
  if (threadIdx.x == 0 && blockIdx.x == 0) {
    ez[src_idx] = value;
  }
}

__global__ void upd_hx_kern(double *__restrict__ hx, const double *__restrict__ ez, int nx, int ny,
                            double chy) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  if (i < nx && j < ny - 1) {
    int idx = IDX2D(i, j, nx);
    hx[idx] -= chy * (ez[IDX2D(i, j + 1, nx)] - ez[idx]);
  }
}

__global__ void upd_hy_kern(double *__restrict__ hy, const double *__restrict__ ez, int nx, int ny,
                            double chx) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  if (i < nx - 1 && j < ny) {
    int idx = IDX2D(i, j, nx - 1);
    hy[idx] += chx * (ez[IDX2D(i + 1, j, nx)] - ez[IDX2D(i, j, nx)]);
  }
}
__global__ void upd_ez_kern(double *__restrict__ ez, const double *__restrict__ hx, const double *__restrict__ hy,
                            int nx, int ny, double cex, double cey) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  if (i > 0 && i < nx - 1 && j > 0 && j < ny - 1) {
    int idx = IDX2D(i, j, nx);
    double dhy = hy[IDX2D(i, j, nx - 1)] - hy[IDX2D(i - 1, j, nx - 1)];
    double dhx = hx[idx] - hx[IDX2D(i, j - 1, nx)];
    ez[idx] += cex * dhy - cey * dhx;
  }
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

    // precomputed constant over iterations
    double chy = sim_params->dt / (sim_params->dy * phys_params->mu),
           chx = sim_params->dt / (sim_params->dx * phys_params->mu);
    double cex = sim_params->dt / (sim_params->dx * phys_params->eps),
           cey = sim_params->dt / (sim_params->dy * phys_params->eps);
    size_t src_idx = IDX2D(sim_params->nx/2,sim_params->ny / 2,sim_params->nx);
    dim3 blockDim(32, 16);  // 256 threads per block
    dim3 gridDim_hx((sim_params->nx + blockDim.x - 1) / blockDim.x,
                    (sim_params->ny - 1 + blockDim.y - 1) / blockDim.y);
    dim3 gridDim_hy((sim_params->nx - 1 + blockDim.x - 1) / blockDim.x,
                 (sim_params->ny + blockDim.y - 1) / blockDim.y);
    dim3 gridDim_ez((sim_params->nx + blockDim.x - 1) / blockDim.x,
                    (sim_params->ny + blockDim.y - 1) / blockDim.y);
  double start = GET_TIME();
  for (int n = 0; n < sim_params->nt; n++) {
    if (n && (n % (sim_params->nt / 10)) == 0) {
      double time_sofar = GET_TIME() - start;
      double eta = (sim_params->nt - n) * time_sofar / n;
      printf("Computing time step %d/%d (ETA: %g seconds)     \r", n,
             sim_params->nt, eta);
      fflush(stdout);
    }
    upd_hx_kern<<<gridDim_hx, blockDim>>>(hx.values, ez.values, sim_params->nx,
                                          sim_params->ny, chy);
    // Update hy on GPU
    upd_hy_kern<<<gridDim_hy, blockDim>>>(hy.values, ez.values, sim_params->nx,
                                          sim_params->ny, chx);
    // Update ez on GPU
    upd_ez_kern<<<gridDim_ez, blockDim>>>(ez.values, hx.values, hy.values,
                                          sim_params->nx, sim_params->ny, cex,
                                          cey);
    //CUDA guarantees in-order execution within the same stream

    // impose source
    switch (problem_id) {
    case 1:
    case 2:
{

  double t = n * sim_params->dt;
  double source_val = sin(2.0 * M_PI * 2.4e9 * t);
  //faster version using stream
  apply_source_kern<<<1, 1>>>(ez.values, src_idx, source_val);
  break;
}
      // sinusoidal excitation at 2.4 GHz in the middle of the domain
      // source_val = sin(2. * M_PI * 2.4e9 * t);
      // cudaMemcpy(&ez.values[src_idx], &source_val, sizeof(double),
      //           cudaMemcpyHostToDevice);
      break;
    default:
      printf("Error: unknown source\n");
      break;
    }

    // output step data in VTK format
    if (sim_params->sampling_rate && !(n % sim_params->sampling_rate)) {
      struct data ez_host;
      ez_host.name = ez.name;
      ez_host.nx = ez.nx;
      ez_host.ny = ez.ny;
      ez_host.dx = ez.dx;
      ez_host.dy = ez.dy;
      ez_host.values = (double *)malloc(ez.nx * ez.ny * sizeof(double));
      cudaMemcpy(ez_host.values, ez.values, ez.nx * ez.ny * sizeof(double),
                 cudaMemcpyDeviceToHost);
      write_data_vtk(&ez_host, n, 0);
      free(ez_host.values);
    }
  }

  // write VTK manifest, linking to individual step data files
  write_manifest_vtk("ez", sim_params->dt, sim_params->nt,
                     sim_params->sampling_rate, 1);
  //write_manifest_vtk("hx", sim_params->dt, sim_params->nt, sampling_rate, 1);
  //write_manifest_vtk("hy", sim_params->dt, sim_params->nt, sampling_rate, 1);

  double time = GET_TIME() - start;
  printf("\nDone: %g seconds (%g MUpdates/s)\n", time,
         1.e-6 * (double)sim_params->nx * (double)sim_params->ny *
             (double)sim_params->nt / time);


  cudaFree(ez.values);
  cudaFree(hx.values);
  cudaFree(hy.values);

  return EXIT_SUCCESS;
}