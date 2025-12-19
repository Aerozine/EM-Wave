#include "data.h"
#include "main.h"
#include "solver.h"

#include <cuda_runtime.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>


// Constant memory for physics parameters
__constant__ double c_chx, c_chy, c_chz, c_cex, c_cey, c_cez;

#define GETK(data, nx, ny, i, j, k)                                            \
  ((data)[(nx) * (ny) * (k) + (nx) * (j) + (i)])
#define SETK(data, nx, ny, i, j, k, val)                                       \
  ((data)[(nx) * (ny) * (k) + (nx) * (j) + (i)] = (val))

__global__ void upd_ex_kern(double *ex, const double *hy, const double *hz,
                            int nx, int ny, int nz) {
  __shared__ double s_hz[34][6][6]; // 32+2 for halo, 4+2, 4+2
  __shared__ double s_hy[34][6][6];

  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;

  int tx = threadIdx.x;
  int ty = threadIdx.y;
  int tz = threadIdx.z;

  // Load into shared memory with halos
  if (i < nx - 1 && j < ny - 1 && k < nz - 1) {
    s_hz[tx][ty][tz] = GETK(hz, nx - 1, ny - 1, i, j, k);
    s_hy[tx][ty][tz] = GETK(hy, nx - 1, ny, i, j, k);

    // Load halo regions
    if (ty == 0 && j > 0) {
      s_hz[tx][ty + 5][tz] = GETK(hz, nx - 1, ny - 1, i, j - 1, k);
    }
    if (tz == 0 && k > 0) {
      s_hy[tx][ty][tz + 5] = GETK(hy, nx - 1, ny, i, j, k - 1);
    }
  }
  __syncthreads();

  if (!(i < nx - 1 && j > 0 && j < ny - 1 && k > 0 && k < nz - 1))
    return;

  double ex_ijk;
  if (ty > 0 && tz > 0) {
    // Use shared memory
    ex_ijk = GETK(ex, nx, ny, i, j, k) +
             c_cey * (s_hz[tx][ty][tz] - s_hz[tx][ty - 1][tz]) -
             c_cez * (s_hy[tx][ty][tz] - s_hy[tx][ty][tz - 1]);
  } else {
    // Fall back to global memory for boundary cases
    ex_ijk = GETK(ex, nx, ny, i, j, k) +
             c_cey * (GETK(hz, nx - 1, ny - 1, i, j, k) -
                      GETK(hz, nx - 1, ny - 1, i, j - 1, k)) -
             c_cez * (GETK(hy, nx - 1, ny, i, j, k) -
                      GETK(hy, nx - 1, ny, i, j, k - 1));
  }
  SETK(ex, nx, ny, i, j, k, ex_ijk);
}

__global__ void upd_ey_kern(double *ey, const double *hx, const double *hz,
                            int nx, int ny, int nz) {
  __shared__ double s_hx[34][6][6];
  __shared__ double s_hz[34][6][6];

  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;

  int tx = threadIdx.x;
  int ty = threadIdx.y;
  int tz = threadIdx.z;

  if (i < nx && j < ny - 1 && k < nz - 1) {
    if (i < nx - 1)
      s_hz[tx][ty][tz] = GETK(hz, nx - 1, ny - 1, i, j, k);
    s_hx[tx][ty][tz] = GETK(hx, nx, ny - 1, i, j, k);

    if (tx == 0 && i > 0) {
      s_hz[tx + 33][ty][tz] = GETK(hz, nx - 1, ny - 1, i - 1, j, k);
    }
    if (tz == 0 && k > 0) {
      s_hx[tx][ty][tz + 5] = GETK(hx, nx, ny - 1, i, j, k - 1);
    }
  }
  __syncthreads();

  if (!(i > 0 && i < nx - 1 && j < ny - 1 && k > 0 && k < nz - 1))
    return;

  double ey_ijk;
  if (tx > 0 && tz > 0 && i < nx - 1) {
    ey_ijk = GETK(ey, nx, ny, i, j, k) +
             c_cez * (s_hx[tx][ty][tz] - s_hx[tx][ty][tz - 1]) -
             c_cex * (s_hz[tx][ty][tz] - s_hz[tx - 1][ty][tz]);
  } else {
    ey_ijk = GETK(ey, nx, ny, i, j, k) +
             c_cez * (GETK(hx, nx, ny - 1, i, j, k) -
                      GETK(hx, nx, ny - 1, i, j, k - 1)) -
             c_cex * (GETK(hz, nx - 1, ny - 1, i, j, k) -
                      GETK(hz, nx - 1, ny - 1, i - 1, j, k));
  }
  SETK(ey, nx, ny, i, j, k, ey_ijk);
}

__global__ void upd_ez_kern(double *ez, const double *hx, const double *hy,
                            int nx, int ny, int nz) {
  __shared__ double s_hx[34][6][6];
  __shared__ double s_hy[34][6][6];

  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;

  int tx = threadIdx.x;
  int ty = threadIdx.y;
  int tz = threadIdx.z;

  if (i < nx && j < ny && k < nz - 1) {
    if (i < nx - 1)
      s_hy[tx][ty][tz] = GETK(hy, nx - 1, ny, i, j, k);
    if (j < ny - 1)
      s_hx[tx][ty][tz] = GETK(hx, nx, ny - 1, i, j, k);

    if (tx == 0 && i > 0 && i < nx) {
      s_hy[tx + 33][ty][tz] = GETK(hy, nx - 1, ny, i - 1, j, k);
    }
    if (ty == 0 && j > 0 && j < ny) {
      s_hx[tx][ty + 5][tz] = GETK(hx, nx, ny - 1, i, j - 1, k);
    }
  }
  __syncthreads();

  if (!(i > 0 && i < nx - 1 && j > 0 && j < ny - 1 && k < nz - 1))
    return;

  double ez_ijk;
  if (tx > 0 && ty > 0 && i < nx - 1 && j < ny - 1) {
    ez_ijk = GETK(ez, nx, ny, i, j, k) +
             c_cex * (s_hy[tx][ty][tz] - s_hy[tx - 1][ty][tz]) -
             c_cey * (s_hx[tx][ty][tz] - s_hx[tx][ty - 1][tz]);
  } else {
    ez_ijk = GETK(ez, nx, ny, i, j, k) +
             c_cex * (GETK(hy, nx - 1, ny, i, j, k) -
                      GETK(hy, nx - 1, ny, i - 1, j, k)) -
             c_cey * (GETK(hx, nx, ny - 1, i, j, k) -
                      GETK(hx, nx, ny - 1, i, j - 1, k));
  }
  SETK(ez, nx, ny, i, j, k, ez_ijk);
}

__global__ void upd_hx_kern(double *hx, const double *ey, const double *ez,
                            int nx, int ny, int nz) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  if (!(i < nx && j < ny - 1 && k < nz - 1))
    return;
  double hx_ijk =
      GETK(hx, nx, ny - 1, i, j, k) +
      c_chz * (GETK(ey, nx, ny, i, j, k + 1) - GETK(ey, nx, ny, i, j, k)) -
      c_chy * (GETK(ez, nx, ny, i, j + 1, k) - GETK(ez, nx, ny, i, j, k));
  SETK(hx, nx, ny - 1, i, j, k, hx_ijk);
}

__global__ void upd_hy_kern(double *hy, const double *ex, const double *ez,
                            int nx, int ny, int nz) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;
  if (!(i < nx - 1 && j < ny && k < nz - 1))
    return;
  double hy_ijk =
      GETK(hy, nx - 1, ny, i, j, k) +
      c_chx * (GETK(ez, nx, ny, i + 1, j, k) - GETK(ez, nx, ny, i, j, k)) -
      c_chz * (GETK(ex, nx, ny, i, j, k + 1) - GETK(ex, nx, ny, i, j, k));
  SETK(hy, nx - 1, ny, i, j, k, hy_ijk);
}

__global__ void upd_hz_kern(double *hz, const double *ex, const double *ey,
                            int nx, int ny, int nz) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  int k = blockIdx.z * blockDim.z + threadIdx.z;

  if (!(i < nx - 1 && j < ny - 1 && k < nz))
    return;
  double hz_ijk =
      GETK(hz, nx - 1, ny - 1, i, j, k) +
      c_chy * (GETK(ex, nx, ny, i, j + 1, k) - GETK(ex, nx, ny, i, j, k)) -
      c_chx * (GETK(ey, nx, ny, i + 1, j, k) - GETK(ey, nx, ny, i, j, k));
  SETK(hz, nx - 1, ny - 1, i, j, k, hz_ijk);
}

// Kernel to apply source on GPU
__global__ void apply_source_kern(double *ez, int src_idx, double source_val) {
  if (threadIdx.x == 0 && blockIdx.x == 0) {
    ez[src_idx] = source_val;
  }
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

  // Precompute constants and copy to constant memory
  double chx = sim_params->dt / (sim_params->dx * phys_params->mu);
  double chy = sim_params->dt / (sim_params->dy * phys_params->mu);
  double chz = sim_params->dt / (sim_params->dz * phys_params->mu);
  double cex = sim_params->dt / (sim_params->dx * phys_params->eps);
  double cey = sim_params->dt / (sim_params->dy * phys_params->eps);
  double cez = sim_params->dt / (sim_params->dz * phys_params->eps);

  CUDA_CHECK(cudaMemcpyToSymbol(c_chx, &chx, sizeof(double)));
  CUDA_CHECK(cudaMemcpyToSymbol(c_chy, &chy, sizeof(double)));
  CUDA_CHECK(cudaMemcpyToSymbol(c_chz, &chz, sizeof(double)));
  CUDA_CHECK(cudaMemcpyToSymbol(c_cex, &cex, sizeof(double)));
  CUDA_CHECK(cudaMemcpyToSymbol(c_cey, &cey, sizeof(double)));
  CUDA_CHECK(cudaMemcpyToSymbol(c_cez, &cez, sizeof(double)));

  // Optimized block dimensions: 32x4x4 = 512 threads
  dim3 blockDim(32, 4, 4);

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

  // Create streams for parallel execution and I/O
  cudaStream_t streamHx, streamHy, streamHz;
  cudaStream_t streamEx, streamEy, streamEz;
  cudaStream_t copyStream;

  CUDA_CHECK(cudaStreamCreate(&streamHx));
  CUDA_CHECK(cudaStreamCreate(&streamHy));
  CUDA_CHECK(cudaStreamCreate(&streamHz));
  CUDA_CHECK(cudaStreamCreate(&streamEx));
  CUDA_CHECK(cudaStreamCreate(&streamEy));
  CUDA_CHECK(cudaStreamCreate(&streamEz));
  CUDA_CHECK(cudaStreamCreate(&copyStream));

  // Allocate pinned host memory for async transfers
  double *h_ex_pinned = NULL, *h_ey_pinned = NULL, *h_ez_pinned = NULL;
  if (sim_params->sampling_rate) {
    CUDA_CHECK(
        cudaMallocHost(&h_ex_pinned, ex.nx * ex.ny * ex.nz * sizeof(double)));
    CUDA_CHECK(
        cudaMallocHost(&h_ey_pinned, ey.nx * ey.ny * ey.nz * sizeof(double)));
    CUDA_CHECK(
        cudaMallocHost(&h_ez_pinned, ez.nx * ez.ny * ez.nz * sizeof(double)));
  }

  // Precompute source index
  int src_idx = (sim_params->nz / 2) * sim_params->nx * sim_params->ny +
                (sim_params->ny / 2) * sim_params->nx + (sim_params->nx / 2);

  double start = GET_TIME();

  for (int n = 0; n < sim_params->nt; n++) {
    if (n && (n % (sim_params->nt / 10)) == 0) {
      double time_sofar = GET_TIME() - start;
      double eta = (sim_params->nt - n) * time_sofar / n;
      printf("Computing time step %d/%d (ETA: %g seconds)     \r", n,
             sim_params->nt, eta);
      fflush(stdout);
    }

    // Update H fields in parallel on separate streams
    upd_hx_kern<<<gridDim_hx, blockDim, 0, streamHx>>>(
        hx.values, ey.values, ez.values, sim_params->nx, sim_params->ny,
        sim_params->nz);

    upd_hy_kern<<<gridDim_hy, blockDim, 0, streamHy>>>(
        hy.values, ex.values, ez.values, sim_params->nx, sim_params->ny,
        sim_params->nz);

    upd_hz_kern<<<gridDim_hz, blockDim, 0, streamHz>>>(
        hz.values, ex.values, ey.values, sim_params->nx, sim_params->ny,
        sim_params->nz);

    // Synchronize H-field streams before updating E-fields
    CUDA_CHECK(cudaStreamSynchronize(streamHx));
    CUDA_CHECK(cudaStreamSynchronize(streamHy));
    CUDA_CHECK(cudaStreamSynchronize(streamHz));

    // Update E fields in parallel on separate streams
    upd_ex_kern<<<gridDim_e, blockDim, 0, streamEx>>>(
        ex.values, hy.values, hz.values, sim_params->nx, sim_params->ny,
        sim_params->nz);

    upd_ey_kern<<<gridDim_e, blockDim, 0, streamEy>>>(
        ey.values, hx.values, hz.values, sim_params->nx, sim_params->ny,
        sim_params->nz);

    upd_ez_kern<<<gridDim_e, blockDim, 0, streamEz>>>(
        ez.values, hx.values, hy.values, sim_params->nx, sim_params->ny,
        sim_params->nz);

    // Synchronize E-field streams before applying source
    CUDA_CHECK(cudaStreamSynchronize(streamEx));
    CUDA_CHECK(cudaStreamSynchronize(streamEy));
    CUDA_CHECK(cudaStreamSynchronize(streamEz));

    // Impose source on GPU
    double t = n * sim_params->dt;
    switch (problem_id) {
    case 1:
    case 2: {
      double source_val = sin(2. * M_PI * 2.4e9 * t);
      apply_source_kern<<<1, 1, 0, streamEz>>>(ez.values, src_idx, source_val);
      break;
    }
    default:
      printf("Error: unknown source\n");
      break;
    }

    // Async output with separate stream
    if (sim_params->sampling_rate && !(n % sim_params->sampling_rate)) {
      // Wait for previous copy to finish before starting new one
      CUDA_CHECK(cudaStreamSynchronize(copyStream));

      // Ensure E-field updates are complete before copying
      CUDA_CHECK(cudaStreamSynchronize(streamEx));
      CUDA_CHECK(cudaStreamSynchronize(streamEy));
      CUDA_CHECK(cudaStreamSynchronize(streamEz));

      // Async copy to pinned memory
      CUDA_CHECK(cudaMemcpyAsync(h_ex_pinned, ex.values,
                                 ex.nx * ex.ny * ex.nz * sizeof(double),
                                 cudaMemcpyDeviceToHost, copyStream));
      CUDA_CHECK(cudaMemcpyAsync(h_ey_pinned, ey.values,
                                 ey.nx * ey.ny * ey.nz * sizeof(double),
                                 cudaMemcpyDeviceToHost, copyStream));
      CUDA_CHECK(cudaMemcpyAsync(h_ez_pinned, ez.values,
                                 ez.nx * ez.ny * ez.nz * sizeof(double),
                                 cudaMemcpyDeviceToHost, copyStream));

      // Synchronize copy stream before file I/O
      CUDA_CHECK(cudaStreamSynchronize(copyStream));

      // Write to files (could be further optimized with threading)
      struct data ex_host = ex;
      ex_host.values = h_ex_pinned;
      write_data_vtk(&ex_host, n, 0);

      struct data ey_host = ey;
      ey_host.values = h_ey_pinned;
      write_data_vtk(&ey_host, n, 0);

      struct data ez_host = ez;
      ez_host.values = h_ez_pinned;
      write_data_vtk(&ez_host, n, 0);
    }
  }

  // Synchronize all streams before cleanup
  CUDA_CHECK(cudaStreamSynchronize(streamHx));
  CUDA_CHECK(cudaStreamSynchronize(streamHy));
  CUDA_CHECK(cudaStreamSynchronize(streamHz));
  CUDA_CHECK(cudaStreamSynchronize(streamEx));
  CUDA_CHECK(cudaStreamSynchronize(streamEy));
  CUDA_CHECK(cudaStreamSynchronize(streamEz));
  CUDA_CHECK(cudaStreamSynchronize(copyStream));

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

  // Cleanup
  if (sim_params->sampling_rate) {
    cudaFreeHost(h_ex_pinned);
    cudaFreeHost(h_ey_pinned);
    cudaFreeHost(h_ez_pinned);
  }

  cudaStreamDestroy(streamHx);
  cudaStreamDestroy(streamHy);
  cudaStreamDestroy(streamHz);
  cudaStreamDestroy(streamEx);
  cudaStreamDestroy(streamEy);
  cudaStreamDestroy(streamEz);
  cudaStreamDestroy(copyStream);

  free_data(&ex);
  free_data(&ey);
  free_data(&ez);
  free_data(&hx);
  free_data(&hy);
  free_data(&hz);

  return EXIT_SUCCESS;
}