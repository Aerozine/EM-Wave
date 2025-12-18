#include "solver.h"
#include "data.h"
#include "main.h"
#include "params.h"

#include <math.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct area {
  int ndim;
  int *start;
  int *end;
  double *origin;
};

void free_area(struct area *object) {
  free(object->start);
  free(object->end);
  free(object->origin);
  free(object);
}

// Function that sets all parameters relating to the computation area of the
// current process. It returns the area, and sets the neighbours in mpi_params.
struct area *get_area(struct SimulationParams *sim_params,
                      struct MpiParams *mpi_params) {

  // Allocation of memory
  struct area *current_area = malloc(sizeof(struct area));
  if (!current_area)
    return NULL;

  current_area->ndim = sim_params->ndim;
  current_area->start = NULL;
  current_area->end = NULL;
  current_area->origin = NULL;
  current_area->start = malloc(sizeof(int) * current_area->ndim);
  current_area->end = malloc(sizeof(int) * current_area->ndim);
  current_area->origin = malloc(sizeof(double) * current_area->ndim);
  if (!current_area->start || !current_area->end) {
    free(current_area->start);
    free(current_area->end);
    free(current_area->origin);
    free(current_area);
    return NULL;
  }

  // Setting the parameters
  if (!mpi_params->use_mpi) {
    for (int i = 0; i < sim_params->ndim; i++) {
      current_area->start[i] = 0;
      current_area->end[i] = sim_params->size_of_space[i];
      current_area->origin[i] = 0.;
    }
  } else {
    // Second allocation of memory
    mpi_params->procs_per_dim = malloc(sizeof(int) * sim_params->ndim);
    mpi_params->periods = malloc(sizeof(int) * sim_params->ndim);
    mpi_params->coords = malloc(sizeof(int) * sim_params->ndim);
    mpi_params->neighbours =
        malloc(sizeof(struct neighbour) * 2 * sim_params->ndim);
    if (!mpi_params->procs_per_dim || !mpi_params->periods ||
        !mpi_params->coords || !mpi_params->neighbours) {
      free(mpi_params->procs_per_dim);
      free(mpi_params->periods);
      free(mpi_params->coords);
      free(mpi_params->neighbours);
      free_area(current_area);
      return NULL;
    }

    for (int i = 0; i < sim_params->ndim; i++) {
      mpi_params->procs_per_dim[i] = 0;
      mpi_params->periods[i] = 0;
      mpi_params->coords[i] = 0;
    }

    // Getting the repartition from MPI
    MPI_Dims_create(mpi_params->num_ranks, sim_params->ndim,
                    mpi_params->procs_per_dim);
    MPI_Cart_create(MPI_COMM_WORLD, sim_params->ndim, mpi_params->procs_per_dim,
                    mpi_params->periods, 0, &(mpi_params->cart_comm));
    MPI_Cart_coords(mpi_params->cart_comm, mpi_params->rank, sim_params->ndim,
                    mpi_params->coords);

    for (int i = 0; i < sim_params->ndim; i++) {
      // Getting start and end of domain
      int base_size =
          sim_params->size_of_space[i] / mpi_params->procs_per_dim[i];
      int remainder =
          sim_params->size_of_space[i] % mpi_params->procs_per_dim[i];

      int start = mpi_params->coords[i] * base_size +
                  (mpi_params->coords[i] < remainder ? mpi_params->coords[i]
                                                     : remainder);
      int size = base_size + (mpi_params->coords[i] < remainder ? 1 : 0);
      int end = start + size - 1;

      current_area->start[i] = start;
      current_area->end[i] = end;

      // Getting neighbours
      int negative, positive;
      MPI_Cart_shift(mpi_params->cart_comm, i, 1, &negative, &positive);

      // Note that if there's no neighbour, the rank is set to -1
      if (negative != MPI_PROC_NULL) {
        mpi_params->neighbours[2 * i].rank = negative;
        mpi_params->neighbours[2 * i].pos = 2 * i;
      } else {
        mpi_params->neighbours[2 * i].rank = -1;
      }
      if (positive != MPI_PROC_NULL) {
        mpi_params->neighbours[2 * i + 1].rank = positive;
        mpi_params->neighbours[2 * i + 1].pos = 2 * i + 1;
      } else {
        mpi_params->neighbours[2 * i + 1].rank = -1;
      }
    }
  }

  return current_area;
}

void free_all_solve_pointers(int nb_neighbours, struct area *proc_area,
                             struct data *hx, struct data *hy, struct data *hz,
                             struct data *ex, struct data *ey, struct data *ez,
                             MPI_Request *send_requests,
                             MPI_Request *recv_requests, float **sent_data,
                             float **received_data, bool *received_neighbour,
                             double *origins) {
  free_area(proc_area);
  free_data(hx);
  free_data(hy);
  free_data(hz);
  free_data(ex);
  free_data(ey);
  free_data(ez);
  free(send_requests);
  free(recv_requests);

  for (int i = 0; i < nb_neighbours; i++) {
    free(sent_data[i]);
    free(received_data[i]);
  }
  free(sent_data);
  free(received_data);
  free(received_neighbour);
  free(origins);
}

inline void h_loop(struct data *hx, struct data *hy, struct data *hz,
                   struct data *ex, struct data *ey, struct data *ez,
                   struct SimulationParams *sim_params,
                   struct PhysicalParams *phys_params,
                   struct MpiParams *mpi_params, float **received_data,
                   bool *received_neighbour, MPI_Request *requests) {
  // Main loop
  float chx1 = sim_params->steps[sim_params->ndim] /
               (sim_params->steps[2] * phys_params->mu);
  float chx2 = sim_params->steps[sim_params->ndim] /
               (sim_params->steps[1] * phys_params->mu);
  float chy1 = sim_params->steps[sim_params->ndim] /
               (sim_params->steps[0] * phys_params->mu);
  float chy2 = sim_params->steps[sim_params->ndim] /
               (phys_params->mu * sim_params->steps[2]);
  float chz1 = sim_params->steps[sim_params->ndim] /
               (phys_params->mu * sim_params->steps[1]);
  float chz2 = sim_params->steps[sim_params->ndim] /
               (phys_params->mu * sim_params->steps[0]);
  for (int k = 0; k < mpi_params->sizes[2] - 1; k++) {
    for (int j = 0; j < mpi_params->sizes[1] - 1; j++) {
      for (int i = 0; i < mpi_params->sizes[0] - 1; i++) {
        float hx_ij = GET(hx, i, j, k) +
                      chx1 * (GET(ey, i, j, k + 1) - GET(ey, i, j, k)) -
                      chx2 * (GET(ez, i, j + 1, k) - GET(ez, i, j, k));
        SET(hx, i, j, k, hx_ij);

        float hy_ij = GET(hy, i, j, k) +
                      chy1 * (GET(ez, i + 1, j, k) - GET(ez, i, j, k)) -
                      chy2 * (GET(ex, i, j, k + 1) - GET(ex, i, j, k));
        SET(hy, i, j, k, hy_ij);

        float hz_ij = GET(hz, i, j, k) +
                      chz1 * (GET(ex, i, j + 1, k) - GET(ex, i, j, k)) -
                      chz2 * (GET(ey, i + 1, j, k) - GET(ey, i, j, k));
        SET(hz, i, j, k, hz_ij);
      }
    }
  }

  // hx loop for i = mpi_params->sizes[0] - 1
  int i = mpi_params->sizes[0] - 1;
  for (int k = 0; k < mpi_params->sizes[2] - 1; k++) {
    for (int j = 0; j < mpi_params->sizes[1] - 1; j++) {
      float hx_ij = GET(hx, i, j, k) +
                    chx1 * (GET(ey, i, j, k + 1) - GET(ey, i, j, k)) -
                    chx2 * (GET(ez, i, j + 1, k) - GET(ez, i, j, k));
      SET(hx, i, j, k, hx_ij);
    }
  }

  // hy loop for j = mpi_params->sizes[1] - 1
  int j = mpi_params->sizes[1] - 1;
  for (int k = 0; k < mpi_params->sizes[2] - 1; k++) {
    for (int i = 0; i < mpi_params->sizes[0] - 1; i++) {
      float hy_ij = GET(hy, i, j, k) +
                    chy1 * (GET(ez, i + 1, j, k) - GET(ez, i, j, k)) -
                    chy2 * (GET(ex, i, j, k + 1) - GET(ex, i, j, k));
      SET(hy, i, j, k, hy_ij);
    }
  }

  // hz loop for k = mpi_params->sizes[2] - 1
  int k = mpi_params->sizes[2] - 1;
  for (int j = 0; j < mpi_params->sizes[1] - 1; j++) {
    for (int i = 0; i < mpi_params->sizes[0] - 1; i++) {
      float hz_ij = GET(hz, i, j, k) +
                    chz1 * (GET(ex, i, j + 1, k) - GET(ex, i, j, k)) -
                    chz2 * (GET(ey, i + 1, j, k) - GET(ey, i, j, k));
      SET(hz, i, j, k, hz_ij);
    }
  }

  // Neighbours updating
  MPI_Wait(&(requests[2 * X_END]), MPI_STATUS_IGNORE);
  MPI_Wait(&(requests[2 * X_END + 1]), MPI_STATUS_IGNORE);
  MPI_Wait(&(requests[2 * Y_END]), MPI_STATUS_IGNORE);
  MPI_Wait(&(requests[2 * Y_END + 1]), MPI_STATUS_IGNORE);
  MPI_Wait(&(requests[2 * Z_END]), MPI_STATUS_IGNORE);
  MPI_Wait(&(requests[2 * Z_END + 1]), MPI_STATUS_IGNORE);

  // X_END main condition
  if (received_neighbour[X_END]) {
    int i = mpi_params->sizes[0] - 1;
    for (int k = 0; k < mpi_params->sizes[2]; k++) {
      for (int j = 0; j < mpi_params->sizes[1]; j++) {
        if (k != mpi_params->sizes[2] - 1) {
          float hy_ij =
              GET(hy, i, j, k) +
              chy1 * (received_data[X_END]
                                   [mpi_params->send_array_sizes[X_END >> 1] +
                                    k * ez->ny + j] -
                      GET(ez, i, j, k)) -
              chy2 * (GET(ex, i, j, k + 1) - GET(ex, i, j, k));
          SET(hy, i, j, k, hy_ij);
        }

        if (j != mpi_params->sizes[1] - 1) {
          float hz_ij =
              GET(hz, i, j, k) +
              chz1 * (GET(ex, i, j + 1, k) - GET(ex, i, j, k)) -
              chz2 * (received_data[X_END][k * ey->ny + j] - GET(ey, i, j, k));
          SET(hz, i, j, k, hz_ij);
        }
      }
    }

    if (received_neighbour[Y_END]) {
      int j = mpi_params->sizes[1] - 1;
      for (int k = 0; k < mpi_params->sizes[2]; k++) {
        float hz_ij =
            GET(hz, i, j, k) +
            chz1 * (received_data[Y_END][k * ey->nx + i] - GET(ex, i, j, k)) -
            chz2 * (received_data[X_END][k * ey->ny + j] - GET(ey, i, j, k));

        SET(hz, i, j, k, hz_ij);
      }
    }

    if (received_neighbour[Z_END]) {
      int k = mpi_params->sizes[2] - 1;
      for (int j = 0; j < mpi_params->sizes[1]; j++) {
        float hy_ij =
            GET(hy, i, j, k) +
            chy1 *
                (received_data[X_END][mpi_params->send_array_sizes[X_END >> 1] +
                                      k * ez->ny + j] -
                 GET(ez, i, j, k)) -
            chy2 * (received_data[Z_END][j * ex->nx + i] - GET(ex, i, j, k));
        SET(hy, i, j, k, hy_ij);
      }
    }
  }

  // Y_END main condition
  if (received_neighbour[Y_END]) {
    int j = mpi_params->sizes[1] - 1;
    for (int k = 0; k < mpi_params->sizes[2]; k++) {
      for (int i = 0; i < mpi_params->sizes[0]; i++) {
        if (k != mpi_params->sizes[2] - 1) {
          float hx_ij =
              GET(hx, i, j, k) +
              chx1 * (GET(ey, i, j, k + 1) - GET(ey, i, j, k)) -
              chx2 * (received_data[Y_END]
                                   [mpi_params->send_array_sizes[Y_END >> 1] +
                                    k * ez->nx + i] -
                      GET(ez, i, j, k));
          SET(hx, i, j, k, hx_ij);
        }

        if (i != mpi_params->sizes[0] - 1) {
          float hz_ij =
              GET(hz, i, j, k) +
              chz1 * (received_data[Y_END][k * ex->nx + i] - GET(ex, i, j, k)) -
              chz2 * (GET(ey, i + 1, j, k) - GET(ey, i, j, k));
          SET(hz, i, j, k, hz_ij);
        }
      }
    }

    if (received_neighbour[Z_END]) {
      int k = mpi_params->sizes[2] - 1;
      for (int i = 0; i < mpi_params->sizes[0]; i++) {
        float hx_ij =
            GET(hx, i, j, k) +
            chx1 *
                (received_data[Z_END][mpi_params->send_array_sizes[Z_END >> 1] +
                                      j * ey->nx + i] -
                 GET(ey, i, j, k)) -
            chx2 *
                (received_data[Y_END][mpi_params->send_array_sizes[Y_END >> 1] +
                                      k * ez->nx + i] -
                 GET(ez, i, j, k));
        SET(hx, i, j, k, hx_ij);
      }
    }
  }

  // Main Z_END condition
  if (received_neighbour[Z_END]) {
    int k = mpi_params->sizes[2] - 1;
    for (int j = 0; j < mpi_params->sizes[1]; j++) {
      for (int i = 0; i < mpi_params->sizes[0]; i++) {
        if (j != mpi_params->sizes[1] - 1) {
          float hx_ij =
              GET(hx, i, j, k) +
              chx1 * (received_data[Z_END]
                                   [mpi_params->send_array_sizes[Z_END >> 1] +
                                    j * ey->nx + i] -
                      GET(ey, i, j, k)) -
              chx2 * (GET(ez, i, j + 1, k) - GET(ez, i, j, k));
          SET(hx, i, j, k, hx_ij);
        }

        if (i != mpi_params->sizes[0] - 1) {
          float hy_ij =
              GET(hy, i, j, k) +
              chy1 * (GET(ez, i + 1, j, k) - GET(ez, i, j, k)) -
              chy2 * (received_data[Z_END][j * ex->nx + i] - GET(ex, i, j, k));
          SET(hy, i, j, k, hy_ij);
        }
      }
    }
  }
}

inline void e_loop(struct data *hx, struct data *hy, struct data *hz,
                   struct data *ex, struct data *ey, struct data *ez,
                   struct SimulationParams *sim_params,
                   struct PhysicalParams *phys_params,
                   struct MpiParams *mpi_params, float **received_data,
                   bool *received_neighbour, MPI_Request *requests) {
  float cex1 = sim_params->steps[sim_params->ndim] /
               (phys_params->eps * sim_params->steps[1]);
  float cex2 = sim_params->steps[sim_params->ndim] /
               (phys_params->eps * sim_params->steps[2]);
  float cey1 = sim_params->steps[sim_params->ndim] /
               (phys_params->eps * sim_params->steps[2]);
  float cey2 = sim_params->steps[sim_params->ndim] /
               (phys_params->eps * sim_params->steps[0]);
  float cez1 = sim_params->steps[sim_params->ndim] /
               (sim_params->steps[0] * phys_params->eps),
        cez2 = sim_params->steps[sim_params->ndim] /
               (sim_params->steps[1] * phys_params->eps);

  for (int k = 1; k < mpi_params->sizes[2]; k++) {
    for (int j = 1; j < mpi_params->sizes[1]; j++) {
      for (int i = 1; i < mpi_params->sizes[0]; i++) {
        float ex_ij = GET(ex, i, j, k) +
                      cex1 * (GET(hz, i, j, k) - GET(hz, i, j - 1, k)) -
                      cex2 * (GET(hy, i, j, k) - GET(hy, i, j, k - 1));
        SET(ex, i, j, k, ex_ij);

        float ey_ij = GET(ey, i, j, k) +
                      cey1 * (GET(hx, i, j, k) - GET(hx, i, j, k - 1)) -
                      cey2 * (GET(hz, i, j, k) - GET(hz, i - 1, j, k));
        SET(ey, i, j, k, ey_ij);

        float ez_ij = GET(ez, i, j, k) +
                      cez1 * (GET(hy, i, j, k) - GET(hy, i - 1, j, k)) -
                      cez2 * (GET(hx, i, j, k) - GET(hx, i, j - 1, k));
        SET(ez, i, j, k, ez_ij);
      }
    }
  }

  // ex loop for i = 0
  int i = 0;
  for (int k = 1; k < mpi_params->sizes[2]; k++) {
    for (int j = 1; j < mpi_params->sizes[1]; j++) {
      float ex_ij = GET(ex, i, j, k) +
                    cex1 * (GET(hz, i, j, k) - GET(hz, i, j - 1, k)) -
                    cex2 * (GET(hy, i, j, k) - GET(hy, i, j, k - 1));
      SET(ex, i, j, k, ex_ij);
    }
  }

  // ey loop for j = 0
  int j = 0;
  for (int k = 1; k < mpi_params->sizes[2]; k++) {
    for (int i = 1; i < mpi_params->sizes[0]; i++) {
      float ey_ij = GET(ey, i, j, k) +
                    cey1 * (GET(hx, i, j, k) - GET(hx, i, j, k - 1)) -
                    cey2 * (GET(hz, i, j, k) - GET(hz, i - 1, j, k));
      SET(ey, i, j, k, ey_ij);
    }
  }

  // ez loop for k = 0
  int k = 0;
  for (int j = 1; j < mpi_params->sizes[1]; j++) {
    for (int i = 1; i < mpi_params->sizes[0]; i++) {
      float ez_ij = GET(ez, i, j, k) +
                    cez1 * (GET(hy, i, j, k) - GET(hy, i - 1, j, k)) -
                    cez2 * (GET(hx, i, j, k) - GET(hx, i, j - 1, k));
      SET(ez, i, j, k, ez_ij);
    }
  }

  MPI_Wait(&(requests[2 * X_START]), MPI_STATUS_IGNORE);
  MPI_Wait(&(requests[2 * X_START + 1]), MPI_STATUS_IGNORE);
  MPI_Wait(&(requests[2 * Y_START]), MPI_STATUS_IGNORE);
  MPI_Wait(&(requests[2 * Y_START + 1]), MPI_STATUS_IGNORE);
  MPI_Wait(&(requests[2 * Z_START]), MPI_STATUS_IGNORE);
  MPI_Wait(&(requests[2 * Z_START + 1]), MPI_STATUS_IGNORE);

  if (received_neighbour[X_START]) {
    int i = 0;
    for (int k = 0; k < mpi_params->sizes[2]; k++) {
      for (int j = 0; j < mpi_params->sizes[1]; j++) {
        if (k != 0) {
          float ey_ij =
              GET(ey, i, j, k) +
              cey1 * (GET(hx, i, j, k) - GET(hx, i, j, k - 1)) -
              cey2 * (GET(hz, i, j, k) -
                      received_data[X_START]
                                   [mpi_params->send_array_sizes[X_START >> 1] +
                                    k * hz->ny + j]);
          SET(ey, i, j, k, ey_ij);
        }

        if (j != 0) {
          float ez_ij = GET(ez, i, j, k) +
                        cez1 * (GET(hy, i, j, k) -
                                received_data[X_START][k * hy->ny + j]) -
                        cez2 * (GET(hx, i, j, k) - GET(hx, i, j - 1, k));
          SET(ez, i, j, k, ez_ij);
        }
      }
    }

    if (received_neighbour[Y_START]) {
      int j = 0;
      for (int k = 0; k < mpi_params->sizes[2]; k++) {
        float ez_ij =
            GET(ez, i, j, k) +
            cez1 * (GET(hy, i, j, k) - received_data[X_START][k * hy->ny + j]) -
            cez2 * (GET(hx, i, j, k) - received_data[Y_START][k * hx->nx + i]);
        SET(ez, i, j, k, ez_ij);
      }
    }

    if (received_neighbour[Z_START]) {
      int k = 0;
      for (int j = 0; j < mpi_params->sizes[1]; j++) {
        float ey_ij =
            GET(ey, i, j, k) +
            cey1 * (GET(hx, i, j, k) - received_data[Z_START][j * hx->nx + i]) -
            cey2 * (GET(hz, i, j, k) -
                    received_data[X_START]
                                 [mpi_params->send_array_sizes[X_START >> 1] +
                                  k * hz->ny + j]);
        SET(ey, i, j, k, ey_ij);
      }
    }
  }

  if (received_neighbour[Y_START]) {
    int j = 0;
    for (int k = 0; k < mpi_params->sizes[2]; k++) {
      for (int i = 0; i < mpi_params->sizes[0]; i++) {
        if (k != 0) {
          float ex_ij =
              GET(ex, i, j, k) +
              cex1 * (GET(hz, i, j, k) -
                      received_data[Y_START]
                                   [mpi_params->send_array_sizes[Y_START >> 1] +
                                    k * hz->nx + i]) -
              cex2 * (GET(hy, i, j, k) - GET(hy, i, j, k - 1));
          SET(ex, i, j, k, ex_ij);
        }

        if (i != 0) {
          float ez_ij = GET(ez, i, j, k) +
                        cez1 * (GET(hy, i, j, k) - GET(hy, i - 1, j, k)) -
                        cez2 * (GET(hx, i, j, k) -
                                received_data[Y_START][hx->nx * k + i]);
          SET(ez, i, j, k, ez_ij);
        }
      }
    }

    if (received_neighbour[Z_START]) {
      int k = 0;
      for (int i = 0; i < mpi_params->sizes[0]; i++) {
        float ex_ij =
            GET(ex, i, j, k) +
            cex1 * (GET(hz, i, j, k) -
                    received_data[Y_START]
                                 [mpi_params->send_array_sizes[Y_START >> 1] +
                                  k * hz->nx + i]) -
            cex2 * (GET(hy, i, j, k) -
                    received_data[Z_START]
                                 [mpi_params->send_array_sizes[Z_START >> 1] +
                                  hy->nx * j + i]);
        SET(ex, i, j, k, ex_ij);
      }
    }
  }

  if (received_neighbour[Z_START]) {
    int k = 0;
    for (int j = 0; j < mpi_params->sizes[1]; j++) {
      for (int i = 0; i < mpi_params->sizes[0]; i++) {
        if (j != 0) {
          float ex_ij =
              GET(ex, i, j, k) +
              cex1 * (GET(hz, i, j, k) - GET(hz, i, j - 1, k)) -
              cex2 * (GET(hy, i, j, k) -
                      received_data[Z_START]
                                   [mpi_params->send_array_sizes[Z_START >> 1] +
                                    j * hy->nx + i]);
          SET(ex, i, j, k, ex_ij);
        }

        if (i != 0) {
          float ey_ij = GET(ey, i, j, k) +
                        cey1 * (GET(hx, i, j, k) -
                                received_data[Z_START][j * hx->nx + i]) -
                        cey2 * (GET(hz, i, j, k) - GET(hz, i - 1, j, k));
          SET(ey, i, j, k, ey_ij);
        }
      }
    }
  }
}

static inline void send_data(struct data *dat, MPI_Request *requests,
                             struct MpiParams *mpi_params, float **sent_data) {
  const char *name = dat->name;

  if (strcmp(name, "ex") == 0) {
    if (mpi_params->neighbours[Y_START].rank != -1) {
      for (int k = 0; k < mpi_params->sizes[2]; k++)
        for (int i = 0; i < mpi_params->sizes[0]; i++)
          sent_data[Y_START][k * dat->nx + i] = GET(dat, i, 0, k);

      requests[2 * Y_START] = MPI_REQUEST_NULL;

      MPI_Isend(sent_data[Y_START], mpi_params->send_array_sizes[Y_START >> 1],
                MPI_FLOAT, mpi_params->neighbours[Y_START].rank, 2 * Y_START,
                mpi_params->cart_comm, &(requests[2 * Y_START]));
    }

    if (mpi_params->neighbours[Z_START].rank != -1) {
      for (int j = 0; j < mpi_params->sizes[1]; j++)
        for (int i = 0; i < mpi_params->sizes[0]; i++)
          sent_data[Z_START][j * dat->nx + i] = GET(dat, i, j, 0);

      requests[2 * Z_START] = MPI_REQUEST_NULL;

      MPI_Isend(sent_data[Z_START], mpi_params->send_array_sizes[Z_START >> 1],
                MPI_FLOAT, mpi_params->neighbours[Z_START].rank, 2 * Z_START,
                mpi_params->cart_comm, &(requests[2 * Z_START]));
    }
  }

  if (strcmp(name, "ey") == 0) {
    if (mpi_params->neighbours[X_START].rank != -1) {
      for (int k = 0; k < mpi_params->sizes[2]; k++) {
        for (int j = 0; j < mpi_params->sizes[1]; j++) {
          sent_data[X_START][k * dat->ny + j] = GET(dat, 0, j, k);
        }
      }

      requests[2 * X_START] = MPI_REQUEST_NULL;

      MPI_Isend(sent_data[X_START], mpi_params->send_array_sizes[X_START >> 1],
                MPI_FLOAT, mpi_params->neighbours[X_START].rank, 2 * X_START,
                mpi_params->cart_comm, &(requests[2 * X_START]));
    }

    if (mpi_params->neighbours[Z_START].rank != -1) {
      for (int j = 0; j < mpi_params->sizes[1]; j++)
        for (int i = 0; i < mpi_params->sizes[0]; i++)
          sent_data[Z_START][mpi_params->send_array_sizes[Z_START >> 1] +
                             j * dat->nx + i] = GET(dat, i, j, 0);

      requests[2 * Z_START + 1] = MPI_REQUEST_NULL;

      MPI_Isend(sent_data[Z_START] + mpi_params->send_array_sizes[Z_START >> 1],
                mpi_params->send_array_sizes[Z_START >> 1], MPI_FLOAT,
                mpi_params->neighbours[Z_START].rank, 2 * Z_START + 1,
                mpi_params->cart_comm, &(requests[2 * Z_START + 1]));
    }
  }

  if (strcmp(name, "ez") == 0) {
    if (mpi_params->neighbours[X_START].rank != -1) {
      for (int k = 0; k < mpi_params->sizes[2]; k++)
        for (int j = 0; j < mpi_params->sizes[1]; j++)
          sent_data[X_START][mpi_params->send_array_sizes[X_START >> 1] +
                             k * dat->ny + j] = GET(dat, 0, j, k);

      requests[2 * X_START + 1] = MPI_REQUEST_NULL;

      MPI_Isend(sent_data[X_START] + mpi_params->send_array_sizes[X_START >> 1],
                mpi_params->send_array_sizes[X_START >> 1], MPI_FLOAT,
                mpi_params->neighbours[X_START].rank, 2 * X_START + 1,
                mpi_params->cart_comm, &(requests[2 * X_START + 1]));
    }

    if (mpi_params->neighbours[Y_START].rank != -1) {
      for (int k = 0; k < mpi_params->sizes[2]; k++)
        for (int i = 0; i < mpi_params->sizes[0]; i++)
          sent_data[Y_START][mpi_params->send_array_sizes[Y_START >> 1] +
                             k * dat->nx + i] = GET(dat, i, 0, k);

      requests[2 * Y_START + 1] = MPI_REQUEST_NULL;

      MPI_Isend(sent_data[Y_START] + mpi_params->send_array_sizes[Y_START >> 1],
                mpi_params->send_array_sizes[Y_START >> 1], MPI_FLOAT,
                mpi_params->neighbours[Y_START].rank, 2 * Y_START + 1,
                mpi_params->cart_comm, &(requests[2 * Y_START + 1]));
    }
  }

  if (strcmp(name, "hx") == 0) {
    if (mpi_params->neighbours[Y_END].rank != -1) {
      for (int k = 0; k < mpi_params->sizes[2]; k++)
        for (int i = 0; i < mpi_params->sizes[0]; i++)
          sent_data[Y_END][k * dat->nx + i] =
              GET(dat, i, mpi_params->sizes[Y_END >> 1] - 1, k);

      requests[2 * Y_END] = MPI_REQUEST_NULL;

      MPI_Isend(sent_data[Y_END], mpi_params->send_array_sizes[Y_END >> 1],
                MPI_FLOAT, mpi_params->neighbours[Y_END].rank, 2 * Y_END,
                mpi_params->cart_comm, &(requests[2 * Y_END]));
    }

    if (mpi_params->neighbours[Z_END].rank != -1) {
      for (int j = 0; j < mpi_params->sizes[1]; j++)
        for (int i = 0; i < mpi_params->sizes[0]; i++)
          sent_data[Z_END][j * dat->nx + i] =
              GET(dat, i, j, mpi_params->sizes[Z_END >> 1] - 1);

      requests[2 * Z_END] = MPI_REQUEST_NULL;

      MPI_Isend(sent_data[Z_END], mpi_params->send_array_sizes[Z_END >> 1],
                MPI_FLOAT, mpi_params->neighbours[Z_END].rank, 2 * Z_END,
                mpi_params->cart_comm, &(requests[2 * Z_END]));
    }
  }

  if (strcmp(name, "hy") == 0) {
    if (mpi_params->neighbours[X_END].rank != -1) {
      for (int k = 0; k < mpi_params->sizes[2]; k++)
        for (int j = 0; j < mpi_params->sizes[1]; j++)
          sent_data[X_END][k * dat->ny + j] =
              GET(dat, mpi_params->sizes[0] - 1, j, k);

      requests[2 * X_END] = MPI_REQUEST_NULL;

      MPI_Isend(sent_data[X_END], mpi_params->send_array_sizes[X_END >> 1],
                MPI_FLOAT, mpi_params->neighbours[X_END].rank, 2 * X_END,
                mpi_params->cart_comm, &(requests[2 * X_END]));
    }

    if (mpi_params->neighbours[Z_END].rank != -1) {
      for (int j = 0; j < mpi_params->sizes[1]; j++)
        for (int i = 0; i < mpi_params->sizes[0]; i++)
          sent_data[Z_END][mpi_params->send_array_sizes[Z_END >> 1] +
                           j * dat->nx + i] =
              GET(dat, i, j, mpi_params->sizes[2] - 1);

      requests[2 * Z_END + 1] = MPI_REQUEST_NULL;

      MPI_Isend(sent_data[Z_END] + mpi_params->send_array_sizes[Z_END >> 1],
                mpi_params->send_array_sizes[Z_END >> 1], MPI_FLOAT,
                mpi_params->neighbours[Z_END].rank, 2 * Z_END + 1,
                mpi_params->cart_comm, &(requests[2 * Z_END + 1]));
    }
  }

  if (strcmp(name, "hz") == 0) {
    if (mpi_params->neighbours[X_END].rank != -1) {
      for (int k = 0; k < mpi_params->sizes[2]; k++)
        for (int j = 0; j < mpi_params->sizes[1]; j++)
          sent_data[X_END][mpi_params->send_array_sizes[X_END >> 1] +
                           k * dat->ny + j] =
              GET(dat, mpi_params->sizes[0] - 1, j, k);

      requests[2 * X_END + 1] = MPI_REQUEST_NULL;

      MPI_Isend(sent_data[X_END] + mpi_params->send_array_sizes[X_END >> 1],
                mpi_params->send_array_sizes[X_END >> 1], MPI_FLOAT,
                mpi_params->neighbours[X_END].rank, 2 * X_END + 1,
                mpi_params->cart_comm, &(requests[2 * X_END + 1]));
    }

    if (mpi_params->neighbours[Y_END].rank != -1) {
      for (int k = 0; k < mpi_params->sizes[2]; k++)
        for (int i = 0; i < mpi_params->sizes[0]; i++)
          sent_data[Y_END][mpi_params->send_array_sizes[Y_END >> 1] +
                           k * dat->nx + i] =
              GET(dat, i, mpi_params->sizes[1] - 1, k);

      requests[2 * Y_END + 1] = MPI_REQUEST_NULL;

      MPI_Isend(sent_data[Y_END] + mpi_params->send_array_sizes[Y_END >> 1],
                mpi_params->send_array_sizes[Y_END >> 1], MPI_FLOAT,
                mpi_params->neighbours[Y_END].rank, 2 * Y_END + 1,
                mpi_params->cart_comm, &(requests[2 * Y_END + 1]));
    }
  }
}

inline void receive_data(struct data *dat, float **received_data,
                         bool *received_neighbour, MPI_Request *requests,
                         struct MpiParams *mpi_params) {
  const char *name = dat->name;

  if (strcmp(name, "ex") == 0) {
    if (mpi_params->neighbours[Y_END].rank != -1) {
      requests[2 * Y_END] = MPI_REQUEST_NULL;
      MPI_Irecv(received_data[Y_END], mpi_params->send_array_sizes[Y_END >> 1],
                MPI_FLOAT, mpi_params->neighbours[Y_END].rank, 2 * Y_START,
                mpi_params->cart_comm, &(requests[2 * Y_END]));
      received_neighbour[Y_END] = true;
    }

    if (mpi_params->neighbours[Z_END].rank != -1) {
      requests[2 * Z_END] = MPI_REQUEST_NULL;
      MPI_Irecv(received_data[Z_END], mpi_params->send_array_sizes[Z_END >> 1],
                MPI_FLOAT, mpi_params->neighbours[Z_END].rank, 2 * Z_START,
                mpi_params->cart_comm, &(requests[2 * Z_END]));
      received_neighbour[Z_END] = true;
    }
  }

  if (strcmp(name, "ey") == 0) {
    if (mpi_params->neighbours[X_END].rank != -1) {
      requests[2 * X_END] = MPI_REQUEST_NULL;
      MPI_Irecv(received_data[X_END], mpi_params->send_array_sizes[X_END >> 1],
                MPI_FLOAT, mpi_params->neighbours[X_END].rank, 2 * X_START,
                mpi_params->cart_comm, &(requests[2 * X_END]));
      received_neighbour[X_END] = true;
    }

    if (mpi_params->neighbours[Z_END].rank != -1) {
      requests[2 * Z_END + 1] = MPI_REQUEST_NULL;
      MPI_Irecv(received_data[Z_END] + mpi_params->send_array_sizes[Z_END >> 1],
                mpi_params->send_array_sizes[Z_END >> 1], MPI_FLOAT,
                mpi_params->neighbours[Z_END].rank, 2 * Z_START + 1,
                mpi_params->cart_comm, &(requests[2 * Z_END + 1]));
      received_neighbour[Z_END] = true;
    }
  }

  if (strcmp(name, "ez") == 0) {
    if (mpi_params->neighbours[Y_END].rank != -1) {
      requests[2 * Y_END + 1] = MPI_REQUEST_NULL;
      MPI_Irecv(received_data[Y_END] + mpi_params->send_array_sizes[Y_END >> 1],
                mpi_params->send_array_sizes[Y_END >> 1], MPI_FLOAT,
                mpi_params->neighbours[Y_END].rank, 2 * Y_START + 1,
                mpi_params->cart_comm, &(requests[2 * Y_END + 1]));
      received_neighbour[Y_END] = true;
    }

    if (mpi_params->neighbours[X_END].rank != -1) {
      requests[2 * X_END + 1] = MPI_REQUEST_NULL;
      MPI_Irecv(received_data[X_END] + mpi_params->send_array_sizes[X_END >> 1],
                mpi_params->send_array_sizes[X_END >> 1], MPI_FLOAT,
                mpi_params->neighbours[X_END].rank, 2 * X_START + 1,
                mpi_params->cart_comm, &(requests[2 * X_END + 1]));
      received_neighbour[X_END] = true;
    }
  }

  if (strcmp(name, "hx") == 0) {
    if (mpi_params->neighbours[Y_START].rank != -1) {
      requests[2 * Y_START] = MPI_REQUEST_NULL;
      MPI_Irecv(received_data[Y_START],
                mpi_params->send_array_sizes[Y_START >> 1], MPI_FLOAT,
                mpi_params->neighbours[Y_START].rank, 2 * Y_END,
                mpi_params->cart_comm, &(requests[2 * Y_START]));
      received_neighbour[Y_START] = true;
    }

    if (mpi_params->neighbours[Z_START].rank != -1) {
      requests[2 * Z_START] = MPI_REQUEST_NULL;
      MPI_Irecv(received_data[Z_START],
                mpi_params->send_array_sizes[Z_START >> 1], MPI_FLOAT,
                mpi_params->neighbours[Z_START].rank, 2 * Z_END,
                mpi_params->cart_comm, &(requests[2 * Z_START]));
      received_neighbour[Z_START] = true;
    }
  }

  if (strcmp(name, "hy") == 0) {
    if (mpi_params->neighbours[X_START].rank != -1) {
      requests[2 * X_START] = MPI_REQUEST_NULL;
      MPI_Irecv(received_data[X_START],
                mpi_params->send_array_sizes[X_START >> 1], MPI_FLOAT,
                mpi_params->neighbours[X_START].rank, 2 * X_END,
                mpi_params->cart_comm, &(requests[2 * X_START]));
      received_neighbour[X_START] = true;
    }

    if (mpi_params->neighbours[Z_START].rank != -1) {
      requests[2 * Z_START + 1] = MPI_REQUEST_NULL;
      MPI_Irecv(received_data[Z_START] +
                    mpi_params->send_array_sizes[Z_START >> 1],
                mpi_params->send_array_sizes[Z_START >> 1], MPI_FLOAT,
                mpi_params->neighbours[Z_START].rank, 2 * Z_END + 1,
                mpi_params->cart_comm, &(requests[2 * Z_START + 1]));
      received_neighbour[Z_START] = true;
    }
  }

  if (strcmp(name, "hz") == 0) {
    if (mpi_params->neighbours[X_START].rank != -1) {
      requests[2 * X_START + 1] = MPI_REQUEST_NULL;
      MPI_Irecv(received_data[X_START] +
                    mpi_params->send_array_sizes[X_START >> 1],
                mpi_params->send_array_sizes[X_START >> 1], MPI_FLOAT,
                mpi_params->neighbours[X_START].rank, 2 * X_END + 1,
                mpi_params->cart_comm, &(requests[2 * X_START + 1]));
      received_neighbour[X_START] = true;
    }

    if (mpi_params->neighbours[Y_START].rank != -1) {
      requests[2 * Y_START + 1] = MPI_REQUEST_NULL;
      MPI_Irecv(received_data[Y_START] +
                    mpi_params->send_array_sizes[Y_START >> 1],
                mpi_params->send_array_sizes[Y_START >> 1], MPI_FLOAT,
                mpi_params->neighbours[Y_START].rank, 2 * Y_END + 1,
                mpi_params->cart_comm, &(requests[2 * Y_START + 1]));
      received_neighbour[Y_START] = true;
    }
  }
}

int solve(struct SimulationParams *sim_params,
          struct PhysicalParams *phys_params, struct PerformanceData *perf_data,
          struct MpiParams *mpi_params) {
  int nb_neighbours = sim_params->ndim * 2;

  struct area *proc_area = get_area(sim_params, mpi_params);
  for (int i = 0; i < sim_params->ndim; i++) {
    mpi_params->sizes[i] = proc_area->end[i] - proc_area->start[i] + 1;
  }
  mpi_params->send_sizes[0][0] = mpi_params->sizes[1];
  mpi_params->send_sizes[0][1] = mpi_params->sizes[2];
  mpi_params->send_sizes[1][0] = mpi_params->sizes[0];
  mpi_params->send_sizes[1][1] = mpi_params->sizes[2];
  mpi_params->send_sizes[2][0] = mpi_params->sizes[0];
  mpi_params->send_sizes[2][1] = mpi_params->sizes[1];
  for (int i = 0; i < sim_params->ndim; i++)
    mpi_params->send_array_sizes[i] =
        mpi_params->send_sizes[i][0] * mpi_params->send_sizes[i][1];

  DEBUG_PRINT(
      "Proc %d, nx = %d (%d to %d), ny = %d (%d to %d), nz = %d (%d to %d)\n",
      mpi_params->rank, mpi_params->sizes[0], proc_area->start[0],
      proc_area->end[0], mpi_params->sizes[1], proc_area->start[1],
      proc_area->end[1], mpi_params->sizes[2], proc_area->start[2],
      proc_area->end[2]);
  double *origins = NULL;
  origins = malloc(sizeof(double) * sim_params->ndim);
  if (!origins) {
    printf("Not enough memory\n");
    return EXIT_FAILURE;
  }

  for (int i = 0; i < sim_params->ndim; i++) {
    origins[i] = proc_area->start[i] * sim_params->steps[i];
    // TODO: fix when more than 2 procs in a given dir of space
    if (origins[i] != 0)
      origins[i] -= 1 * sim_params->steps[i];
  }

  // Initializing data
  struct data ex, ey, ez, hx, hy, hz;
  float val = 0.;
  if (init_data(&ex, "ex", mpi_params->sizes[0], mpi_params->sizes[1],
                mpi_params->sizes[2], sim_params->steps[0],
                sim_params->steps[1], sim_params->steps[2], origins[0],
                origins[1], origins[2], val) ||
      init_data(&ey, "ey", mpi_params->sizes[0], mpi_params->sizes[1],
                mpi_params->sizes[2], sim_params->steps[0],
                sim_params->steps[1], sim_params->steps[2], origins[0],
                origins[1], origins[2], val) ||
      init_data(&ez, "ez", mpi_params->sizes[0], mpi_params->sizes[1],
                mpi_params->sizes[2], sim_params->steps[0],
                sim_params->steps[1], sim_params->steps[2], origins[0],
                origins[1], origins[2], val) ||
      init_data(&hx, "hx", mpi_params->sizes[0], mpi_params->sizes[1],
                mpi_params->sizes[2], sim_params->steps[0],
                sim_params->steps[1], sim_params->steps[2], origins[0],
                origins[1], origins[2], val) ||
      init_data(&hy, "hy", mpi_params->sizes[0], mpi_params->sizes[1],
                mpi_params->sizes[2], sim_params->steps[0],
                sim_params->steps[1], sim_params->steps[2], origins[0],
                origins[1], origins[2], val) ||
      init_data(&hz, "hz", mpi_params->sizes[0], mpi_params->sizes[1],
                mpi_params->sizes[2], sim_params->steps[0],
                sim_params->steps[1], sim_params->steps[2], origins[0],
                origins[1], origins[2], val)) {
    printf("Error: could not allocate data\n");
    return EXIT_FAILURE;
  }

  // Setting up MPI requests variables
  // Note the initialization to NULL to be able to free even if not malloc yet
  MPI_Request *send_requests = NULL;
  MPI_Request *recv_requests = NULL;
  float **sent_data = NULL;
  float **received_data = NULL;
  bool *received_neighbour = NULL;
  received_neighbour = malloc(sizeof(bool) * nb_neighbours);
  send_requests = malloc(sizeof(MPI_Request) * nb_neighbours * 2);
  recv_requests = malloc(sizeof(MPI_Request) * nb_neighbours * 2);
  sent_data = malloc(sizeof(float *) * nb_neighbours);
  received_data = malloc(sizeof(float *) * nb_neighbours);

  if (!send_requests || !recv_requests || !sent_data || !received_data ||
      !received_neighbour) {
    printf("Error: allocation problem in initial MPI structs init\n");
    free_all_solve_pointers(nb_neighbours, proc_area, &hx, &hy, &hz, &ex, &ey,
                            &ez, send_requests, recv_requests, sent_data,
                            received_data, received_neighbour, origins);
    return EXIT_FAILURE;
  }

  for (int i = 0; i < nb_neighbours; i++) {
    received_neighbour[i] = false;
    sent_data[i] = NULL;
    received_data[i] = NULL;
  }
  for (int i = 0; i < nb_neighbours; i++) {
    int index = i >> 1;

    // ROW MAJOR
    // 2 fields per array. First one is the one with the lowest index (x, y,
    // z)
    sent_data[i] =
        calloc(2 * mpi_params->send_array_sizes[index], sizeof(float));
    received_data[i] =
        calloc(2 * mpi_params->send_array_sizes[index], sizeof(float));

    if (!sent_data[i] || !received_data[i]) {
      printf("Error: allocation problem in second MPI structs init\n");
      free_all_solve_pointers(nb_neighbours, proc_area, &hx, &hy, &hz, &ex, &ey,
                              &ez, send_requests, recv_requests, sent_data,
                              received_data, received_neighbour, origins);
      return EXIT_FAILURE;
    }
  }

  // Initializing the arrays
  for (int i = 0; i < nb_neighbours * 2; i++) {
    recv_requests[i] = MPI_REQUEST_NULL;
    send_requests[i] = MPI_REQUEST_NULL;
  }

  double start = GET_TIME();

  // Time loop
  for (int t = 0; t < sim_params->size_of_space[sim_params->ndim]; t++) {
    if (mpi_params->rank == 0 && t &&
        (t % (sim_params->size_of_space[sim_params->ndim] / 10)) == 0) {
      double time_sofar = GET_TIME() - start;
      double eta =
          (sim_params->size_of_space[sim_params->ndim] - t) * time_sofar / t;
      printf("Computing time step %d/%d (ETA: %g seconds)     \r", t,
             sim_params->size_of_space[sim_params->ndim], eta);
      fflush(stdout);
    }

    // Send E field, then receive it for the H loop
    receive_data(&ex, received_data, received_neighbour, recv_requests,
                 mpi_params);
    receive_data(&ey, received_data, received_neighbour, recv_requests,
                 mpi_params);
    receive_data(&ez, received_data, received_neighbour, recv_requests,
                 mpi_params);

    send_data(&ex, send_requests, mpi_params, sent_data);
    send_data(&ey, send_requests, mpi_params, sent_data);
    send_data(&ez, send_requests, mpi_params, sent_data);

    h_loop(&hx, &hy, &hz, &ex, &ey, &ez, sim_params, phys_params, mpi_params,
           received_data, received_neighbour, recv_requests);

    // Send H field
    send_data(&hx, send_requests, mpi_params, sent_data);
    send_data(&hy, send_requests, mpi_params, sent_data);
    send_data(&hz, send_requests, mpi_params, sent_data);

    // Receive H field for the E loop
    receive_data(&hx, received_data, received_neighbour, recv_requests,
                 mpi_params);
    receive_data(&hy, received_data, received_neighbour, recv_requests,
                 mpi_params);
    receive_data(&hz, received_data, received_neighbour, recv_requests,
                 mpi_params);

    e_loop(&hx, &hy, &hz, &ex, &ey, &ez, sim_params, phys_params, mpi_params,
           received_data, received_neighbour, recv_requests);

    // To make sure loops are a bit coordinated
    MPI_Waitall(nb_neighbours, send_requests, MPI_STATUSES_IGNORE);

    // impose source
    int source_x = sim_params->size_of_space[0] / 2;
    int source_y = sim_params->size_of_space[1] / 2;
    int source_z = sim_params->size_of_space[2] / 2;
    double n = t * sim_params->steps[sim_params->ndim];
    if ((proc_area->start[0] <= source_x && proc_area->end[0] >= source_x) &&
        (proc_area->start[1] <= source_y && proc_area->end[1] >= source_y) &&
        (proc_area->start[2] <= source_z && proc_area->end[2] >= source_z)) {
      // DEBUG_PRINT("Process %d applying source at %d, %d\n",
      //             mpi_params->rank + 1, source_x, source_y);
      switch (sim_params->problem_id) {
      case 1:
      case 2:
      case 3:
      case 4:
        // sinusoidal excitation at 2.4 GHz in the middle of the domain
        int idx = source_x - proc_area->start[0],
            idy = source_y - proc_area->start[1],
            idz = source_z - proc_area->start[2];
        float source = sin(2. * M_PI * 2.4e9 * n);
        SET(&ex, idx, idy, idz, source);
        SET(&ey, idx, idy, idz, source);
        SET(&ez, idx, idy, idz, source);
        break;
      default:
        printf("Error: unknown source\n");
        break;
      }
    }

    // output step data in VTK format
    if (sim_params->sampling_rate && !(t % sim_params->sampling_rate)) {
      write_data_vtk(&ex, t, mpi_params->rank);
      write_data_vtk(&ey, t, mpi_params->rank);
      write_data_vtk(&ez, t, mpi_params->rank);
    }
  }

  if (sim_params->sampling_rate && mpi_params->rank == 0) {
    write_manifest_vtk("ex", sim_params->steps[sim_params->ndim],
                       sim_params->size_of_space[sim_params->ndim],
                       sim_params->sampling_rate, mpi_params->num_ranks);
    write_manifest_vtk("ey", sim_params->steps[sim_params->ndim],
                       sim_params->size_of_space[sim_params->ndim],
                       sim_params->sampling_rate, mpi_params->num_ranks);
    write_manifest_vtk("ez", sim_params->steps[sim_params->ndim],
                       sim_params->size_of_space[sim_params->ndim],
                       sim_params->sampling_rate, mpi_params->num_ranks);
  }
  double time = GET_TIME() - start;

  double MUps_per_sec =
      1.e-6 *
      ((double)(mpi_params->sizes[0]) * (double)(mpi_params->sizes[1]) *
       (double)(mpi_params->sizes[2]) *
       (double)sim_params->size_of_space[sim_params->ndim]) /
      time;
  perf_data->time = time;
  perf_data->MUps_per_sec = MUps_per_sec;

  free_all_solve_pointers(nb_neighbours, proc_area, &hx, &hy, &hz, &ex, &ey,
                          &ez, send_requests, recv_requests, sent_data,
                          received_data, received_neighbour, origins);

  printf("\n");

  return EXIT_SUCCESS;
}
