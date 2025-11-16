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
                             struct data *hx, struct data *hy, struct data *ez,
                             MPI_Request *send_requests,
                             MPI_Request *recv_requests, float **sent_data,
                             float **received_data, bool *received_neighbour,
                             double *origins) {
  free_area(proc_area);
  free_data(hx);
  free_data(hy);
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

inline void hx_loop(struct data *hx, struct data *ez,
                    struct SimulationParams *sim_params,
                    struct PhysicalParams *phys_params,
                    struct MpiParams *mpi_params, float **received_data,
                    bool *received_neighbour, MPI_Request *requests) {
  float chy = sim_params->steps[sim_params->ndim] /
              (sim_params->steps[1] * phys_params->mu);
  for (int j = 0; j < mpi_params->sizes[1] - 1; j++) {
    for (int i = 0; i < mpi_params->sizes[0]; i++) {
      float hx_ij = GET(hx, i, j) - chy * (GET(ez, i, j + 1) - GET(ez, i, j));
      SET(hx, i, j, hx_ij);
    }
  }

  MPI_Wait(&(requests[Y_END]), MPI_STATUS_IGNORE);
  if (received_neighbour[Y_END]) {
    for (int i = 0; i < mpi_params->sizes[0]; i++) {
      float hx_ij = GET(hx, i, mpi_params->sizes[1] - 1) -
                    chy * (received_data[Y_END][i] -
                           GET(ez, i, mpi_params->sizes[1] - 1));
      SET(hx, i, mpi_params->sizes[1] - 1, hx_ij);
    }
  } // Don't update if no neighbours.
}

inline void hy_loop(struct data *hy, struct data *ez,
                    struct SimulationParams *sim_params,
                    struct PhysicalParams *phys_params,
                    struct MpiParams *mpi_params, float **received_data,
                    bool *received_neighbour, MPI_Request *requests) {
  float chx = sim_params->steps[sim_params->ndim] /
              (sim_params->steps[0] * phys_params->mu);
  for (int j = 0; j < mpi_params->sizes[1]; j++) {
    for (int i = 0; i < mpi_params->sizes[0] - 1; i++) {
      float hy_ij = GET(hy, i, j) + chx * (GET(ez, i + 1, j) - GET(ez, i, j));
      SET(hy, i, j, hy_ij);
    }
  }

  MPI_Wait(&(requests[X_END]), MPI_STATUS_IGNORE);
  if (received_neighbour[X_END]) {
    for (int j = 0; j < mpi_params->sizes[1]; j++) {
      float hy_ij = GET(hy, mpi_params->sizes[0] - 1, j) +
                    chx * (received_data[X_END][j] -
                           GET(ez, mpi_params->sizes[0] - 1, j));
      SET(hy, mpi_params->sizes[0] - 1, j, hy_ij);
    }
  }
}

inline void ez_loop(struct data *hx, struct data *hy, struct data *ez,
                    struct SimulationParams *sim_params,
                    struct PhysicalParams *phys_params,
                    struct MpiParams *mpi_params, float **received_data,
                    bool *received_neighbour, MPI_Request *requests) {
  float cex = sim_params->steps[sim_params->ndim] /
              (sim_params->steps[0] * phys_params->eps),
        cey = sim_params->steps[sim_params->ndim] /
              (sim_params->steps[1] * phys_params->eps);
  for (int j = 1; j < mpi_params->sizes[1]; j++) {
    for (int i = 1; i < mpi_params->sizes[0]; i++) {
      float ez_ij = GET(ez, i, j) + cex * (GET(hy, i, j) - GET(hy, i - 1, j)) -
                    cey * (GET(hx, i, j) - GET(hx, i, j - 1));
      SET(ez, i, j, ez_ij);
    }
  }

  MPI_Wait(&(requests[X_START]), MPI_STATUS_IGNORE);
  if (received_neighbour[X_START]) {
    for (int j = 1; j < mpi_params->sizes[1]; j++) {
      float ez_ij = GET(ez, 0, j) +
                    cex * (GET(hy, 0, j) - received_data[X_START][j]) -
                    cey * (GET(hx, 0, j) - GET(hx, 0, j - 1));
      SET(ez, 0, j, ez_ij);
    }

    MPI_Wait(&(requests[Y_START]), MPI_STATUS_IGNORE);
    if (received_neighbour[Y_START]) {
      float ez_ij = GET(ez, 0, 0) +
                    cex * (GET(hy, 0, 0) - received_data[X_START][0]) -
                    cey * (GET(hx, 0, 0) - received_data[Y_START][0]);
      SET(ez, 0, 0, ez_ij);
    }
  }

  MPI_Wait(&(requests[Y_START]), MPI_STATUS_IGNORE);
  if (received_neighbour[Y_START]) {
    for (int i = 1; i < mpi_params->sizes[0]; i++) {
      float ez_ij = GET(ez, i, 0) + cex * (GET(hy, i, 0) - GET(hy, i - 1, 0)) -
                    cey * (GET(hx, i, 0) - received_data[Y_START][i]);
      SET(ez, i, 0, ez_ij);
    }
  }
}

inline void send_data(struct data *dat, MPI_Request *requests,
                      struct MpiParams *mpi_params, float **sent_data) {
  const char *name = dat->name;

  if (strcmp(name, "ez") == 0) {
    // Copy to good sent_data
    if (mpi_params->neighbours[X_START].rank != -1) {
      for (int i = 0; i < mpi_params->send_sizes[0]; i++)
        sent_data[X_START][i] = GET(dat, 0, i);

      requests[X_START] = MPI_REQUEST_NULL;
      MPI_Isend(sent_data[X_START], mpi_params->send_sizes[0], MPI_FLOAT,
                mpi_params->neighbours[X_START].rank, X_START,
                mpi_params->cart_comm, &(requests[X_START]));
    }

    if (mpi_params->neighbours[Y_START].rank != -1) {
      for (int i = 0; i < mpi_params->send_sizes[1]; i++)
        sent_data[Y_START][i] = GET(dat, i, 0);

      requests[Y_START] = MPI_REQUEST_NULL;
      MPI_Isend(sent_data[Y_START], mpi_params->send_sizes[1], MPI_FLOAT,
                mpi_params->neighbours[Y_START].rank, Y_START,
                mpi_params->cart_comm, &(requests[Y_START]));
    }

    // Send data to good neighbours
  } else if (strcmp(name, "hx") == 0 &&
             mpi_params->neighbours[Y_END].rank != -1) {
    for (int i = 0; i < mpi_params->send_sizes[1]; i++)
      sent_data[Y_END][i] = GET(dat, i, mpi_params->sizes[1] - 1);

    requests[Y_END] = MPI_REQUEST_NULL;

    MPI_Isend(sent_data[Y_END], mpi_params->send_sizes[1], MPI_FLOAT,
              mpi_params->neighbours[Y_END].rank, Y_END, mpi_params->cart_comm,
              &(requests[Y_END]));
  } else if (strcmp(name, "hy") == 0 &&
             mpi_params->neighbours[X_END].rank != -1) {
    for (int i = 0; i < mpi_params->send_sizes[0]; i++)
      sent_data[X_END][i] = GET(dat, mpi_params->sizes[0] - 1, i);

    requests[X_END] = MPI_REQUEST_NULL;

    MPI_Isend(sent_data[X_END], mpi_params->send_sizes[0], MPI_FLOAT,
              mpi_params->neighbours[X_END].rank, X_END, mpi_params->cart_comm,
              &(requests[X_END]));
  }
}

inline void receive_data(struct data *dat, float **received_data,
                         bool *received_neighbour, MPI_Request *requests,
                         struct MpiParams *mpi_params) {
  const char *name = dat->name;

  if (strcmp(name, "ez") == 0) {
    if (mpi_params->neighbours[X_END].rank != -1) {
      requests[X_END] = MPI_REQUEST_NULL;
      MPI_Irecv(received_data[X_END], mpi_params->send_sizes[0], MPI_FLOAT,
                mpi_params->neighbours[X_END].rank, X_START,
                mpi_params->cart_comm, &(requests[X_END]));
      received_neighbour[X_END] = true;
    }

    if (mpi_params->neighbours[Y_END].rank != -1) {
      requests[Y_END] = MPI_REQUEST_NULL;
      MPI_Irecv(received_data[Y_END], mpi_params->send_sizes[1], MPI_FLOAT,
                mpi_params->neighbours[Y_END].rank, Y_START,
                mpi_params->cart_comm, &(requests[Y_END]));
      received_neighbour[Y_END] = true;
    }
  } else if (strcmp(name, "hx") == 0 &&
             mpi_params->neighbours[Y_START].rank != -1) {
    requests[Y_START] = MPI_REQUEST_NULL;
    MPI_Irecv(received_data[Y_START], mpi_params->send_sizes[1], MPI_FLOAT,
              mpi_params->neighbours[Y_START].rank, Y_END,
              mpi_params->cart_comm, &(requests[Y_START]));
    received_neighbour[Y_START] = true;
  } else if (strcmp(name, "hy") == 0 &&
             mpi_params->neighbours[X_START].rank != -1) {
    requests[X_START] = MPI_REQUEST_NULL;
    MPI_Irecv(received_data[X_START], mpi_params->send_sizes[0], MPI_FLOAT,
              mpi_params->neighbours[X_START].rank, X_END,
              mpi_params->cart_comm, &(requests[X_START]));
    received_neighbour[X_START] = true;
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
  mpi_params->send_sizes[0] = mpi_params->sizes[1];
  mpi_params->send_sizes[1] = mpi_params->sizes[0];

  DEBUG_PRINT("Proc %d, nx = %d (%d to %d), ny = %d (%d to %d)\n",
              mpi_params->rank, mpi_params->sizes[0], proc_area->start[0],
              proc_area->end[0], mpi_params->sizes[1], proc_area->start[1],
              proc_area->end[1]);
  double *origins = NULL;
  origins = malloc(sizeof(double) * sim_params->ndim);
  if (!origins) {
    printf("Not enough memory\n");
    return EXIT_FAILURE;
  }

  for (int i = 0; i < sim_params->ndim; i++) {
    origins[i] = proc_area->start[i] * sim_params->steps[i];
    if (origins[i] != 0)
      origins[i] -= 1 * sim_params->steps[i];
  }

  struct data ez, hx, hy;
  float val = 0.;
  // nx, ny
  if (init_data(&ez, "ez", mpi_params->sizes[0], mpi_params->sizes[1],
                sim_params->steps[0], sim_params->steps[1], origins[0],
                origins[1], val) ||
      init_data(&hx, "hx", mpi_params->sizes[0], mpi_params->sizes[1],
                sim_params->steps[0], sim_params->steps[1], origins[0],
                origins[1], val) ||
      init_data(&hy, "hy", mpi_params->sizes[0], mpi_params->sizes[1],
                sim_params->steps[0], sim_params->steps[1], origins[0],
                origins[1], val)) {
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
  send_requests = malloc(sizeof(MPI_Request) * nb_neighbours);
  recv_requests = malloc(sizeof(MPI_Request) * nb_neighbours);
  sent_data = malloc(sizeof(float *) * nb_neighbours);
  received_data = malloc(sizeof(float *) * nb_neighbours);

  if (!send_requests || !recv_requests || !sent_data || !received_data ||
      !received_neighbour) {
    printf("Error: allocation problem in initial MPI structs init\n");
    free_all_solve_pointers(nb_neighbours, proc_area, &hx, &hy, &ez,
                            send_requests, recv_requests, sent_data,
                            received_data, received_neighbour, origins);
    return EXIT_FAILURE;
  }

  for (int i = 0; i < nb_neighbours; i++) {
    received_neighbour[i] = false;
    sent_data[i] = NULL;
    received_data[i] = NULL;
  }
  for (int i = 0; i < nb_neighbours; i++) {
    // an x line is of length ny, inversly for an y line
    // * 3 for hx, hy, ez
    int size = (i < 2) ? mpi_params->send_sizes[0] : mpi_params->send_sizes[1];

    sent_data[i] = malloc(sizeof(float) * size);
    received_data[i] = malloc(sizeof(float) * size);

    if (!sent_data[i] || !received_data[i]) {
      printf("Error: allocation problem in second MPI structs init\n");
      free_all_solve_pointers(nb_neighbours, proc_area, &hx, &hy, &ez,
                              send_requests, recv_requests, sent_data,
                              received_data, received_neighbour, origins);
      return EXIT_FAILURE;
    }
  }

  // Initializing the arrays
  for (int i = 0; i < nb_neighbours; i++) {
    recv_requests[i] = MPI_REQUEST_NULL;
    send_requests[i] = MPI_REQUEST_NULL;

    int size = (i < 2) ? mpi_params->send_sizes[0] : mpi_params->send_sizes[1];
    for (int j = 0; j < size; j++) {
      sent_data[i][j] = 0.;
      received_data[i][j] = 0.;
    }
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

    send_data(&ez, send_requests, mpi_params, sent_data);
    receive_data(&ez, received_data, received_neighbour, recv_requests,
                 mpi_params);

    hx_loop(&hx, &ez, sim_params, phys_params, mpi_params, received_data,
            received_neighbour, recv_requests);

    send_data(&hx, send_requests, mpi_params, sent_data);

    hy_loop(&hy, &ez, sim_params, phys_params, mpi_params, received_data,
            received_neighbour, recv_requests);

    send_data(&hy, send_requests, mpi_params, sent_data);

    receive_data(&hx, received_data, received_neighbour, recv_requests,
                 mpi_params);
    receive_data(&hy, received_data, received_neighbour, recv_requests,
                 mpi_params);

    ez_loop(&hx, &hy, &ez, sim_params, phys_params, mpi_params, received_data,
            received_neighbour, recv_requests);

    MPI_Waitall(nb_neighbours, send_requests, MPI_STATUSES_IGNORE);

    // impose source
    int source_x = sim_params->size_of_space[0] / 2;
    int source_y = sim_params->size_of_space[1] / 2;
    double n = t * sim_params->steps[sim_params->ndim];
    if ((proc_area->start[0] <= source_x && proc_area->end[0] >= source_x) &&
        (proc_area->start[1] <= source_y && proc_area->end[1] >= source_y)) {
      // DEBUG_PRINT("Process %d applying source at %d, %d\n",
      //             mpi_params->rank + 1, source_x, source_y);
      switch (sim_params->problem_id) {
      case 1:
      case 2:
        // sinusoidal excitation at 2.4 GHz in the middle of the domain
        SET(&ez, source_x - proc_area->start[0], source_y - proc_area->start[1],
            sin(2. * M_PI * 2.4e9 * n));
        break;
      case 3:
      case 4:
        SET(&ez, source_x - proc_area->start[0], source_y - proc_area->start[1],
            sin(2. * M_PI * 2.4e9 * n));
        break;
      default:
        printf("Error: unknown source\n");
        break;
      }
    }

    // output step data in VTK format
    if (sim_params->sampling_rate && !(t % sim_params->sampling_rate)) {
      write_data_vtk(&ez, t, mpi_params->rank);
      // write_data_vtk(&hx, t, mpi_params->rank);
      // write_data_vtk(&hy, t, mpi_params->rank);
    }
  }

  if (mpi_params->rank == 0) {
    write_manifest_vtk("ez", sim_params->steps[sim_params->ndim],
                       sim_params->size_of_space[sim_params->ndim],
                       sim_params->sampling_rate, mpi_params->num_ranks);
    // write_manifest_vtk("hx", sim_params->steps[sim_params->ndim],
    //                    sim_params->size_of_space[sim_params->ndim],
    //                    sim_params->sampling_rate, mpi_params->num_ranks);
    // write_manifest_vtk("hy", sim_params->steps[sim_params->ndim],
    //                    sim_params->size_of_space[sim_params->ndim],
    //                    sim_params->sampling_rate, mpi_params->num_ranks);
  }
  double time = GET_TIME() - start;

  double MUps_per_sec =
      1.e-6 *
      ((double)(mpi_params->sizes[0]) * (double)(mpi_params->sizes[1]) *
       (double)sim_params->size_of_space[sim_params->ndim]) /
      time;
  perf_data->time = time;
  perf_data->MUps_per_sec = MUps_per_sec;

  free_all_solve_pointers(nb_neighbours, proc_area, &hx, &hy, &ez,
                          send_requests, recv_requests, sent_data,
                          received_data, received_neighbour, origins);

  printf("\n");

  return EXIT_SUCCESS;
}
