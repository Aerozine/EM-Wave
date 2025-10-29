#include "solver.h"
#include "data.h"
#include "main.h"
#include "params.h"

#include <math.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

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
      int end = start + size;

      current_area->start[i] = start;
      current_area->end[i] = end;

      // Getting neighbours
      int negative, positive;
      MPI_Cart_shift(mpi_params->cart_comm, i, 1, &negative, &positive);
      printf("Proc %d, dim %d, neg %d, pos %d\n", mpi_params->rank + 1, i,
             negative, positive);

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
                             MPI_Request *recv_requests, double **sent_data,
                             double **received_data, int *sizes) {
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
  free(sizes);
}

void send_data_new(int nb_neighbours, double **sent_data, int *sizes,
                   struct data *hx, struct data *hy, struct data *ez,
                   struct MpiParams *mpi_params, MPI_Request *requests) {
  for (int i = 0; i < nb_neighbours; i++)
    requests[i] = MPI_REQUEST_NULL;

  // First, we set the sent_data

  // XSTART
  for (int i = 0; i < sizes[X_START] - 1; i++) {
    // hx
    sent_data[X_START][i] = GET(hx, 0, i);
    // hy
    sent_data[X_START][sizes[X_START] + i] = GET(hy, 0, i);
    // ez
    sent_data[X_START][sizes[X_START] * 2 + i] = GET(ez, 0, i);
  }

  // XEND
  for (int i = 0; i < sizes[X_END] - 1; i++) {
    // hx
    sent_data[X_END][i] = GET(hx, sizes[Y_START] - 1, i);
    // hy
    sent_data[X_END][sizes[X_END] + i] = GET(hy, sizes[Y_START] - 2, i);
    // ez
    sent_data[X_END][sizes[X_END] * 2 + i] = GET(ez, sizes[Y_START] - 1, i);
  }

  // YSTART
  for (int i = 0; i < sizes[Y_START] - 1; i++) {
    // hx
    sent_data[Y_START][i] = GET(hx, i, 0);
    // hy
    sent_data[Y_START][sizes[Y_START] + i] = GET(hy, i, 0);
    // ez
    sent_data[Y_START][sizes[Y_START] * 2 + i] = GET(ez, i, 0);
  }

  // YEND
  for (int i = 0; i < sizes[Y_END] - 1; i++) {
    // hx
    sent_data[Y_END][i] = GET(hx, i, sizes[X_START] - 2);
    // hy
    sent_data[Y_END][sizes[Y_END] + i] = GET(hy, i, sizes[X_START] - 1);
    // ez
    sent_data[Y_END][sizes[Y_END] * 2 + i] = GET(ez, i, sizes[X_START] - 1);
  }

  // Now, we can send all of that
  for (int i = 0; i < nb_neighbours; i++) {
    struct neighbour current_neighbour = mpi_params->neighbours[i];
    if (current_neighbour.rank != -1) {
      switch (current_neighbour.pos) {
      case X_START:
        MPI_Isend(sent_data[X_START], sizes[X_START] * 3, MPI_DOUBLE,
                  current_neighbour.rank, X_START, mpi_params->cart_comm,
                  &(requests[X_START]));
        break;
      case X_END:
        MPI_Isend(sent_data[X_END], sizes[X_END] * 3, MPI_DOUBLE,
                  current_neighbour.rank, X_END, mpi_params->cart_comm,
                  &(requests[X_END]));
        break;
      case Y_START:
        MPI_Isend(sent_data[Y_START], sizes[Y_START] * 3, MPI_DOUBLE,
                  current_neighbour.rank, Y_START, mpi_params->cart_comm,
                  &(requests[Y_START]));
        break;
      case Y_END:
        MPI_Isend(sent_data[Y_END], sizes[Y_END] * 3, MPI_DOUBLE,
                  current_neighbour.rank, Y_END, mpi_params->cart_comm,
                  &(requests[Y_END]));
        break;
      default:
        break;
      }
    }
  }
}

int receive_data_new(int nb_neighbours, double **received_data, int *sizes,
                     struct data *hx, struct data *hy, struct data *ez,
                     struct MpiParams *mpi_params, MPI_Request *requests) {
  for (int i = 0; i < nb_neighbours; i++)
    requests[i] = MPI_REQUEST_NULL;

  bool *received_neighbour = malloc(sizeof(bool) * nb_neighbours);
  if (!received_neighbour) {
    return EXIT_FAILURE;
  }
  for (int i = 0; i < nb_neighbours; i++)
    received_neighbour[i] = false;

  // Iterate over all neighbours
  for (int i = 0; i < nb_neighbours; i++) {
    struct neighbour current_neighbour = mpi_params->neighbours[i];
    if (current_neighbour.rank != -1) {
      switch (current_neighbour.pos) {
      case X_START:
        MPI_Irecv(received_data[X_START], sizes[X_START] * 3, MPI_DOUBLE,
                  current_neighbour.rank, X_END, mpi_params->cart_comm,
                  &(requests[X_START]));
        received_neighbour[X_START] = true;
        break;
      case X_END:
        MPI_Irecv(received_data[X_END], sizes[X_END] * 3, MPI_DOUBLE,
                  current_neighbour.rank, X_START, mpi_params->cart_comm,
                  &(requests[X_END]));
        received_neighbour[X_END] = true;
        break;
      case Y_START:
        MPI_Irecv(received_data[Y_START], sizes[Y_START] * 3, MPI_DOUBLE,
                  current_neighbour.rank, Y_END, mpi_params->cart_comm,
                  &(requests[Y_START]));
        received_neighbour[Y_START] = true;
        break;
      case Y_END:
        MPI_Irecv(received_data[Y_END], sizes[Y_END] * 3, MPI_DOUBLE,
                  current_neighbour.rank, Y_START, mpi_params->cart_comm,
                  &(requests[Y_END]));
        received_neighbour[Y_END] = true;
        break;
      default:
        break;
      }
    }
  }

  for (int i = 0; i < nb_neighbours; i++) {
    if (received_neighbour[i]) {
      switch (i) {
      case X_START:
        for (int j = 0; j < sizes[i] - 1; j++) {
          SET(hx, 0, j, received_data[i][j]);
          SET(hy, 0, j, received_data[i][sizes[i] + j]);
          SET(ez, 0, j, received_data[i][2 * sizes[i] + j]);
        }
        break;
      case X_END:
        for (int j = 0; j < sizes[i] - 1; j++) {
          SET(hx, sizes[Y_START] - 1, j, received_data[i][j]);
          SET(hy, sizes[Y_START] - 2, j, received_data[i][sizes[i] + j]);
          SET(ez, sizes[Y_START] - 1, j, received_data[i][2 * sizes[i] + j]);
        }
        break;
      case Y_START:
        for (int j = 0; j < sizes[i] - 1; j++) {
          SET(hx, j, 0, received_data[i][j]);
          SET(hy, j, 0, received_data[i][sizes[i] + j]);
          SET(ez, j, 0, received_data[i][2 * sizes[i] + j]);
        }
        break;
      case Y_END:
        for (int j = 0; j < sizes[i] - 1; j++) {
          SET(hx, j, sizes[X_START] - 2, received_data[i][j]);
          SET(hy, j, sizes[X_START] - 1, received_data[i][sizes[i] + j]);
          SET(ez, j, sizes[X_START] - 1, received_data[i][2 * sizes[i] + j]);
        }
        break;
      default:
        break;
      }
    }
  }

  free(received_neighbour);
  return EXIT_SUCCESS;
}

int solve(struct SimulationParams *sim_params,
          struct PhysicalParams *phys_params, struct PerformanceData *perf_data,
          struct MpiParams *mpi_params) {
  int nb_neighbours = sim_params->ndim * 2;

  struct area *proc_area = get_area(sim_params, mpi_params);
  int nx = proc_area->end[0] - proc_area->start[0];
  int ny = proc_area->end[1] - proc_area->start[1];

  struct data ez, hx, hy;
  double val = 0.;
  if (init_data(&ez, "ez", nx, ny, sim_params->steps[0], sim_params->steps[1],
                proc_area->start[0] * sim_params->steps[0],
                proc_area->start[1] * sim_params->steps[1], val) ||
      init_data(&hx, "hx", nx, ny - 1, sim_params->steps[0],
                sim_params->steps[1],
                proc_area->start[0] * sim_params->steps[0],
                proc_area->start[1] * sim_params->steps[1], val) ||
      init_data(&hy, "hy", nx - 1, ny, sim_params->steps[0],
                sim_params->steps[1],
                proc_area->start[0] * sim_params->steps[0],
                proc_area->start[1] * sim_params->steps[1], val)) {
    printf("Error: could not allocate data\n");
    return EXIT_FAILURE;
  }

  double start = GET_TIME();

  DEBUG_PRINT("Starting MPI structs allocation\n");

  // Setting up MPI requests variables
  // Note the initialization to NULL to be able to free even if not malloc yet
  MPI_Request *send_requests = NULL;
  MPI_Request *recv_requests = NULL;
  double **sent_data = NULL;
  double **received_data = NULL;
  int *sizes = NULL;
  send_requests = malloc(sizeof(MPI_Request) * nb_neighbours);
  recv_requests = malloc(sizeof(MPI_Request) * nb_neighbours);
  sent_data = malloc(sizeof(double *) * nb_neighbours);
  received_data = malloc(sizeof(double *) * nb_neighbours);
  sizes = malloc(sizeof(int) * nb_neighbours);

  if (!send_requests || !recv_requests || !sent_data || !received_data ||
      !sizes) {
    printf("Error: allocation problem in initial MPI structs init\n");
    free_all_solve_pointers(nb_neighbours, proc_area, &hx, &hy, &ez,
                            send_requests, recv_requests, sent_data,
                            received_data, sizes);
    return EXIT_FAILURE;
  }

  for (int i = 0; i < nb_neighbours; i++) {
    sent_data[i] = NULL;
    received_data[i] = NULL;
  }
  for (int i = 0; i < nb_neighbours; i++) {
    // an x line is of length ny, inversly for an y line
    sizes[i] = (i < 2) ? ny : nx;
    // * 3 for hx, hy, ez
    sent_data[i] = malloc(sizeof(double) * sizes[i] * 3);
    received_data[i] = malloc(sizeof(double) * sizes[i] * 3);

    if (!sent_data[i] || !received_data[i]) {
      printf("Error: allocation problem in second MPI structs init\n");
      free_all_solve_pointers(nb_neighbours, proc_area, &hx, &hy, &ez,
                              send_requests, recv_requests, sent_data,
                              received_data, sizes);
      return EXIT_FAILURE;
    }
  }

  DEBUG_PRINT("Entering time loop\n");
  // Time loop
  for (int t = 0; t < sim_params->size_of_space[sim_params->ndim]; t++) {
    // Sending the data
    send_data_new(nb_neighbours, sent_data, sizes, &hx, &hy, &ez, mpi_params,
                  send_requests);
    receive_data_new(nb_neighbours, received_data, sizes, &hx, &hy, &ez,
                     mpi_params, recv_requests);

    // Receiving the data

    MPI_Waitall(nb_neighbours, send_requests, MPI_STATUSES_IGNORE);
    MPI_Waitall(nb_neighbours, recv_requests, MPI_STATUSES_IGNORE);

    if (t && (t % (sim_params->size_of_space[sim_params->ndim] / 10)) == 0) {
      double time_sofar = GET_TIME() - start;
      double eta =
          (sim_params->size_of_space[sim_params->ndim] - t) * time_sofar / t;
      printf("Computing time step %d/%d (ETA: %g seconds)     \r", t,
             sim_params->size_of_space[sim_params->ndim], eta);
      fflush(stdout);
    }

    double chy = sim_params->steps[sim_params->ndim] /
                 (sim_params->steps[1] * phys_params->mu);
    for (int j = 0; j < ny - 1; j++) {
      for (int i = 0; i < nx; i++) {
        double hx_ij =
            GET(&hx, i, j) - chy * (GET(&ez, i, j + 1) - GET(&ez, i, j));
        SET(&hx, i, j, hx_ij);
      }
    }

    // hy loop
    // chx = dt / (dx * mu)
    double chx = sim_params->steps[sim_params->ndim] /
                 (sim_params->steps[0] * phys_params->mu);
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx - 1; i++) {
        double hy_ij =
            GET(&hy, i, j) + chx * (GET(&ez, i + 1, j) - GET(&ez, i, j));
        SET(&hy, i, j, hy_ij);
      }
    }

    // update ez
    // cex = dt / (dx * epsilon)
    // cey = dt / (dy * epsilon)
    double cex = sim_params->steps[sim_params->ndim] /
                 (sim_params->steps[0] * phys_params->eps),
           cey = sim_params->steps[sim_params->ndim] /
                 (sim_params->steps[1] * phys_params->eps);
    for (int j = 1; j < ny - 1; j++) {
      for (int i = 1; i < nx - 1; i++) {
        double ez_ij = GET(&ez, i, j) +
                       cex * (GET(&hy, i, j) - GET(&hy, i - 1, j)) -
                       cey * (GET(&hx, i, j) - GET(&hx, i, j - 1));
        SET(&ez, i, j, ez_ij);
      }
    }

    // impose source
    int source_x = sim_params->size_of_space[0] / 2 - 2;
    int source_y = sim_params->size_of_space[1] / 2 - 2;
    double n = t * sim_params->steps[sim_params->ndim];
    if ((proc_area->start[0] <= source_x && proc_area->end[0] > source_x) &&
        (proc_area->start[1] <= source_y && proc_area->end[1] > source_y)) {
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
    }
  }
  if (mpi_params->rank == 0)
    write_manifest_vtk("ez", sim_params->steps[sim_params->ndim],
                       sim_params->size_of_space[sim_params->ndim],
                       sim_params->sampling_rate, mpi_params->num_ranks);

  double time = GET_TIME() - start;
  double MUps_per_sec = 1.e-6 * (double)nx * (double)ny *
                        (double)sim_params->size_of_space[sim_params->ndim] /
                        time;
  perf_data->time = time;
  perf_data->MUps_per_sec = MUps_per_sec;

  free_all_solve_pointers(nb_neighbours, proc_area, &hx, &hy, &ez,
                          send_requests, recv_requests, sent_data,
                          received_data, sizes);

  printf("\n");

  return EXIT_SUCCESS;
}
