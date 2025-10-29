#include "solver.h"
#include "data.h"
#include "main.h"
#include "params.h"

#include <math.h>
#include <mpi.h>
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

// Yes, this function is too big.
// But it only does one thing, which is apply the finite difference scheme.
int solve(struct SimulationParams *sim_params,
          struct PhysicalParams *phys_params, struct PerformanceData *perf_data,
          struct MpiParams *mpi_params) {

  // Getting the relevant area
  struct area *computation_area = get_area(sim_params, mpi_params);
  int nx = computation_area->end[0] - computation_area->start[0];
  int ny = computation_area->end[1] - computation_area->start[1];
  DEBUG_PRINT("Process %d out of %d, has nx = %d and ny = %d\n\tstart_x = %d, "
              "start_y = %d, end_x = %d, end_y = %d\n",
              mpi_params->rank + 1, mpi_params->num_ranks, nx, ny,
              computation_area->start[0], computation_area->start[1],
              computation_area->end[0], computation_area->end[1]);

  // Initialization of data structures
  struct data ez, hx, hy;
  if (init_data(&ez, "ez", nx, ny, sim_params->steps[0], sim_params->steps[1],
                computation_area->start[0] * sim_params->steps[0],
                computation_area->start[1] * sim_params->steps[1], 0.) ||
      init_data(&hx, "hx", nx, ny - 1, sim_params->steps[0],
                sim_params->steps[1],
                computation_area->start[0] * sim_params->steps[0],
                computation_area->start[1] * sim_params->steps[1], 0.) ||
      init_data(&hy, "hy", nx - 1, ny, sim_params->steps[0],
                sim_params->steps[1],
                computation_area->start[0] * sim_params->steps[0],
                computation_area->start[1] * sim_params->steps[1], 0.)) {
    printf("Error: could not allocate data\n");
    return EXIT_FAILURE;
  }

  // Start of the algorithm
  double start = GET_TIME();
  DEBUG_PRINT("Starting the computation on process %d out of %d\n",
              mpi_params->rank + 1, mpi_params->num_ranks);

  // Setting up some things for the time loop
  MPI_Request requests[4];
  double *temp_x_start = NULL, *temp_x_end = NULL;
  double **received_data = NULL;
  int *counts = NULL;
  temp_x_start = malloc(sizeof(double) * ny);
  temp_x_end = malloc(sizeof(double) * ny);
  received_data = malloc(sizeof(double *) * 2 * sim_params->ndim);
  counts = malloc(sizeof(int) * 2 * sim_params->ndim);
  if (!temp_x_start || !temp_x_end || !received_data || !counts) {
    free(temp_x_start);
    free(temp_x_end);
    free(received_data);
    free(counts);
    free_area(computation_area);
    return EXIT_FAILURE;
  }

  for (int i = 0; i < 2 * sim_params->ndim; i++) {
    counts[i] = computation_area->end[i / 2] - computation_area->start[i / 2];
    received_data[i] = malloc(sizeof(double) * counts[i]);
    if (!received_data[i]) {
      for (int j = 0; j < i; j++) {
        free(received_data[j]);
      }
      free(received_data);
      free(temp_x_end);
      free(temp_x_start);
      free_area(computation_area);
      return EXIT_FAILURE;
    }
  }

  // Time loop
  for (int n = 0; n < sim_params->size_of_space[sim_params->ndim]; n++) {

    // Sending part
    // tags : 0 to xstart, 1 to xend, 2 to ystart, 3 to yend
    // First, we setup the buffers for X_START and X_END
    for (int i = 0; i < ny; i++) {
      temp_x_start[i] = GET(&ez, 0, i);
      temp_x_end[i] = GET(&hy, nx - 1, i);
    }

    // Then, we wend the data to the neighbours
    for (int i = 0; i < 2 * sim_params->ndim; i++) {
      struct neighbour current_neighbour = mpi_params->neighbours[i];
      if (current_neighbour.rank != -1) {
        switch (current_neighbour.pos) {
        case X_START:
          DEBUG_PRINT("Process %d sending to %d, tag %d\n",
                      mpi_params->rank + 1, current_neighbour.rank + 1, 0);
          MPI_Isend(temp_x_start, ny, MPI_DOUBLE, current_neighbour.rank, 0,
                    mpi_params->cart_comm, &requests[0]);
          break;
        case Y_START:
          DEBUG_PRINT("Process %d sending to %d, tag %d\n",
                      mpi_params->rank + 1, current_neighbour.rank + 1, 2);
          MPI_Isend(&(ez.values), nx, MPI_DOUBLE, current_neighbour.rank, 2,
                    mpi_params->cart_comm, &requests[1]);
          break;
        case X_END:
          DEBUG_PRINT("Process %d sending to %d, tag %d\n",
                      mpi_params->rank + 1, current_neighbour.rank + 1, 1);
          MPI_Isend(temp_x_end, ny, MPI_DOUBLE, current_neighbour.rank, 1,
                    mpi_params->cart_comm, &requests[2]);
          break;
        case Y_END:
          DEBUG_PRINT("Process %d sending to %d, tag %d\n",
                      mpi_params->rank + 1, current_neighbour.rank + 1, 3);
          MPI_Isend(&(hx.values), nx, MPI_DOUBLE, current_neighbour.rank, 3,
                    mpi_params->cart_comm, &requests[3]);
          break;
        default:
          break;
        }
      }
    }

    // Now, we get into the computation
    if (n && (n % (sim_params->size_of_space[sim_params->ndim] / 10)) == 0) {
      double time_sofar = GET_TIME() - start;
      double eta =
          (sim_params->size_of_space[sim_params->ndim] - n) * time_sofar / n;
      printf("Computing time step %d/%d (ETA: %g seconds)     \r", n,
             sim_params->size_of_space[sim_params->ndim], eta);
      fflush(stdout);
    }

    // hx loop
    // chy = dt / (dy * mu)
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
    int source_x = sim_params->size_of_space[0] / 2;
    int source_y = sim_params->size_of_space[1] / 2;
    double t = n * sim_params->steps[sim_params->ndim];
    if ((computation_area->start[0] <= source_x &&
         computation_area->end[0] > source_x) &&
        (computation_area->start[1] <= source_y &&
         computation_area->end[1] > source_y)) {
      // DEBUG_PRINT("Process %d applying source at %d, %d\n",
      //             mpi_params->rank + 1, source_x, source_y);
      switch (sim_params->problem_id) {
      case 1:
      case 2:
        // sinusoidal excitation at 2.4 GHz in the middle of the domain
        SET(&ez, source_x, source_y, sin(2. * M_PI * 2.4e9 * t));
        break;
      case 3:
        SET(&ez, source_x, source_y, sin(2. * M_PI * 2.4e9 * t));
        break;
      default:
        printf("Error: unknown source\n");
        break;
      }
    }

    // Now, we receive the data from the neighbours
    for (int i = 0; i < 2 * sim_params->ndim; i++) {
      struct neighbour current_neighbour = mpi_params->neighbours[i];
      if (current_neighbour.rank != -1) {
        int tag = (current_neighbour.pos % 2) ? current_neighbour.pos - 1
                                              : current_neighbour.pos + 1;
        int size = (tag < 2) ? ny : nx;
        DEBUG_PRINT("Process %d waiting for data from %d, tag %d\n",
                    mpi_params->rank + 1, current_neighbour.rank + 1, tag);
        MPI_Recv(received_data[i], size, MPI_DOUBLE, current_neighbour.rank,
                 tag, mpi_params->cart_comm, MPI_STATUS_IGNORE);
      }
    }

    // Update what needs to be updated
    for (int i = 0; i < 2 * sim_params->ndim; i++) {
      struct neighbour current_neighbour = mpi_params->neighbours[i];
      if (current_neighbour.rank != -1) {
        switch (current_neighbour.pos) {
        case Y_END:
          for (int j = 0; j < nx; j++) {
            double hx_y_end_i =
                GET(&hx, j, ny - 1) -
                chy * (received_data[Y_END][j] - GET(&ez, j, ny - 2));
            SET(&hx, nx - 1, j, hx_y_end_i);
          }
          break;
        default:
          break;
        }
      }
    }

    // output step data in VTK format
    if (sim_params->sampling_rate && !(n % sim_params->sampling_rate)) {
      write_data_vtk(&ez, n, mpi_params->rank);
      // write_data_vtk(&hx, n, 0);
      // write_data_vtk(&hy, n, 0);
    }
  }

  // write VTK manifest, linking to individual step data files
  if (mpi_params->rank == 0)
    write_manifest_vtk("ez", sim_params->steps[sim_params->ndim],
                       sim_params->size_of_space[sim_params->ndim],
                       sim_params->sampling_rate, mpi_params->num_ranks);
  // write_manifest_vtk("hx", sim_params->steps[sim_params->ndim],
  // sim_params->size_of_space[sim_params->ndim], sampling_rate, 1);
  // write_manifest_vtk("hy", sim_params->steps[sim_params->ndim],
  // sim_params->size_of_space[sim_params->ndim], sampling_rate, 1);

  double time = GET_TIME() - start;
  double MUps_per_sec = 1.e-6 * (double)nx * (double)ny *
                        (double)sim_params->size_of_space[sim_params->ndim] /
                        time;

  perf_data->time = time;
  perf_data->MUps_per_sec = MUps_per_sec;

  free_data(&ez);
  free_data(&hx);
  free_data(&hy);
  free_area(computation_area);

  free(temp_x_start);
  free(temp_x_end);

  printf("\n");

  return EXIT_SUCCESS;
}
