#ifndef __PARAMS__
#define __PARAMS__

#include <stdbool.h>
#include <mpi.h>

enum NeighbourPosition { X_START, X_END, Y_START, Y_END, Z_START, Z_END };

struct neighbour {
    int rank;
    enum NeighbourPosition pos;
};

struct SimulationParams {
    int sampling_rate, problem_id, ndim;
    int *size_of_space;
    double *steps;
};

struct PerformanceData {
    double time, MUps_per_sec;
};

struct PhysicalParams {
    double eps, mu;
};

struct MpiParams {
    bool use_mpi;
    int rank, num_ranks;
    int *procs_per_dim, *periods, *coords;
    struct neighbour *neighbours;
    MPI_Comm cart_comm;
    int *sizes;
    int **send_sizes;
    int *send_array_sizes;
};

int init_params(struct SimulationParams *sim_params,
                struct PerformanceData *perf_data,
                struct PhysicalParams *phys_params,
                struct MpiParams *mpi_params);

int set_params(struct SimulationParams *sim_params);

void free_params(struct SimulationParams *sim_params,
                 struct MpiParams *mpi_params);

#endif