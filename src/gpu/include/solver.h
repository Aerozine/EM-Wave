#ifndef __SOLVER__
#define __SOLVER__

#include "main.h"
#include <time.h>

#define GET_TIME() ((real)clock() / CLOCKS_PER_SEC) // cpu time

int solve(struct SimulationParams *sim_params,
          struct PhysicalParams *phys_params, int problem_id);

// from cuda doc to help dbg
#define CUDA_CHECK(call)                                                       \
  do {                                                                         \
    cudaError_t err = call;                                                    \
    if (err != cudaSuccess) {                                                  \
      fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(err),     \
              __FILE__, __LINE__);                                             \
      exit(EXIT_FAILURE);                                                      \
    }                                                                          \
  } while (0)

#endif
