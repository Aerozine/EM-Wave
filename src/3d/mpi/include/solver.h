#ifndef __SOLVER__
#define __SOLVER__

#include <time.h>

#include "params.h"

#if defined(_OPENMP)
#include <omp.h>
#define GET_TIME() (omp_get_wtime()) // wall time
#else
#define GET_TIME() ((double)clock() / CLOCKS_PER_SEC) // cpu time
#endif

int solve(struct SimulationParams *sim_params,
          struct PhysicalParams *phys_params, struct PerformanceData *perf_data,
          struct MpiParams *mpi_params);

#endif
