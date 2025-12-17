#ifndef __SOLVER__
#define __SOLVER__

#include <time.h>

#include "params.h"

#define GET_TIME() ((double)clock() / CLOCKS_PER_SEC) // cpu time

int solve(struct SimulationParams *sim_params,
          struct PhysicalParams *phys_params, struct PerformanceData *perf_data);

#endif
