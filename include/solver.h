#ifndef __SOLVER__
#define __SOLVER__

#include "main.h"
#include <time.h>

#define GET_TIME() ((double)clock() / CLOCKS_PER_SEC) // cpu time

int solve(struct SimulationParams *sim_params,
          struct PhysicalParams *phys_params, int problem_id,
          struct PerformanceData *perf_data);

#endif
