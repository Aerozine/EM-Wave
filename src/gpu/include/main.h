#ifndef __MAIN_HPC__
#define __MAIN_HPC__

struct SimulationParams {
  real dx, dy, dt;
  int nx, ny, nt, sampling_rate;
};

struct PhysicalParams {
  real eps, mu;
};

#ifdef DEBUG
#define DEBUG_PRINT(...) fprintf(stderr, __VA_ARGS__);
#else
#define DEBUG_PRINT(...) (void(0));
#endif

#endif
