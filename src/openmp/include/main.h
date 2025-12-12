#ifndef __MAIN_HPC__
#define __MAIN_HPC__

struct SimulationParams {
  double dx, dy, dt;
  int nx, ny, nt, sampling_rate;
};

struct PhysicalParams {
  double eps, mu;
};

#endif
