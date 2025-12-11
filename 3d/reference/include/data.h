#ifndef __DATA__
#define __DATA__

#ifndef RES_FOLDER
#define RES_FOLDER "data/"
#endif

#include <stdint.h>

#define GET(data, i, j, k)                                                     \
  ((data)->values[((data)->nx * (data)->ny) * (k) + (data)->nx * (j) + (i)])
#define SET(data, i, j, k, val)                                                \
  ((data)->values[((data)->nx * (data)->ny) * (k) + (data)->nx * (j) + (i)] =    \
       (val))


struct data {
  const char *name;
  int nx, ny, nz;
  double dx, dy, dz;
  float *values;
};

int init_data(struct data *data, const char *name, int nx, int ny, int nz, double dx,
              double dy, double dz, float val);
void free_data(struct data *data);
int write_data_vtk(struct data *data, int step, int rank);
int write_manifest_vtk(const char *name, double dt, int nt, int sampling_rate,
                       int numranks);

#endif
