#ifndef __DATA__
#define __DATA__

#ifndef RES_FOLDER
#define RES_FOLDER "data/"
#endif

#include <stdint.h>

#define GET(data, i, j, k)                                                     \
  ((data)->values[((data)->nxny) * (k) + (data)->nx * (j) + (i)])
#define SET(data, i, j, k, val)                                                \
  ((data)->values[((data)->nxny) * (k) + (data)->nx * (j) + (i)] = (val))

struct data {
  const char *name;
  int nx, ny, nz, nxny;
  float dx, dy, dz, origin_x, origin_y, origin_z;
  float *values;
};

int init_data(struct data *data, const char *name, int nx, int ny, int nz,
              float dx, float dy, float dz, float origin_x, float origin_y,
              float origin_z, float val);
void free_data(struct data *data);
int write_data_vtk(struct data *data, int step, int rank);
int write_manifest_vtk(const char *name, double dt, int nt, int sampling_rate,
                       int numranks);

#endif
