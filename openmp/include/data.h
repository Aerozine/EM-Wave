#ifndef __DATA__
#define __DATA__

#ifndef RES_FOLDER
#define RES_FOLDER "data/"
#endif

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define GET(data, i, j) ((data)->values[(data)->nx * (j) + (i)])
#define SET(data, i, j, val) ((data)->values[(data)->nx * (j) + (i)] = (val))

struct data {
  const char *name;
  int nx, ny;
  double dx, dy;
  float *values;
};

int init_data(struct data *data, const char *name, int nx, int ny, double dx,
              double dy, float val);
void free_data(struct data *data);
int write_data_vtk(struct data *data, int step, int rank);
int write_manifest_vtk(const char *name, double dt, int nt, int sampling_rate,
                       int numranks);

#endif
