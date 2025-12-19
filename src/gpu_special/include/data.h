#ifndef __DATA__
#define __DATA__
#ifndef RES_FOLDER
#define RES_FOLDER "data/"
#endif
#define PRECISION 1 // Use 0 for float, 1 for double
#if PRECISION == 0
typedef float real;
#else
typedef double real;
#endif
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// Row-major indexing for coalesced access
#define PIDX(i, j, pitch) ((j) * ((pitch) / sizeof(real)) + (i))

struct data {
  const char *name;
  int nx, ny;
  real dx, dy, *values;
  size_t pitch;
};
int init_data(struct data *data, const char *name, int nx, int ny, real dx,
              real dy, real val);
void free_data(struct data *data);
int write_data_vtk(struct data *data, int step, int rank);
int write_manifest_vtk(const char *name, real dt, int nt, int sampling_rate,
                       int numranks);

#endif
