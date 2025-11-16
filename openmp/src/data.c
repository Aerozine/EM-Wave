#include "data.h"

int init_data(struct data *data, const char *name, int nx, int ny, double dx,
              double dy, float val) {
  data->name = name;
  data->nx = nx;
  data->ny = ny;
  data->dx = dx;
  data->dy = dy;
  data->values = (float *)malloc(nx * ny * sizeof(float));
  if (!data->values) {
    printf("Error: Could not allocate data\n");
    return 1;
  }
  for (int i = 0; i < nx * ny; i++)
    data->values[i] = val;
  return 0;
}

void free_data(struct data *data) { free(data->values); }

int write_data_vtk(struct data *data, int step, int rank) {
  char out[512];
  if (strlen(data->name) > 256) {
    printf("Error: data name too long for output VTK file\n");
    return 1;
  }
  sprintf(out, RES_FOLDER "%s_rank%d_%d.vti", data->name, rank, step);

  FILE *fp = fopen(out, "wb");
  if (!fp) {
    printf("Error: Could not open output VTK file '%s'\n", out);
    return 1;
  }

  uint64_t num_points = data->nx * data->ny;
  uint64_t num_bytes = num_points * sizeof(float);

  fprintf(fp,
          "<?xml version=\"1.0\"?>\n"
          "<VTKFile"
          " type=\"ImageData\""
          " version=\"1.0\""
          " byte_order=\"LittleEndian\""
          " header_type=\"UInt64\""
          ">\n"
          "  <ImageData"
          " WholeExtent=\"0 %d 0 %d 0 %d\""
          " Spacing=\"%lf %lf %lf\""
          " Origin=\"%lf %lf %lf\""
          ">\n"
          "    <Piece Extent=\"0 %d 0 %d 0 %d\">\n"
          "      <PointData Scalars=\"scalar_data\">\n"
          "        <DataArray"
          " type=\"Float32\""
          " Name=\"%s\""
          " format=\"appended\""
          " offset=\"0\""
          ">\n"
          "        </DataArray>\n"
          "      </PointData>\n"
          "    </Piece>\n"
          "  </ImageData>\n"
          "  <AppendedData encoding=\"raw\">\n_",
          data->nx - 1, data->ny - 1, 0, data->dx, data->dy, 0., 0., 0., 0.,
          data->nx - 1, data->ny - 1, 0, data->name);

  fwrite(&num_bytes, sizeof(uint64_t), 1, fp);
  fwrite(data->values, sizeof(float), num_points, fp);

  fprintf(fp, "  </AppendedData>\n"
              "</VTKFile>\n");

  fclose(fp);

  return 0;
}

int write_manifest_vtk(const char *name, double dt, int nt, int sampling_rate,
                       int numranks) {
  char out[512];
  if (strlen(name) > 256) {
    printf("Error: name too long for Paraview manifest file\n");
    return 1;
  }
  sprintf(out, "%s.pvd", name);

  FILE *fp = fopen(out, "wb");
  if (!fp) {
    printf("Error: Could not open output VTK manifest file '%s'\n", out);
    return 1;
  }

  fprintf(fp, "<VTKFile"
              " type=\"Collection\""
              " version=\"0.1\""
              " byte_order=\"LittleEndian\">\n"
              "  <Collection>\n");

  for (int n = 0; n < nt; n++) {
    if (sampling_rate && !(n % sampling_rate)) {
      double t = n * dt;
      for (int rank = 0; rank < numranks; rank++) {
        fprintf(fp,
                "    <DataSet"
                " timestep=\"%g\""
                " part=\"%d\""
                " file='" RES_FOLDER "%s_rank%d_%d.vti'/>\n",
                t, rank, name, rank, n);
      }
    }
  }

  fprintf(fp, "  </Collection>\n"
              "</VTKFile>\n");
  fclose(fp);
  return 0;
}
