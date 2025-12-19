# Propagation of EM waves

![](doc/emag_waves.png)

This is a project realised for the _High Performance Scientific Computing_ class
at the University of Liège.

The goal is to simulate, using a finite difference scheme, the propagation of
electromagnetic waves.

This implementation exists in several versions:

-   a serial one,
-   an openmp one,
-   a MPI one,
-   a CUDA one.

Both 2D and 3D simulations can be run.

## Building and running

The code can be compiled with the following commands :

```bash
DIM3=<dim> BUILD="<target>" make -j$(nprocs)
```

The `<dim>` option can be either `0` or `1`, where the former corresponds to building for a 2D simulation, and the latter is for a 3D one.

The possible targets are:

-   `reference`, for the serial build,
-   `openmp`, for the openmp build,
-   `mpi`, for the MPI build,
-   `gpu`, for the CUDA build.

All of those options are available for both the 2D and 3D cases. For the 2D case, there is an additional target, `stability`, corresponding to the stability analysis. To build it, a Python environment with numpy and matplotlib is required.

Running the simulation depends on the build. For:

-   `reference`, run `./hpc_project <id>`, where `<id>` is the problem id,
-   `openmp`, run `./hpc_project <id> [size] [threads]`, where `[size]` is the size of the sides of the spatial grid (which is a square), and `[threads]` is the number of threads to use,
-   `mpi`, run `mpirun -n <procs> hpc_project <id> [size]`, where `<procs>` is the number of processes to run.

Those commands should of course be adapted to the cluster.
