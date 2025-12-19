# MPI - Algorithm

This will explain the use of MPI in this project.

## main.c

First, MPI is initialized, as well as a `MpiParams` structure. It will contain the parameters relating to the MPI side. It then calls `solve`.

## solve.c

### `solve`

The entry point here is the function `solve`. The first thing it does is get the computation area for the particular process, via the function `get_area`.

Afterwards, it computes the number of nodes according to each direction of space, and stores them in `mpi_params->sizes`, as well as the size of the arrays to send for each direction of space (stored in `mpi_params->send_sizes`). It also determines the origin of the computation domain in the total axes (stored in `origins`), to be able to write the VTI files afterwards.

The data structures `ez`, `hx`, and `hy` are then created. They all are $nx * ny$, to make sending and receiving data easier.

Then, all of the necessary arrays for sending and receiving data are created. The send and receive MPI requests are in `send_requests` and `recv_requests`. The sent data is stored in `sent_data`. The first level is the neighbour (represented by NeighbourPosition), and the second is the data. The data is an array of size $nx$ if the neighbour is in direction $y$, and $ny$ in direction $x$. The arrays are stored such that the array corresponding to a certain neighbour only stores the data needed by/from that neighbour. Thus, the data arrays are of size $nx$ or $ny$. An array of boolean `received_neighbour` is used to track if data has been received from a given neighbour or not.

The data is sent and received as soon as the computation is done, and the waits are right before the data is needed. It is all asynchronous, to allow for better performance.

The function then prints some information on the advancement of the computation. Afterwards, it calls the three big functions: `hx_loop`, `hy_loop`, and `ez_loop`. Those are also defined as `inline` functions, to not loose too much performance.

We need to exchange data between the `hx`, `hy` loops and `ez` loop. This is because `ez` uses the recently computed magnetic field, but, if no data is exchanged, it'll use the old ones at the borders, causing errors and instability in the scheme. Thus, we at least need to exchange the magnetic fields before computing `ez`.

### `get_area`

The function `get_area` does two main things: determine the relevant computational area, and get the neighbouring processes in each direction of space.

To get the position in the decomposition of the domain, MPI functions are used, namely `MPI_Dims_create`, `MPI_Cart_create`, and `MPI_Cart_coords`. The number of processes in each direction of space is stored in `mpi_params->procs_per_dim`. The communication object is `mpi_params->cart_comm`, and the coordinates of the current process are in `mpi_params->coords`.

The start and end of the processe's computation domain are then determined, taking into account a possible remainder in the division. The start and end coordinates in each direction of space get stored in `area->start` and `area->end`. They are both inclusive.

The neighbours are then determined, using `MPI_Cart_shift`. They are stored in `mpi_params->neighbours`, as a `neighbour` structure. The position of neighbours is tracked with a `NeighbourPosition` enum. For direction $i$ of space, the neighbour at the start is denoted with $2 * i$, and the one at the end with $2 * i + 1$. It is also possible to use `X_START`, `X_END`, `Y_START`, and `Y_END` for more clarity. The rank of the neighbour is then stored at `mpi_params->neighbours[<position>].rank`.

### Send data & receive data

Sends & receives data from given struct.

Ez :

-   send to Y_START
-   receive from Y_END
-   idem for X

Hx:

-   send to Y_END
-   receive from Y_START

Hy:

-   send to X_END
-   receive from X_START

For each neighbour, you just have to send one data. => send, receive optimized.

### `hx_loop`

This function first computes a mathematical parameter for the scheme. The main loop does the computation on $[0, nx[ \times [0, ny - 1[$ as the scheme necessitates a value of `ez` present the following column ($y$ direction). It is then possible to compute the $ny - 1$ column only if data has been received from `Y_END`. In that case, it makes one additional loop on $x$ to fill out that column.

### `hy_loop`

This does basically the same thing as `hx_loop`, except with $x$ and $y$ inverted. The relevant neighbour is then `X_END`.

### `ez_loop`

First, the two physical parameters are computed. Afterwards, the main computation double loop is executed. It starts at $1$ for each coordinate, as the scheme needs the $h_x$ and $h_y$ of the previous nodes in both directions of space.

If `X_START` is received, then a loop on $y$ can be done for $x=0$. It must still start at $y = 1$. If `Y_START` has also been received, then $(0, 0)$ can be updated. There is another loop for `Y_START`.
