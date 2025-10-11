#set page(paper: "a4")
#set par(first-line-indent: (amount: 10pt, all: true), justify: true)

#align(center,
text(16pt)[*Question 2 - Analysis of the arithmetic intensity*]
)

#v(10pt)

= Introduction

This analysis will limit itself to the computation of the solution, and more specifically to the looping parts. The justification of such approximation is that those loops are the most demanding part of the algorithm by far.

= General structure

The algorithm consists of a top-level #raw("for") loop, which is responsible for the time evolution. The number of iterations here is equal to the number of time steps $n_t$.

The algorithm then consists in three times the same structure, which is a double #raw("for") loop. They are responsible for updating $H_x$, $H_y$, and $E_z$ respectively. For each of those structures, a coefficient is computed before entering the double loop.

After those loops, the source is imposed. This part does not contain any loop, as it only consists in setting $E_z$ to a certain value. This constitutes the end of the time loop.

= Analysis of a double #raw("for") loop of $H$

As the computation scheme is a little different for $E_z$, the loops for $H_x$ and $H_y$ will first be considered. They are extremely similar.

Before the loops start, a coefficient $(Delta t) / (mu Delta q)$ is computed, where $q = x$ or $y$. This operation consists in a multiplication and a division, with access to three variables from data structures. Even as this coefficient does not change between time iterations, recomputing it each time can be beneficial, as it is possible that the compiler makes sure that it is kept in the cache for the duration of the two #raw("for") loops. This could only be verified by inspecting the assembly code. In conclusion, this operation consists in $2$ flops, and $2$ memory access, as two of the needed variables are part of the same structure.

Let us now consider the two #raw("for") loops. The number of iterations of the loop linked to $H_q$ is the number of grid points along $q$, hence $n_q$. If $r$ is the other direction of space (excluding $z$), then the remaining loop does $n_r - 1$ iterations. The $-1$ comes from the fact that forward derivation is used, and that the $H$ and $E$ grids are staggered. The total number of iterations is then $n_q (n_r - 1)$. If $n_r$ is big - which is a reasonable approximation, as a supercomputer is used here - then it can be considered that $n_q n_r$ iterations are done.

Finally, let us consider the operation done at each iteration. Three different values are needed : $H_q (i, j)$, $E_z (i + 1, j)$, and $E_z (i, j)$. We can assume that the two $E_z$ values are obtained with the same cache line, making that only one memory access. Thus, if those values are not in cache, two memory accesses must be performed.

The algorithm performs one subtraction, one multiplication, and finally one addition. It then sets the new value. This does not cause any more memory access, as the variable set is the $H_q (i, j)$ value needed for the calculation. Three flops are then needed at each iteration.

To determine whether the memory accesses or the flops are the most demanding operation, one must take into account the memory bandwidth and the processing power. On a machine with $2.8$ TFLOPS/s and $200$ GB/S of memory bandwidth, one can compare the ratio between the two with the 3/2 ratio that theoretically saturates both. Here, the ratio obtained is 14. It can be concluded that a big performance limit comes from the low memory bandwidth.

This calculation isn't right, as less memory accesses are needed, given the cache. On NIC5, the L1 data cache is $32$ KB, the L2 is $512$ KB, and the L3 is 32MB per 8 cores. The L3 cache is thus 128MB for 32 cores. This can theoretically store $128 * 1024^2 / 8 = 16 med 777 med 216$ doubles. If this is equally allocated between the two data accesses, it reduces the need for. 

= Questions

- Does the 200 GB/S of memory bandwidth take into account the cache ? Does it need to, as it is the limiting factor ? Or do we not take into account optimizations (requesting data a bit before needing it) ?
- How to take into account the cache line ?
- Can I just compute the ratio between the flops and memory accesses and then generalize it ?
