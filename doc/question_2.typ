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

Finally, let us consider the operation done at each iteration. Three different values are needed : $H_q (i, j)$, $E_z (i + 1, j)$, and $E_z (i, j)$. Thus, three variables must be accessed at each iteration. One can argue that $E_z (i + 1, j)$ is the $E_z (i, j)$ of the next iteration. It can then be kept in cache, reducing the number of accesses to the main memory. It can be concluded that two memory accesses are necessary at each iteration.

The algorithm performs one subtraction, one multiplication, and finally one addition. It then sets the new value. This does not cause any more memory access, as the variable set is the $H_q (i, j)$ value needed for the calculation. Three flops are then needed at each iteration.

Finally, it can be seen that, for each iteration, $3$ floating point operations are executed, and two memory accesses are performed.

= Application to a real system

Let us consider a computer which as the following specifications : $2.8$ TFLOPS/s of computing power, and a memory bandwidth of $200$ GB/s. Converting into flops, we see that it can execute $2.8 dot.op 10^(12)$ floating point operations each second. The number of bytes accessible each second is equal to $2 dot.op 10^(11)$.

Let us now consider that the variables are stored as doubles. This means that they each take up $8$ bytes. We then conclude that $(2 dot.op 10^(11)) / 8 = 25 dot.op 10^9$ variables can be accessed from memory each second. If we compute the ratio between the flops/s and variables accessed per second, we get $(2.8 dot.op 10^(12)) / (25 dot.op 10^9) = 112$. However, it was determined at the previous section that the ratio allowing for the saturation of both computing power and memory bandwidth was $3/2$. We then conclude that the main limitation here is the memory bandwidth.

If we now assume that the variables are stored as floats and not doubles, we can calculate that $50 dot 10^9$ variables can be accessed each second. The ratio then becomes equal to $56$. We see that there is a still a great limitation coming from the memory bandwidth, but an improvement in performance of two can be achieved.

= Implications of the cache

In the sections above, the cache was only mentioned as being able to store a value of $E_z$ between two iterations. In fact, as the memory is accessed in chunks of size equal to the cache line, the cache will store quite a bit more than that. However, this will not influence the calculations above.

Let us first consider the case where only the $E_z$ value is stored in cache. Such a behaviour could represent cache trashing. This will heavily impact the performance, not because of the memory bandwidth, but because of the latency inherent to memory accesses. The CPU will spend its time waiting for data to arrive, do a ridiculously low number of operations on it, then wait again.

This is not representative of the real world. Two things, when taken together, will heavily improve performance : the cache, and memory prefetching. The CPU will only operate on data found in cache, and hardware (and software, as we compile with #raw("gcc -O3")) prefetchers will make sure that the needed data is available on time. This will of course only mitigate the memory access latency. The number of variables available from memory each second is constant and determined by the memory bandwidth. The CPU will still end up waiting for data to arrive, but the cause won't be latency.
