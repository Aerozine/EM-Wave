#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

omp_strong_close = np.loadtxt("omp_strong_close.txt", unpack=True)
omp_strong_spread = np.loadtxt("omp_strong_spread.txt", unpack=True)
omp_weak_close = np.loadtxt("omp_weak_close.txt", unpack=True)
omp_weak_spread = np.loadtxt("omp_weak_spread.txt", unpack=True)

mpi_strong_close = np.loadtxt("mpi_strong_close.txt", unpack=True)
mpi_strong_spread = np.loadtxt("mpi_strong_spread.txt", unpack=True)
mpi_weak_close = np.loadtxt("mpi_weak_close.txt", unpack=True)
mpi_weak_spread = np.loadtxt("mpi_weak_spread.txt", unpack=True)

base_weak = 2000
sizes_weak = base_weak * np.array([1, 2, 4, 8])
cores_weak = np.array([1, 4, 16, 64])

size_strong = 20000
cores_strong = np.array([1, 2, 4, 8, 16, 32, 64])


# =============
# OMP STRONG per core
# =============

figname = "graphs/OMP_strong_per_core"

fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(cores_strong, omp_strong_close / cores_strong, "-o", label="Close placement")
ax.plot(
    cores_strong, omp_strong_spread / cores_strong, "-o", label="Spread out placement"
)
ax.set_xlabel("Number of cores", fontsize=16)
ax.set_xscale("log")
ax.grid()
ax.set_xticks(cores_strong)
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_ylabel("MUpdates/seconds/core", fontsize=16)
ax.set_title("OpenMP - Strong scalability", fontsize=16)
ax.legend(fontsize=16)
ax.tick_params(axis="both", which="major", labelsize=12)
ax.tick_params(axis="both", which="minor", labelsize=12)
plt.tight_layout()
plt.savefig(figname + ".png")
plt.savefig(figname + ".svg")

# =============
# OMP STRONG
# =============

figname = "graphs/OMP_strong"

fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(cores_strong, omp_strong_close, "-o", label="Close placement")
ax.plot(cores_strong, omp_strong_spread, "-o", label="Spread out placement")
ax.set_xlabel("Number of cores", fontsize=16)
ax.set_xscale("log")
ax.grid()
ax.set_xticks(cores_strong)
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_ylabel("MUpdates/seconds", fontsize=16)
ax.set_title("OpenMP - Strong scalability", fontsize=16)
ax.legend(fontsize=16)
ax.tick_params(axis="both", which="major", labelsize=12)
ax.tick_params(axis="both", which="minor", labelsize=12)
plt.tight_layout()
plt.savefig(figname + ".png")
plt.savefig(figname + ".svg")

# =============
# OMP STRONG speedup
# =============

figname = "graphs/OMP_strong_speedup"

fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(
    cores_strong, omp_strong_close / omp_strong_close[0], "-o", label="Close placement"
)
ax.plot(
    cores_strong,
    omp_strong_spread / omp_strong_spread[0],
    "-o",
    label="Spread out placement",
)
ax.set_xlabel("Number of cores", fontsize=16)
ax.set_xscale("log")
ax.grid()
ax.set_xticks(cores_strong)
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_ylabel("Speedup", fontsize=16)
ax.set_title("Speedup considering strong scalability with OpenMP", fontsize=16)
ax.legend(fontsize=16)
ax.tick_params(axis="both", which="major", labelsize=12)
ax.tick_params(axis="both", which="minor", labelsize=12)
plt.tight_layout()
plt.savefig(figname + ".png")
plt.savefig(figname + ".svg")

# =============
# OMP STRONG efficiency
# =============

figname = "graphs/OMP_strong_efficiency"

fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(
    cores_strong,
    omp_strong_close / omp_strong_close[0] / cores_strong,
    "-o",
    label="Close placement",
)
ax.plot(
    cores_strong,
    omp_strong_spread / omp_strong_spread[0] / cores_strong,
    "-o",
    label="Spread out placement",
)
ax.set_xlabel("Number of cores", fontsize=16)
ax.set_xscale("log")
ax.grid()
ax.set_xticks(cores_strong)
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_ylabel("Efficiency", fontsize=16)
ax.set_title("Efficiency considering strong scalability with OpenMP", fontsize=16)
ax.legend(fontsize=16)
ax.tick_params(axis="both", which="major", labelsize=12)
ax.tick_params(axis="both", which="minor", labelsize=12)
plt.tight_layout()
plt.savefig(figname + ".png")
plt.savefig(figname + ".svg")

# =============
# OMP weak speedup
# =============

figname = "graphs/OMP_weak_speedup"

fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(cores_weak, omp_weak_close / omp_weak_close[0], "-o", label="Close placement")
ax.plot(
    cores_weak,
    omp_weak_spread / omp_weak_spread[0],
    "-o",
    label="Spread out placement",
)
ax.set_xlabel("Number of cores", fontsize=16)
ax.set_xscale("log")
ax.grid()
ax.set_xticks(cores_weak)
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_ylabel("Speedup", fontsize=16)
ax.set_title("Speedup considering weak scalability with OpenMP", fontsize=16)
ax.legend(fontsize=16)
ax.tick_params(axis="both", which="major", labelsize=12)
ax.tick_params(axis="both", which="minor", labelsize=12)
plt.tight_layout()
plt.savefig(figname + ".png")
plt.savefig(figname + ".svg")

# =============
# OMP WEAK efficiency
# =============

figname = "graphs/OMP_weak_efficiency"

fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(
    cores_weak,
    omp_weak_close / omp_weak_close[0] / cores_weak,
    "-o",
    label="Close placement",
)
ax.plot(
    cores_weak,
    omp_weak_spread / omp_weak_spread[0] / cores_weak,
    "-o",
    label="Spread out placement",
)
ax.set_xlabel("Number of cores", fontsize=16)
ax.set_xscale("log")
ax.grid()
ax.set_xticks(cores_weak)
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_ylabel("Efficiency", fontsize=16)
ax.set_title("Efficiency considering weak scalability with OpenMP", fontsize=16)
ax.legend(fontsize=16)
ax.tick_params(axis="both", which="major", labelsize=12)
ax.tick_params(axis="both", which="minor", labelsize=12)
plt.tight_layout()
plt.savefig(figname + ".png")
plt.savefig(figname + ".svg")

# =============
# OMP WEAK per cores
# =============

figname = "graphs/OMP_weak"

fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(cores_weak, omp_weak_close / cores_weak, "-o", label="Weak - close")
ax.plot(cores_weak, omp_weak_spread / cores_weak, "-o", label="Weak - spread")
ax.set_xlabel("Number of cores", fontsize=16)
ax.set_xscale("log")
ax.grid()
ax.set_xticks(cores_weak)
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_ylabel("MUpdates/seconds/core", fontsize=16)
ax.set_title("OpenMP - Weak scalability", fontsize=16)
ax.legend(fontsize=16)
ax.tick_params(axis="both", which="major", labelsize=12)
ax.tick_params(axis="both", which="minor", labelsize=12)
plt.tight_layout()
plt.savefig(figname + ".png")
plt.savefig(figname + ".svg")


# =============
# MPI STRONG per core
# =============

figname = "graphs/MPI_strong_per_core"

fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(cores_strong, mpi_strong_close / cores_strong, "-o", label="Close placement")
ax.plot(
    cores_strong, mpi_strong_spread / cores_strong, "-o", label="Spread out placement"
)
ax.set_xlabel("Number of cores", fontsize=16)
ax.set_xscale("log")
ax.grid()
ax.set_xticks(cores_strong)
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_ylabel("MUpdates/seconds/core", fontsize=16)
ax.set_title("MPI - Strong scalability", fontsize=16)
ax.legend(fontsize=16)
ax.tick_params(axis="both", which="major", labelsize=12)
ax.tick_params(axis="both", which="minor", labelsize=12)
plt.tight_layout()
plt.savefig(figname + ".png")
plt.savefig(figname + ".svg")

# =============
# MPI STRONG
# =============

figname = "graphs/MPI_strong"

fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(cores_strong, mpi_strong_close, "-o", label="Close placement")
ax.plot(cores_strong, mpi_strong_spread, "-o", label="Spread out placement")
ax.set_xlabel("Number of cores", fontsize=16)
ax.set_xscale("log")
ax.grid()
ax.set_xticks(cores_strong)
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_ylabel("MUpdates/seconds", fontsize=16)
ax.set_title("MPI - Strong scalability", fontsize=16)
ax.legend(fontsize=16)
ax.tick_params(axis="both", which="major", labelsize=12)
ax.tick_params(axis="both", which="minor", labelsize=12)
plt.tight_layout()
plt.savefig(figname + ".png")
plt.savefig(figname + ".svg")

# =============
# MPI STRONG speedup
# =============

figname = "graphs/MPI_strong_speedup"

fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(
    cores_strong, mpi_strong_close / mpi_strong_close[0], "-o", label="Close placement"
)
ax.plot(
    cores_strong,
    mpi_strong_spread / mpi_strong_spread[0],
    "-o",
    label="Spread out placement",
)
ax.set_xlabel("Number of cores", fontsize=16)
ax.set_xscale("log")
ax.grid()
ax.set_xticks(cores_strong)
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_ylabel("Speedup", fontsize=16)
ax.set_title("Speedup considering strong scalability with MPI", fontsize=16)
ax.legend(fontsize=16)
ax.tick_params(axis="both", which="major", labelsize=12)
ax.tick_params(axis="both", which="minor", labelsize=12)
plt.tight_layout()
plt.savefig(figname + ".png")
plt.savefig(figname + ".svg")

# =============
# MPI STRONG efficiency
# =============

figname = "graphs/MPI_strong_efficiency"

fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(
    cores_strong,
    mpi_strong_close / mpi_strong_close[0] / cores_strong,
    "-o",
    label="Close placement",
)
ax.plot(
    cores_strong,
    mpi_strong_spread / mpi_strong_spread[0] / cores_strong,
    "-o",
    label="Spread out placement",
)
ax.set_xlabel("Number of cores", fontsize=16)
ax.set_xscale("log")
ax.grid()
ax.set_xticks(cores_strong)
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_ylabel("Efficiency", fontsize=16)
ax.set_title("Efficiency considering strong scalability with MPI", fontsize=16)
ax.legend(fontsize=16)
ax.tick_params(axis="both", which="major", labelsize=12)
ax.tick_params(axis="both", which="minor", labelsize=12)
plt.tight_layout()
plt.savefig(figname + ".png")
plt.savefig(figname + ".svg")

# =============
# MPI weak speedup
# =============

figname = "graphs/MPI_weak_speedup"

fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(cores_weak, mpi_weak_close / mpi_weak_close[0], "-o", label="Close placement")
ax.plot(
    cores_weak,
    mpi_weak_spread / mpi_weak_spread[0],
    "-o",
    label="Spread out placement",
)
ax.set_xlabel("Number of cores", fontsize=16)
ax.set_xscale("log")
ax.grid()
ax.set_xticks(cores_weak)
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_ylabel("Speedup", fontsize=16)
ax.set_title("Speedup considering weak scalability with MPI", fontsize=16)
ax.legend(fontsize=16)
ax.tick_params(axis="both", which="major", labelsize=12)
ax.tick_params(axis="both", which="minor", labelsize=12)
plt.tight_layout()
plt.savefig(figname + ".png")
plt.savefig(figname + ".svg")

# =============
# MPI WEAK efficiency
# =============

figname = "graphs/MPI_weak_efficiency"

fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(
    cores_weak,
    mpi_weak_close / mpi_weak_close[0] / cores_weak,
    "-o",
    label="Close placement",
)
ax.plot(
    cores_weak,
    mpi_weak_spread / mpi_weak_spread[0] / cores_weak,
    "-o",
    label="Spread out placement",
)
ax.set_xlabel("Number of cores", fontsize=16)
ax.set_xscale("log")
ax.grid()
ax.set_xticks(cores_weak)
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_ylabel("Efficiency", fontsize=16)
ax.set_title("Efficiency considering weak scalability with MPI", fontsize=16)
ax.legend(fontsize=16)
ax.tick_params(axis="both", which="major", labelsize=12)
ax.tick_params(axis="both", which="minor", labelsize=12)
plt.tight_layout()
plt.savefig(figname + ".png")
plt.savefig(figname + ".svg")

# =============
# MPI WEAK per cores
# =============

figname = "graphs/MPI_weak"

fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(cores_weak, mpi_weak_close / cores_weak, "-o", label="Weak - close")
ax.plot(cores_weak, mpi_weak_spread / cores_weak, "-o", label="Weak - spread")
ax.set_xlabel("Number of cores", fontsize=16)
ax.set_xscale("log")
ax.grid()
ax.set_xticks(cores_weak)
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_ylabel("MUpdates/seconds/core", fontsize=16)
ax.set_title("MPI - Weak scalability", fontsize=16)
ax.legend(fontsize=16)
ax.tick_params(axis="both", which="major", labelsize=12)
ax.tick_params(axis="both", which="minor", labelsize=12)
plt.tight_layout()
plt.savefig(figname + ".png")
plt.savefig(figname + ".svg")


# =============
# Comparison STRONG
# =============

figname = "graphs/comparison_strong"

fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(cores_strong, omp_strong_close, "-o", label="OMP - close")
ax.plot(cores_strong, mpi_strong_close, "-o", label="MPI - close")
ax.plot(cores_strong, omp_strong_spread, "-o", label="OMP - spread")
ax.plot(cores_strong, mpi_strong_spread, "-o", label="MPI - spread")
ax.set_xlabel("Number of cores", fontsize=16)
ax.set_xscale("log")
ax.grid()
ax.set_xticks(cores_strong)
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_ylabel("MUpdates/seconds", fontsize=16)
ax.set_title("Comparison of the scalability with OMP and MPI", fontsize=16)
ax.legend(fontsize=16)
ax.tick_params(axis="both", which="major", labelsize=12)
ax.tick_params(axis="both", which="minor", labelsize=12)
plt.tight_layout()
plt.savefig(figname + ".png")
plt.savefig(figname + ".svg")
