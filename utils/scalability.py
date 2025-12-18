#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

omp_strong_close = np.loadtxt("omp_strong_close.txt", unpack=True)
omp_strong_spread = np.loadtxt("omp_strong_spread.txt", unpack=True)
omp_weak_close = np.loadtxt("omp_weak_close.txt", unpack=True)
omp_weak_spread = np.loadtxt("omp_weak_spread.txt", unpack=True)

base_weak = 2000
sizes_weak = base_weak * np.array([1, 2, 4, 8])
cores_weak = np.array([1, 4, 16, 64])

size_strong = 20000
cores_strong = np.array([1, 2, 4, 8, 16, 32, 64])


# =============
# OMP STRONG
# =============

figname = "graphs/OMP_strong"

fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(cores_strong, omp_strong_close, "-o", label="Strong - close")
ax.plot(cores_strong, omp_strong_spread, "-o", label="Strong - spread")
ax.set_xlabel("Number of cores", fontsize=16)
# ax.set_xscale("log")
ax.grid()
ax.set_xticks(cores_strong)
ax.set_ylabel("MUpdates/seconds", fontsize=16)
ax.legend(fontsize=16)
ax.tick_params(axis="both", which="major", labelsize=12)
ax.tick_params(axis="both", which="minor", labelsize=12)
plt.tight_layout()
plt.savefig(figname + ".png")
plt.savefig(figname + ".svg")

# =============
# OMP WEAK
# =============

figname = "graphs/OMP_weak"

fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(cores_weak, omp_weak_close / cores_weak, "-o", label="Weak - close")
ax.plot(cores_weak, omp_weak_spread / cores_weak, "-o", label="Weak - spread")
ax.set_xlabel("Number of cores", fontsize=16)
# ax.set_xscale("log")
ax.grid()
ax.set_xticks(cores_weak)
ax.set_ylabel("MUpdates/seconds/cores", fontsize=16)
ax.legend(fontsize=16)
ax.tick_params(axis="both", which="major", labelsize=12)
ax.tick_params(axis="both", which="minor", labelsize=12)
plt.tight_layout()
plt.savefig(figname + ".png")
plt.savefig(figname + ".svg")
