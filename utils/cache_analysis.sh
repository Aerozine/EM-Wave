#!/bin/bash

perf stat -e cache-references,cache-misses,L1-dcache-load-misses,L1-dcache-loads,LLC-load-misses,LLC-loads $@