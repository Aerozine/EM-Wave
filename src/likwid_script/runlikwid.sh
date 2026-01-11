#!/bin/bash

# -------- configuration --------
CORE=0
BIN=./hpc_project
ARGS="1"
# --------------------------------

groups=(
  L2CACHE
  FLOPS_DP
  NUMA
  ICACHE
  CACHE
  MEMREAD
  DIVIDE
  L3
  DATA
  BRANCH
  L2
  FLOPS_SP
  ENERGY
  L3CACHE
  MEM
  CLOCK
  TLB
  MEMWRITE
  CPI
)

for g in "${groups[@]}"; do
  echo "Running LIKWID group: $g"
  likwid-perfctr -f -C $CORE -g $g $BIN $ARGS > reference/result_$g
done
