#!/bin/bash
# Peak memory of the scenarios for the terms added after the calibration
# (several histories, census trajectories, repeated records, reliability),
# two runs each in a fresh Julia process:  ./run.sh > measured.tsv
cd "$(dirname "$0")"
for f in l_hist4 l_growth l_records l_reliability; do
  for rep in 1 2; do
    julia --startup-file=no -t 8 --project=../../.. ../measure_later.jl $f.json 2>&1 | tail -1
  done
done
