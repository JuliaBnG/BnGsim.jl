#!/bin/bash
# Run time by phase (0.37 time model): every scenario of scenarios/,
# later/, phases/ and check/ and the examples, or the scenarios given as
# arguments, two processes each (each process runs replicates 1 and 2,
# the first cold), at 8 threads. Run on an otherwise idle host from this
# directory:  ./run_phases.sh > phases.tsv
cd "$(dirname "$0")"
files=("$@")
[ ${#files[@]} -eq 0 ] && files=(scenarios/*.json later/*.json phases/*.json check/*.json ../../examples/*.json)
for f in "${files[@]}"; do
  for run in 1 2; do
    julia --startup-file=no -t 8 --project=../.. measure_phases.jl "$f" 2>&1 | grep -v '^\s*$' | sed "s|^|$run\t|"
  done
done
