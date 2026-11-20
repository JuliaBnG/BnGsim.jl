#!/bin/bash
# Peak memory and run time of the calibration (h*, p*, gblup*, ssgblup*) and
# held-out (v*) scenarios, two runs each, each in a fresh Julia process.
# Run on an otherwise idle host from this directory:  ./run.sh > measured.tsv
cd "$(dirname "$0")/scenarios"
for f in h500_500 h1000_1000 h2000_2000 h2000_4000 h3000_3000 p1000_5 p2000_10 p4000_10 \
         gblup_2000 gblup_4000 gblup_6000 gblup_8000 ssgblup_2000 ssgblup_4000 ssgblup_6000 ssgblup_8000 \
         v1 v2 v3 v4 v5 v6; do
  for rep in 1 2; do
    julia --startup-file=no -t 8 --project=../../.. ../measure.jl $f.json 2>&1 | tail -1
  done
done
