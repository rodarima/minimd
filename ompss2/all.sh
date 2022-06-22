#!/bin/sh

for testcase in data/*; do
  casename=$(basename "$testcase")
  nnodes=$(echo "$casename" | awk -F. '{print gensub("n", "", "g", $1)}')

  echo $casename
  ./run.sh "$casename" "$nnodes" 2>&1 | grep "^time"
done
