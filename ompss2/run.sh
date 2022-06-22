#!/bin/sh

inputcase=$1
nnodes=$2

if [ -z "$inputcase" ]; then
  inputcase="n2.b2"
fi

if [ -z "$nnodes" ]; then
  nnodes="2"
fi

rm core.miniMD.xeon0*

export I_MPI_DEBUG=5
export I_MPI_PMI_LIBRARY=/nix/store/1clca213zr12c4y096b905zry99gdc4c-slurm-16.05.8.1/lib/libpmi.so

# To run in xeon we need a matching slurm
export PATH=/nix/store/1clca213zr12c4y096b905zry99gdc4c-slurm-16.05.8.1/bin/srun:$PATH

echo "running: srun -N 2 -W 1 -K1 ./miniMD -i data/$inputcase/input.conf -r data/$inputcase/ref"
srun -N $nnodes -W 1 -K1 strace -ff -o strace ./miniMD -i data/$inputcase/input.conf -r data/$inputcase/ref
