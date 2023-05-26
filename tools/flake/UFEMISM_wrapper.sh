#!/usr/bin/env sh
nproc=1
if [[ ! -z "$1" ]]; then nproc=$1; fi
mpirun -n $nproc $( dirname -- "${BASH_SOURCE[0]}" )/UFEMISM_program $2
