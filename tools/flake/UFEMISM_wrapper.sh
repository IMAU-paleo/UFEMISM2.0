#!/usr/bin/env sh
mpirun -n $1 $( dirname -- "${BASH_SOURCE[0]}" )/UFEMISM_program $2
