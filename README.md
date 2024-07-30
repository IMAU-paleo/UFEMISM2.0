# UFEMISM2.0
Version 2.0 of the Utrecht FinitE voluMe Ice-Sheet Model

## Quick Start

This assumes a Linux or MacOS system with a working fortran compile chain
and installed PeTSC, NetCDF and MPI (e.g. OpenMPI) libraries.

Example for bash:

```bash
  # 1. Copy compilation and run files from templates/ to UFEMISM2.0/
  cp templates/compile_all.sh .
  cp templates/run_UFEMISM.sh .

  # 2. Modify src/Makefile_include_local.txt to your local settings and
  #    compilation preferences    

  # 3. Modify src/Makefile so it points to Makefile_include_local.txt

  # 4. Compile the model
  ./compile_all.sh

  # 5. Run a test simulation
  ./run_UFEMISM.sh
```

## Running with nix

If you have access to a linux or macos machine with the nix package manager. UFEMISM can be run with:

```
nix --extra-experimental-features flakes --extra-experimental-features nix-command run github:IMAU-paleo/UFEMISM2.0/main?dir=tools/flake 2 <config_file>
```

In this case `2` specifies the number of processors. This is a bit different then the normal syntax. Because
normally UFEMISM is run with `mpirun -n 2` prefixed instead of `2` being part of the UFEMISM command line. 

`main` can be replaced with any tag, branch, or commit hash, to run the code from that point in the repository.
