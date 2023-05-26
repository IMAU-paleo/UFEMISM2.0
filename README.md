# UFEMISM2.0
Version 2.0 of the Utrecht FinitE voluMe Ice-Sheet Model

## Running with nix

If you have access to a linux machine with the nix package manager. UFEMISM can be run with:

```
nix run github:IMAU-paleo/UFEMISM2.0?dir=tools/flake 2 <config_file>
```

where `2` specifies the number of processors. This is a bit different then the normal syntax. Because
normally UFEMISM is run with `mpirun -n 2` prefixed instead of `2` being part of the UFEMISM command line. 
