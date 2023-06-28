# UFEMISM2.0
Version 2.0 of the Utrecht FinitE voluMe Ice-Sheet Model

## Running with nix

If you have access to a linux or macos machine with the nix package manager. UFEMISM can be run with:

```
nix --extra-experimental-features flakes --extra-experimental-features nix-command run github:IMAU-paleo/UFEMISM2.0/main?dir=tools/flake 2 <config_file>
```

In this case `2` specifies the number of processors. This is a bit different then the normal syntax. Because
normally UFEMISM is run with `mpirun -n 2` prefixed instead of `2` being part of the UFEMISM command line. 

`main` can be replaced with any tag, branch, or commit hash, to run the code from that point in the repository.
