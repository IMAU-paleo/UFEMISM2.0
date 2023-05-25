{
  inputs.nixpkgs.url = "github:nixos/nixpkgs/nixos-22.11";

  outputs = { self, nixpkgs }:
  let
    pkgs = nixpkgs.legacyPackages.x86_64-linux;
    mypython = pkgs.python310;
  in with pkgs; with mypython.pkgs; {
    devShell.x86_64-linux = mkShell { 
      buildInputs = [
          xarray
          netcdf4
          matplotlib
          gcc
          gfortran
          openmpi
          netcdf
          netcdffortran
          petsc
          lapack
          pkg-config
          valgrind
          gdb
          octaveFull
          octavePackages.netcdf
          ncview
      ];
    };
    packages.x86_64-linux.UFEMISM = 
    stdenv.mkDerivation {
       name = "UFEMISM";
       src = self;
       sourceRoot = "source/src";
       buildInputs = [ gfortran
                       openmpi
                       netcdf4
                       netcdffortran
                       petsc
                       lapack
                       pkg-config ];
       installPhase = "
         mkdir -p $out/bin; 
         install -t $out/bin UFEMISM_program";
    };


    # nix run <loc>#UFEMISM_program
    apps.x86_64-linux.UFEMISM_program = {
    type = "app";
    program = "${self.packages.x86_64-linux.UFEMISM}/bin/UFEMISM_program";
    };


    # Default nix run
    apps.x86_64-linux.default = self.apps.x86_64-linux.UFEMISM_program;

  };

}
