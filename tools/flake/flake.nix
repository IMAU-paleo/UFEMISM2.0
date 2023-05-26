{
  description = "distributed, mesh based, ice sheet simulation software";
  inputs.nixpkgs.url = "github:nixos/nixpkgs/nixos-unstable";
  inputs.flake-utils.url = "github:numtide/flake-utils";

  outputs = { self, nixpkgs, flake-utils }:
  flake-utils.lib.eachDefaultSystem (system:
    let
      pkgs = nixpkgs.legacyPackages.${system};
      mypython = pkgs.python310;
    in with pkgs; with mypython.pkgs; {
      devShell = mkShell { 
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
        
      packages = rec {
        UFEMISM = 
        stdenv.mkDerivation {
           name = "UFEMISM";
           src = ../..;
           buildInputs = [ pkgs.gfortran
                           pkgs.makeWrapper
                           pkgs.openmpi
                           pkgs.netcdf
                           pkgs.netcdffortran
                           pkgs.hdf5
                           pkgs.curl
                           pkgs.petsc
                           pkgs.lapack
                           pkgs.pkg-config ]
                           ++ lib.optionals (stdenv.isDarwin) [ darwin.CF ];
           buildPhase = "
             cp  tools/flake/Makefile_include_nix.txt src/Makefile_include_local.txt
             cd src
             make
           ";

           installPhase = "
             mkdir -p $out/bin; 
             install -t $out/bin UFEMISM_program
             install -t $out/bin ../tools/flake/UFEMISM_wrapper.sh";
           postFixup =
             let
               dependency_path = pkgs.lib.makeBinPath [ openmpi ];
             in
             ''
               wrapProgram "$out/bin/UFEMISM_wrapper.sh" --prefix PATH : "${dependency_path}"
             '';
        };
        default = UFEMISM;
      };



      apps = rec {
        # nix run <loc>#UFEMISM_program
        UFEMISM_program = {
          type = "app";
          program = "${self.packages.${system}.UFEMISM}/bin/UFEMISM_program";
        };
        # nix run <loc>#UFEMISM_wrapper
        UFEMISM_wrapper = {
          type = "app";
          program = "${self.packages.${system}.UFEMISM}/bin/UFEMISM_wrapper.sh";
        };
        # Default nix run
        default = UFEMISM_wrapper;
      };

    });
}
