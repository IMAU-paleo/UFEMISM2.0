#! /bin/csh -f

# Compile all the UPSY code.
#
# Usage: ./compile_UPSY.csh  [VERSION]  [SELECTION]
#
#   [VERSION]: dev, perf
#
#     dev : developer's version; extra compiler flags (see Makefile_include, COMPILER_FLAGS_CHECK),
#           plus run-time assertions (i.e. DO_ASSERTIONS = true) and resource tracking (i.e.
#           DO_RESOURCE_TRACKING = true). Useful for tracking coding errors, but slow.
#
#     perf: performance version; all of the above disabled. Coding errors will much more often
#           simply result in segmentation faults, but runs much faster.
#
#   [SELECTION]: changed, clean
#
#     changed: (re)compile only changed modules. Works 99% of the time, fails when you change the
#              definition of a derived type without recompiling all of the modules thst use that
#              type. In that case, better to just do a clean compilation.
#
#     clean: recompile all modules. Always works, but slower.
#

# Safety
if ($#argv != 2) goto usage

echo ""

#  Confirm user's compilation choices
set version   = $argv[1]
if ($version == 'dev') then
  echo "dev: compiling UPSY, developers' version"
else if ($version == 'perf') then
  echo "perf: compiling UPSY, performance version"
else
  goto usage
endif

set selection = $argv[2]
if ($selection == 'changed') then
  echo "changed: (re)compiling changed modules only"
else if ($selection == 'clean') then
  echo "clean: recompiling all modules"
else
  goto usage
endif

echo ""

# If no build directory exists, create it
if (! -d build) mkdir build

# For a "clean" build, remove all build files first
if ($selection == 'clean') rm -rf build/*

# For a "changed" build, remove only the CMake cache file
if ($selection == 'changed') rm -f build/CMakeCache.txt

# Use CMake to build UPSY, with Ninja to determine module dependencies;
# use different compiler flags for the development/performance build
cd build

if ($version == 'dev') then

  cmake -G Ninja -DPETSC_DIR=`brew --prefix petsc` \
    -DDO_ASSERTIONS=ON \
    -DDO_RESOURCE_TRACKING=ON \
    -DEXTRA_Fortran_FLAGS="\
      -fdiagnostics-color=always;\
      -Og;\
      -Wall;\
      -ffree-line-length-none;\
      -cpp;\
      -Werror=implicit-interface;\
      -fimplicit-none;\
      -g;\
      -march=native;\
      -fcheck=all;\
      -fbacktrace;\
      -finit-real=nan;\
      -finit-integer=-42;\
      -finit-character=33" ..

else if ($version == 'perf') then

  cmake -G Ninja -DPETSC_DIR=`brew --prefix petsc` \
    -DDO_ASSERTIONS=OFF \
    -DDO_RESOURCE_TRACKING=OFF \
    -DEXTRA_Fortran_FLAGS="\
      -fdiagnostics-color=always;\
      -O3;\
      -Wall;\
      -ffree-line-length-none;\
      -cpp;\
      -fimplicit-none;\
      -g;\
      -march=native" ..

endif

ninja -v
cd ..

# Copy compiled program
if ($version == 'dev') then

  rm -f UPSY_unit_test_program_dev
  mv build/UPSY_unit_test_program UPSY_unit_test_program_dev
  rm -f UPSY_unit_test_program
  cp UPSY_unit_test_program_dev UPSY_unit_test_program

else if ($version == 'perf') then

  rm -f UPSY_unit_test_program_perf
  mv build/UPSY_unit_test_program UPSY_unit_test_program_perf
  rm -f UPSY_unit_test_program
  cp UPSY_unit_test_program_perf UPSY_unit_test_program

endif

exit 0

usage:

echo ""
echo "Usage: ./compile_UPSY.csh  [VERSION]  [SELECTION]"
echo ""
echo "  [VERSION]: dev, perf"
echo ""
echo "    dev : developer's version; extra compiler flags (see Makefile_include, COMPILER_FLAGS_CHECK),"
echo "          plus run-time assertions (i.e. DO_ASSERTIONS = true) and resource tracking (i.e."
echo "          DO_RESOURCE_TRACKING = true). Useful for tracking coding errors, but slow."
echo ""
echo "    perf: performance version; all of the above disabled. Coding errors will much more often"
echo "          simply result in segmentation faults, but runs much faster."
echo ""
echo "  [SELECTION]: changed, clean"
echo ""
echo "    changed: (re)compile only changed modules. Works 99% of the time, fails when you change the"
echo "             definition of a derived type without recompiling all of the modules thst use that"
echo "             type. In that case, better to just do a clean compilation."
echo ""
echo "    clean: recompile all modules. Always works, but slower."
echo ""

exit 1